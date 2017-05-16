import json
import os
from math import floor, sqrt
from datetime import datetime, timezone

import pandas as pd

from .user import Instrument, Site
from .util import MixedModelFit, weighted_mean

__author__ = "Eric Dose :: New Mexico Mira Project, Albuquerque"

AN_TOP_DIRECTORY = 'J:/Astro/Images/C14'

#######
#
#  Usage:
#    precheck(an_rel_directory='20170509', instrument_name='Borea')


def precheck(an_rel_directory=None, instrument_name='Borea'):
    """
    First check on new FITS files, at beginning of processing. Prints exceptions and summary stats.
    :param an_rel_directory: [string]
    :param instrument_name: [string
    :return: [None]
    """
    if an_rel_directory is None:
        print('Fie! you must provide an_rel_directory.')
        return
    # Get list of all files in directory, and list of all FITS files. Compare & report non-FITS.
    # For each FITS file: (1) Verify object name and FITS name prefix match, (2) verify platesolved,
    #    (3) verify platesolved focal length close to that specified in Instrument object,
    #    (4) verify FOV file present.
    # Write out errors & warnings. Write out summary stats. Log.


class SkyModel:
    def __init__(self, an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None, filter=None,
                 instrument_name='Borea', site_name='DSW',
                 max_cat_mag_error=0.01, max_inst_mag_sigma=0.03, max_color_vi=+2.5,
                 saturation_adu=None,
                 fit_sky_bias=True, fit_vignette=True, fit_xy=False,
                 fit_transform=False, fit_extinction=True, do_plots=True):
        """
        Constructs a sky model using mixed-model regression on df_master.
           Normally used by make_model()
        :param an_top_directory: e.g., 'J:\Astro\Images\C14' [string]
        :param an_rel_directory: e.g., '20170504'. The dir 'Photometry' is subdir of this. [string]
        :param filter: name of filter to which this model applies [string, e.g., 'V' or 'R']
        :param instrument_name: name of Instrument, e.g., 'Borea' [string; name of Instrument obj]
        :param site_name: name of observing site, e.g., 'DSW' [string; name of Site object]
        :param max_cat_mag_error: maximum catalog error allowed to stars in model [float]
        :param max_inst_mag_sigma: max instrument magnitude error allowed star observations [float]
        :param max_color_vi: maximum V-I color allowed to stars in model [float]
        :param saturation_adu: ccd ADUs that constitute saturation [float; from Instrument if None] 
        :param fit_sky_bias: True to fit sky bias term [bool]
        :param fit_vignette: True to fit vignette term (dist^2 from ccd center) [bool]
        :param fit_xy: True to fit X and Y gradient terms [bool]
        :param fit_transform: True to fit transform terms; else use values from Instument obj [bool]
        :param fit_extinction: True to fit extinction terms; else use values from Site obj [bool]
        Parameter 'fit_star_id' is not included in this version (would lead to crossed RE vars).
        """

        self.an_top_directory = an_top_directory
        self.an_rel_directory = an_rel_directory
        self.filter = filter
        self.instrument_name = instrument_name
        self.site_name = site_name
        self.max_cat_mag_error = max_cat_mag_error
        self.max_inst_mag_sigma = max_inst_mag_sigma
        self.max_color_vi = max_color_vi
        if saturation_adu is not None:
            self.saturation_adu = saturation_adu
        else:
            instrument = Instrument(self.instrument_name)
            self.saturation_adu = instrument.camera['saturation_adu']
        self.fit_sky_bias = fit_sky_bias
        self.fit_vignette = fit_vignette
        self.fit_xy = fit_xy
        self.fit_transform = fit_transform
        self.fit_extinction = fit_extinction
        self.dep_var_name = 'InstMag_with_offsets'
        self.df_model = None    # data to/from regression, one row per input pt [pandas DataFrame]
        self.mm_fit = None      # placeholder [MixedModelFit object]
        self.df_star = None     # one row per unique model star [pandas DataFrame]
        self.extinction = None  # scalar result, placeholder
        self.transform = None   # "
        self.vignette = None    # "
        self.x = None           # "
        self.y = None           # "
        self.sky_bias = None    # "
        self.converged = False  # "
        self.n_obs = None       # "
        self.n_images = None    # "
        self.sigma = None       # "
        self.df_image = None    # one row per image, placeholder

        # Rows from df_master, as curated by user in file omit.txt:
        df, warning_lines = apply_omit_txt(self.an_top_directory, self.an_rel_directory)

        # Remove rows for several causes of ineligibility to help form a sky model:
        df = df[df['Filter'] == self.filter]
        df = df[df['StarType'] == 'Comp']
        df = df[df['CatMag'].notnull()]
        df = df[df['Airmass'].notnull()]
        df = df[df['InstMagSigma'] <= self.max_inst_mag_sigma]
        df = df[df['MaxADU_Ur'].notnull()]
        df = df[df['MaxADU_Ur'] <= self.saturation_adu]
        df = df[df['CI'].notnull()]
        df = df[df['CI'] <= self.max_color_vi]
        df = df[df['CatMagError'].notnull()]
        df = df[df['CatMagError'] <= self.max_cat_mag_error]
        self.df_model = df

        self._prep_and_do_regression()
        self._build_output()
        if do_plots:
            self.plots()
        # self.to_json_file()  # GIVE UP on JSON -- it can't handle DataFrames and Series.
        write_stare_comps_txt_stub(self.an_top_directory, self.an_rel_directory)

    # @classmethod
    # def from_json(cls, an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None, filter=None):
    #     """
    #     Alternate constructor, reads from JSON file previously written.
    #        Normally used by predict(), which requires final (immutable) Models.
    #     :param an_top_directory: path to an_rel_folder [str]
    #     :param an_rel_directory: folder for this instrument on this Astronight [string]
    #     :param filter: the filter to which this model applies [string, e.g., 'V' or 'R']
    #     :return: newly constructed Model object [class Model]
    #     """
    #     json_fullpath = os.path.join(an_top_directory, an_rel_directory, "model-" + filter + ".txt")
    #     with open(json_fullpath, 'r') as f:
    #         d = json.load(f)
    #     # Series and DataFrames (pandas) were stored as dictionaries, so convert them back:
    #
    #
    #
    #     # Distribute dictionary elements to their original object attributes:
    #     return cls(an_top_directory=an_top_directory, an_rel_directory=an_rel_directory,
    #                filter=filter,
    #                instrument_name=d['instrument_name'],
    #                max_inst_mag_sigma=d['max_inst_mag_sigma'],
    #                max_cat_mag_error=d['max_cat_mag_error'],
    #                max_color_vi=+d['max_color_vi'], saturated_adu=d['saturation_adu'],
    #                fit_skyBias=d['fit_sky_bias'],
    #                fit_vignette=d['fit_vignette'], fit_xy=d['fit_xy'],
    #                fit_transform=d['fit_transform'], fit_extinction=d['fit_extinction'])
    #
    # def to_json_file(self):
    #     """
    #     Writes (most of) the current object to a JSON file.
    #        Does NOT include df_master (which is huge, but can be found in the same directory).
    #        Also does not include the statsmodels::MixedModelLm which cannot be serialized.
    #     :return: True if file successfully written, else False.
    #     """
    #     json_fullpath = os.path.join(self.an_top_directory, self.an_rel_directory, 'Photometry',
    #                                  "model-" + self.filter + ".json")
    #     json_dict = vars(self).copy()
    #
    #     # Convert pandas DataFrames to dictionaries without json-illegal int64, etc.
    #     json_dict['df_model'] = convert_pandas_to_json_compatible_dict(json_dict['df_model'])
    #     json_dict['df_star'] = convert_pandas_to_json_compatible_dict(json_dict['df_star'])
    #     # Convert pandas Series to dictionaries without json-illegal int64, etc.
    #     json_dict['mm_fit'].fitted_values = \
    #         convert_pandas_to_json_compatible_dict(json_dict['mm_fit'].fitted_values)
    #     json_dict['mm_fit'].group_values = \
    #         convert_pandas_to_json_compatible_dict(json_dict['mm_fit'].group_values)
    #     json_dict['mm_fit'].residuals = \
    #         convert_pandas_to_json_compatible_dict(json_dict['mm_fit'].residuals)
    #
    #     with open(json_fullpath, 'w') as f:
    #         json.dump(json_dict, f, indent=4)

    def _prep_and_do_regression(self):

        # Initiate dependent-variable offset, which will aggregate all such offset terms:
        dep_var_offset = self.df_model['CatMag'].copy()  # *copy* CatMag, or it will be damaged

        # Build fixed-effect (x) variable list and construct dep-var offset:
        fixed_effect_var_list = []
        if self.fit_transform:
            fixed_effect_var_list.append('CI')
        else:
            instrument = Instrument(self.instrument_name)
            transform_vi = instrument.filters[self.filter]['transform']['V-I']
            dep_var_offset += transform_vi * self.df_model['CI']
        if self.fit_extinction:
            fixed_effect_var_list.append('Airmass')
        else:
            site = Site(self.site_name)
            extinction = site.extinction[self.filter]
            dep_var_offset += extinction * self.df_model['Airmass']
        if self.fit_sky_bias:
            if sum([x != 0 for x in self.df_model['SkyBias']]) > int(len(self.df_model) / 2):
                fixed_effect_var_list.append('SkyBias')
        if self.fit_vignette:
            fixed_effect_var_list.append('Vignette')
        if self.fit_xy:
            fixed_effect_var_list.extend(['X1024', 'Y1024'])

        # Build 'random-effect' variable:
        random_effect_var_name = 'FITSfile'  # cirrus effect is per-image

        # Build dependent (y) variable:
        self.df_model[self.dep_var_name] = self.df_model['InstMag'] - dep_var_offset

        # Execute regression:
        self.mm_fit = MixedModelFit(data=self.df_model, dep_var=self.dep_var_name,
                                    fixed_vars=fixed_effect_var_list,
                                    group_var=random_effect_var_name)

    def _build_output(self):
        """
        Builds appropriate output attributes for external use. 
        :return: None
        """
        # Add 1/obs regression data as new df_model columns:
        self.df_model['FittedValue'] = self.mm_fit.df_observations['FittedValue']
        self.df_model['Residual'] = self.mm_fit.df_observations['Residual']
        dep_var_offset = self.df_model['InstMag'] - self.df_model[self.dep_var_name]
        self.df_model['FittedInstMag'] = self.df_model['FittedValue'] + dep_var_offset

        # Build df_star (star ID and count only):
        self.df_star = self.df_model[['Serial', 'ModelStarID']].groupby('ModelStarID').count()
        self.df_star['ModelStarID'] = self.df_star.index

        # Build df_image (from Mixed Model random effect), 1 row per FITS file,
        #    index = FITSfile, columns = FITSfile, JD_mid, Value:
        df = self.mm_fit.df_random_effects.copy()
        df = df.rename(columns={'GroupName': 'FITSfile', 'GroupValue': 'Value'})
        df_xref = self.df_model[['FITSfile', 'JD_mid']].drop_duplicates()
        df = pd.merge(df, df_xref, on='FITSfile', how='left', sort=False).sort_values(by='JD_mid')
        self.df_image = df.copy()

        # Extract and store scalar results:
        if self.fit_transform:
            self.transform = self.mm_fit.df_fixed_effects.loc['CI', 'Value']  # .loc(row, col)
        else:
            instrument = Instrument(self.instrument_name)
            self.transform = instrument.filters[self.filter]['transform']['V-I']

        if self.fit_extinction:
            self.extinction = self.mm_fit.df_fixed_effects.loc['Airmass', 'Value']
        else:
            site = Site(self.site_name)
            self.extinction = site.extinction[self.filter]

        self.vignette = self.mm_fit.df_fixed_effects.loc['Vignette', 'Value'] \
            if self.fit_vignette is True else 0.0
        self.x = self.mm_fit.df_fixed_effects.loc['X1024', 'Value'] if self.fit_xy is True else 0.0
        self.y = self.mm_fit.df_fixed_effects['Y1024', 'Value'] if self.fit_xy is True else 0.0
        self.sky_bias = self.mm_fit.df_fixed_effects.loc['SkyBias', 'Value'] \
            if self.fit_sky_bias is True else 0.0
        self.converged = self.mm_fit.converged
        self.n_obs = len(self.df_model)
        self.n_images = len(self.df_model['FITSfile'].drop_duplicates())
        self.sigma = self.mm_fit.sigma
        print('\n', len(self.df_model), ' observations --> sigma=',
              round((1000.0 * self.sigma), 1), ' mMag')

    def plots(self):
        import numpy as np
        import matplotlib.pyplot as plt

        # Setup for all Figures.
        obs_is_std = [name.startswith('Std_') for name in self.df_model['FOV']]
        obs_point_colors = ['darkgreen' if obs_is_std[i] is True else 'black'
                            for i, x in enumerate(obs_is_std)]
        image_is_std = [name.startswith('Std_') for name in self.df_image['FITSfile']]
        image_point_colors = ['darkgreen' if image_is_std[i] is True else 'black'
                              for i, x in enumerate(image_is_std)]
        jd_floor = floor(min(self.df_model['JD_mid']))
        obs_jd_fract = self.df_model['JD_mid']-jd_floor
        xlabel_jd = 'JD(mid)-' + str(jd_floor)
        obs_residuals_mmag = self.df_model['Residual'] * 1000.0

        # FIGURE 1 (Q-Q plot):
        from scipy.stats import norm
        df_y = self.df_model.copy()[['Residual', 'Serial']]
        df_y['Residual'] *= 1000.0
        df_y['Colors'] = obs_point_colors  # keep colors attached to correct residuals when sorted
        df_y = df_y.sort_values(by='Residual')
        n = len(df_y)
        t_values = [norm.ppf((k-0.5)/n) for k in range(1, n+1)]

        # Construct Q-Q plot:
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(12, 8))  # (width, height) in "inches"
        ax.grid(True, color='lightgray', zorder=-1000)
        ax.set_title('Q-Q plot of Residuals: ' +
                     self.an_rel_directory + '    ' + self.filter +
                     ' filter', color='darkblue', fontsize=20, weight='bold')
        ax.set_xlabel('t (sigma.residuals = ' + str(round(1000.0 * self.sigma, 1)) + ' mMag)')
        ax.set_ylabel('Residual (mMag)')
        ax.scatter(x=t_values, y=df_y['Residual'], alpha=0.6, color=df_y['Colors'], zorder=+1000)

        # Label potential outliers:
        mean_y = df_y['Residual'].mean()
        std_y = df_y['Residual'].std()
        z_score_y = (df_y['Residual'] - mean_y) / std_y
        df_y['T'] = t_values
        df_to_label = df_y[abs(z_score_y) >= 2.0]
        for x, y, label in zip(df_to_label['T'], df_to_label['Residual'], df_to_label['Serial']):
            ax.annotate(label, xy=(x, y), xytext=(4, -4),
                        textcoords='offset points', ha='left', va='top', rotation=-40)

        # Add reference line:
        x_low = 1.10 * min(df_y['T'])
        x_high = 1.10 * max(df_y['T'])
        y_low = x_low * std_y
        y_high = x_high * std_y
        ax.plot([x_low, x_high], [y_low, y_high], color='gray', zorder=-100, linewidth=1)

        # Add annotation: number of observations:
        fig.text(x=0.5, y=0.87,
                 s=str(len(self.df_model)) + ' observations in model.',
                 verticalalignment='top', horizontalalignment='center',
                 fontsize=12)
        plt.show()

        # FIGURE 2 (multiplot): Set up plot grid and style parameters:
        fig, axes = plt.subplots(ncols=4, nrows=3, figsize=(16, 10))  # (width, height) in "inches"

        def make_labels(ax, title, xlabel, ylabel, zero_line=True):
            ax.set_title(title, y=0.89)
            ax.set_xlabel(xlabel, labelpad=-27)
            ax.set_ylabel(ylabel, labelpad=-8)
            if zero_line is True:
                ax.axhline(y=0, color='lightgray', linewidth=1, zorder=-100)

        # Cirrus Plot (one point per image):
        ax = axes[0, 0]
        make_labels(ax, 'Image Cirrus Plot', xlabel_jd, 'mMag')
        ax.scatter(x=self.df_image['JD_mid']-jd_floor, y=self.df_image['Value'] * 1000.0,
                   alpha=0.6, color=image_point_colors)

        # Sky background vs JD_mid:
        ax = axes[0, 1]
        make_labels(ax, 'Sky background vs JD_mid', xlabel_jd, 'Sky ADU',
                    zero_line=False)
        ax.scatter(x=obs_jd_fract, y=self.df_model['SkyADU'],
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs Instrument Magnitude:
        ax = axes[0, 2]
        make_labels(ax, 'Residuals vs Instrument Mag', 'Instrument Mag', 'mMag')
        ax.scatter(x=self.df_model['InstMag'], y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs Max ADU:
        ax = axes[0, 3]
        xlabel_text = 'Max ADU, uncalibrated [log scale]'
        make_labels(ax, 'Residuals vs Max ADU', xlabel_text, 'mMag')
        ax.set_xlabel(xlabel_text, labelpad=-30)
        ax.set_xscale('log')
        x_scale_min = min(1000.0, 0.9 * min(self.df_model['MaxADU_Ur']))
        x_scale_max = 1.1 * max(self.df_model['MaxADU_Ur'])
        ax.set_xlim(x_scale_min, x_scale_max)
        # ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
        ax.scatter(x=self.df_model['MaxADU_Ur'], y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs Sky background:
        ax = axes[1, 0]
        make_labels(ax, 'Residuals vs Sky background', 'Sky ADUs', 'mMag')
        ax.scatter(x=self.df_model['SkyADU'], y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs JD:
        ax = axes[1, 1]
        make_labels(ax, 'Residuals vs JD', xlabel_jd, 'mMag')
        ax.scatter(x=obs_jd_fract, y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs Color Index:
        ax = axes[1, 2]
        make_labels(ax, 'Residuals vs Color Index', 'Color Index (V-I)', 'mMag')
        ax.scatter(x=self.df_model['CI'], y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs Airmass:
        ax = axes[1, 3]
        make_labels(ax, 'Residuals vs Airmass', 'Airmass', 'mMag')
        ax.scatter(x=self.df_model['Airmass'], y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs Exposure Time:
        ax = axes[2, 0]
        make_labels(ax, 'Residuals vs Exposure Time', 'seconds', 'mMag')
        ax.scatter(x=self.df_model['Exposure'], y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs Vignette:
        ax = axes[2, 1]
        make_labels(ax, 'Residuals vs Vignette', 'pixels from CCD center', 'mMag')
        ax.scatter(x=1024*np.sqrt(self.df_model['Vignette']), y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs X:
        ax = axes[2, 2]
        make_labels(ax, 'Residuals vs X', 'X pixels from CCD center', 'mMag')
        ax.scatter(x=1024*self.df_model['X1024'], y=self.df_model['Residual'] * 1000.0, alpha=0.6,
                   color=obs_point_colors)

        # Residuals vs Y:
        ax = axes[2, 3]
        make_labels(ax, 'Residuals vs Y', 'Y pixels from CCD center', 'mMag')
        ax.scatter(x=1024*self.df_model['Y1024'], y=self.df_model['Residual'] * 1000.0, alpha=0.6,
                   color=obs_point_colors)

        # Finish the figure, and show the entire plot:
        fig.tight_layout(rect=(0, 0, 1, 0.925))
        fig.subplots_adjust(left=0.06, bottom=0.06, right=0.94, top=0.85, wspace=0.25, hspace=0.25)
        fig.suptitle(self.an_rel_directory +
                     '                   ' + self.filter + ' filter                ' +
                     '{:%Y-%m-%d     %H:%M  utc}'.format(datetime.now(timezone.utc)),
                     color='darkblue', fontsize=20, weight='bold')
        plt.show()

    def predict(self, predict_input):
        """
        Uses current model to predict best star magnitudes for *observed* InstMag and other inputs.
        :param predict_input: data, all needed input columns, including groups [pandas DataFrame]
        :return: a dependent variable prediction for each input row [pandas Series of floats,
            with index = index of predict_input]
        """
        # First, verify that input is a DataFrame and that all needed columns are present.
        #    Names must be same as in model.
        if not isinstance(predict_input, pd.DataFrame):
            print('>>>>> SkyModel.predict(predict_input): predict_input is not a pandas DataFrame.')
            return None
        required_input_columns = ['Serial', 'FITSfile', 'InstMag', 'CI', 'Airmass']
        if self.fit_sky_bias:
            required_input_columns.append('SkyBias')
        if self.fit_vignette:
            required_input_columns.append('Vignette')
        if self.fit_xy:
            required_input_columns.extend(['X1024', 'Y1024'])
        all_present = all([name in predict_input.columns for name in required_input_columns])
        if not all_present:
            print('>>>>> SkyModel.predict(): at least one column missing.')
            print('      Current model requires these columns:')
            print('         ' + ', '.join(required_input_columns))
            return None

        # Make copy of dataframe; add bogus CatMag column (required by model):
        df_for_predict = predict_input.copy()
        bogus_cat_mag = 0.0
        df_for_predict['CatMag'] = bogus_cat_mag  # totally bogus local value, to be reversed later

        # Execute MixedModelFit.predict(), giving Intercept + bogus CatMag + FEs + REs (pd.Series):
        raw_predictions = self.mm_fit.predict(df_for_predict)

        # Now, the tricky part: estimating best magnitudes for unknowns/targets
        #   (effectively get best CatMag per star).
        # eq A: predict() = Intercept + bogus CatMag + FEs + REs (per MixedModelFit.predict()).
        # eq B: obs InstMag - offsets ~~ Intercept + best CatMag + FEs + REs (as in regression).
        # so (eq B - eq A) -> best CatMag ~~ obs InstMag - offsets - predict() + bogus CatMag
        # We still need to get: offsets (just as done in regression), then estimate best catMag:

        # Compute dependent-variable offsets for unknown stars:
        dep_var_offsets = pd.Series(len(df_for_predict) * [0.0], index=raw_predictions.index)
        if self.fit_transform is False:
            dep_var_offsets += self.transform * df_for_predict['CI']
        if self.fit_extinction is False:
            dep_var_offsets += self.extinction * df_for_predict['Airmass']

        # Extract best CatMag d'un seul coup, per (eq B - eq A), above:
        predicted_star_mags = \
            df_for_predict['InstMag'] - dep_var_offsets - raw_predictions + bogus_cat_mag
        return predicted_star_mags


class PredictionSet:
    def __init__(self, an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None,
                 instrument_name='Borea', site_name='DSW',
                 max_inst_mag_sigma=0.05, skymodel_list=None):
        """
        Constructs a prediction set, i.e., a set of best estimates of comp, check, and target star
        magnitudes, ready for marking up (curating) and reporting to the AAVSO (for example). 
        :param an_top_directory: e.g., 'J:\Astro\Images\C14' [string]
        :param an_rel_directory: e.g., '20170504'. The dir 'Photometry' is subdir of this. [string]
        :param instrument_name: name of Instrument, e.g., 'Borea' [string; name of Instrument obj]
        :param site_name: name of observing site, e.g., 'DSW' [string; name of Site object]
        :param max_inst_mag_sigma: max instrument magnitude error allowed star observations [float]
        :param skymodel_list: SkyModel objects ready to be used [list of SkyModel objs]
        """
        self.an_top_directory = an_top_directory
        self.an_rel_directory = an_rel_directory
        self.instrument_name = instrument_name
        self.instrument = Instrument(self.instrument_name)
        self.site_name = site_name
        self.site = Site(self.site_name)
        self.max_inst_mag_sigma = max_inst_mag_sigma
        # ---- we'll deal with transform preferences in a later version.
        #      for now, we'll stick with V-I for all filters
        # if transform_preferences is None:
        #     transform_preferences = self.get_transform_preferences(self.instrument)
        # else:
        #     self.transform_preferences = transform_preferences
        self.saturation_adu = self.instrument.camera['saturation_adu']
        self.skymodels = dict((skymodel.filter, skymodel) for skymodel in skymodel_list)
        self.images_with_targets = None

        df_eligible_obs, warning_lines = apply_omit_txt(self.an_top_directory,
                                                        self.an_rel_directory)

        df_curated_comp_obs = self.curate_stare_comps(df_eligible_obs)

        df_estimates_comps = self.estimate_comps(df_curated_comp_obs)

        df_cirrus_effect = self.estimate_cirrus_effect(df_estimates_comps)

        df_estimates_checks_targets = self.estimate_checks_targets()

        self.compute_all_errors()

        write_report_map_stub()

        self._to_json()  # maybe.


    def curate_stare_comps(self, df_eligible_obs):
        """
        Using user's stare_comps.txt in this an_rel_directory, 
           limit comps applied to stare observations.
        :param df_eligible_obs: dataframe of obs (typically all obs eligible to this pt) [DataFrame]
        :return: dataframe of all observations remaining eligible after this curation [DataFrame]
        """
        fullpath = os.path.join(self.an_top_directory, self.an_rel_directory,
                                'Photometry', 'stare_comps.txt')
        if not os.path.exists(fullpath):
            write_omit_txt_stub(self.an_top_directory, self.an_rel_directory)
            return df_eligible_obs.copy()  # no change or warnings, since stare_comps.txt was absent

        with open(fullpath) as f:
            lines = f.readlines()
        lines = [line.split(";")[0] for line in lines]  # remove all comments
        lines = [line.strip() for line in lines]  # remove lead/trail blanks
        lines = [line for line in lines if line != '']  # remove empty lines
        df_curated_obs = df_eligible_obs.copy()  # default in case no lines in stare_comps.txt
        warning_lines = []

        for line in lines:
            warning_line = None
            if line.startswith("#COMPS"):
                parms, warning_line = get_line_parms(line, "#OBS", 3, None)
                if parms is not None:
                    fov_name = parms[0]
                    filter_name = parms[1]
                    comp_ids = ' '.join(parms[2:]).split()
                    rows_to_remove = (df_eligible_obs['StarType'] == 'Comp') and \
                        (df_eligible_obs['FOV'] == fov_name) and \
                        (df_eligible_obs['Filter'] == filter_name) and \
                        (df_eligible_obs['StarID'] not in comp_ids)
                    df_curated_obs = df_eligible_obs[~ rows_to_remove]
                else:
                    warning_lines.append(warning_line)
            else:
                warning_line = 'Directive not understood: \'' + line + '\'.'
                warning_lines.append(warning_line)
        return df_curated_obs, warning_lines

    def estimate_comps(self, df_eligible_obs):
        """
        Get raw, predicted mag estimates for ALL eligible comp-star observations, in all filters.
        :param df_eligible_obs: data for all eligible star observations [pandas DataFrame]  
        :return: dataframe of 
        """
        df = df_eligible_obs.copy()
        # Select rows, then remove rows that suffer various causes of ineligibility:
        df = df[df['StarType'] == 'Comp']
        df = df[df['MaxADU_Ur'].notnull()]
        df = df[df['MaxADU_Ur'] <= self.saturation_adu]
        df = df[df['InstMagSigma'] <= self.max_inst_mag_sigma]
        df = df[df['CatMag'].notnull()]
        df = df[df['CI'].notnull()]
        df = df[df['Airmass'].notnull()]
        df_eligible_comp_obs = df

        images_with_eligible_comps = df_eligible_comp_obs['FITSfile'].drop_duplicates()
        self.images_with_targets = \
            (df_eligible_obs[df_eligible_obs['StarType'].lower() == 'target'])['FITSfile'] \
            .drop_duplicates()  # eliminates standard FOVs etc
        images_with_targets_and_comps = [im for im in self.images_with_targets
                                         if im in images_with_eligible_comps]

        # Get best magnitude estimates for all eligible comps:
        df_list = []
        for filter_name, skymodel in self.skymodels.items():
            comp_obs_to_include = (df_eligible_comp_obs['Filter'] == filter_name) and \
                                   (df_eligible_comp_obs['FITSfile'] in images_with_targets)
            df_predict_input = df_eligible_comp_obs[comp_obs_to_include]
            predict_output = skymodel.predict(df_predict_input)

            df_estimates_this_model = \
                df_predict_input[['Serial', 'ModelStarID', 'FITSfile', 'StarID', 'Chart',
                                  'Xcentroid', 'Ycentroid', 'InstMag', 'InstMagSigma', 'StarType',
                                  'CatMag', 'CatMagSaved', 'CatMagError', 'Exposure',
                                  'JD_mid', 'Filter', 'Airmass', 'CI', 'SkyBias', 'Vignette']]
            df_estimates_this_model['EstimatedMag'] = predict_output
            df_list.append(df_estimates_this_model)  # list.append() performs in-place
        df_estimates_comps = pd.concat(df_list, ignore_index=True)
        df_estimates_comps['UseInEnsemble'] = True  # default until falsified
        df_estimates_comps.index = df_estimates_comps['Serial']
        return df_estimates_comps

    def estimate_cirrus_effect(self, df_estimates_comps):
        """
        Generates per-image cirrus effect with more careful inclusion/exclusion criteria than
           were used in the MixedModelFit, esp. in rejecting outlier comp observations.
        :return: data for best per-image cirrus-effects, in magnitudes [pandas DataFrame] 
        """
        df_list = []
        for image in self.images_with_targets:
            df_estimates_comps_this_image = \
                df_estimates_comps[df_estimates_comps['FITSfile'] == image]
            cirrus_effect_from_comps = \
                df_estimates_comps_this_image['EstimatedMag'] - \
                df_estimates_comps_this_image['CatMag']
            sigma2 = df_estimates_comps_this_image['CatMagError']**2 + \
                df_estimates_comps_this_image['InstMagSigma']**2
            least_allowed_sigma = 0.01
            raw_weights = pd.Series([1.0 / max(s2, least_allowed_sigma**2) for s2 in sigma2])
            normalized_weights = raw_weights / sum(raw_weights)
            v2 = sum(normalized_weights**2)
            cirrus_effect_this_image = weighted_mean(cirrus_effect_from_comps, normalized_weights)
            resid2 = (cirrus_effect_from_comps - cirrus_effect_this_image)**2
            cirrus_effect_this_image = df_estimates_comps_this_image.iloc[1]['CatMagError'] \
                if len(df_estimates_comps_this_image) == 1 \
                else sqrt(v2 * sum(normalized_weights * resid2))
            comp_ids_used = df_estimates_comps_this_image['StarID']  # starting point = use all
            num_comps_used = len(df_estimates_comps_this_image)      # "
            num_comps_removed = 0                                    # "

            # Reject this image's worst comp stars and recalculate, if
            #    (we start with at least 4 comp stars) AND (criterion1 >= 16 OR criterion2 >= 20).
            criterion1 = 0  # how many times worse is the worst comp vs avg of other comps
            criterion2 = 0  # square of worst comp's effective t-value, relative to CatMagError
            if len(df_estimates_comps_this_image) >= 4:
                x = normalized_weights * resid2  # a pd.Series; will be >= 0
                criterion1 = max(x) / ((sum(x)-max(x)) / (len(x)-1))
                y = raw_weights * resid2  # a pd.Series; will be >= 0
                criterion2 = max(y)
                selection = pd.Series(len(df_estimates_comps_this_image) * [True])  # is this used??
                if (criterion1 >= 16) or (criterion2 >= 20):
                    c1 = x / ((sum(x)-x) / (len(x)-1))
                    c2 = y
                    score = pd.Series([max(c_1/16.0, c_2/20.0) for (c_1, c_2) in zip(c1, c2)])
                    max_to_remove = floor(len(score) / 4)


def get_df_master(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None):
    """
    Simple utility to read df_master.csv file and return the DataFrame.
    :param an_top_directory: path to an_rel_directory [str]
    :param an_rel_directory: directory for this instrument on this Astronight [string] 
    :return: pandas DataFrame with all comp, check, and target star raw photometric data
         (which dataframe is generated and csv-written by R, as of May 2017).
         The DataFrame index is set to equal column Serial (which was already a kind of index).
    """
    if an_rel_directory is None:
        return None
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'df_master.csv')
    if not os.path.exists(fullpath):
        return None
    df_master = pd.read_csv(fullpath, sep=';')
    df_master.index = df_master['Serial']
    return df_master


def apply_omit_txt(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None):
    """
    Gets df_master and omit.txt and returns filtered DataFrame.
    :param an_top_directory: [string]
    :param an_rel_directory: [string]
    :return: rows of df_master whose omission is NOT requested by user (in omit.txt) [DataFrame] 
    """
    if an_rel_directory is None:
        return None
    df_master = get_df_master(an_top_directory=an_top_directory, an_rel_directory=an_rel_directory)
    if df_master is None:
        return None

    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit.txt')
    if not os.path.exists(fullpath):
        write_omit_txt_stub(an_top_directory, an_rel_directory)
        return df_master.copy(), []  # no change or warnings, since omit.txt absent

    with open(fullpath) as f:
        lines = f.readlines()
    lines = [line.split(";")[0] for line in lines]  # remove all comments
    lines = [line.strip() for line in lines]        # remove lead/trail blanks
    lines = [line for line in lines if line != '']  # remove empty lines
    df_filtered = df_master.copy()  # start with a copy, omit lines per user requests.
    warning_lines = []

    for line in lines:
        warning_line = None
        rows_to_omit = len(df_filtered) * [False]  # default to be overwritten
        if line.startswith('#OBS'):
            parms, warning_line = get_line_parms(line, "#OBS", 2, 2)
            if parms is not None:
                fits_file_name = parms[0] + '.fts'
                star_id = parms[1]
                rows_to_omit = (df_filtered['FITSfile'] == fits_file_name) & \
                               (df_filtered['StarID'] == star_id)  # a pandas Series
            else:
                warning_lines.append(warning_line)
        elif line.startswith('#STAR'):
            parms, warning_line = get_line_parms(line, "#STAR", 2, 3)
            if parms is not None:
                fov = parms[0]
                star_id = parms[1]
                rows_to_omit = (df_filtered['FOV'] == fov) & (df_filtered['StarID'] == star_id)
                if len(parms) == 3:
                    # further restrict omission request
                    filter = parms[2]
                    rows_to_omit = rows_to_omit & (df_filtered['Filter'] == filter)
            else:
                warning_lines.append(warning_line)
        elif line.startswith('#SERIAL'):
            raw_parms, warning_line = get_line_parms(line, "#SERIAL", 1, None)
            if raw_parms is not None:
                parms = ' '.join(raw_parms).split()  # resep. on whitespace;  no empty elements
                serials_to_omit = [int(p) for p in parms]
                rows_to_omit = df_filtered['Serial'].isin(serials_to_omit)
            else:
                warning_lines.append(warning_line)
        elif line.startswith('#IMAGE'):
            parms, warning_line = get_line_parms(line, "#IMAGE", 1, 1)
            if parms is not None:
                image_to_omit = parms[0] + '.fts'
                rows_to_omit = df_filtered['FITSfile'] == image_to_omit
            else:
                warning_lines.append(warning_line)
        elif line.startswith('#JD'):
            parms, warning_line = get_line_parms(line, "#JD", 2, 2)
            if parms is not None:
                jd_floor = floor(min(df_filtered['JD_mid']))
                jd_start = float(parms[0]) + jd_floor
                jd_end = float(parms[1]) + jd_floor
                rows_to_omit = (df_filtered['JD_mid'] >= jd_start) & \
                          (df_filtered['JD_mid'] <= jd_end)
            else:
                warning_lines.append(warning_line)
        else:
            warning_line = 'Directive not understood: \'' + line + '\'.'
            warning_lines.append(warning_line)

        if sum(rows_to_omit) >= 1:
            df_filtered = df_filtered[~ rows_to_omit]  # remove rows as user requested.
        else:
            if warning_line is None:
                warning_lines.append('No rows omitted: \'' + line + '\'.')

    for warning_line in warning_lines:
        print(warning_line)
    print(str(len(df_master) - len(df_filtered)) + ' rows removed via omit.txt. ' +
          str(len(df_filtered)) + ' rows remain.')
    return df_filtered, warning_lines


def write_omit_txt_stub(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None):
    lines = [';----- This is omit.txt for AN directory ' + an_rel_directory,
             ';----- Use this file to omit observations from input to SkyModel (all filters).',
             ';----- Example directive lines:',
             ';',
             ';#OBS    Obj-0000-V, 132 ; to omit star 132 from FITS image Obj-0000-V.fts',
             ';#STAR   FOV, 132, V     ; to omit star 132 from all FITS with FOV '
             'and filter V',
             ';#STAR   FOV, 132        ; to omit star 132 from all FITS with FOV '
             'and ALL filters',
             ';#IMAGE  Obj-0000-V      ; to omit FITS image Obj-0000-V.fts specifically',
             ';#JD     0.72, 1         ; to omit fractional JD from 0.72 through 1',
             ';#SERIAL 123,77 54   6   ; to omit observations by Serial number (many per line OK)',
             ';',
             ';----- Add your directive lines:',
             ';']
    lines = [line + '\n' for line in lines]
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit.txt')
    if not os.path.exists(fullpath):
        with open(fullpath, 'w') as f:
            f.writelines(lines)
        lines_written = len(lines)
    else:
        lines_written = 0
    return lines_written


def write_stare_comps_txt_stub(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None):
    lines = [';----- This is stare_comps.txt for AN directory ' + an_rel_directory,
             ';----- Select comp stars (by FOV, filter, & StarID) from input to predict().',
             ';----- Example directive line:',
             ';',
             ';#COMPS  Obj, V, 132, 133 144    ; to KEEP from FOV \'Obj\': '
             'comp stars \'132\' \'133\' and \'144\' in filter \'V\'',
             ';',
             ';----- Add your directive lines:',
             ';']
    lines = [line + '\n' for line in lines]
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'stare_comps.txt')
    if not os.path.exists(fullpath):
        with open(fullpath, 'w') as f:
            f.writelines(lines)
        lines_written = len(lines)
    else:
        lines_written = 0
    return lines_written


def write_report_map_stub(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None):
    lines = [';----- This is report_map.txt for AN directory ' + an_rel_directory,
             ';----- Use this file to omit and/or combine target observations, for AAVSO report',
             ';----- Example directive line:',
             ';',
             ';#TARGET  GD Cyg    ; to omit this target star altogether from AAVSO report',
             ';#JD      0.65 0.67 ; to omit this JD (fractional) range from AAVSO report',
             ';#SERIAL  34 44,129  32 1202 ; to omit these 5 Serial numbers from AAVSO report',
             ';#COMBINE 80 128    ; to combine (average) these 2 Serial numbers in AAVSO report',
             ';',
             ';----- Add your directive lines:',
             ';']
    lines = [line + '\n' for line in lines]
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'report_map.txt')
    if not os.path.exists(fullpath):
        with open(fullpath, 'w') as f:
            f.writelines(lines)
        lines_written = len(lines)
    else:
        lines_written = 0
    return lines_written



# def convert_pandas_to_json_compatible_dict(pandas_obj):
#         """
#         Python's json package cannot handle int64, etc that pandas insists on putting into
#             its Series and DataFrame constructs. So have them write themselves to JSON text,
#             then have json read the text back in as a dictionary ready for inclusion in an
#             object that json package considers legitimate to write to JSON text. Whew.
#         :return: dict representation of pandas object, ready to include in a JSON file [py dict]
#         """
#         json_text = pandas_obj.to_json()
#         return json.loads(json_text)  # a dictionary without int64, etc


def get_line_parms(line, directive, nparms_min=None, nparms_max=None):
    line_parms = [p.strip() for p in line[(len(directive)):].split(',')]
    valid_num_line_parms = True  # until falsified
    if nparms_min is not None:  # ignore nparms_min if None.
        valid_num_line_parms = valid_num_line_parms & (len(line_parms) >= nparms_min)
    if nparms_max is not None:  # ignore nparms_max if None.
        valid_num_line_parms = valid_num_line_parms & (len(line_parms) <= nparms_max)
    if valid_num_line_parms:
        return line_parms, None
    else:
        return None, 'Line has wrong number of parameters: \'' + line + '\'.'


    # def get_transform_preferences(self, instrument):
    #     # transform_preferences: for each observation filter, one or more filter pairs, in
    #     #    descending order of preference, used to make transform corrections, in the form
    #     #    {'V':['V-I', 'B-V'], 'R':['V-I', 'V-R', 'B-V'], 'I':['V-I', 'R-I', 'V-R', 'B-V']}
    #     #    None -> get from Instrument object [dict of {str:list of strings}]
    #     pass

