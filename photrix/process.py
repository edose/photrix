import json
import os

import pandas as pd

from .user import Instrument, Site
from .util import MixedModelFit

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
                 fit_sky_bias=True, fit_vignette=True, fit_xy=True,
                 fit_transform=False, fit_extinction=True):
        """
        Constructs a sky model using mixed-model regression on df_master.
           Normally used by make_model()
        :param an_top_directory: e.g., 'J:\Astro\Images\C14' [string]
        :param an_rel_directory: e.g., '20170504'. The dir 'Photometry' is subdir of this. [string]
        :param filter: name of filter to which this model applies [string, e.g., 'V' or 'R']
        :param instrument_name: name of Instrument, e.g., 'Borea' [string; name of Instrument obj]
        :param site_name: name of observing site, e.g., 'DSW' [string; name of Site object]
        :param max_cat_mag_error: maximum catalog error allowed to stars in model [float]
        :param max_inst_mag_sigma: maximum instrument magnitude error allowed star observations [float]
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
        self.df_model = None  # data to/from regression, one row per input pt [pandas DataFrame]
        # self.regression_results = None  # placeholder [MixedModelFit object]
        self.df_star = None  # one row per unique model star [pandas DataFrame]
        self.extinction = None  # scalar result
        self.transform = None   # "
        self.vignette = None    # "

        # note: df_master is NOT part of the object being constructed here (much too large).
        df_master = get_df_master(self.an_top_directory, self.an_rel_directory)

        # Remove rows from df_master as specified by user in file omit.txt:
        df, warning_lines = apply_omit_txt(df_master)

        # Remove rows for several ineligibilities:
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
        self._build_output_tables()
        self.make_plots(df_master)
        # self.to_json_file()  # GIVE UP on JSON -- it can't handle DataFrames and Series.
        write_stare_comps_txt_stub(self.an_top_directory, self.an_rel_directory)

    @classmethod
    def from_json(cls, an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None, filter=None):
        """
        Alternate constructor, reads from JSON file previously written.
           Normally used by predict(), which requires final (immutable) Models.
        :param an_top_directory: path to an_rel_folder [str]
        :param an_rel_directory: folder for this instrument on this Astronight [string]
        :param filter: the filter to which this model applies [string, e.g., 'V' or 'R']
        :return: newly constructed Model object [class Model]
        """
        json_fullpath = os.path.join(an_top_directory, an_rel_directory, "model-" + filter + ".txt")
        with open(json_fullpath, 'r') as f:
            d = json.load(f)
        # Series and DataFrames (pandas) were stored as dictionaries, so convert them back:



        # Distribute dictionary elements to their original object attributes:
        return cls(an_top_directory=an_top_directory, an_rel_directory=an_rel_directory,
                   filter=filter,
                   instrument_name=d['instrument_name'],
                   max_inst_mag_sigma=d['max_inst_mag_sigma'],
                   max_cat_mag_error=d['max_cat_mag_error'],
                   max_color_vi=+d['max_color_vi'], saturated_adu=d['saturation_adu'],
                   fit_skyBias=d['fit_sky_bias'],
                   fit_vignette=d['fit_vignette'], fit_xy=d['fit_xy'],
                   fit_transform=d['fit_transform'], fit_extinction=d['fit_extinction'])

    def to_json_file(self):
        """
        Writes (most of) the current object to a JSON file.
           Does NOT include df_master (which is huge, but can be found in the same directory).
           Also does not include the statsmodels::MixedModelLm which cannot be serialized.
        :return: True if file successfully written, else False.
        """
        json_fullpath = os.path.join(self.an_top_directory, self.an_rel_directory, 'Photometry',
                                     "model-" + self.filter + ".json")
        json_dict = vars(self).copy()

        # Convert pandas DataFrames to dictionaries without json-illegal int64, etc.
        json_dict['df_model'] = convert_pandas_to_json_compatible_dict(json_dict['df_model'])
        json_dict['df_star'] = convert_pandas_to_json_compatible_dict(json_dict['df_star'])
        # Convert pandas Series to dictionaries without json-illegal int64, etc.
        json_dict['regression_results'].fitted_values = \
            convert_pandas_to_json_compatible_dict(json_dict['regression_results'].fitted_values)
        json_dict['regression_results'].group_values = \
            convert_pandas_to_json_compatible_dict(json_dict['regression_results'].group_values)
        json_dict['regression_results'].residuals = \
            convert_pandas_to_json_compatible_dict(json_dict['regression_results'].residuals)

        with open(json_fullpath, 'w') as f:
            json.dump(json_dict, f, indent=4)

    def _prep_and_do_regression(self):
        # Build dep variable and new column (whose values may be adjusted below):
        dep_var = 'DepVar'
        self.df_model[dep_var] = self.df_model['InstMag'] - self.df_model['CatMag']

        # Build fixed-effect variable list:
        fixed_effect_vars = []
        if self.fit_transform:
            fixed_effect_vars.append('CI')
        else:
            instrument = Instrument(self.instrument_name)
            transform_vi = instrument.filters[self.filter]['transform']['V-I']
            transform_adjustment = transform_vi * self.df_model['CI']
            self.df_model[dep_var] -= transform_adjustment  # TODO: verify sign is correct.
        if self.fit_extinction:
            fixed_effect_vars.append('Airmass')
        else:
            site = Site(self.site_name)
            extinction = site.extinction[self.filter]
            extinction_adjustment = extinction * self.df_model['Airmass']
            self.df_model[dep_var] -= extinction_adjustment  # TODO: verify sign is correct.
        if self.fit_sky_bias:
            if sum([x != 0 for x in self.df_model['SkyBias']]) > int(len(self.df_model) / 2):
                fixed_effect_vars.append('SkyBias')
        if self.fit_vignette:
            fixed_effect_vars.append('Vignette')
        if self.fit_xy:
            fixed_effect_vars.extend(['X1024', 'Y1024'])

        # Build groups ('random-effect') variable:
        group_var = 'FITSfile'  # cirrus effect is per-image

        # Execute regression:
        self.regression_results = MixedModelFit(data=self.df_model, dep_var=dep_var,
                                                fixed_vars=fixed_effect_vars, group_var=group_var)

    def _build_output_tables(self):
        # Add new columns to df_model:
        self.df_model['Residual'] = self.regression_results.residuals
        self.df_model['FittedValues'] = self.regression_results.fitted_values

        # Build df_star (star ID and count only):
        self.df_star = self.df_model[['Serial', 'ModelStarID']].groupby('ModelStarID').count()

        # Extract and store scalar results:
        # self.transform, .extinction, and .vignette.
        if self.fit_transform:
            self.transform = self.regression_results.coeffs['CI']
        else:
            instrument = Instrument(self.instrument_name)
            self.transform = instrument.filters[self.filter]['transform']['V-I']
        if self.fit_extinction:
            self.extinction = self.regression_results.coeffs['Airmass']
        else:
            site = Site(self.site_name)
            self.extinction = site.extinction[self.filter]
        self.vignette = self.regression_results.coeffs['Vignette'] if self.fit_vignette else 0.0

        print('\n', len(self.df_model), ' observations --> sigma=',
              round((1000.0*self.regression_results.sigma), 1), ' mMag')

    def make_plots(self, df_master=None):
        if df_master is None:
            df_master = get_df_master(self.an_top_directory, self.an_rel_directory)
        pass  # make the make_plots here


def get_df_master(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None):
    """
    Simple utility to read df_master.csv file and return the DataFrame.
    :param an_top_directory: path to an_rel_directory [str]
    :param an_rel_directory: directory for this instrument on this Astronight [string] 
    :return: pandas DataFrame with all comp, check, and target star raw photometric data
         (which dataframe is generated and csv-written by R, as of May 2017).
    """
    if an_rel_directory is None:
        return None
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'df_master.csv')
    if not os.path.exists(fullpath):
        return None
    df_master = pd.read_csv(fullpath, sep=';')
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

    # Nested function (strictly local for convenience):
    def get_line_parms(line, directive, nparms_min, nparms_max):
        line_parms = [p.strip() for p in line[(len(directive)):].split(',')]
        valid_num_line_parms = (len(line_parms) >= nparms_min)
        if nparms_max is not None:  # ignore nparms_max if it is None.
            valid_num_line_parms = valid_num_line_parms and (len(line_parms) <= nparms_max)
        if valid_num_line_parms:
            return line_parms, None
        else:
            return None, 'Line has wrong number of parameters: \'' + line + '\'.'

    for line in lines:
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
            parms, warning_line = get_line_parms(line, "#IMAGE", 2, 2)
            if parms is not None:
                jd_start = float(parms[0])
                jd_end = float(parms[1])
                rows_to_omit = (df_filtered['JD_mid'] >= jd_start) & \
                          (df_filtered['JD_mid'] == jd_end)
            else:
                warning_lines.append(warning_line)
        else:
            warning_lines.append('Directive not understood: \'' + line + '\'.')

        if sum(rows_to_omit) >= 1:
            df_filtered = df_filtered[~ rows_to_omit]  # remove rows as user requested.
        else:
            warning_lines.append('No rows omitted: \'' + line + '\'.')

    for warning_line in warning_lines:
        print(warning_line)
    print('\n_apply_omit_txt() removed ' + str(len(df_master) - len(df_filtered)) + ' rows, ' +
          str(len(df_filtered)) + ' rows remain.')
    return df_filtered, warning_lines


def write_omit_txt_stub(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None):
    lines = [';----- This is omit.txt for AN directory ' + an_rel_directory,
             ';----- Use this file to omit observations from input to SkyModel (all filters).',
             ';----- Example directive lines:',
             ';',
             ';#OBS    Obj-0000-V, 132 ; to omit star 132 from FITS image Obj-0000-V.fts',
             ';#STAR   Obj, 132, V     ; to omit star 132 from all FITS with object Obj'
             'and filter V',
             ';#STAR   Obj, 132        ; to omit star 132 from all FITS with object Obj'
             'and ALL filters',
             ';#IMAGE  Obj-0000-V      ; to omit FITS image Obj-0000-V specifically',
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


def convert_pandas_to_json_compatible_dict(pandas_obj):
        """
        Python's json package cannot handle int64, etc that pandas insists on putting into
            its Series and DataFrame constructs. So have them write themselves to JSON text,
            then have json read the text back in as a dictionary ready for inclusion in an 
            object that json package considers legitimate to write to JSON text. Whew.
        :return: dict representation of pandas object, ready for inclusion in a JSON file [py dict]
        """
        json_text = pandas_obj.to_json()
        return json.loads(json_text)  # a dictionary without int64, etc

