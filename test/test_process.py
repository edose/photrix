import os
import shutil

import numpy as np
import pandas as pd
import pytest

from photrix import process
from photrix.user import Instrument, Site
from photrix.util import MixedModelFit

PHOTRIX_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_TOP_DIRECTORY = os.path.join(PHOTRIX_ROOT_DIRECTORY, "test")


def test__write_omit_txt_stub():
    an_top_directory = TEST_TOP_DIRECTORY
    an_rel_directory = '$an_for_test'
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit.txt')

    # Case: omit.txt does not already exist:
    savepath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit-SAVE.txt')
    shutil.copy2(fullpath, savepath)  # make a copy to restore later.
    if os.path.exists(fullpath):
        os.remove(fullpath)
    assert os.path.exists(fullpath) is False
    lines_written = process.write_omit_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
                                                an_rel_directory=an_rel_directory)
    assert os.path.exists(fullpath) is True
    with open(fullpath, 'r') as f:
        lines = f.readlines()
    assert len(lines) == lines_written == 13
    assert all([line.startswith(';') for line in lines])
    if os.path.exists(fullpath):
        os.remove(fullpath)
    shutil.move(savepath, fullpath)  # restore saved copy.

    # Case: omit.txt does already exist (written just above):
    lines_written = process.write_omit_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
                                                an_rel_directory=an_rel_directory)
    assert os.path.exists(fullpath) is True
    assert lines_written == 0
    with open(fullpath, 'r') as f:
        lines = f.readlines()
    assert all([line.startswith(';') for line in lines[0:13]])  # in case addl lines have been added


def test__write_stare_comps_txt_stub():
    an_top_directory = TEST_TOP_DIRECTORY
    an_rel_directory = '$an_for_test'
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'stare_comps.txt')

    # Case: stare_comps.txt does not already exist:
    savepath = os.path.join(an_top_directory, an_rel_directory,
                            'Photometry', 'stare_comps-SAVE.txt')
    shutil.copy2(fullpath, savepath)  # make a copy to restore later.
    if os.path.exists(fullpath):
        os.remove(fullpath)
    assert os.path.exists(fullpath) is False
    lines_written = process.write_stare_comps_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
                                                       an_rel_directory=an_rel_directory)
    assert os.path.exists(fullpath) is True
    assert lines_written == 8
    with open(fullpath, 'r') as f:
        lines = f.readlines()
    assert all([line.startswith(';') for line in lines])
    if os.path.exists(fullpath):
        os.remove(fullpath)
    shutil.move(savepath, fullpath)  # restore saved copy.

    # Case: stare_comps.txt does already exist (as written just above):
    lines_written = process.write_stare_comps_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
                                                       an_rel_directory=an_rel_directory)
    assert os.path.exists(fullpath) is True
    assert lines_written == 0
    with open(fullpath, 'r') as f:
        lines = f.readlines()
    assert all([line.startswith(';') for line in lines[0:7]])  # in case addl lines have been added


def test_get_df_master():
    an_top_directory = TEST_TOP_DIRECTORY
    an_rel_directory = '$an_for_test'
    df_master = process.get_df_master(an_top_directory, an_rel_directory)
    assert len(df_master) == 2285
    assert len(df_master.columns) == 35
    assert all([col in df_master.columns for col in ['Serial', 'FITSfile', 'Filter']])
    assert list(df_master.index) == list(df_master['Serial'])


def test_apply_omit_txt():
    an_top_directory = TEST_TOP_DIRECTORY
    an_rel_directory = '$an_for_test'

    # Nested function
    def do_apply_omit_lines(directive_lines):
        fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit.txt')
        # Write new omit.txt with test directive lines:
        _overwrite_omit_txt(an_top_directory=TEST_TOP_DIRECTORY,
                            an_rel_directory=an_rel_directory,
                            directive_lines=directive_lines)
        # Make & return output data:
        df_filtered, warning_lines = process.apply_omit_txt(an_top_directory=TEST_TOP_DIRECTORY,
                                                            an_rel_directory=an_rel_directory)
        return df_filtered, warning_lines

    # Case: #OBS directives:
    df_master = process.get_df_master(an_top_directory=TEST_TOP_DIRECTORY,
                                      an_rel_directory=an_rel_directory)  # fresh copy
    df_filtered, warning_lines = \
        do_apply_omit_lines(['#OBS   Std_SA107-0001-V, 123  ;  one existing star obs',
                             '#OBS   RT Oph-0002-R, 152     ;  another existing star obs',
                             '#OBS   RZ Mon-0001-V, 16      ;  one non-existant star obs (warning)',
                             '#CRAZY_DIRECTIVE  XXX,999     ;  raise warning'
                           ])
    assert set(df_filtered.columns) == set(df_master.columns)
    assert set(df_master.Serial) - set(df_filtered.Serial) == set([671, 1797])
    assert len(warning_lines) == 2
    assert warning_lines[0].find('#OBS   RZ Mon-0001-V') != -1  # -1 is returned if a is not in b
    assert warning_lines[1].find('#CRAZY_DIRECTIVE ') != -1  # -1 is returned if a is not in b

    # Case: #STAR directives:
    df_master = process.get_df_master(an_top_directory=TEST_TOP_DIRECTORY,
                                      an_rel_directory=an_rel_directory)  # fresh copy
    df_filtered, warning_lines = \
        do_apply_omit_lines(['#STAR  Std_SA107, 123, V  ;  one existing star, V-filter only',
                             '#STAR  RU Lyn, 136        ;  one existing star, all filters',
                             '#STAR  RZ Mon, 16         ;  one non-existant star (warning)'
                           ])
    assert set(df_filtered.columns) == set(df_master.columns)
    assert set(df_master.Serial) - set(df_filtered.Serial) == set([671] +
                                                                  [709, 722, 735, 747] + [])
    assert len(warning_lines) == 1
    assert warning_lines[0].find('#STAR  RZ Mon, 16') != -1  # -1 is returned if a is not in b

    # Case: #IMAGE directives:
    df_master = process.get_df_master(an_top_directory=TEST_TOP_DIRECTORY,
                                      an_rel_directory=an_rel_directory)  # fresh copy
    df_filtered, warning_lines = \
        do_apply_omit_lines(['#IMAGE  Std_SA107-0002-R  ;  one existing image',
                             '#IMAGE  RU Lyn-0999-X     ;  non-existent image (warning)'
                           ])
    assert set(df_filtered.columns) == set(df_master.columns)
    assert set(df_master.Serial) - set(df_filtered.Serial) == set(list(range(683, 696)) + [])
    assert len(warning_lines) == 1
    assert warning_lines[0].find('#IMAGE  RU Lyn-0999-X') != -1  # -1 is returned if a is not in b

    # Case: #JD directives:
    df_master = process.get_df_master(an_top_directory=TEST_TOP_DIRECTORY,
                                      an_rel_directory=an_rel_directory)  # fresh copy
    df_filtered, warning_lines = \
        do_apply_omit_lines(['#JD 0.668, 0.6695   ;  two images',
                             '#JD 0.2, 0.3        ;  no images (warning)'])
    assert set(df_filtered.columns) == set(df_master.columns)
    assert set(df_master.Serial) - set(df_filtered.Serial) == set(list(range(489, 512)) + [])
    assert len(warning_lines) == 1
    assert warning_lines[0].find('#JD 0.2, 0.3') != -1  # -1 is returned if a is not in b

    # Case: #SERIAL directives:
    df_master = process.get_df_master(an_top_directory=TEST_TOP_DIRECTORY,
                                      an_rel_directory=an_rel_directory)  # fresh copy
    df_filtered, warning_lines = \
        do_apply_omit_lines(['#SERIAL 99999                            ;  no images (warning)',
                             '#SERIAL  12, 14,15,4444 677,, 5 777 765,2000, 14  ; all ex 4444 OK'])
    assert set(df_filtered.columns) == set(df_master.columns)
    assert set(df_master.Serial) - set(df_filtered.Serial) == \
           set([12, 14, 15, 677, 5, 777, 765, 2000, 14])  # 4444 excluded; not an actual Serial #
    assert len(warning_lines) == 1
    assert warning_lines[0].find('#SERIAL 99999') != -1  # -1 is returned if a is not in b


def test_class_skymodel():
    an_top_directory = TEST_TOP_DIRECTORY
    an_rel_directory = '$an_for_test'
    test_filter = 'V'
    directive_lines = ['#SERIAL  348 203 1884 678 182 177 1653 1880 ;  V outliers',
                       '#IMAGE   QZ Aql-0001-V  ;  crazy cirrus term',
                       '#SERIAL  352 690  ;  R outliers',
                       '#SERIAL  703 875 193  ;  I outliers']  # in actual R processing, for comp.
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit.txt')

    # Write new omit.txt with test directive lines:
    process.write_omit_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
                                an_rel_directory=an_rel_directory)
    with open(fullpath) as f:
        lines = f.readlines()
        directive_lines = [line + '\n' for line in directive_lines]
    lines.extend(directive_lines)
    with open(fullpath, 'w') as f:
        f.writelines(lines)

    modelV = process.SkyModel(an_top_directory=TEST_TOP_DIRECTORY,
                              an_rel_directory=an_rel_directory, filter=test_filter,
                              fit_extinction=False, do_plots=False)

    # Test attributes from inputs:
    assert modelV.an_top_directory == TEST_TOP_DIRECTORY
    assert modelV.an_rel_directory == an_rel_directory
    assert modelV.filter == test_filter
    assert modelV.instrument_name == 'Borea'  # default
    assert modelV.site_name == 'DSW'  # default
    assert modelV.max_cat_mag_error == 0.01  # default
    assert modelV.max_inst_mag_sigma == 0.03  # default
    assert modelV.max_color_vi == +2.5  # default
    assert modelV.saturation_adu == Instrument(modelV.instrument_name).camera['saturation_adu']
    assert modelV.fit_sky_bias is True
    assert modelV.fit_vignette is True
    assert modelV.fit_xy is False
    assert modelV.fit_transform is False
    assert modelV.fit_extinction is False

    # Test results attributes:
    assert modelV.converged is True
    assert modelV.n_obs == 96
    assert len(modelV.df_model) == modelV.n_obs
    assert modelV.n_images == 18
    assert len(modelV.df_model['FITSfile'].drop_duplicates()) == modelV.n_images
    assert isinstance(modelV.mm_fit, MixedModelFit)
    assert len(modelV.mm_fit.df_fixed_effects) == 3
    assert modelV.transform ==\
           (Instrument(modelV.instrument_name)).filters[modelV.filter]['transform']['V-I']
    assert modelV.extinction == Site(modelV.site_name).extinction[modelV.filter]
    assert modelV.vignette == pytest.approx(-0.00578, abs=0.00005)
    assert modelV.x == 0
    assert modelV.y == 0
    assert modelV.sky_bias == pytest.approx(0.6630, abs=0.0001)
    assert modelV.sigma == pytest.approx(0.0136, abs=0.0005)

    # Test SkyModel.predict_fixed_only():
    # n_rows = 3
    # df_input = modelV.df_model[:n_rows]  # first rows
    df_input = pd.DataFrame({'Serial': [9997, 9998, 9999],
                             'SkyBias': [0.55, 0.9, 0.5454],
                             'Vignette': [0.322, 0, 1],
                             'CI': [0.577, 2.2, 0.12],
                             'Airmass': [1.57, 1.0, 2.1],
                             'FITSfile': ['BG Gem-0001-V.fts', 'BG Gem-0001-V.fts',
                                          'Std_SA35-0001-V.fts'],
                             'InstMag': [-7.68043698, -10.7139893, -6.500945076]},
                            index=['9997a', '9998a', '9999a'])
    expected_star_mags = [12.34, 9.2, 13.333]  # ideal CatMags
    mag_predictions_fixed_only = modelV.predict_fixed_only(df_input)
    random_effect_values = modelV.df_image.loc[df_input['FITSfile'], 'Value']
    # Remember: we SUBTRACT random effects (because original fit was
    #    InstMag ~ CatMag + Random Effects + offsets + fixed effects:
    mag_predictions = mag_predictions_fixed_only.values - random_effect_values.values
    assert list(mag_predictions) == pytest.approx(expected_star_mags)


def test_curate_stare_comps():
    an_top_directory = TEST_TOP_DIRECTORY
    an_rel_directory = '$an_for_test'

    # Nested function
    def do_apply_stare_comps_lines(directive_lines):
        filename = 'stare_comps.txt'
        fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', filename)

        # Write new stare_comps.txt with test directive lines:
        process.write_stare_comps_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
                                           an_rel_directory=an_rel_directory)
        with open(fullpath) as f:
            lines = f.readlines()
            directive_lines = [line + '\n' for line in directive_lines]
        lines.extend(directive_lines)
        with open(fullpath, 'w') as f:
            f.writelines(lines)

        # Make output dataframe (apply_omit_txt():
        df_eligible_obs, _ = process.apply_omit_txt(an_top_directory=TEST_TOP_DIRECTORY,
                                                    an_rel_directory=an_rel_directory)
        df_curated_obs, warning_lines = process.curate_stare_comps(
            an_top_directory=TEST_TOP_DIRECTORY,
            an_rel_directory=an_rel_directory,
            df_in=df_eligible_obs)
        return df_eligible_obs, df_curated_obs, warning_lines

    # Case: #COMPS directive:
    #   (These are NOT the comps originally removed in processing this AN;
    #       they are set to make a better test.)
    df_eligible_obs, df_curated_obs, warning_lines = do_apply_stare_comps_lines(
        ['#COMPS  V1023 Her , V , 117,120,123  ;  one FOV, keep 3 comps',
         '#CRAZY_DIRECTIVE  XXX,999     ;  raise warning'
         ])
    rows_expected_removed = (df_eligible_obs['FOV'] == 'V1023 Her') &\
                            (df_eligible_obs['StarType'] == 'Comp') &\
                            (df_eligible_obs['Filter'] == 'V') &\
                            (~ df_eligible_obs['StarID'].isin(['117', '120', '123']))
    assert sum(rows_expected_removed) == 90  # a check before proceeding
    serials_expected_removed = set((df_eligible_obs[rows_expected_removed])['Serial'])
    serials_actually_removed = set(df_eligible_obs['Serial']) - set(df_curated_obs['Serial'])
    assert serials_actually_removed == serials_expected_removed
    starids_removed = df_eligible_obs.loc[list(serials_actually_removed), 'StarID']
    assert all(starids_removed.isin(['111']))
    assert len(warning_lines) == 1
    assert warning_lines[0].find('Directive not understood') != -1  # -1 returned if a is not in b
    assert warning_lines[0].find('#CRAZY_DIRECTIVE ') != -1  # -1 returned if a is not in b
    # end of test_curate_stare_comps().


def test_class_predictionset():
    an_top_directory = TEST_TOP_DIRECTORY
    an_rel_directory = '$an_for_test'

    # Ensure omit.txt and stare_comps.txt are set up before we start (no backups)
    #    with directives used in original R processing of 20170504:
    _overwrite_omit_txt(an_top_directory=an_top_directory, an_rel_directory=an_rel_directory)
    _overwrite_stare_comps_txt(an_top_directory=an_top_directory, an_rel_directory=an_rel_directory)

    # Construct skymodel objects:
    skymodel_list = []
    for test_filter in ['V', 'R', 'I']:
        skymodel_this_filter = process.SkyModel(an_top_directory=an_top_directory,
                                                an_rel_directory=an_rel_directory,
                                                filter=test_filter,
                                                fit_extinction=False, do_plots=False)
        skymodel_list.append(skymodel_this_filter)

    ps = process.PredictionSet(an_top_directory=an_top_directory,
                               an_rel_directory='$an_for_test',
                               instrument_name='Borea',
                               site_name='DSW',
                               max_inst_mag_sigma=0.05,
                               skymodel_list=skymodel_list)

    # Test basic attributes:
    assert ps.an_top_directory == an_top_directory
    assert ps.an_rel_directory == an_rel_directory
    assert ps.instrument_name == 'Borea'
    assert ps.site_name == 'DSW'
    assert ps.max_inst_mag_sigma == 0.05
    assert ps.saturation_adu == \
        Instrument(instrument_name=ps.instrument_name).camera['saturation_adu']
    assert len(ps.images_with_targets_and_comps) == 191  # matches R::images_with_comps

    # Test df_all_eligible_obs (match R::df_filtered of line 39/30, version 1.1.4):
    assert ps.df_all_eligible_obs.shape == (2259, 35)
    assert len(ps.df_all_eligible_obs['ModelStarID'].drop_duplicates()) == 464
    assert len(ps.df_all_eligible_obs['StarID'].drop_duplicates()) == 154
    assert len(ps.df_all_eligible_obs['FITSfile'].drop_duplicates()) == 217
    assert (ps.df_all_eligible_obs['JD_mid'] - 2457878.0).mean() == \
        pytest.approx(0.7671787, 0.000001)

    # Test df_all_curated_obs (match R::df_filtered of line 40/41, version 1.1.4,
    #    these will not have changed from df_all_eligible_obs, as stare_comps removed no comps):
    assert ps.df_all_curated_obs.shape == ps.df_all_eligible_obs.shape
    assert len(ps.df_all_curated_obs['ModelStarID'].drop_duplicates()) ==\
        len(ps.df_all_eligible_obs['ModelStarID'].drop_duplicates())
    assert len(ps.df_all_curated_obs['StarID'].drop_duplicates()) ==\
               len(ps.df_all_eligible_obs['StarID'].drop_duplicates())
    assert len(ps.df_all_curated_obs['FITSfile'].drop_duplicates()) == \
               len(ps.df_all_eligible_obs['FITSfile'].drop_duplicates())
    assert (ps.df_all_curated_obs['JD_mid'] - 2457878.0).mean() == \
               (ps.df_all_eligible_obs['JD_mid'] - 2457878.0).mean()

    # Test df_comp_mags (match R::df_estimates_comps of line 80/81, version 1.1.4):
    # In R: df_estimates_comps includes only images with targets.
    # In py/photrix: ps.df_comp_mags includes ALL images with eligible comps incl Std FOVs etc.
    # Test R-equivalent dataframe r_df_estimates_comps:
    rows_to_keep = [ff in ps.images_with_targets_and_comps
                    for ff in ps.df_comp_mags['FITSfile'].values]
    r_df_estimates_comps = (ps.df_comp_mags.copy()).loc[rows_to_keep, :]
    assert r_df_estimates_comps.shape == (885, 21)
    assert len(r_df_estimates_comps['ModelStarID'].drop_duplicates()) == 271
    assert len(r_df_estimates_comps['StarID'].drop_duplicates()) == 79
    assert all(r_df_estimates_comps['StarType'] == 'Comp')
    assert any(np.isnan(r_df_estimates_comps['CatMag'])) is False
    assert r_df_estimates_comps['CatMag'].mean() == pytest.approx(11.74600113, abs=0.000001)
    assert (r_df_estimates_comps['JD_mid'] - 2457878.0).mean() == \
               pytest.approx(0.7606676, 0.000001)

    # Test actual dataframe ps.df_comp_mags:
    assert ps.df_comp_mags.shape == (1126, 21)
    assert len(ps.df_comp_mags['ModelStarID'].drop_duplicates()) == 310
    assert len(ps.df_comp_mags['StarID'].drop_duplicates()) == 101
    assert any(np.isnan(ps.df_comp_mags['CatMag'])) is False
    assert ps.df_comp_mags['CatMag'].mean() == pytest.approx(11.795893, abs=0.000001)
    assert (ps.df_comp_mags['JD_mid'] - 2457878.0).mean() == \
               pytest.approx(0.755980, 0.000001)

    images_r_df_estimates_comps = set(r_df_estimates_comps['FITSfile'])
    images_df_comp_mags = set(ps.df_comp_mags['FITSfile'])
    assert images_r_df_estimates_comps == set(ps.images_with_targets_and_comps)
    assert images_df_comp_mags == set(ps.images_with_eligible_comps)

    # Rigorous test of df_comp_mags via calculation from scratch:
    df = ps.df_comp_mags
    df_comp = df.loc[df['ModelStarID'] == 'RZ Mon_145']
    this_skymodel = [s_l for s_l in skymodel_list if s_l.filter =='V'][0]
    fe = this_skymodel.mm_fit.df_fixed_effects.Value
    inst_mag = df_comp['InstMag'].iloc[0]
    dep_var_offsets = this_skymodel.transform * df_comp['CI'].iloc[0] +\
        this_skymodel.extinction * df_comp['Airmass'].iloc[0]
    raw_mm_prediction = fe.Intercept + \
        fe.SkyBias * df_comp['SkyBias'].iloc[0] +\
        fe.Vignette * df_comp['Vignette'].iloc[0]
    estimated_mag_predicted = inst_mag - dep_var_offsets - raw_mm_prediction
    # The next line does NOT include cirrus/image effect:
    assert estimated_mag_predicted == pytest.approx(df_comp['EstimatedMag'].iloc[0], abs=0.000001)

    # Test df_cirrus_effect (SUPERSET of (not =) R::df_cirrus_effect line 164/165, version 1.2.0),
    #     that is, now includes images without targets (e.g., Std FOVs):
    assert ps.df_cirrus_effect.shape == (212, 8)
    assert set(ps.df_cirrus_effect.columns) == set(['Image', 'CirrusEffect', 'CirrusSigma',
                                                    'Criterion1', 'Criterion2', 'NumCompsUsed',
                                                    'CompIDsUsed', 'NumCompsRemoved'])
    assert ps.df_cirrus_effect['NumCompsRemoved'].sum() == 24
    assert ps.df_cirrus_effect.loc['SS Gem-0003-I.fts', 'CompIDsUsed'] == '104,110,113,95'
    assert ps.df_cirrus_effect.loc['Std_SA107-0003-I.fts', 'NumCompsUsed'] == 11
    assert ps.df_cirrus_effect.loc['Std_SA107-0003-I.fts', 'NumCompsRemoved'] == 1
    # Test the SUBSET of df_cirrus_effect with only images w/ targets (matching R:df_cirrus_effect):
    rows_with_targets = [(im in ps.images_with_targets_and_comps)
                         for im in ps.df_cirrus_effect['Image']]
    df_cirrus_effect_with_targets = ps.df_cirrus_effect[rows_with_targets]
    assert df_cirrus_effect_with_targets.shape == (191, 8)
    assert set(df_cirrus_effect_with_targets.columns) == \
           set(['Image', 'CirrusEffect', 'CirrusSigma',
                'Criterion1', 'Criterion2', 'NumCompsUsed',
                'CompIDsUsed', 'NumCompsRemoved'])
    assert df_cirrus_effect_with_targets['NumCompsRemoved'].sum() == 17
    assert df_cirrus_effect_with_targets.loc['SS Gem-0003-I.fts', 'CompIDsUsed'] == '104,110,113,95'
    assert 'Std...' not in df_cirrus_effect_with_targets['Image']



# ---------------  INTERNAL TEST-HELPER FUNCTIONS ----------------------------------------------

def _overwrite_omit_txt(an_top_directory, an_rel_directory, directive_lines=None):
    # No backup. Just write it:
    header = [';----- This is omit.txt for AN directory ' + an_rel_directory,
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
    if directive_lines is None:
        # The original 20170504 lines:
        directive_lines_to_write = ['#SERIAL  348 203 1884 678 182 177 1653 1880 ;  V outliers',
                                    '#IMAGE   QZ Aql-0001-V  ;  crazy cirrus term',
                                    '#SERIAL  352 690  ;  R outliers',
                                    '#SERIAL  703 875 193  ;  I outliers']
    else:
        directive_lines_to_write = directive_lines
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit.txt')
    all_text = [line + '\n' for line in (header + directive_lines_to_write)]
    with open(fullpath, 'w') as f:
        f.writelines(all_text)


def _overwrite_stare_comps_txt(an_top_directory, an_rel_directory, directive_lines=None):
    # No backup. Just write it:
    #    (These ARE the comps originally kept in processing this AN. I think there are no
    #        comps removed (kept them all), but that's OK since a restrictive set was tested
    #        in test_stare_comps() above.)
    header = [';----- This is stare_comps.txt for AN directory ' + an_rel_directory,
              ';----- Select comp stars (by FOV, filter, & StarID) from input to predict_fixed_only().',
              ';----- Example directive line:',
              ';',
              ';#COMPS  Obj, V, 132, 133 144    ; to KEEP from FOV \'Obj\': '
              'comp stars \'132\' \'133\' and \'144\' in filter \'V\'',
              ';',
              ';----- Add your directive lines:',
              ';']
    if directive_lines is None:
        # The original directive lines in AN20170504.
        directive_lines_to_write = ['#COMPS  V1023 Her , V , 117,120,111',
                                    '#COMPS  V1023 Her , I , 117,120,111']
    else:
        directive_lines_to_write = directive_lines
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'stare_comps.txt')

    all_text = [line + '\n' for line in (header + directive_lines_to_write)]
    with open(fullpath, 'w') as f:
        f.writelines(all_text)
