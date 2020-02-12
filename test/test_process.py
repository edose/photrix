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


def test__get_line_parms():
    assert process._get_line_parms('#SERIAL  12, 34 44 , 42  ;  this is a comment',
                                   '#SERIAL', True, 1, None) == (['12', '34', '44', '42'], None)
    assert process._get_line_parms('#SERIAL  ST Tri-0000-V  ;  this is a comment',
                                   '#SERIAL', False, 1, None) == (['ST Tri-0000-V'], None)
    assert process._get_line_parms('#JD  0.25 0.5 ;  this is a comment',
                                   '#JD', True, 2, 2) == (['0.25', '0.5'], None)
    parms = process._get_line_parms('#JD  0.3  ;  this is a comment', '#JD', True, 2, 2)
    assert parms[0] is None
    assert 'wrong number of parameters: ' in parms[1]
    parms = process._get_line_parms('#JX  0.3  ;  this is a comment', '#JD', True, 2, 2)
    assert parms[0] is None
    assert 'does not begin with correct directive: ' in parms[1]


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
    lines_written = process._write_omit_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
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
    lines_written = process._write_omit_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
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
    lines_written = process._write_stare_comps_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
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
    lines_written = process._write_stare_comps_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
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
    assert len(df_master.columns) == 36
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
        df_filtered, warning_lines = process._apply_omit_txt(an_top_directory=TEST_TOP_DIRECTORY,
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
    # TODO: add tests for fixed-effect log_adu (parallel to FE sky_bias).
    an_top_directory = TEST_TOP_DIRECTORY
    an_rel_directory = '$an_for_test'
    test_filter = 'V'
    directive_lines = ['#SERIAL  348 203 1884 678 182 177 1653 1880 ;  V outliers',
                       '#IMAGE   QZ Aql-0001-V  ;  crazy cirrus term',
                       '#SERIAL  352 690  ;  R outliers',
                       '#SERIAL  703 875 193  ;  I outliers']  # in actual R processing, for comp.
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit.txt')

    # Write new omit.txt with test directive lines:
    process._write_omit_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
                                 an_rel_directory=an_rel_directory)
    with open(fullpath) as f:
        lines = f.readlines()
        directive_lines = [line + '\n' for line in directive_lines]
    lines.extend(directive_lines)
    with open(fullpath, 'w') as f:
        f.writelines(lines)

    # Case 1: Model WITHOUT log_adu (CCD nonlinearity) term.
    modelV = process.SkyModel(an_top_directory=TEST_TOP_DIRECTORY,
                              an_rel_directory=an_rel_directory, filter=test_filter,
                              fit_extinction=False, fit_log_adu=False, do_plots=False)
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
    assert modelV.fit_log_adu is False
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
        (Instrument(modelV.instrument_name)).transform(modelV.filter, 'V-I')
    assert modelV.extinction == Site(modelV.site_name).extinction[modelV.filter]
    assert modelV.vignette == pytest.approx(-0.00603, abs=0.0001)  # changed
    assert modelV.x == 0
    assert modelV.y == 0
    assert modelV.sky_bias == pytest.approx(0.6671, abs=0.001)
    assert modelV.log_adu == 0
    assert modelV.sigma == pytest.approx(0.0143, abs=0.001)
    # Test SkyModel._predict_fixed_only():
    df_input = pd.DataFrame({'Serial': [9997, 9998, 9999],
                             'SkyBias': [0.55, 0.9, 0.5454],
                             'Vignette': [0.322, 0, 1],
                             'CI': [0.577, 2.2, 0.12],
                             'Airmass': [1.57, 1.0, 2.1],
                             'FITSfile': ['BG Gem-0001-V.fts', 'BG Gem-0001-V.fts',
                                          'Std_SA35-0001-V.fts'],
                             'InstMag': [-7.68043698, -10.7139893, -6.500945076]},
                            index=['9997a', '9998a', '9999a'])
    expected_star_mags = [12.3553, 9.2302, 13.2887]  # ideal CatMags
    mag_predictions_fixed_only = modelV._predict_fixed_only(df_input)
    random_effect_values = modelV.df_image.loc[df_input['FITSfile'], 'Value']
    # Remember: we SUBTRACT random effects (because original fit was
    #    InstMag ~ CatMag + Random Effects + offsets + fixed effects:
    mag_predictions = mag_predictions_fixed_only.values - random_effect_values.values
    assert list(mag_predictions) == pytest.approx(expected_star_mags, abs=0.0005)

    # Case 2: Model WITH log_adu (CCD nonlinearity) term.
    modelV = process.SkyModel(an_top_directory=TEST_TOP_DIRECTORY,
                              an_rel_directory=an_rel_directory, filter=test_filter,
                              fit_extinction=False, fit_log_adu=True, do_plots=False)
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
    assert modelV.fit_log_adu is True
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
    assert len(modelV.mm_fit.df_fixed_effects) == 4
    assert modelV.transform ==\
        (Instrument(modelV.instrument_name)).transform(modelV.filter, 'V-I')
    assert modelV.extinction == Site(modelV.site_name).extinction[modelV.filter]
    assert modelV.vignette == pytest.approx(-0.00050, abs=0.0001)
    assert modelV.x == 0
    assert modelV.y == 0
    assert modelV.sky_bias == pytest.approx(0.5500, abs=0.001)
    assert modelV.log_adu == pytest.approx(-0.0284, abs=0.001)
    assert modelV.sigma == pytest.approx(0.0135, abs=0.001)
    # Test SkyModel._predict_fixed_only():
    df_input = pd.DataFrame({'Serial': [9997, 9998, 9999],
                             'SkyBias': [0.55, 0.9, 0.5454],
                             'LogADU': [3.1, 3.5, 3.22],
                             'Vignette': [0.322, 0, 1],
                             'CI': [0.577, 2.2, 0.12],
                             'Airmass': [1.57, 1.0, 2.1],
                             'FITSfile': ['BG Gem-0001-V.fts', 'BG Gem-0001-V.fts',
                                          'Std_SA35-0001-V.fts'],
                             'InstMag': [-7.68043698, -10.7139893, -6.500945076]},
                            index=['9997a', '9998a', '9999a'])
    expected_star_mags = [12.3921, 9.3211, 13.3284]  # ideal CatMags
    mag_predictions_fixed_only = modelV._predict_fixed_only(df_input)
    random_effect_values = modelV.df_image.loc[df_input['FITSfile'], 'Value']
    # Remember: we SUBTRACT random effects (because original fit was
    #    InstMag ~ CatMag + Random Effects + offsets + fixed effects:
    mag_predictions = mag_predictions_fixed_only.values - random_effect_values.values
    assert list(mag_predictions) == pytest.approx(expected_star_mags, abs=0.001)


def test_curate_stare_comps():
    an_top_directory = TEST_TOP_DIRECTORY
    an_rel_directory = '$an_for_test'

    # Nested function for testing:  ----------------------------------------------
    def _do_apply_stare_comps_lines(directive_lines):
        filename = 'stare_comps.txt'
        fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', filename)

        # Write new stare_comps.txt with test directive lines:
        os.remove(fullpath)
        process._write_stare_comps_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
                                            an_rel_directory=an_rel_directory)
        with open(fullpath) as f:
            lines = f.readlines()
            directive_lines = [line + '\n' for line in directive_lines]
        lines.extend(directive_lines)
        with open(fullpath, 'w') as f:
            f.writelines(lines)

        # Make output dataframe (_apply_omit_txt():
        df_eligible_obs, _ = process._apply_omit_txt(an_top_directory=TEST_TOP_DIRECTORY,
                                                     an_rel_directory=an_rel_directory)
        df_curated_obs, warning_lines = process._curate_stare_comps(
            an_top_directory=TEST_TOP_DIRECTORY,
            an_rel_directory=an_rel_directory,
            df_in=df_eligible_obs)
        return df_eligible_obs, df_curated_obs, warning_lines
    # -------------------------------------------------------------------------------

    # Case: #COMPS directive:
    #   (These are NOT the comps originally removed in processing this AN;
    #       they are set to make a better test.)
    df_eligible_obs, df_curated_obs, warning_lines = _do_apply_stare_comps_lines(
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


def test_solve_for_real_ci():
    # Start with ideal values, make observed values, then test that fn recovers ideal Color Index:
    ideal_mags = {'V': 12.5, 'I': 8.8}
    ci_filters = ['V', 'I']
    transforms = {'I': -0.044, 'V': 0.025}
    ideal_ci = ideal_mags['V'] - ideal_mags['I']
    untransformed_mags = {cif: ideal_mags[cif] + transforms[cif] * ideal_ci for cif in ci_filters}
    real_ci = process._solve_for_real_ci(untransformed_mags, ci_filters, transforms)
    assert real_ci == pytest.approx(ideal_ci, abs=0.000001)


def test_extract_ci_points():
    # Case 1:
    df_star_id = pd.DataFrame({'Serial': [2, 4, 6, 8],
                               'ModelStarID': 'AU XXX_111',
                               'Filter': ['V', 'R', 'I', 'I'],
                               'JD_num': [0.66, 0.67, 0.68, 0.69],
                               'CI': 0.0,
                               'UntransformedMag': [12.8, 10.2, 8.3, 8.2]})
    ci_filters = ['V', 'I']
    transforms = {'V': 0.025, 'R': -0.08, 'I': 0.052, 'XX': 1.0}
    df_result = process._extract_ci_points(df_star_id=df_star_id, ci_filters=ci_filters,
                                           transforms=transforms)
    assert len(df_result) == 1
    assert df_result.loc[0, 'JD_num'] == \
           (df_star_id.loc[0, 'JD_num'] + df_star_id.loc[2, 'JD_num']) / 2.0
    assert df_result.loc[0, 'CI'] == pytest.approx((12.8-8.3) / (1.0 + 0.025 - 0.052), abs=0.000001)

    # Case 2:
    df_star_id = pd.DataFrame({'Serial': np.arange(1, 9),
                               'ModelStarID': 'AU XXX_111',
                               'Filter': ['V', 'R', 'I', 'I', 'V', 'I', 'V', 'R'],
                               'JD_num': 0.66 + np.arange(8)*0.01,
                               'CI': 0.0,
                               'UntransformedMag': [12.8, 10.2, 8.3, 8.2, 12.6, 8.1, 12.5, 10.11]})
    ci_filters = ['V', 'I']
    transforms = {'V': 0.025, 'R': -0.08, 'I': 0.052, 'XX': 1.0}
    df_result = process._extract_ci_points(df_star_id=df_star_id, ci_filters=ci_filters,
                                           transforms=transforms)
    assert len(df_result) == 4
    assert list(df_result['JD_num']) == pytest.approx([0.67, 0.695, 0.705, 0.715], abs=0.00001)
    untransformed_ci = [12.8-8.3, 12.6-8.2, 12.6-8.1, 12.5-8.1]
    real_ci = [uci / (1.0 + 0.025 - 0.052) for uci in untransformed_ci]
    assert list(df_result['CI']) == pytest.approx(real_ci, abs=0.000001)


def test_impute_target_ci():
    # Make a df for each interpolation case, concatenate them into one df, and test function:
    df_1_ci_point = pd.DataFrame({'Serial': 1 + np.arange(4),
                                  'ModelStarID': 'AU_1_PT',
                                  'StarType': 'Target',
                                  'Filter': ['V', 'R', 'I', 'I'],
                                  'JD_mid': 2457800.35 + 0.01*np.arange(4),
                                  'CI': 0.0,
                                  'UntransformedMag': [12.8, 10.2, 8.3, 8.2]})
    df_2_ci_points = pd.DataFrame({'Serial': 10 + np.arange(6),
                                   'ModelStarID': 'AU_2_PTS',
                                   'StarType': 'Target',
                                   'Filter': ['V', 'R', 'I', 'I', 'V', 'V'],
                                   'JD_mid': 2457800.45 + 0.01*np.arange(6),
                                   'CI': 0.0,
                                   'UntransformedMag': [12.9, 10.5, 8.7, 8.9, 12.6, 12.7]})
    df_3_ci_points = pd.DataFrame({'Serial': 20 + np.arange(8),
                                   'ModelStarID': 'AU_3_PTS',
                                   'StarType': 'Check',
                                   'Filter': ['V', 'R', 'I', 'I', 'V', 'V', 'R', 'I'],
                                   'JD_mid': 2457800.45 + 0.01*np.arange(8),
                                   'CI': 0.0,
                                   'UntransformedMag': [12.5, 10.5, 8.5, 8.4, 12.4, 12.3, 10, 8.8]})
    df_4_ci_points = pd.DataFrame({'Serial': 40 + np.arange(10),
                                   'ModelStarID': 'AU_4_PTS',
                                   'StarType': 'Check',
                                   'Filter': ['V', 'R', 'I', 'I', 'V', 'V', 'R', 'I', 'R', 'V'],
                                   'JD_mid': 2457800.55 + 0.005*np.arange(10),
                                   'CI': 0.0,
                                   'UntransformedMag': [12.8, 10.2, 8.3, 8.2, 12.2, 12.3, 10.2,
                                                        8.8, 10.2, 12.5]})
    df_7_ci_points = pd.DataFrame({'Serial': 60 + np.arange(11),
                                   'ModelStarID': 'AU_7_PTS',
                                   'StarType': 'Check',
                                   'Filter': ['V', 'I', 'V', 'V', 'I', 'V', 'V', 'I',
                                              'V', 'V', 'I'],
                                   'JD_mid': 2457800.65 + 0.005*np.arange(11),
                                   'CI': 0.0,
                                   'UntransformedMag': [12.6, 8.3, 12.2, 12.8, 8.3, 12.2,
                                                        12.8, 8.3, 12.2, 12.8, 8.4]})
    df_no_ci_points = pd.DataFrame({'Serial': 80 + np.arange(4),
                                    'ModelStarID': 'AU_NO_PTS',
                                    'StarType': 'Target',
                                    'Filter': ['V', 'R', 'V', 'V'],
                                    'JD_mid': 2457800.75 + 0.02*np.arange(4),
                                    'CI': 0.0,
                                    'UntransformedMag': [12.8, 10.2, 8.3, 8.2]})
    df_predictions_checks_targets = pd.concat([df_1_ci_point, df_2_ci_points, df_3_ci_points,
                                              df_4_ci_points, df_7_ci_points, df_no_ci_points])
    df_predictions_checks_targets.index = df_predictions_checks_targets['Serial']
    ci_filters = ['V', 'I']
    transforms = {'V': 0.025, 'R': -0.08, 'I': 0.052, 'XX': 1.0}
    df_updated = process._impute_target_ci(df_predictions_checks_targets, ci_filters, transforms)
    factor = 1 / (1 + transforms['V'] - transforms['I'])  # for solving for real CI from obs CI

    # Test case: 1 CI point:
    df_1 = df_updated[df_updated['ModelStarID'] == 'AU_1_PT']  # fn output for this star
    # Verify no change
    assert list(df_1['Serial']) == list(1 + np.arange(4))
    assert list(df_1['StarType']) == len(df_1) * ['Target']
    assert list(df_1['JD_mid']) == list(df_1_ci_point['JD_mid'])
    assert list(df_1['UntransformedMag']) == list(df_1_ci_point['UntransformedMag'])
    ci_points = factor * (df_1_ci_point.iloc[0]['UntransformedMag'] -
                            df_1_ci_point.iloc[2]['UntransformedMag'])
    # Verify Color Index values:
    assert list(df_1['CI']) == len(df_1) * [pytest.approx(ci_points, abs=0.0001)]

    # Test case: 2 CI points:
    df_2 = df_updated[df_updated['ModelStarID'] == 'AU_2_PTS']  # fn output for this star
    assert list(df_2['Serial']) == list(10 + np.arange(6))
    assert list(df_2['StarType']) == len(df_2) * ['Target']
    assert list(df_2['JD_mid']) == list(df_2_ci_points['JD_mid'])
    assert list(df_2['UntransformedMag']) == list(df_2_ci_points['UntransformedMag'])
    ci_points = [factor * (df_2_ci_points.iloc[i]['UntransformedMag'] -
                           df_2_ci_points.iloc[j]['UntransformedMag'])
                 for (i, j) in [(0, 2), (4, 3)]]
    assert df_2.iloc[0]['CI'] == pytest.approx(ci_points[0], abs=0.0001)
    assert df_2.iloc[1]['CI'] == pytest.approx(ci_points[0], abs=0.0001)
    assert df_2.iloc[2]['CI'] == pytest.approx(4.1110, abs=0.0001)
    assert df_2.iloc[3]['CI'] == pytest.approx(3.9054, abs=0.0001)
    assert df_2.iloc[4]['CI'] == pytest.approx(ci_points[1], abs=0.0001)
    assert df_2.iloc[5]['CI'] == pytest.approx(ci_points[1], abs=0.0001)

    # Test case: 3 CI points:
    df_3 = df_updated[df_updated['ModelStarID'] == 'AU_3_PTS']  # fn output for this star
    assert list(df_3['Serial']) == list(20 + np.arange(8))
    assert list(df_3['StarType']) == len(df_3) * ['Check']
    assert list(df_3['JD_mid']) == list(df_3_ci_points['JD_mid'])
    assert list(df_3['UntransformedMag']) == list(df_3_ci_points['UntransformedMag'])
    assert df_3.iloc[0]['CI'] == df_3.iloc[1]['CI'] == pytest.approx(4.19664, abs=0.0001)
    assert df_3.iloc[2]['CI'] == pytest.approx(4.09387, abs=0.0001)
    assert df_3.iloc[3]['CI'] == pytest.approx(3.9911, abs=0.0001)
    assert df_3.iloc[4]['CI'] == pytest.approx(3.8883, abs=0.0001)
    assert df_3.iloc[5]['CI'] == pytest.approx(3.7855, abs=0.0001)
    assert df_3.iloc[6]['CI'] == df_3.iloc[7]['CI'] == pytest.approx(3.68277, abs=0.0001)

    # Test case: 4 CI points (make a nice spline):
    df_4 = df_updated[df_updated['ModelStarID'] == 'AU_4_PTS']  # fn output for this star
    assert list(df_4['Serial']) == list(40 + np.arange(10))
    assert list(df_4['StarType']) == len(df_4) * ['Check']
    assert list(df_4['JD_mid']) == list(df_4_ci_points['JD_mid'])
    assert list(df_4['UntransformedMag']) == list(df_4_ci_points['UntransformedMag'])
    assert df_4.iloc[0]['CI'] == df_4.iloc[1]['CI'] == pytest.approx(4.62487, abs=0.0001)
    assert df_4.iloc[2]['CI'] == pytest.approx(4.47805, abs=0.0001)
    assert df_4.iloc[3]['CI'] == pytest.approx(4.24314, abs=0.0001)
    assert df_4.iloc[4]['CI'] == pytest.approx(3.97886, abs=0.0001)
    assert df_4.iloc[5]['CI'] == pytest.approx(3.74394, abs=0.0001)
    assert df_4.iloc[6]['CI'] == pytest.approx(3.59712, abs=0.0001)
    assert df_4.iloc[7]['CI'] == pytest.approx(3.59712, abs=0.0001)
    assert df_4.iloc[8]['CI'] == df_4.iloc[9]['CI'] == pytest.approx(3.80267, abs=0.0001)

    # Test case: 7 CI points:
    df_7 = df_updated[df_updated['ModelStarID'] == 'AU_7_PTS']  # fn output for this star
    assert list(df_7['Serial']) == list(60 + np.arange(11))
    assert list(df_7['StarType']) == len(df_7) * ['Check']
    assert list(df_7['JD_mid']) == list(df_7_ci_points['JD_mid'])
    assert list(df_7['UntransformedMag']) == list(df_7_ci_points['UntransformedMag'])
    assert df_7.iloc[0]['CI'] == pytest.approx(4.41340, abs=0.0001)
    assert df_7.iloc[2]['CI'] == pytest.approx(4.19340, abs=0.0001)
    assert df_7.iloc[5]['CI'] == pytest.approx(4.08470, abs=0.0001)
    assert df_7.iloc[8]['CI'] == pytest.approx(3.68003, abs=0.0001)
    assert df_7.iloc[10]['CI'] == pytest.approx(4.52080, abs=0.0001)

    # Test case: NO CI points at all:
    df_no = df_updated[df_updated['ModelStarID'] == 'AU_NO_PTS']  # fn output for this star
    assert list(df_no['Serial']) == list(80 + np.arange(4))
    assert list(df_no['StarType']) == len(df_no) * ['Target']
    assert list(df_no['JD_mid']) == list(df_no_ci_points['JD_mid'])
    assert list(df_no['UntransformedMag']) == list(df_no_ci_points['UntransformedMag'])
    assert all(np.isnan(df_no['CI']))


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

    # Test df_all_eligible_obs (added col LogADU in R ver 1.2.1):
    assert ps.df_all_eligible_obs.shape == (2259, 36)
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

    # Test df_comp_mags (match R::df_estimates_comps of line 80/81, version 1.1.4),
    #    added column LogADU in R ver 1.2.1.
    # In R: df_estimates_comps includes only images with targets.
    # In py/photrix: ps.df_comp_mags includes ALL images with eligible comps incl Std FOVs etc.
    # Test R-equivalent dataframe r_df_estimates_comps:
    rows_to_keep = [ff in ps.images_with_targets_and_comps
                    for ff in ps.df_comp_mags['FITSfile'].values]
    r_df_estimates_comps = (ps.df_comp_mags.copy()).loc[rows_to_keep, :]
    assert r_df_estimates_comps.shape == (885, 22)
    assert len(r_df_estimates_comps['ModelStarID'].drop_duplicates()) == 271
    assert len(r_df_estimates_comps['StarID'].drop_duplicates()) == 79
    assert all(r_df_estimates_comps['StarType'] == 'Comp')
    assert any(np.isnan(r_df_estimates_comps['CatMag'])) is False
    assert r_df_estimates_comps['CatMag'].mean() == pytest.approx(11.74600113, abs=0.000001)
    assert (r_df_estimates_comps['JD_mid'] - 2457878.0).mean() == \
               pytest.approx(0.7606676, 0.000001)

    # Test actual dataframe ps.df_comp_mags (nb: added column LogADU in R ver 1.2.1):
    assert ps.df_comp_mags.shape == (1126, 22)
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
    this_skymodel = [s_l for s_l in skymodel_list if s_l.filter == 'V'][0]
    fe = this_skymodel.mm_fit.df_fixed_effects.Value
    inst_mag = df_comp['InstMag'].iloc[0]
    dep_var_offsets = this_skymodel.transform * df_comp['CI'].iloc[0] +\
        this_skymodel.extinction * df_comp['Airmass'].iloc[0]
    raw_mm_prediction = fe.Intercept + \
        fe.SkyBias * df_comp['SkyBias'].iloc[0] +\
        fe.Vignette * df_comp['Vignette'].iloc[0] +\
        fe.LogADU * df_comp['LogADU'].iloc[0]
    estimated_mag_predicted = inst_mag - dep_var_offsets - raw_mm_prediction
    # The next line does NOT include cirrus/image effect:
    assert estimated_mag_predicted == pytest.approx(df_comp['EstimatedMag'].iloc[0], abs=0.000001)

    # Test df_cirrus_effect (SUPERSET of (not =) R::df_cirrus_effect line 164/165, version 1.2.0),
    #     that is, now includes images without targets (e.g., Std FOVs):
    assert ps.df_cirrus_effect.shape == (212, 8)
    assert set(ps.df_cirrus_effect.columns) == set(['Image', 'CirrusEffect', 'CirrusSigma',
                                                    'Criterion1', 'Criterion2', 'NumCompsUsed',
                                                    'CompIDsUsed', 'NumCompsRemoved'])
    assert ps.df_cirrus_effect['NumCompsRemoved'].sum() == 13
    assert ps.df_cirrus_effect.loc['SS Gem-0003-I.fts', 'CompIDsUsed'] == '104,110,113,95'
    assert ps.df_cirrus_effect.loc['Std_SA32-0009-I.fts', 'NumCompsUsed'] == 14
    assert ps.df_cirrus_effect.loc['Std_SA32-0009-I.fts', 'NumCompsRemoved'] == 1
    # Test the SUBSET of df_cirrus_effect with only images w/ targets (matching R:df_cirrus_effect):
    rows_with_targets = [(im in ps.images_with_targets_and_comps)
                         for im in ps.df_cirrus_effect['Image']]
    df_cirrus_effect_with_targets = ps.df_cirrus_effect[rows_with_targets]
    assert df_cirrus_effect_with_targets.shape == (191, 8)
    assert set(df_cirrus_effect_with_targets.columns) == \
           set(['Image', 'CirrusEffect', 'CirrusSigma',
                'Criterion1', 'Criterion2', 'NumCompsUsed',
                'CompIDsUsed', 'NumCompsRemoved'])
    assert df_cirrus_effect_with_targets['NumCompsRemoved'].sum() == 11
    assert df_cirrus_effect_with_targets.loc['SS Gem-0003-I.fts', 'CompIDsUsed'] == '104,110,113,95'
    assert 'Std...' not in df_cirrus_effect_with_targets['Image']

    # Test result df_transformed from _compute_transformed_mags():
    expected_columns = set(['Serial', 'ModelStarID', 'FITSfile', 'StarID', 'Chart',
                            'Xcentroid', 'Ycentroid', 'InstMag', 'InstMagSigma', 'StarType',
                            'CatMag', 'CatMagError', 'Exposure', 'JD_mid', 'Filter',
                            'Airmass', 'CI', 'SkyBias', 'Vignette', 'LogADU',
                            'UseInEnsemble', 'CirrusEffect', 'CirrusSigma', 'CompIDsUsed',
                            'Image', 'NumCompsRemoved', 'NumCompsUsed', 'JD_num',
                            'TransformedMag', 'ModelSigma', 'TotalSigma', 'FOV',
                            'MaxADU_Ur', 'FWHM', 'SkyADU', 'SkySigma'])
    assert set(ps.df_transformed.columns) == expected_columns
    assert ps.df_transformed.shape == (532, 36)
    assert list(ps.df_transformed['Serial'].iloc[[0, 10, 500]]) == [441, 332, 1589]
    assert ps.df_transformed['TransformedMag'].sum() == pytest.approx(6408.9, abs=1)  # changed
    assert ps.df_transformed['TotalSigma'].sum() == pytest.approx(10.103, abs=0.01)


def test_stare_comps():
    df_test = pd.DataFrame({'Serial': range(10),
                            'FOV': ['A', 'B', 'A', 'A', 'B', 'A', 'B', 'A', 'B', 'A'],
                            'StarID': ['Star1', 'Star1', 'Star1', 'Star2', 'Star2', 'Star2', 'Ono',
                                       'Star1', 'Star2', 'Star2'],
                            'Filter': ['V', 'X', 'V', 'V', 'V', 'X', 'V', 'V', 'V', 'V'],
                            'CompIDsUsed': ['12,23,34,45', '12,23,34', '12,34,45', '12,23,34,45',
                                            '12,23,45', '12,23,34', '12,45', '12,23,34',
                                            '12,23,34,45', '34']})
    df_test.index = df_test['Serial']

    result_A1 = process.get_stare_comps(df_test, fov='A', star_id='Star1', this_filter='V')
    assert len(result_A1) == 5
    assert result_A1[0].startswith('EDIT file ')
    assert result_A1[4].strip() == '4 comps -> 2 images qualify  -->   #COMPS A, V, 12, 23, 34, 45'

    result_A2 = process.get_stare_comps(df_test, fov='A', star_id='Star2', this_filter='V')
    assert len(result_A2) == 5
    assert result_A2[1].strip() == '1 comps -> 2 images qualify  -->   #COMPS A, V, 34'
    assert result_A2[3].strip() == '3 comps -> 1 images qualify  -->   #COMPS A, V, 12, 23, 34'

    # Case: only one image (won't happen for real stares, but an edge case):
    result_B1 = process.get_stare_comps(df_test, fov='B', star_id='Star1', this_filter='V')
    assert len(result_B1) == 2
    assert result_B1[1].strip() == '>>> One or zero qualifying images in dataframe.'


def test_transform_model():
    an_top_directory = TEST_TOP_DIRECTORY
    an_rel_directory = '$an_for_test'
    df_master = process.get_df_master(an_top_directory=an_top_directory,
                                      an_rel_directory=an_rel_directory)

    # Case: one FOV with one image:
    tm = process.TransformModel(an_top_directory=an_top_directory,
                                an_rel_directory=an_rel_directory,
                                filter='V', ci_type='V-I', fovs_to_include='RZ Mon',
                                instrument_name='Borea', site_name='DSW')
    assert tm.is_valid
    assert len(tm.image_list) == 1
    assert len(tm.fitted_values) == 12
    assert list(tm.fitted_values[:2]) == pytest.approx([-20.2926, -20.2896], abs=0.0001)
    assert tm.param_values['Intercept'] == pytest.approx(-20.1782, abs=0.0001)
    assert set(tm.param_values.index) == set(['Intercept', 'SkyBias', 'LogADU', 'CI'])
    assert len(tm.residuals) == 12
    assert list(tm.residuals[:2]) == pytest.approx([-0.00568, -0.00614], abs=0.00001)
    assert tm.transform_value == pytest.approx(-0.09770, abs=0.0001)
    assert tm.sigma == pytest.approx(0.02151, abs=0.0001)
    assert tm.image_effect is None

    # Case: one FOV with multiple images:
    tm = process.TransformModel(an_top_directory=an_top_directory,
                                an_rel_directory=an_rel_directory,
                                filter='V', ci_type='V-I', fovs_to_include='Std_SA32',
                                instrument_name='Borea', site_name='DSW')
    assert tm.is_valid
    assert len(tm.image_list) == 4
    assert all([n.startswith('Std_SA32-00') for n in tm.image_list])
    assert all([n.endswith('-V.fts') for n in tm.image_list])
    assert tm.param_values['Intercept'] == pytest.approx(-20.488572, abs=0.0001)

    # Case: list of FOVs:
    tm = process.TransformModel(an_top_directory=an_top_directory,
                                an_rel_directory=an_rel_directory,
                                filter='V', ci_type='V-I',
                                fovs_to_include=['Std_SA32', 'RZ Mon'],
                                instrument_name='Borea', site_name='DSW')
    assert tm.is_valid
    assert len(tm.image_list) == 5
    assert tm.param_values['CI'] == pytest.approx(-0.027827, abs=0.0001)

    # Case: fovs_to_include = "Standards" (test data selection only):
    tm = process.TransformModel(an_top_directory=an_top_directory,
                                an_rel_directory=an_rel_directory,
                                filter='V', ci_type='V-I', fovs_to_include='Standards',
                                instrument_name='Borea', site_name='DSW')
    assert tm.is_valid
    assert len(tm.image_list) == 7
    assert all([n.startswith('Std_') for n in tm.image_list])
    assert all([n.endswith('-V.fts') for n in tm.image_list])
    assert tm.param_values['CI'] == pytest.approx(-0.029374, abs=0.0001)

    # Case: fovs_to_include = "All" (test data selection only):
    tm = process.TransformModel(an_top_directory=an_top_directory,
                                an_rel_directory=an_rel_directory,
                                filter='V', ci_type='V-I', fovs_to_include='All',
                                instrument_name='Borea', site_name='DSW')
    assert tm.is_valid
    all_images = df_master['FITSfile'].drop_duplicates().tolist()
    all_fovs = df_master['FOV'].drop_duplicates().tolist()
    assert len(tm.image_list) == 119
    assert set(tm.image_list) <= set(all_images)
    assert len(tm.fov_list) == 26
    assert set(tm.fov_list) <= set(all_fovs)
    assert tm.param_values['CI'] == pytest.approx(-0.047128, abs=0.0001)
    assert tm.transform_value == tm.param_values['CI']
    assert tm.transform_sigma == pytest.approx(0.00606, abs=0.0001)
    assert tm.sigma == pytest.approx(0.0215, abs=0.0001)

    # Case: fovs_to_include = "All" w/ different color index (test data selection only):
    tm = process.TransformModel(an_top_directory=an_top_directory,
                                an_rel_directory=an_rel_directory,
                                filter='R', ci_type='R-I', fovs_to_include='All',
                                instrument_name='Borea', site_name='DSW')
    assert tm.is_valid
    assert len(tm.image_list) == 18
    assert tm.transform_value == pytest.approx(+0.077440, abs=0.0001)


def test_predictionset_aavso_report():
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
    df_report = ps.aavso_report(write_file=True, return_df=True)

    assert df_report.shape == (291, 17)


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
              ';----- Select comp stars (by FOV, filter, & StarID) from input to _predict_fixed_only().',
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


def test_make_df_master():
    # For now, do this in external directories. Later, set it up in photrix test subdirectories.
    # For now, we do not construct df_master.csv files here, we construct them separately,
    #     then only compare them here.

    # Read in dataframe as constructed by R software (2015-7):
    df_r = pd.read_csv('C:/Astro/Images/Borea Photrix/20170710-R/Photometry/df_master.csv', sep=';',)
    df_r.index = df_r['Serial']
    # Execute test function (photrix, python), and read dataframe back in:
    # process.make_df_master('C:/Astro/Images/Borea Photrix', '20170710-py', ask_user=False)
    df_py = pd.read_csv('C:/Astro/Images/Borea Photrix/20170710-py/Photometry/df_master.csv', sep=';')
    df_py.index = df_py['Serial']

    assert set(df_py.columns) == set(df_r.columns)  # columns must be same, differing order OK.
    assert len(df_py) == len(df_r)
    assert df_py['Serial'].tolist() == list(range(1, 1 + len(df_py)))

    # Check column dtypes are OK:
    for col in df_py.columns:
        col_types_same = (df_py[col].dtype == df_r[col].dtype)
        col_type_ok = col_types_same or\
                      (str(df_py[col].dtype).startswith('float') and
                       str(df_r[col].dtype).startswith('int'))
        if not col_type_ok:
            print('py column', col, 'is', str(df_py[col].dtype), 'but r has', df_r[col].dtype)
        assert col_type_ok

    # Check that object and integer column contents are exactly equal,
    #     and that float column contents are approximately equal.
    print()
    assert [col[:3] in ['obj', 'str', 'int', 'flo'] for col in df_py.columns]
    py_cols = [col for col in df_py.columns
               if str(df_py[col].dtype).lower()[:3] in ['obj', 'str', 'int']] +\
              [col for col in df_py.columns
               if str(df_py[col].dtype).lower()[:3] == 'flo']  # exact cols first, then floats.
    for col in py_cols:
        print('starting col', col)
        py_type = str(df_py[col].dtype).lower()[:3]
        py_values = list(df_py[col])
        r_values = list(df_r[col])
        if col == 'UTC_start':
            assert all([py.split('+')[0] == r.replace('T', ' ')
                        for (py, r) in zip(py_values, r_values)])  # py is ISO 8601, R ~ different.
        elif col == 'SkyADU':
            mean_shift = sum(py_values)/len(py_values) - sum(r_values)/len(r_values)
            assert mean_shift == pytest.approx(0, abs=2)  # seems to be the shift.
        elif col in ['SkySigma', 'FWHM']:
            pct_abs_diff = [abs(py - r) / ((py + r) / 2) * 100.0
                            for (py, r) in zip(py_values, r_values)
                            if py > 0 and r > 0]
            mean_pct_abs_diff = sum(pct_abs_diff) / len(pct_abs_diff)
            assert mean_pct_abs_diff < 5.0
        elif col in ['Vignette', 'X1024', 'Y1024']:
            abs_diff = [abs(py-r) for (py, r) in zip(py_values, r_values)]
            mean_abs_diff = sum(abs_diff)/len(abs_diff)
            assert mean_abs_diff < 0.004
        elif col in ['Xcentroid', 'Ycentroid']:
            abs_diff = [abs(py - r) for (py, r) in zip(py_values, r_values)]
            mean_abs_diff = sum(abs_diff) / len(abs_diff)
            assert mean_abs_diff < 0.2
        elif col == 'SkyBias':
            pass  # probably going to delete this term before too long, so don't bother testing.
        elif col == 'InstMag':
            stdev_diff_bright = (pd.Series([py - r
                                            for (py, r) in zip(py_values, r_values)
                                            if py < -6])).std()
            assert stdev_diff_bright < 0.015
        elif col == 'InstMagSigma':
            max_diff = (pd.Series([abs(py - r)
                                     for (py, r) in zip(py_values, r_values)
                                     if (py > 0) and (py < 0.03)])).max()
            assert max_diff < 0.015
        elif col == 'CatMag':
            both_nan = [np.isnan(py) and np.isnan(r) for (py, r) in zip(py_values, r_values)]
            equal = [py == r for (py, r) in zip(py_values, r_values)]
            cat_mag_match = [b or e for (b,e) in zip(both_nan, equal)]
            assert all(cat_mag_match)
        elif col == 'CatMagError':
            both_nan = [np.isnan(py) and np.isnan(r) for (py, r) in zip(py_values, r_values)]
            equal = [py == r for (py, r) in zip(py_values, r_values)]
            r_nan_py_max = []
            for (py, r, f, star) in zip(py_values, r_values, df_py['Filter'], df_py['ModelStarID']):
                this_tf = False
                if np.isnan(r) and not np.isnan(py):
                    this_max = max(df_py.loc[(df_py['Filter'] == f)
                                             & (df_py['ModelStarID'] == star)
                                             & (df_py['CatMagError'] is not None),
                                             'CatMagError'])
                    this_tf = (py == this_max)
                r_nan_py_max.append(this_tf)
            cat_mag_match = [(b or e or rp) for (b, e, rp) in zip(both_nan, equal, r_nan_py_max)]
            assert all(cat_mag_match)
        elif col == 'CI':
            py_nan = set(np.isnan(py_values))
            r_nan = set(np.isnan(r_values))
            assert py_nan == r_nan
            max_diff = pd.Series([abs(py - r) for (py, r) in zip(py_values, r_values)]).max()
            assert max_diff < 0.001
        elif py_type in ['obj', 'str', 'int']:
            if not py_values == r_values:
                print('col', col, 'differ')
            assert py_values == r_values
        else:
            # Here, both are presumed to be floats:
            py_values = [float(v) for v in py_values]
            r_values = [float(v) for v in r_values]
            if not py_values == pytest.approx(r_values):
                print('col', col, 'differ')
            assert py_values == pytest.approx(r_values)



