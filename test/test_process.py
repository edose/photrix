import os
import shutil

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
    def do_apply_omit_txt(directive_lines):
        fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit.txt')

        omit_txt_backed_up = _backup_omit_txt(an_top_directory, an_rel_directory)

        # Write new omit.txt with test directive lines:
        process.write_omit_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
                                    an_rel_directory=an_rel_directory)
        with open(fullpath) as f:
            lines = f.readlines()
            directive_lines = [line + '\n' for line in directive_lines]
        lines.extend(directive_lines)
        with open(fullpath, 'w') as f:
            f.writelines(lines)

        # Make output data:
        df_filtered, warning_lines = process.apply_omit_txt(an_top_directory=TEST_TOP_DIRECTORY,
                                                            an_rel_directory=an_rel_directory)

        if omit_txt_backed_up:
            _restore_omit_txt(an_top_directory, an_rel_directory)
        if not os.path.exists(fullpath):
            process.write_omit_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
                                        an_rel_directory=an_rel_directory)
        return df_filtered, warning_lines

    # Case: #OBS directives:
    df_master = process.get_df_master(an_top_directory=TEST_TOP_DIRECTORY,
                                      an_rel_directory=an_rel_directory)  # fresh copy
    df_filtered, warning_lines = \
        do_apply_omit_txt(['#OBS   Std_SA107-0001-V, 123  ;  one existing star obs',
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
        do_apply_omit_txt(['#STAR  Std_SA107, 123, V  ;  one existing star, V-filter only',
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
        do_apply_omit_txt(['#IMAGE  Std_SA107-0002-R  ;  one existing image',
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
        do_apply_omit_txt(['#JD 0.668, 0.6695   ;  two images',
                           '#JD 0.2, 0.3        ;  no images (warning)'])
    assert set(df_filtered.columns) == set(df_master.columns)
    assert set(df_master.Serial) - set(df_filtered.Serial) == set(list(range(489, 512)) + [])
    assert len(warning_lines) == 1
    assert warning_lines[0].find('#JD 0.2, 0.3') != -1  # -1 is returned if a is not in b

    # Case: #SERIAL directives:
    df_master = process.get_df_master(an_top_directory=TEST_TOP_DIRECTORY,
                                      an_rel_directory=an_rel_directory)  # fresh copy
    df_filtered, warning_lines = \
        do_apply_omit_txt(['#SERIAL 99999                            ;  no images (warning)',
                           '#SERIAL  12, 14,15,4444 677,, 5 777 765,2000, 14  ; all but 4444 OK '])
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
                       '#SERIAL  703 875 193  ;  I outliers'] # in actual R processing, for comp.
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit.txt')

    omit_txt_backed_up = _backup_omit_txt(an_top_directory, an_rel_directory)

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

    # Restore omit.txt if it was backed up, else write stub:
    if omit_txt_backed_up:
        _restore_omit_txt(an_top_directory, an_rel_directory)
    if not os.path.exists(fullpath):
        process.write_omit_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
                                    an_rel_directory=an_rel_directory)

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

    # Test SkyModel.predict():
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
                            index=['9997a', '9998a', '9999a'])  # predict() has no access to CatMag
    expected_star_mags = [12.34, 9.2, 13.333]  # ideal CatMags
    mag_predictions = modelV.predict(df_input)
    assert list(mag_predictions) == pytest.approx(expected_star_mags)


def _backup_omit_txt(an_top_directory, an_rel_directory):
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit.txt')
    backup_path = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit-SAVE.txt')
    omit_txt_backed_up = False
    if os.path.exists(fullpath):
        if os.path.exists(backup_path):
            os.remove(backup_path)
        shutil.copy2(fullpath, backup_path)  # make a copy to restore later.
        omit_txt_backed_up = True
        os.remove(fullpath)
    return omit_txt_backed_up


def _restore_omit_txt(an_top_directory, an_rel_directory):
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit.txt')
    backup_path = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit-SAVE.txt')

    if os.path.exists(backup_path):
        if os.path.exists(fullpath):
            os.remove(fullpath)
        shutil.move(backup_path, fullpath)  # restore saved copy.