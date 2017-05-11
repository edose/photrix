import os
import shutil

import pandas as pd
import pytest

from photrix import process

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


def test_apply_omit_txt():
    an_top_directory = TEST_TOP_DIRECTORY
    an_rel_directory = '$an_for_test'

    # Nested function
    def do_apply_omit_txt(directive_lines):

        fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit.txt')

        # Backup omit.txt:
        backup_path = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit-SAVE.txt')

        omit_txt_backed_up = False
        if os.path.exists(fullpath):
            if os.path.exists(backup_path):
                os.remove(backup_path)
            shutil.copy2(fullpath, backup_path)  # make a copy to restore later.
            omit_txt_backed_up = True
            os.remove(fullpath)
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

        # Restore omit.txt if it was backed up, else write stub:
        if omit_txt_backed_up and os.path.exists(backup_path):
            if os.path.exists(fullpath):
                os.remove(fullpath)
            shutil.move(backup_path, fullpath)  # restore saved copy.
        if not os.path.exists(fullpath):
            process.write_omit_txt_stub(an_top_directory=TEST_TOP_DIRECTORY,
                                        an_rel_directory=an_rel_directory)

        return df_filtered, warning_lines

    df_master = process.get_df_master(an_top_directory=TEST_TOP_DIRECTORY,
                                      an_rel_directory=an_rel_directory)

    # Case: #OBS directives:
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


def test_class_skymodel():
    an_top_directory = TEST_TOP_DIRECTORY
    an_rel_directory = '$an_for_test'
    modelV = process.SkyModel(an_top_directory=TEST_TOP_DIRECTORY,
                              an_rel_directory=an_rel_directory, filter='V')

    # Test elementary attributes:
    assert modelV.an_top_directory == TEST_TOP_DIRECTORY

