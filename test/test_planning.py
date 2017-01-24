import os
import pandas as pd
import pytest

from photrix import planning
# from photrix import util
from photrix.fov import Fov
from photrix.user import Astronight
# from photrix.web import get_aavso_webobs_raw_table

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"

PHOTRIX_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_FOV_DIRECTORY = os.path.join(PHOTRIX_ROOT_DIRECTORY, "test", "$fovs_for_test")
TEST_DATA_DIRECTORY = os.path.join(PHOTRIX_ROOT_DIRECTORY, "test", "$data_for_test")
#
# def test_make_df_fov():
#     print()
#     df_fov = planning.make_df_fov(fov_directory=TEST_FOV_DIRECTORY)
#     assert df_fov.shape == (7, 6)
#     assert list(df_fov.columns) == \
#            ["fov_name", "main_target", "fov_priority", "obs_style", "ra", "dec"]
#     # print(df_fov)
#
#
# def test_df_fov_by_fov_priority():
#     df_all = planning.make_df_fov(fov_directory=TEST_FOV_DIRECTORY)
#     # Test normal case (including standard fovs).
#     min_priority = 3
#     df = planning.df_fov_by_fov_priority(df_all, min_priority)
#     assert df.shape == (6, 6)
#     assert set(df["fov_name"]) == \
#         {"AU Aur Modified", "ASAS J162540_19123", "Std_SA100", "ST Tri", "NSV 14581", "AU Aur"}
#     priority_ok = df["fov_priority"] >= min_priority
#     is_standard_fov = df["obs_style"] == "Standard"
#     assert all(priority_ok | is_standard_fov)
#
#     # Test exclude standard fovs.
#     min_priority = 3
#     df = planning.df_fov_by_fov_priority(df_all, min_priority, include_std_fovs=False)
#     assert df.shape == (5, 6)
#     assert set(df["fov_name"]) == \
#         {"AU Aur Modified", "ASAS J162540_19123", "ST Tri", "NSV 14581", "AU Aur"}
#
#     # Test default min_priority (which should return all fovs).
#     df = planning.df_fov_by_fov_priority(df_all)  # absent min_fov_priority
#     assert df.shape == df_all.shape
#     assert set(df["fov_name"]) == set(df_all["fov_name"])  # returns original df
#
#
# def test_df_fov_by_obs_style():
#     df_all = planning.make_df_fov(fov_directory=TEST_FOV_DIRECTORY)
#     # Test normal case (list).
#     df = planning.df_fov_by_obs_styles(df_all, ["Standard", "LPV"])
#     assert df.shape == (3, 6)
#     assert set(df["fov_name"]) == {"AU Aur Modified", "Std_SA100", "AU Aur"}
#
#     # Test normal case (single string, not as list).
#     df = planning.df_fov_by_obs_styles(df_all, "LPV")
#     assert df.shape == (2, 6)
#     assert set(df["fov_name"]) == {"AU Aur Modified", "AU Aur"}
#
#     # Test absent list of obs styles (thus include all fovs).
#     df = planning.df_fov_by_obs_styles(df_all)
#     assert df.shape == df_all.shape
#     assert set(df["fov_name"]) == set(df_all["fov_name"])  # returns original df
#


def test_df_fov_night():
    df_all = planning.make_df_fov(fov_directory=TEST_FOV_DIRECTORY)
    # print(df_all, "\n\n")
    df = planning.df_fov_night(df_all, an_string='20160919', site_name='BDO_Kansas',
                               min_moon_degrees=0, remove_unobservables=False)
    # print(df)


def test_aavso_webobs():
    # This function tests only class AavsoWebobs, that is,
    #    it assumes correct operation of photrix.web.get_aavso_webobs_raw_table().

    # Test NON-STARE case 1.
    # Use dataframe stored as .csv file, rather than using web now:
    data_fullpath = os.path.join(TEST_DATA_DIRECTORY, "ST_Tri_192.csv")
    df192 = pd.read_csv(data_fullpath, index_col=0)
    aw = planning.AavsoWebobs(dataframe=df192, filters='V', stare=False)
    assert len(aw.table) == len(df192)
    assert aw.star_id.lower() == 'ST Tri'.lower()
    assert aw.filters == ['V']
    aw_filters_lower = [f.lower() for f in aw.filters]
    subtable_expected = aw.table[[f.lower() in aw_filters_lower for f in aw.table['filter']]]
    assert list(aw.subtable['jd']) == list(subtable_expected['jd'])
    assert list(aw.subtable['filter']) == list(subtable_expected['filter'])
    assert aw.star_id == df192['target_name'].iloc[0]
    assert aw.stare is False
    irow_latest_jd = list(aw.subtable['jd']).index(max(aw.subtable['jd']))
    assert aw.latest_jd == subtable_expected['jd'].iloc[irow_latest_jd]
    assert aw.latest_mag == subtable_expected['mag'].iloc[irow_latest_jd]
    assert aw.latest_mag_filter == subtable_expected['filter'].iloc[irow_latest_jd]

    # NON-STARE case 2: no data for given filter:
    aw = planning.AavsoWebobs(dataframe=df192, filters='not a filter', stare=False)
    assert aw.stare is False
    assert aw.latest_jd is None
    assert aw.latest_mag is None
    assert aw.latest_mag_filter is None

    # NON-STARE case 3: two filters, both with data:
    aw = planning.AavsoWebobs(dataframe=df192, filters=['V', 'R'], stare=False)
    assert aw.stare is False
    assert aw.latest_jd == pytest.approx(2457373.766, abs=0.001)
    assert aw.latest_mag == pytest.approx(14.027, abs=0.001)
    assert aw.latest_mag_filter == 'V'
    aw = planning.AavsoWebobs(dataframe=df192, filters=['R', 'V'], stare=False)
    assert aw.stare is False
    assert aw.latest_jd == pytest.approx(2457373.766, abs=0.001)
    assert aw.latest_mag == pytest.approx(14.027, abs=0.001)
    assert aw.latest_mag_filter == 'V'

    # NON-STARE case 4: two filters, one of which has no data:
    aw = planning.AavsoWebobs(dataframe=df192, filters=['not a filter', 'R'], stare=False)
    assert aw.stare is False
    assert aw.latest_jd == pytest.approx(2455190.9, abs=0.001)
    assert aw.latest_mag == pytest.approx(13.578, abs=0.001)  # 13.578
    assert aw.latest_mag_filter == 'R'
    aw = planning.AavsoWebobs(dataframe=df192, filters=['R', 'not a filter'], stare=False)
    assert aw.stare is False
    assert aw.latest_jd == pytest.approx(2455190.9, abs=0.001)
    assert aw.latest_mag == pytest.approx(13.578, abs=0.001)  # 13.578
    assert aw.latest_mag_filter == 'R'

    # Test STARE cases, using df192 (but not aw) from above section:
    # STARE case 1: use df192 as is, V filter.
    aw = planning.AavsoWebobs(dataframe=df192, filters='V', stare=True)
    assert aw.stare is True
    assert aw.latest_jd == pytest.approx(2457373.766, abs=0.001)
    assert aw.latest_mag == pytest.approx(14.027, abs=0.001)
    assert aw.latest_mag_filter == 'V'

    # STARE case 2: truncated dataframe, giving different latest observation:
    df192_trunc = df192[10:]  # drop latest 10 observations
    aw = planning.AavsoWebobs(dataframe=df192_trunc, filters='V', stare=True)
    assert aw.stare is True
    assert aw.latest_jd == pytest.approx(2457373.738, abs=0.001)
    assert aw.latest_mag == pytest.approx(14.013, abs=0.001)
    assert aw.latest_mag_filter == 'V'

    # STARE case 3: full dataframe, different filter (which is used only in later obs):
    aw = planning.AavsoWebobs(dataframe=df192, filters='I', stare=True)
    assert aw.stare is True
    assert aw.latest_jd == pytest.approx(2457373.769, abs=0.001)
    assert aw.latest_mag == pytest.approx(13.531, abs=0.001)
    assert aw.latest_mag_filter == 'I'

    # STARE case 4: full dataframe, different filter (which is used only in earlier obs):
    aw = planning.AavsoWebobs(dataframe=df192, filters='R', stare=True)
    assert aw.stare is True
    assert aw.latest_jd == pytest.approx(2455190.9, abs=0.001)
    assert aw.latest_mag == pytest.approx(13.578, abs=0.001)  # 13.578
    assert aw.latest_mag_filter == 'R'

    # STARE case 5: three filters, two of which have valid stares:
    df192_trunc = df192[-30:]  # keep oldest 30 observations
    aw = planning.AavsoWebobs(dataframe=df192_trunc, filters=['V', 'R', 'not a filter'],
                              stare=True)
    assert aw.stare is True
    assert aw.latest_jd == pytest.approx(2455190.892, abs=0.001)
    assert aw.latest_mag == pytest.approx(13.867, abs=0.001)
    assert aw.latest_mag_filter == 'V'
    df192_trunc = df192[-29:]  # keep oldest 29 observations
    aw = planning.AavsoWebobs(dataframe=df192_trunc, filters=['V', 'R'], stare=True)
    assert aw.stare is True
    assert aw.latest_jd == pytest.approx(2455190.892, abs=0.001)
    assert aw.latest_mag == pytest.approx(13.619, abs=0.001)
    assert aw.latest_mag_filter == 'R'

    # STARE case 6: a filter with no data (gives empty set):
    aw = planning.AavsoWebobs(dataframe=df192, filters='not a filter', stare=True)
    assert aw.stare is True
    assert aw.latest_jd is None
    assert aw.latest_mag is None
    assert aw.latest_mag_filter is None

    # STARE case 7: two filters, only one has enough data for a valid stare.
    df192_trunc = df192[:30]  # keep latest 30 observations
    aw = planning.AavsoWebobs(dataframe=df192_trunc, filters=['V', 'I'], stare=True)
    assert aw.stare is True
    assert aw.latest_jd == pytest.approx(2457373.766, abs=0.001)
    assert aw.latest_mag == pytest.approx(14.027, abs=0.001)
    assert aw.latest_mag_filter == 'V'


def test_get_an_priority():
    fov = Fov('ST Tri')
    an = Astronight('20160301', 'BDO_Kansas')
    an_priority = planning.get_an_priority(fov=fov, an=an, stare=False)

