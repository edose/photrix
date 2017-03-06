import os
import pandas as pd
import pytest
from collections import namedtuple
import random

from photrix import planning
# from photrix import util
from photrix.fov import Fov
from photrix.user import Astronight
# from photrix.web import get_aavso_webobs_raw_table

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"

PHOTRIX_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_FOV_DIRECTORY = os.path.join(PHOTRIX_ROOT_DIRECTORY, "test", "$fovs_for_test")
TEST_DATA_DIRECTORY = os.path.join(PHOTRIX_ROOT_DIRECTORY, "test", "$data_for_test")


def test_make_df_fov():
    df_fov = planning.make_df_fov(fov_directory=TEST_FOV_DIRECTORY)
    required_columns = ["fov_name", "main_target", "fov_priority",
                                    "obs_style", "ra", "dec", 'period',
                                    'target_type', 'max_exposure', 'radec']
    assert list(df_fov.columns) == required_columns  # in order
    required_fov_names = set(['ASAS J162540_19123', 'AU Aur Modified', 'AU Aur', 'AU Aur BURN',
                              'HAT_P_3', 'NSV 14581', 'ST Tri', 'Std_SA100'])
    assert set(df_fov.fov_name) == required_fov_names  # not necessarily in order
    assert df_fov.shape == (len(required_fov_names), len(required_columns))
    assert df_fov['fov_name'].is_monotonic_increasing


def test_filter_df_fov_by_obs_styles():
    df_all = planning.make_df_fov(fov_directory=TEST_FOV_DIRECTORY)
    # Test normal case (as a list).
    df = planning.filter_df_fov_by_obs_styles(df_all, ["Standard", "LPV"])
    assert list(df.columns) == list(df_all.columns)
    assert set(df["fov_name"]) == set(["AU Aur Modified", "Std_SA100", "AU Aur"])
    assert df.shape == (3, len(df_all.columns))

    # Test normal case (as a single string, not as list):
    df = planning.filter_df_fov_by_obs_styles(df_all, "LPV")
    assert list(df.columns) == list(df_all.columns)
    assert set(df["fov_name"]) == set(["AU Aur Modified", "AU Aur"])
    assert df.shape == (2, len(df_all.columns))

    # Test empty obs_style list, returns input df_fov:
    df = planning.filter_df_fov_by_obs_styles(df_all, [])
    assert list(df.columns) == list(df_all.columns)
    assert set(df['fov_name']) == set(df_all['fov_name'])
    assert df.shape == (len(df_all), len(df_all.columns))

    # Test default (no action, returns input df_fov):
    df = planning.filter_df_fov_by_obs_styles(df_all)
    assert list(df.columns) == list(df_all.columns)
    assert set(df["fov_name"]) == set(df_all["fov_name"])  # returns original df
    assert df.shape == (len(df_all), len(df_all.columns))


def test_filter_df_fov_by_fov_priority():
    df_all = planning.make_df_fov(fov_directory=TEST_FOV_DIRECTORY)
    # Test normal case (including standard fovs).
    min_priority = 3
    df = planning.filter_df_fov_by_fov_priority(df_all, min_priority)
    assert df.shape == (7, len(df_all.columns))
    assert list(df.columns) == list(df_all.columns)
    assert set(df["fov_name"]) == \
        {"AU Aur BURN", "AU Aur Modified", "ASAS J162540_19123",
         "Std_SA100", "ST Tri", "NSV 14581", "AU Aur"}
    priority_ok = df["fov_priority"] >= min_priority
    is_standard_fov = (df["obs_style"].str.lower() == "Standard".lower())
    assert all(priority_ok | is_standard_fov)

    # Test exclude standard fovs.
    min_priority = 3
    df = planning.filter_df_fov_by_fov_priority(df_all, min_priority, include_std_fovs=False)
    assert df.shape == (6, len(df_all.columns))
    assert list(df.columns) == list(df_all.columns)
    assert set(df["fov_name"]) == \
        {"AU Aur BURN", "AU Aur Modified", "ASAS J162540_19123", "ST Tri", "NSV 14581", "AU Aur"}
    priority_ok = df["fov_priority"] >= min_priority
    is_standard_fov = (df["obs_style"].str.lower() == "Standard".lower())
    assert all(priority_ok)
    assert all(~ is_standard_fov)

    # Test default min_priority (which should return all fovs).
    df = planning.filter_df_fov_by_fov_priority(df_all)  # absent min_fov_priority
    assert df.shape == df_all.shape
    assert list(df.columns) == list(df_all.columns)
    assert set(df["fov_name"]) == set(df_all["fov_name"])  # returns original df


def test_filter_df_fov_available():
    df_all = planning.make_df_fov(fov_directory=TEST_FOV_DIRECTORY)
    # Test normal case, no moon effect.
    df = planning.complete_df_fov_an(df_all, an_string='20160919', site_name='BDO_Kansas',
                                     min_moon_degrees=0, remove_unobservables=False)
    moon_deg_expected = [143, 45, 45, 45, 118, 74, 25, 88]
    seconds_expected = [8693, 22905, 22905, 22905, 4323, 35846, 28034, 0]
    # print('\n', df)
    assert all([list(df['moon_deg'])[i] == pytest.approx(moon_deg_expected[i], abs=1)
                for i in range(len(df))])
    assert all([list(df['seconds'])[i] == pytest.approx(seconds_expected[i], abs=60)
                for i in range(len(df))])

    # Test normal case, big moon effect.
    df = planning.complete_df_fov_an(df_all, an_string='20160919', site_name='BDO_Kansas',
                                     min_moon_degrees=80, remove_unobservables=False)
    moon_deg_expected = [143, 45, 45, 45, 118, 74, 25, 88]
    seconds_expected = [8693, 0, 0, 0, 4323, 4776, 0, 0]
    assert all([list(df['moon_deg'])[i] == pytest.approx(moon_deg_expected[i], abs=1)
                for i in range(len(df))])
    assert all([list(df['seconds'])[i] == pytest.approx(seconds_expected[i], abs=60)
                for i in range(len(df))])
    # print(df)

    # Test case: 6 months before previous case (for opposite fov availabilities):
    df_all = planning.make_df_fov(fov_directory=TEST_FOV_DIRECTORY)
    # Test normal case, no moon effect.
    df = planning.complete_df_fov_an(df_all, an_string='20160319', site_name='BDO_Kansas',
                                     min_moon_degrees=65, remove_unobservables=True)
    moon_deg_expected = [97, 70, 70, 70, 62, 86, 96]
    seconds_expected = [19589, 15848, 15878, 15878, 2165, 36198, 5120]
    assert all([list(df['moon_deg'])[i] == pytest.approx(moon_deg_expected[i], abs=1)
                for i in range(len(df))])
    assert all([list(df['seconds'])[i] == pytest.approx(seconds_expected[i], abs=60)
                for i in range(len(df))])
    # print('\n', df)

    # Test case: near-new moon (moon 3% phase, no factor at all whatever targets' sky position):
    df_all = planning.make_df_fov(fov_directory=TEST_FOV_DIRECTORY)
    # Test normal case, no moon effect.
    df = planning.complete_df_fov_an(df_all, an_string='20161226', site_name='BDO_Kansas',
                                     min_moon_degrees=120, remove_unobservables=True)
    moon_deg_expected = [37, 148, 148, 148, 78, 109, 146, 118]
    seconds_expected = [4122, 40384, 40384, 40384, 19071, 45262, 29656, 27029]
    assert all([list(df['moon_deg'])[i] == pytest.approx(moon_deg_expected[i], abs=1)
                for i in range(len(df))])
    assert all([list(df['seconds'])[i] == pytest.approx(seconds_expected[i], abs=60)
                for i in range(len(df))])
    # print('\n', df)

    # Test case: moon very near at least one FOV:
    df_all = planning.make_df_fov(fov_directory=TEST_FOV_DIRECTORY)
    # Test normal case, no moon effect.
    df = planning.complete_df_fov_an(df_all, an_string='20170113', site_name='BDO_Kansas',
                                     min_moon_degrees=70, remove_unobservables=True)
    # print('\n', df)
    moon_deg_expected = [100, 65, 65, 65, 83, 90]
    seconds_expected = [8423, 2674, 2674, 2674, 44439, 24531]
    assert all([list(df['moon_deg'])[i] == pytest.approx(moon_deg_expected[i], abs=1)
                for i in range(len(df))])
    assert all([list(df['seconds'])[i] == pytest.approx(seconds_expected[i], abs=60)
                for i in range(len(df))])


def test_reorder_actions():
    Plan = namedtuple('Plan', ['plan_id', 'action_list'])
    Action = namedtuple('Action', ['action_type', 'action_data'])
    action_list_1_ordered = [Action(action_type='Plan', action_data=3),
                             Action(action_type='Chill', action_data=2),
                             Action(action_type='Waituntil', action_data='x'),
                             Action(action_type='Quitat', action_data=14.3),
                             Action(action_type='AFINTERVAL', action_data=90),
                             Action(action_type='autofocus', action_data=None),
                             Action(action_type='fov', action_data='AU Aur'),
                             Action(action_type='Stare', action_data='ST Tri'),
                             Action(action_type='Burn', action_data='Burner'),
                             Action(action_type='Comment', action_data='this is a comment'),
                             Action(action_type='Chain', action_data='AN_B.txt')]
    plan_1_ordered = Plan(plan_id='111', action_list=action_list_1_ordered)
    action_list_2_ordered = [Action(action_type='Plan', action_data=3),
                             Action(action_type='Quitat', action_data=14.3),
                             Action(action_type='AFINTERVAL', action_data=90),
                             Action(action_type='Comment', action_data='this is a comment'),
                             Action(action_type='fov', action_data='AU Aur'),
                             Action(action_type='flats', action_data=None),
                             Action(action_type='SHUTDOWN', action_data=None)]
    plan_2_ordered = Plan(plan_id='22222', action_list=action_list_2_ordered)
    plan_list_ordered = [plan_1_ordered, plan_2_ordered]

    n_trials = 20
    for i_trial in range(n_trials):
            # Can't use random.sample directly: keep all fov, stare, and burn in original order.
            indices_ordered = [0, 1, 2, 3, 4, 'observations', 10]
            indices_disordered = random.sample(indices_ordered, len(indices_ordered))
            index_obs = indices_disordered.index('observations')
            indices_disordered[index_obs:index_obs+1] = [5, 6, 7, 8, 9]
            list_1_disordered = [action_list_1_ordered[i] for i in indices_disordered]

            indices_ordered = [0, 1, 2, 'observations', 6]
            indices_disordered = random.sample(indices_ordered, len(indices_ordered))
            index_obs = indices_disordered.index('observations')
            indices_disordered[index_obs:index_obs+1] = [3, 4, 5]
            list_2_disordered = [action_list_2_ordered[i] for i in indices_disordered]

            plan_list_disordered = [Plan(plan_id=plan_list_ordered[0].plan_id,
                                         action_list=list_1_disordered),
                                    Plan(plan_id=plan_list_ordered[1].plan_id,
                                         action_list=list_2_disordered)]
            plan_list_reordered = planning.reorder_actions(plan_list_disordered)
            assert plan_list_reordered[0].plan_id == plan_list_ordered[0].plan_id
            assert plan_list_reordered[0].action_list == plan_list_reordered[0].action_list
            assert plan_list_reordered[1].plan_id == plan_list_ordered[1].plan_id
            assert plan_list_reordered[1].action_list == plan_list_reordered[1].action_list
            # print('loop ' + str(i_trial) + ' ok.')


def test_class_aavso_webobs():
    # This function tests only class AavsoWebobs, that is,
    #    it assumes correct operation of photrix.web.get_aavso_webobs_raw_table().

    # Use dataframe stored as .csv file, rather than using web now:
    data_fullpath = os.path.join(TEST_DATA_DIRECTORY, "ST_Tri_192.csv")
    df192 = pd.read_csv(data_fullpath, index_col=0)
    aw = planning.AavsoWebobs(dataframe=df192)
    assert len(aw.table) == len(df192)
    assert aw.star_id == df192['target_name'].iloc[0]


def test_class_local_obs_cache():
    # First, make a small df_fov_avail.
    df_fov = planning.make_df_fov()
    df_fov_avail = planning.complete_df_fov_an(df_fov, an_string='20170127', site_name='BDO_Kansas',
                                               min_moon_degrees=60, remove_unobservables=True)
    df_fov = df_fov.iloc[0:20]


