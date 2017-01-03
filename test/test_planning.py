import os
import pandas as pd
from photrix import planning
from photrix import fov
from photrix import util

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"

PHOTRIX_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_FOV_DIRECTORY = os.path.join(PHOTRIX_ROOT_DIRECTORY, "test", "$fovs_for_test")


def test_make_df_fov():
    print()
    df_fov = planning.make_df_fov(fov_directory=TEST_FOV_DIRECTORY)
    assert df_fov.shape == (7, 6)
    assert list(df_fov.columns) == \
           ["fov_name", "main_target", "fov_priority", "obs_style", "ra", "dec"]
    # print(df_fov)

def test_df_fov_by_fov_priority():
    df_all = planning.make_df_fov(fov_directory=TEST_FOV_DIRECTORY)
    # Test normal case (including standard fovs).
    min_priority = 3
    df = planning.df_fov_by_fov_priority(df_all, min_priority)
    assert df.shape == (6, 6)
    assert set(df["fov_name"]) == \
        {"AU Aur Modified", "ASAS J162540_19123", "Std_SA100", "ST Tri", "NSV 14581", "AU Aur"}
    priority_ok = df["fov_priority"] >= min_priority
    is_standard_fov = df["obs_style"] == "Standard"
    assert all(priority_ok | is_standard_fov)

    # Test exclude standard fovs.
    min_priority = 3
    df = planning.df_fov_by_fov_priority(df_all, min_priority, include_std_fovs=False)
    assert df.shape == (5, 6)
    assert set(df["fov_name"]) == \
        {"AU Aur Modified", "ASAS J162540_19123", "ST Tri", "NSV 14581", "AU Aur"}

    # Test absent priority (thus include all fovs).
    df = planning.df_fov_by_fov_priority(df_all)  # absent min_fov_priority
    assert df.shape == df_all.shape
    assert set(df["fov_name"]) == set(df_all["fov_name"])  # returns original df


def test_df_fov_by_obs_style():
    df_all = planning.make_df_fov(fov_directory=TEST_FOV_DIRECTORY)
    # Test normal case (list).
    df = planning.df_fov_by_obs_styles(df_all, ["Standard", "LPV"])
    assert df.shape == (3, 6)
    assert set(df["fov_name"]) == {"AU Aur Modified", "Std_SA100", "AU Aur"}

    # Test normal case (single string, not as list).
    df = planning.df_fov_by_obs_styles(df_all, "LPV")
    assert df.shape == (2, 6)
    assert set(df["fov_name"]) == {"AU Aur Modified", "AU Aur"}

    # Test absent list of obs styles (thus include all fovs).
    df = planning.df_fov_by_obs_styles(df_all)
    assert df.shape == df_all.shape
    assert set(df["fov_name"]) == set(df_all["fov_name"])  # returns original df


def test_df_fov_by_observable():
    # TODO: probably need to write several fns in util module.
    assert 1==1  # placeholder to prevent test failure for now.
