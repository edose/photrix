import os
import pandas as pd
from collections import defaultdict, namedtuple

# from photrix.util import *
from photrix.fov import make_fov_dict
from photrix.user import Astronight

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"

FOV_DIRECTORY = "C:/Dev/Photometry/FOV/"
FOV_OBSERVING_STYLES = ["Standard", "Stare", "Monitor", "LPV"]
MIN_AVAILABLE_SECONDS = 900
FITS_DIRECTORY = "J:/Astro/Images"


def make_df_fov(fov_directory=FOV_DIRECTORY, fov_names_selected=None):
    """
    Takes fov_directory, name list --> returns basic fov data frame df_fov.
    """
    fov_dict = make_fov_dict(fov_directory, fov_names_selected)
    fov_names = list(fov_dict.keys())
    df_fov = pd.DataFrame({'fov_name': fov_names})  # 1 column ('fov_name') only.
    df_fov['main_target'] = [fov_dict[name].main_target for name in fov_names]  # add col
    df_fov['fov_priority'] = [fov_dict[name].priority for name in fov_names]  # add col
    df_fov['obs_style'] = [fov_dict[name].observing_style for name in fov_names]  # add col
    df_fov['ra'] = [fov_dict[name].ra for name in fov_names]  # add col
    df_fov['dec'] = [fov_dict[name].dec for name in fov_names]  # add col
    return df_fov


def df_fov_by_fov_priority(df_fov, min_fov_priority=None, include_std_fovs=True):
    """
    Returns fov filtered by minimum fov_priority.
    Optionally includes all standard FOVs"""
    if min_fov_priority is None:
        return df_fov
    fov_priority_ok = df_fov["fov_priority"] >= min_fov_priority
    if include_std_fovs:
        is_standard_fov = df_fov["obs_style"] == "Standard"
        return df_fov[fov_priority_ok | is_standard_fov]
    else:
        return df_fov[fov_priority_ok]


def df_fov_by_obs_styles(df_fov, obs_style_list=None):
    if obs_style_list is None:
        return df_fov
    if isinstance(obs_style_list, str):
        obs_style_list = [obs_style_list]
    if len(obs_style_list) <= 0:
        return df_fov
    obs_style_list_lower = [style.lower() for style in obs_style_list]
    return df_fov[[style.lower() in obs_style_list_lower for style in df_fov.obs_style]]


def df_fov_by_observable(df_fov, an_string=None, site_name="BDO_Kansas", min_moon_degrees=45):
    if an_string is None or site_name == "":
        return df_fov
    an = Astronight(an_string, site_name)
    site = an.site
    df_fov.assign(avail_sec=lambda df: an.ts_observable(site.min_altitude, df.ra, df.dec).seconds)
    if min_moon_degrees > 0:
        pass  # TODO: calc moon distance at mid availability, only when moon is up (tricky).
    return df_fov[df_fov.avail_sec > MIN_AVAILABLE_SECONDS  &
                  df_fov.moon_degrees > min_moon_degrees]

# def make_monitor_df(fov_directory=FOV_DIRECTORY, fov_names_selected=None,
#                     an_string="20160101", site_name="BDO_Kansas"):
#     # TODO: Rewrite as calls to several granular functions.
#     fov_dict = make_fov_dict(fov_directory=fov_directory, fov_names_selected=fov_names_selected)
#
#     fov_dict = fov_dict_by_obs_style(fov_dict, ['monitor', 'lpv'])  # only those of this obs style
#     fov_dict = fov_dict_min_priority(fov_dict, 1)  # only those with at least minimum priority
#     fov_dict = fov_dict_available(fov_dict, an_string, site_name)  # only those available this night
#
#
#
#
#
#
#     # Eliminate fovs that are too recently observed by me.
#
#
#
#     max_gap_days = max([fov.gap_score_days[2] for fov in fov_dict.keys()])
#     local_obs_age_dict = get_local_obs_age_dict(fov_dict=fov_dict, max_days=max_gap_days)
#     fov_dict = {name: fov for (name, fov) in fov_dict.items()
#                 if fov.calc_priority_score(local_obs_age_dict[name]) > 0}
#
#     # Eliminate fovs that are too recently observed by AAVSO.
#     aavso_obs_age_dict = get_aavso_obs_age_dict(fov_dict=fov_dict, max_days=max_days)
#     fov_dict = {name: fov for (name, fov) in fov_dict.items()
#                 if fov.calc_priority_score(aavso_obs_age_dict[name]) > 0}
#
#     # Construct and return the data frame.
#     out_dict = dict()
#     for fov_name in fov_dict.keys():
#         avail = available_dict[fov_name]  # is a Timespan objectb
#         local_age = local_obs_age_dict[fov_name]
#         aavso_age = aavso_obs_age_dict[fov_name]
#         if local_age is None:
#             if aavso_age is None:
#                 age = None
#             else:
#                 age = aavso_age
#         else:
#             if aavso_age is None:
#                 age = local_age
#             else:
#                 age = max(local_age, aavso_age)
#         out_dict[fov_name] = (fov_name, avail.start, avail.end, avail.seconds, age)
#
#     fov_df = pd.DataFrame()
#
#     return fov_df
#
#
# def get_aavso_reports(report_dir=None, earliest_an=None):
#     report_dict = {}
#     for root, dirs, files in os.walk('J:/Astro/Images/C14/'):
#         if root.endswith("Photometry"):
#             report = [file for file in files if file.startswith("AAVSO")]
#             if len(report) >= 1:
#                 report_fullpath = os.path.join(root, report[0])
#                 with open(report_fullpath) as report_file:
#                     lines = report_file.readlines()
#
#
#     Report = namedtuple('Report', ['JD', 'lines'])
#
#
#
#
# def get_local_obs_age_dict(fov_dict=None, report_dir=None, target_an=None, limit_days=366):
#     # TODO: finish writing get_local_obs_age_dict()
#     """
#     report_dir: directory in which all relevant AAVSO reports reside, as "J:/Astro/2016/Photometry".
#     target_an: target astronight from which to count days, as "20151216".
#     limit_days: days into the past to look up old AAVSO reports.
#     Returns dict of (fov_name, days_since_last_local_obs).
#     """
#     if report_dir is not None and limit_days >= 1:
#         fov_age_dict = {name: None for name in fov_dict.keys()}  # empty dict to start
#         #  TODO: get report_list <- [report_text] for every eligible AAVSO report, latest to earliest.
#
#         for report_text in report_list:
#             #  TODO: get jd_dict <- {fov_name: latest jd_obs} for each main target in this AAVSO report.
#
#
#
#
#
#         for an_dir in dir_list:
#             an_dict = defaultdict(lambda: None)
#             #  read AAVSO report, fill an_dict with target: latest JD
#             for fov_name, fov in fov_dict.items():
#                 an_age = an_dict[fov.main_target]
#                 if an_age is not None:
#                     dict_age = fov_age_dict[fov_name]
#                     if dict_age is not None:
#                         if an_age < dict_age:
#                             fov_age_dict[fov_name] = an_age
#                     else:
#                         fov_age_dict[fov_name] = an_age
#     return fov_age_dict
#
#
# def get_aavso_obs_age_dict(fov_dict=None, max_days=None):
#     # TODO: finish writing get_aavso_obs_age_dict()
#     fov_age_dict = defaultdict(lambda: None)
#     for fov_name, fov in fov_dict.items():
#         # x is some function to get JD or age of most recent aavso observation
#         # may be more sophisticated, to weight visible, V, and other observations differently.
#         # watch out for age vs JD...must be some reference JD for any age calculations (AN middark?)
#         fov_age_dict[fov_name] = x(fov.main_target)
#     return fov_age_dict
