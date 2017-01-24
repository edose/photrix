import os
import os.path
import numpy as np
import pandas as pd
from datetime import datetime, timedelta, timezone
from collections import defaultdict, namedtuple

from photrix.util import RaDec
from photrix.fov import make_fov_dict
from photrix.user import Astronight
from photrix.web import get_aavso_webobs_raw_table

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"

FOV_DIRECTORY = "C:/Dev/Photometry/FOV/"
FOV_OBSERVING_STYLES = ["Standard", "Stare", "Monitor", "LPV"]
MIN_AVAILABLE_SECONDS = 900
MIN_MOON_DEGREES_DEFAULT = 45
FITS_DIRECTORY = "J:/Astro/Images"

PHOTRIX_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCAL_OBS_CACHE_FULLPATH = os.path.join(PHOTRIX_ROOT_DIRECTORY, "local_obs_cache.csv")

AAVSO_WEBOBS_ROWS_TO_GET = 100
STARE_ROWS_MINIMUM = 10


def make_df_fov(fov_directory=FOV_DIRECTORY, fov_names_selected=None):
    """
    Returns basic fov data frame, by reading FOV files (all or selected) within given directory.
    :param fov_directory:
    :param fov_names_selected: default = all FOV files within given directory.
    :return: basic data frame with columns: fov_name, main_target, fov_priority, obs_style,
         ra, dec.
    """
    fov_dict = make_fov_dict(fov_directory, fov_names_selected)
    fov_names = list(fov_dict.keys())
    df_fov = pd.DataFrame({'fov_name': fov_names})  # 1 column ('fov_name') only.
    # Add columns.
    df_fov['main_target'] = [fov_dict[name].main_target for name in fov_names]
    df_fov['fov_priority'] = [fov_dict[name].priority for name in fov_names]
    df_fov['obs_style'] = [fov_dict[name].observing_style for name in fov_names]
    df_fov['ra'] = [fov_dict[name].ra for name in fov_names]
    df_fov['dec'] = [fov_dict[name].dec for name in fov_names]
    df_fov['period'] = [fov_dict[name].period for name in fov_names]
    df_fov['target_type'] = [fov_dict[name].target_type for name in fov_names]
    df_fov['max_exposure'] = [fov_dict[name].max_exposure for name in fov_names]
    # Construct column 'radec' from 'ra' and 'dec'.
    df_fov['radec'] = RaDec(0, 0)  # dummy value to be replaced (needed to set column object type).
    for ind in df_fov.index:
        ra = df_fov.loc[ind, 'ra']
        dec = df_fov.loc[ind, 'dec']
        df_fov.loc[ind, 'radec'] = RaDec(ra, dec)
    return df_fov


def df_fov_by_obs_styles(df_fov, obs_style_list=None):
    if obs_style_list is None:
        return df_fov
    if isinstance(obs_style_list, str):
        obs_style_list = [obs_style_list]
    if len(obs_style_list) <= 0:
        return df_fov
    obs_style_list_lower = [style.lower() for style in obs_style_list]
    return df_fov[[style.lower() in obs_style_list_lower for style in df_fov.obs_style]]


def df_fov_by_fov_priority(df_fov, min_fov_priority=None, include_std_fovs=True):
    """
    Returns fov filtered by minimum fov_priority.
    Optionally includes all standard FOVs"""
    if min_fov_priority is None:
        return df_fov
    fov_priority_ok = df_fov["fov_priority"] >= min_fov_priority
    if include_std_fovs:
        is_standard_fov = df_fov["obs_style"].str.lower() == "standard"
        return df_fov[fov_priority_ok | is_standard_fov]
    else:
        return df_fov[fov_priority_ok]


def df_fov_night(df_fov, an_string=None, site_name="BDO_Kansas",
                 min_moon_degrees=MIN_MOON_DEGREES_DEFAULT, remove_unobservables=True):
    if an_string is None or site_name == "":
        return df_fov
    an = Astronight(an_string, site_name)
    # print("an dark        >", an.ts_dark)
    # print("an dark_no_moon>", an.ts_dark_no_moon)
    # print("moon phase >", an.moon_phase)
    # Construct columns (specific to night and site) for available obs time, this astronight.
    df_fov['start'] = an.local_middark_utc  # dummy value to be overwritten later.
    df_fov['end'] = df_fov['start']  # "
    df_fov['seconds'] = 0.0  # "
    df_fov['moon_deg'] = 0.0  # "
    for ind in df_fov.index:
        # print('ind', ind)
        ts_obs = an.ts_observable(df_fov.loc[ind, 'radec'], an.site.min_altitude, min_moon_degrees)
        df_fov.loc[ind, 'start'] = ts_obs.start
        df_fov.loc[ind, 'end'] = ts_obs.end
        df_fov.loc[ind, 'seconds'] = ts_obs.seconds
        df_fov.loc[ind, 'moon_deg'] = df_fov.loc[ind, 'radec'].degrees_from(an.moon_radec)
    if remove_unobservables:
        # enough_dark_time = df_fov['seconds'] >= MIN_AVAILABLE_SECONDS
        # far_enough_from_moon = df_fov['moon_deg'] >= min_moon_degrees
        # df_fov = df_fov[enough_dark_time & far_enough_from_moon]
        df_fov = df_fov[(df_fov['seconds'] >= MIN_AVAILABLE_SECONDS) &
                        (df_fov['moon_deg'] >= min_moon_degrees)]
    return df_fov


def available_non_stare(fov_directory=FOV_DIRECTORY, an_string=None, site_name="BDO_Kansas",
                        min_fov_priority=None, min_moon_degrees=MIN_MOON_DEGREES_DEFAULT):
    df = make_df_fov(fov_directory)
    # print(len(df))
    df = df_fov_by_obs_styles(df, ['Standard', 'LPV', 'Monitor'])
    # print(len(df))
    df = df_fov_by_fov_priority(df, 4)
    # print(len(df))
    df = df_fov_night(df, an_string='20170111')
    # print(len(df))
    return df


class AavsoWebobs:
    """
    Holds dataframe and derived data for one star from AAVSO's webobs database.
    Also updates local cache file (to unneeded future calls to webobs).
    Usage: table = AavsoWebobs("AU Aur") [for one obs/night], or
           table = AavsoWebobs("ST Tri", stare=True) [for at least 10 obs/night in filter].
    """
    def __init__(self, star_id=None, dataframe=None, filters='V', stare=False):
        if dataframe is not None:
            self.table = dataframe  # typically for testing only.
            self.star_id = self.table['target_name'].iloc[0]
        else:
            self.table = get_aavso_webobs_raw_table(star_id, AAVSO_WEBOBS_ROWS_TO_GET)  # production case
            self.star_id = star_id

        if isinstance(filters, str):
            filter_list = list()
            filter_list.append(filters)  # python kludge because, e.g., list('xy') gives ['x', 'y']
            filters = filter_list
        valid_filters_lower = [f.lower() for f in filters]
        table_filters_lower = self.table['filter'].str.lower()
        keep_row = [f in valid_filters_lower for f in table_filters_lower]
        self.subtable = self.table[keep_row]
        self.filters = filters  # now a list of strings.
        self.stare = stare      # convenience, a boolean.
        self.latest_jd = None  # default, to be overwritten below.
        self.latest_mag = None  # default, to be overwritten below.
        self.latest_mag_filter = None  # filter in which .latest_mag was taken.

        # Construct derived quantities.
        if not stare:  # if we only need latest single observation
            if len(self.subtable) >= 1:
                self.latest_jd = self.subtable['jd'].iloc[0]
                self.latest_mag = self.subtable['mag'].iloc[0]
                self.latest_mag_filter = self.subtable['filter'].iloc[0]
        else:  # here we need STARE_ROWS_MINIMUM rows, in a single filter, within 1/2 day.
            for this_filter in self.filters:
                table_this_filter = self.subtable[self.subtable['filter'].str.lower() ==
                                                  this_filter.lower()]
                num_tests = len(table_this_filter) - STARE_ROWS_MINIMUM + 1
                if num_tests >= 1:
                    for first_test_irow in range(0, num_tests):
                        test_latest_jd = table_this_filter['jd'].iloc[first_test_irow]
                        if self.latest_jd is None:
                            test_needed = True
                        else:
                            test_needed = test_latest_jd > self.latest_jd
                        if test_needed:
                            earliest_jd = \
                                table_this_filter['jd'].iloc[first_test_irow+STARE_ROWS_MINIMUM-1]
                            if test_latest_jd - earliest_jd <= 0.5:  # all obs w/in 1/2 day range
                                self.latest_jd = table_this_filter['jd'].iloc[first_test_irow]
                                self.latest_mag = table_this_filter['mag'].iloc[first_test_irow]
                                self.latest_mag_filter = \
                                    table_this_filter['filter'].iloc[first_test_irow]


class LocalObsCache:
    def __init__(self, fov_dict=None):
        if os.path.isfile(LOCAL_OBS_CACHE_FULLPATH):
            self.df_cache = pd.read_csv(LOCAL_OBS_CACHE_FULLPATH, index_col=0)
        else:
            #  Create *empty* Dataframe with correct dtypes (incl. timezone-aware datetimes):
            self.df_cache = pd.DataFrame({'fov_name': ['dummy'],
                                          'main_target': ['dummy'],
                                          'obs_style': ['dummy'],
                                          'cache_datetime': [datetime.now(timezone.utc)],
                                          'obs_datetime': [datetime.now(timezone.utc)],
                                          'obs_mag': [0.0],
                                          'obs_mag_filter': ['dummy']},
                                         columns=['fov_name', 'main_target', 'obs_style',
                                                  'cache_datetime', 'obs_datetime',
                                                  'obs_mag', 'obs_mag_filter'])[:0]
            self.df_cache.write_to_csv()
        if fov_dict is None:
            self.fov_dict = make_fov_dict()

    def update_fov_entry(self, fov_name, df_fov, days_tolerance=0.25):
        # If main_target is corrupt, delete all cache lines for that fov.
        this_fov = self.fov_dict.get(fov_name)
        if this_fov is None:
            main_target = None
        else:
            main_target = this_fov.main_target
        if main_target is None:
            rows_to_delete = self.df_cache['fov_name'].str.lower() == this_fov.lower()
            self.df_cache = self.df_cache[rows_to_delete == False]
            return

        # If fov and target names don't match in cache, delete all such fov and target lines.
        rows_with_wrong_target = self.df_cache['fov_name'].str.lower() == this_fov.lower() and\
                                 self.df_cache['main_target'].str.lower() != main_target.lower()
        rows_to_keep = [not row for row in rows_with_wrong_target]
        self.df_cache = self.df_cache[rows_to_keep]
        rows_with_wrong_fov = self.df_cache['main_target'].str.lower() == main_target.lower() and\
                         self.df_cache['fov_name'].str.lower() != this_fov.lower()
        rows_to_keep = [not row for row in rows_with_wrong_fov]
        self.df_cache = self.df_cache[rows_to_keep]

        # Skip fov if not available this Astronight.
        fov_available = fov_name.lower() in df_fov['fov_name'].str.lower()
        if fov_available:
            #  If no df row, add row; if new obs newer than cache line, update row.
            if fov_name.lower() not in self.df_cache['main_target'].str.lower():
                #  Add row to df_cache here.
                pass
            else:
                datetime_now = datetime.now(timezone.utc)
                # datetime_entry =
                #  TODO: we have a problem here...we're not accounting for multiple primary filters.

    def update_all_available_fov_entries(self, df_fov, days_tolerance=0.25):
        #  For each fov in fov_list, call update_fov_entry().
        fov_list = list(set(df_fov['fov_name']))
        for fov_name in fov_list:
            self.update_for_entry(fov_name, df_fov, days_tolerance)

    def write_to_csv(self):
        self.df_cache.to_csv(LOCAL_OBS_CACHE_FULLPATH, index=True)


def get_an_priority(fov, an, stare=False):
    if fov.priority <= 0:
        return 0
    if stare is False:
        filters = 'V'
    else:
        filters = ['V', 'R']
    # TODO: next line really needs to be a more general-purpose function (from multiple sources).
    webobs = AavsoWebobs(fov.main_target, filters=filters, stare=stare)

    jd_latest_obs = webobs.latest_jd
    jd_an = an.local_middark_jd
    gap_days = max(0, jd_an - jd_latest_obs)
    an_priority = fov.calc_priority_score(gap_days)
    return an_priority


def get_local_aavso_reports(report_dir=None, earliest_an=None):
    pass
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



