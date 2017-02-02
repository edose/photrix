import os
import os.path
import numpy as np
import pandas as pd
from datetime import datetime, timezone, timedelta
from collections import namedtuple
import ephem

from photrix.util import RaDec, datetime_utc_from_jd, time_hhmm, datetime_utc_from_hhmm, \
    ra_as_hours, dec_as_hex
from photrix.fov import make_fov_dict, FovError
from photrix.user import Astronight, Instrument
from photrix.web import get_aavso_webobs_raw_table

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"

FOV_DIRECTORY = "C:/Dev/Photometry/FOV/"
FOV_OBSERVING_STYLES = ["Standard", "Stare", "Monitor", "LPV", "Burn"]
MIN_AVAILABLE_SECONDS = 900
MIN_MOON_DEGREES_DEFAULT = 45
FITS_DIRECTORY = "J:/Astro/Images"
DEFAULT_PLAN_DIRECTORY = 'C:/Astro/Plans'

PHOTRIX_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCAL_OBS_CACHE_FULLPATH = os.path.join(PHOTRIX_ROOT_DIRECTORY, "local_obs_cache.csv")

AAVSO_WEBOBS_ROWS_TO_GET = 100
MIN_ROWS_ONE_STARE = 10
MAX_DAYS_ONE_STARE = 0.5

ABSOLUTE_MAX_EXPOSURE_TIME = 600  # seconds


def make_df_fov(fov_directory=FOV_DIRECTORY, fov_names_selected=None):
    """
    Returns new, basic fov data frame, by reading FOV files (all or selected) in given directory.
    :param fov_directory: the directory from which to read FOV files.
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
    df_fov.sort_values(by='fov_name', inplace=True)
    df_fov.reset_index(inplace=True, drop=True)
    return df_fov


def filter_df_fov_by_obs_styles(df_fov, obs_style_list=None):
    """
    Returns df_fov filtered to contain only specified observing styles.
    :param df_fov: input fov dataframe.
    :param obs_style_list: list of observing styles to include (or  string for one style).
       None simply returns the input df_fov.
    :return: filtered df_fov.
    """
    if obs_style_list is None:
        return df_fov
    if isinstance(obs_style_list, str):
        obs_style_list = [obs_style_list]
    if len(obs_style_list) <= 0:
        return df_fov
    obs_style_list_lower = [style.lower() for style in obs_style_list]
    return df_fov[[style.lower() in obs_style_list_lower for style in df_fov.obs_style]]


def filter_df_fov_by_fov_priority(df_fov, min_fov_priority=None, include_std_fovs=True):
    """
    Returns df_fov filtered to contain only fovs with specified minimum fov_priority.
    :param df_fov: input fov dataframe.
    :param min_fov_priority: min fov priority to permit. None simply returns the input df_fov.
    :param include_std_fovs: True to include standard FOVs (even though they have no fov_priority).
    :return: filtered df_fov.
    Optionally includes all standard FOVs (default=include standard fovs).
    """
    if min_fov_priority is None:
        return df_fov
    fov_priority_ok = df_fov["fov_priority"] >= min_fov_priority
    if include_std_fovs:
        is_standard_fov = df_fov["obs_style"].str.lower() == "standard"
        return df_fov[fov_priority_ok | is_standard_fov]
    else:
        return df_fov[fov_priority_ok]


def filter_df_fov_available(df_fov, an_string=None, site_name="BDO_Kansas",
                            min_moon_degrees=MIN_MOON_DEGREES_DEFAULT, remove_unobservables=True):
    if an_string is None or site_name == "":
        return df_fov
    an = Astronight(an_string, site_name)
    # print("\n\nan dark        >", an.ts_dark)
    # print("an dark_no_moon>", an.ts_dark_no_moon)
    # print("moon phase >", an.moon_phase)
    # print("min alt > ", an.site.min_altitude)
    # print("min_moon_degrees > ", min_moon_degrees)
    # Construct columns (specific to night and site) for available obs time, this astronight.
    df_fov['moon_deg'] = 0.0  # "
    df_fov['start'] = an.local_middark_utc  # dummy value to be overwritten later.
    df_fov['end'] = df_fov['start']  # "
    df_fov['seconds'] = 0.0  # "
    df_fov['available'] = ' - '.join(2*[4*' '])  # "

    for ind in df_fov.index:
        ts_obs = an.ts_observable(df_fov.loc[ind, 'radec'], min_alt=an.site.min_altitude,
                                  min_moon_dist=min_moon_degrees)
        df_fov.loc[ind, 'moon_deg'] = df_fov.loc[ind, 'radec'].degrees_from(an.moon_radec)
        df_fov.loc[ind, 'start'] = ts_obs.start
        df_fov.loc[ind, 'end'] = ts_obs.end
        df_fov.loc[ind, 'seconds'] = ts_obs.seconds
        if ts_obs.seconds > 0:
            df_fov.loc[ind, 'available'] = ' - '.join([time_hhmm(ts_obs.start),
                                                       time_hhmm(ts_obs.end)])
    if remove_unobservables:
        # df_fov = df_fov[(df_fov['seconds'] >= MIN_AVAILABLE_SECONDS) &
        #                 (df_fov['moon_deg'] >= min_moon_degrees)]
        df_fov = df_fov[(df_fov['seconds'] >= MIN_AVAILABLE_SECONDS)]
    return df_fov


class LocalObsCache:
    """
    Holds a cache dataframe of most recent relevant observations for ~all FOVs.
    Can hold only one dataframe row per fov (however many filters constitute a previous obs).
    Will go to AAVSO webobs to refresh a database row if fov's main target looks too old.
    Cache dataframe columns are:
        fov_name [string]
        main_target [string]
        obs_style [string]
        cache_datetime: datetime this row was updated [datetime.datetime UTC]
        obs_datetime: datetime of most recent known observation [datetime.datetime UTC]
        obs_mag: magnitude of most recent observation [float]
        obs_mag_filter: filter in which obs_mag was measured [string]
    """
    def __init__(self, fov_dict=None):
        # Read in local cache if it exists, create an empty cache df.
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
            self._write_to_csv()
        if fov_dict is None:
            self.fov_dict = make_fov_dict()
        else:
            self.fov_dict = fov_dict

    def update_fov_entry(self, fov_name, df_fov_avail, update_tolerance=0.25):
        """
        Updates cache's entry for one fov: this class's workhorse method.
        """
        this_fov = self.fov_dict.get(fov_name)  # Fov object for this fov_name.
        if this_fov is None:
            raise FovError
        main_target = this_fov.main_target
        self._curate_df_cache(fov_name, main_target)

        # Skip update of this fov if it's not in the list of fovs available this Astronight.
        fov_available = fov_name.lower() in df_fov_avail['fov_name'].str.lower()
        if not fov_available:
            return

        row_pre_exists = fov_name.lower() not in self.df_cache['fov_name'].str.lower()
        if row_pre_exists:
            now = datetime.now(timezone.utc)
            update_age = (now - self.df_cache.loc[fov_name, 'update_datetime']).total_seconds() \
                         / (24*3600)
            entry_out_of_date = update_age > update_tolerance
            obs_age = (now - self.df_cache.loc[fov_name, 'obs_datetime']).total_seconds() \
                      / (24*3600)
            gap_score_surely_zero = obs_age < this_fov.gap_score_days[0]
        else:
            entry_out_of_date = True
            gap_score_surely_zero = False

        if gap_score_surely_zero:
            return

        if entry_out_of_date:
            recent_observations = AavsoWebobs(star_id=main_target)
            obs_style = this_fov.observing_style
            obs_style_lower = obs_style.lower()
            target_type_lower = this_fov.main_target.lower()
            latest_obs_data = None  # default if no matches.
            if (obs_style_lower, target_type_lower) == ('lpv', 'mira'):
                latest_obs_data = self._latest_single_obs(fov_name, obs_style, recent_observations,
                                                          allow_filters=['V'])
            if (obs_style_lower, target_type_lower) == ('lpv', 'lpv'):
                latest_obs_data = self._latest_single_obs(fov_name, obs_style, recent_observations,
                                                          allow_filters=['V', 'R'])
            if obs_style_lower == 'monitor' and target_type_lower != 'astrometric':
                latest_obs_data = self._latest_single_obs(fov_name, obs_style, recent_observations,
                                                          allow_filters=['V', 'R'])
            if obs_style_lower == 'stare':
                latest_obs_data = self._latest_stare_obs(fov_name, obs_style, recent_observations,
                                                         allow_filters=['V', 'R'])
            if latest_obs_data is not None:
                if row_pre_exists:
                    self.df_cache.loc[fov_name] = latest_obs_data
                else:
                    self.df_cache = self.df_cache.append(latest_obs_data, ignore_index=True)
                    self.df_cache.index = self.df_cache['fov_name']

    def update_all_available_fov_entries(self, df_fov, days_tolerance=0.25):
        #  For each fov in fov_list, call update_fov_entry().
        fov_list = list(set(df_fov['fov_name']))
        for fov_name in fov_list:
            self.update_fov_entry(fov_name, df_fov, days_tolerance)

    @staticmethod
    def _latest_single_obs(fov_name, obs_style, recent_observations, allow_filters):
        allow_filters_lower = [f.lower() for f in allow_filters]
        table_filters_lower = recent_observations.table['filter'].str.lower()
        rows_to_keep = [f.lower() in allow_filters_lower for f in table_filters_lower]
        if sum(rows_to_keep) <= 0:
            return None
        latest_obs = recent_observations.table[rows_to_keep].nlargest(1, 'jd').iloc[0]
        return {'fov_name': [fov_name],
                'main_target': [latest_obs.target_name],
                'obs_style': [obs_style],
                'update_datetime': [datetime.now(timezone.utc)],
                'obs_datetime': [datetime_utc_from_jd(latest_obs.jd)],
                'obs_mag': [latest_obs.mag],
                'obs_mag_filter': [latest_obs.filter]}

    @staticmethod
    def _latest_stare_obs(fov_name, obs_style, recent_observations, allow_filters):
        latest_obs_data = None
        for this_filter in allow_filters:
            obs_already_found_this_filter = False
            this_filter_lower = this_filter.lower()
            table_filters_lower = recent_observations.table['filter'].str.lower()
            rows_to_keep = [f.lower() == this_filter_lower for f in table_filters_lower]
            table_this_filter = recent_observations.table[rows_to_keep].sort_values(by='jd',
                                                                                    ascending=False)
            num_tests = len(table_this_filter) - MIN_ROWS_ONE_STARE + 1
            if num_tests >= 1:
                for first_test_irow in range(0, num_tests):
                    if not obs_already_found_this_filter:
                        test_latest_jd = table_this_filter['jd']\
                            .iloc[first_test_irow]
                        test_earliest_jd = table_this_filter['jd']\
                            .iloc[first_test_irow + MIN_ROWS_ONE_STARE - 1]
                        if test_latest_jd - test_earliest_jd <= MAX_DAYS_ONE_STARE:
                            obs_already_found_this_filter = True
                            replace_latest_obs_data = False  # default
                            if latest_obs_data is None:
                                replace_latest_obs_data = True
                            else:
                                if test_latest_jd > latest_obs_data['jd']:
                                    replace_latest_obs_data = True
                            if replace_latest_obs_data:
                                latest_obs = table_this_filter.iloc[first_test_irow]
                                latest_obs_data = \
                                    {'fov_name': [fov_name],
                                     'main_target': [latest_obs.target_name],
                                     'obs_style': [obs_style],
                                     'update_datetime': [datetime.now(timezone.utc)],
                                     'obs_datetime': [datetime_utc_from_jd(latest_obs.jd)],
                                     'obs_mag': [latest_obs.mag],
                                     'obs_mag_filter': [latest_obs.filter]}
        return latest_obs_data

    def _write_to_csv(self):
        self.df_cache.to_csv(LOCAL_OBS_CACHE_FULLPATH, index=True)

    def _curate_df_cache(self, fov_name, main_target):
        # Curation: If main_target is corrupt, delete all cache lines for that fov.
        if main_target is None:
            rows_to_delete = self.df_cache['fov_name'].str.lower() == fov_name.lower()
            self.df_cache = self.df_cache[rows_to_delete == False]
            return

        # Curation: If fov and target names don't match, delete all such fov and target lines.
        rows_with_wrong_target = \
            self.df_cache['fov_name'].str.lower() == fov_name.lower() and \
            self.df_cache['main_target'].str.lower() != main_target.lower()
        rows_to_keep = [not row for row in rows_with_wrong_target]
        self.df_cache = self.df_cache[rows_to_keep]
        rows_with_wrong_fov = \
            self.df_cache['main_target'].str.lower() == main_target.lower() and \
            self.df_cache['fov_name'].str.lower() != self.this_fov_name.lower()
        rows_to_keep = [not row for row in rows_with_wrong_fov]
        self.df_cache = self.df_cache[rows_to_keep]


class AavsoWebobs:
    """
    Holds dataframe for one star from AAVSO's webobs database.
    Also updates local cache file (to unneeded future calls to webobs).
    For Observation Styles: LPV, Monitor, and Stare; no need for Standard or Burn.
    Usage: table = AavsoWebobs("AU Aur") [for one obs/night], or
           table = AavsoWebobs("ST Tri", stare=True) [for at least 10 obs/night in filter].
    """
    def __init__(self, star_id=None, dataframe=None):
        if dataframe is not None:
            self.table = dataframe  # typically for testing only.
            self.star_id = self.table['target_name'].iloc[0]
        else:
            self.table = get_aavso_webobs_raw_table(star_id,
                                                    AAVSO_WEBOBS_ROWS_TO_GET)  # production case
            self.star_id = star_id


def get_an_priority(fov, an, local_obs_cache):
    if fov.priority <= 0:
        return 0
    obs = local_obs_cache[fov.fov_name]
    if obs is None:
        return 2 * fov.priority  # the maximum possible
    jd_obs = obs['jd']
    days = an.local_middark_jd - jd_obs
    return fov.calc_priority_score(days)


# def get_local_aavso_reports(report_dir=None, earliest_an=None):
#     pass
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
    #     report_dir: directory in which all relevant AAVSO reports reside, as
    #       "J:/Astro/2016/Photometry".
    #     target_an: target astronight from which to count days, as "20151216".
    #     limit_days: days into the past to look up old AAVSO reports.
    #     Returns dict of (fov_name, days_since_last_local_obs).
    #     """
    #     if report_dir is not None and limit_days >= 1:
    #         fov_age_dict = {name: None for name in fov_dict.keys()}  # empty dict to start
    #         #  TODO: get report_list <- [report_text] for every eligible AAVSO report,
    #            latest to earliest.
    #
    #         for report_text in report_list:
    #             #  TODO: get jd_dict
    #                i.e., {fov_name: latest jd_obs} for each main target in AAVSO report.
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


def plan_night(excel_path='c:/24hrs/Scratch Plan.xlsx', site_name='BDO_Kansas',
               instrument_name='Borea', fov_dict=None, an_start_hhmm=None,
               output_directory=DEFAULT_PLAN_DIRECTORY):
    """
    Main user fn to take sketch Excel file and generate Summary and ACP Plan files.
    :param excel_path: full path to Excel file holding all info for one night's observations.
    :param site_name: a Site object for location of observations.
    :param instrument_name: an Instrument object for scope to be used.
    :param fov_dict:
    :param an_start_hhmm:
    :param output_directory:
    :return: Writes out Summary file with dateline, and one or more ACP plan files.
    """

    parsed_list, an = parse_excel(excel_path, site_name)

    if fov_dict is None:
        fov_dict = make_fov_dict()
    instrument = Instrument(instrument_name)
    raw_plan_list = make_raw_plan_list(parsed_list, an)

    reordered_plan_list = reorder_actions(raw_plan_list)

    plan_list_with_durations = insert_raw_durations(reordered_plan_list, fov_dict, instrument)




    # Calculate timeline and action status, add to each action in each plan:
    # Develop projected timeline (starting time for each action in entire night):
    if an_start_hhmm is None:
        an_start_dt = an.ts_dark.start
    else:
        an_start_dt = datetime_utc_from_hhmm(an_start_hhmm, an)
    running_dt = an_start_dt
    for i_plan in range(len(plan_list)):
        this_plan = plan_list[i_plan]  # a Plan namedtuple

        # Scan this plan for #WAITUNTIL or #QUITAT, which alter timelines directly:
        quitat_dt = None
        waituntil_dt = None
        for this_action in this_plan.action_list:  # an Action namedtuple
            if this_action.action_type == 'waituntil':
                if this_action.parsed_action[1] == 'hhmm':
                    waituntil_dt = datetime_utc_from_hhmm(this_action.parsed_action[2], an)
                elif this_action.parsed_action[1] == 'sun_degrees':
                    sun_degrees = this_action.parsed_action[2]
                    site_obs = ephem.Observer()
                    site_obs.lat, site_obs.lon = str(an.site.latitude), str(an.site.longitude)
                    site_obs.elevation = an.site.elevation
                    sun = ephem.Sun(site_obs)
                    site_obs.horizon = str(sun_degrees)
                    waituntil_dt = site_obs.previous_setting(sun, an.local_middark_utc) \
                        .datetime().replace(tzinfo=timezone.utc)
            elif this_action.action_type == 'quitat':
                quitat_dt = datetime_utc_from_hhmm(this_action.parsed_action[1], an)

        # Set starting time for this plan:
        if waituntil_dt is not None:
            if waituntil_dt > running_dt:
                running_dt = waituntil_dt
        # Now construct starting time and completion status for each action in this plan:
        for i_action in range(len(this_plan.action_list)):
            this_action = this_plan.action_list[i_action]  # an Action namedtuple
            action_start_dt = running_dt
            action_expected_end_dt = running_dt + timedelta(seconds=this_action.raw_duration)
            if quitat_dt is None:
                running_dt = action_expected_end_dt
                if action_expected_end_dt > an.ts_dark.end:
                    action_new_status = 'DAWN'
                else:
                    action_new_status = 'ok'
            else:
                if action_expected_end_dt > quitat_dt:
                    if running_dt >= quitat_dt:
                        action_new_status = 'SKIPPED'
                    else:
                        running_dt = quitat_dt
                        action_new_status = 'QUITAT'
                else:
                    running_dt = action_expected_end_dt
                    action_new_status = 'ok'
            new_action = this_action._replace(status=action_new_status, start_utc=action_start_dt)
            this_plan.action_list[i_action] = new_action
        new_plan = this_plan._replace(action_list=this_plan.action_list)
        plan_list[i_plan] = new_plan

    write_acp_plans(plan_list, output_directory)

    write_summary(plan_list, an, output_directory)


def make_raw_plan_list(parsed_list, an):
    # Construct raw master plan_list (with as-yet incomplete actions):
    Plan = namedtuple('Plan', ['plan_id', 'action_list'])
    Action = namedtuple('Action', ['action_type', 'parsed_action', 'an_priority',
                                   'raw_duration', 'status', 'start_utc', 'altitude_deg',
                                   'summary_lines', 'acp_plan_lines'])
    raw_plan_list = []  # will be the master container (a list of Plan namedtuples)
    for parsed_plan in parsed_list:
        print('\n***' + str(parsed_plan[0]))
        for plan_action in parsed_plan[1:]:
            print(str(plan_action))

        plan_id = (parsed_plan[0])[1]
        action_list = []  # will be list of Action namedtuples for this plan only.
        for parsed_action in parsed_plan:
            # Get acp_plan, summary, and raw duration for this action, make an Action namedtuple:
            # action_summary_lines, action_acp_plan_lines, action_raw_duration = \
            #     make_lines_from_one_action(parsed_action, an, fov_dict, instrument)
            this_action = Action(action_type=parsed_action[0],
                                 parsed_action=parsed_action,  # store for later use
                                 an_priority=0.0,  # overwrite later
                                 raw_duration=0.0,  # overwrite later
                                 status='no status',  # overwrite later
                                 start_utc=an.ts_dark.start,  # datetime, overwrite later
                                 altitude_deg=0.0,  # overwrite this later
                                 summary_lines=[],  # overwrite later
                                 acp_plan_lines=[]  # [], overwrite later
                                 )
            action_list.append(this_action)
        this_plan = Plan(plan_id=plan_id,
                         action_list=action_list)  # NB: several action fields are still empty.
        raw_plan_list.append(this_plan)
    return raw_plan_list


def reorder_actions(raw_plan_list):
    """
    Puts actions within each plan in the desired order, returns the updated plan list.
    :param raw_plan_list: the plan list to reorder.
    :return: the reordered plan list
    """
    ideal_action_ordering = [['plan'], ['chill'], ['waituntil'], ['quitat'], ['afinterval'],
                          ['stare', 'fov', 'burn', 'autofocus', 'comment'], ['flats', 'darks'],
                          ['shutdown'], ['chain']]  # sublists are kept in original order
    reordered_plan_list = []
    for plan in raw_plan_list:
        reordered_action_list = []
        for action_order_sublist in ideal_action_ordering:
            for i_action in range(len(plan.action_list)):
                if plan.action_list[i_action].action_type.lower() in action_order_sublist:
                    this_action = plan.action_list[i_action]
                    reordered_action_list.append(this_action)
        num_omitted = len(reordered_action_list) - len(plan.action_list)
        if num_omitted != 0:
            print('*** WARNING: ' + str(num_omitted) + ' actions in plan ' + plan.plan_id +
                  'were omitted during ordering.')
        new_plan = plan._replace(action_list=reordered_action_list)
        reordered_plan_list.append(new_plan)
    return reordered_plan_list


def insert_raw_durations(plan_list, fov_dict, instrument):
    for plan in plan_list:
        for action in plan.action_list:
            this_type = action.action_type.lower()
            if this_type == 'plan':
                action.raw_duration = 60
            elif this_type == 'chill':
                action.raw_duration = 0
            elif this_type == 'waituntil':
                action.raw_duration = None  # special case: no duration but executes a delay
            elif this_type == 'quitat':
                action.raw_duration = None  # special case: no duration but executes a delay
            elif this_type == 'afinterval':
                action.raw_duration = 0  # but prob inserts AUTOFOCUS lines which do have duration
            elif this_type == 'autofocus':
                action.raw_duration = 180  # includes slew & filter wheel changes
            elif this_type == 'comment':
                action.raw_duration = 0
            elif this_type == 'fov':
                # TODO: do these raw durations (prob. copy in from previous code)
                action.raw_duration = 555  # <--------- CODE THIS !!!
            elif this_type == 'stare':
                action.raw_duration = 666  # <--------- CODE THIS !!!
            elif this_type == 'flats':
                action.raw_duration = 777  # <--------- CODE THIS !!!
            elif this_type == 'darks':
                action.raw_duration = 888  # <--------- CODE THIS !!!
            elif this_type == 'chain':
                action.raw_duration = 0  # though the plan chained to will have a duration
            elif this_type == 'shutdown':
                action.raw_duration = 480  # a good guess
            else:
                print('*** WARNING: insert_raw_duration() cannot understand action \'' +
                      action.action_type + 'in plan \'' + plan.plan_id + '\'')
    return plan_list


def write_acp_plans(plan_list, output_directory):
    for this_plan in plan_list:
        # Unpack acp_plan_lines:
        plan_comment_string = this_plan.action_list[0].parsed_action[2]
        if plan_comment_string is None:
            plan_comment_string = ''
        else:
            plan_comment_string = ' ; ' + plan_comment_string
        acp_plan_lines = ['; ACP PLAN ' + this_plan.plan_id + plan_comment_string,
                          '; as generated by photrix at ' +
                          '{:%Y-%m-%d %H:%M  UTC}'.format(datetime.now(timezone.utc))]
        for this_action in this_plan.action_list:
            acp_plan_lines.extend(this_action.acp_plan_lines)
        # Write ACP plan file:
        filename = 'plan_' + this_plan.plan_id + '.txt'
        output_fullpath = os.path.join(output_directory, filename)
        print('PRINT lines for plan ' + this_plan.plan_id + ' to ', output_fullpath)
        with open(output_fullpath, 'w') as this_file:
            this_file.write('\n'.join(acp_plan_lines))


def write_summary(plan_list, an, output_directory):
    # Unpack summary_lines:
    summary_lines = ['SUMMARY for AN ' + an.an_date_string,
                     'as generated by photrix at ' +
                     '{:%Y-%M-%d %H:%M  UTC}'.format(datetime.now(timezone.utc)),
                     an.acp_header_string()]
    for this_plan in plan_list:
        for this_action in this_plan.action_list:
            prefix = this_action.status.rjust(8) + ' ' + time_hhmm(this_action.start_utc) + ' '
            n_blank_prefixes = len(this_action.summary_lines) - 1
            if this_action.action_type.lower() == 'plan':
                line_prefixes = n_blank_prefixes * [len(prefix) * ' '] + [prefix]
            else:
                line_prefixes = [prefix] + n_blank_prefixes * [len(prefix) * ' ']
            prefixed_summary_lines = [line_prefixes[i] + this_action.summary_lines[i]
                                      for i in range(len(this_action.summary_lines))]
            summary_lines.extend(prefixed_summary_lines)
    # Write Summary file:
    output_fullpath = os.path.join(output_directory, 'Summary_' + an.an_date_string + '.txt')
    print('PRINT all summary lines to ', output_fullpath)
    with open(output_fullpath, 'w') as this_file:
        this_file.write('\n'.join(summary_lines))


def parse_excel(excel_path='c:/24hrs/planning.xlsx', site_name='BDO_Kansas'):
    """
    Parses sketch Excel file and returns a list of actions constituting one night's observations.
    :param excel_path: full path to Excel file holding all info for one night's observations.
    :param site_name: a Site object for location of observations.
    :return:
    """
    df = pd.read_excel(excel_path).dropna(axis=0, how='all').dropna(axis=1, how='all')
    nrow = len(df)
    ncol = len(df.columns)
    parsed_list = []  # nested list, one element per ACP plan.
    this_plan_id = ''
    plan_actions = []
    an_date_string = str(df.iloc[0, 0]).strip()
    if 20170101 < int(an_date_string) < 20201231:  # max prob should be related to today's date.
        an = Astronight(an_date_string, site_name)
    else:
        print('>>>>> STOPPING: an_date_string '" + an_date_string + "' SEEMS UNREASONABLE.')
        return
    print('an_date_string: ' + an_date_string)  # TEST

    for irow in range(1, nrow):
        for icol in range(ncol):
            cell = df.iloc[irow, icol]
            if isinstance(cell, str):
                do_this_cell = True
            else:
                do_this_cell = ~np.isnan(cell)
            if do_this_cell:
                # Extract and process substrings from this cell:
                cell_str_as_read = str(cell).strip()
                cell_str_lower = cell_str_as_read.lower()
                split_str = cell_str_as_read.split(';', maxsplit=1)
                command = split_str[0].strip()
                if len(split_str) > 1:
                    comment = split_str[1].rstrip()
                else:
                    comment = None
                print(cell_str_as_read)

                # Determine action type and add action to plan_actions:
                if cell_str_lower.startswith('plan'):
                    # Close previous plan if any, probably with chain to next plan:
                    if len(plan_actions) > 0:
                        parsed_list.append(plan_actions)
                        plan_actions = []
                    # Start next plan:
                    this_plan_id = an_date_string + '_' + command[len('plan'):].strip()
                    plan_actions.append(('Plan', this_plan_id, comment))  # append tuple
                elif cell_str_as_read.startswith(';'):
                    plan_actions.append(('comment', comment))
                elif cell_str_lower.startswith('afinterval'):
                    minutes = command[len('afinterval'):].strip()
                    plan_actions.append(('afinterval', minutes))
                elif cell_str_lower.startswith('autofocus'):
                    plan_actions.append(('autofocus',))
                elif cell_str_lower.startswith('chill'):
                    degrees = command[len('chill'):].strip()
                    plan_actions.append(('chill', degrees))
                elif cell_str_lower.startswith('quitat'):
                    hhmm = command[len('quitat'):].strip().replace(':', '')
                    plan_actions.append(('quitat', hhmm))
                elif cell_str_lower.startswith('waituntil'):
                    value = command[len('waituntil'):].strip().replace(':', '')
                    if float(value) < 0:
                        plan_actions.append(('waituntil', 'sun_degrees', value))
                    else:
                        plan_actions.append(('waituntil', 'hhmm', value))
                elif cell_str_lower.startswith('chain'):
                    next_plan_filename = 'plan_' + an_date_string + '_' + \
                                         command[len('chain'):].strip().upper()
                    if not next_plan_filename.endswith('.txt'):
                        next_plan_filename += '.txt'
                    plan_actions.append(('chain', next_plan_filename))
                elif cell_str_lower.startswith('flats'):
                    if this_plan_id[-2:].lower() != '_z':
                        print('>>>>> WARNING: flats directive encountered but plan_id is' +
                              this_plan_id + ', not the usual "_Z".')
                    plan_actions.append(('flats',))
                elif cell_str_lower.startswith('burn'):
                    value = command[len('burn'):].strip()
                    this_fov_name, ra_string, dec_string = tuple(value.rsplit(maxsplit=2))
                    plan_actions.append(('burn', this_fov_name.strip(),
                                         ra_string.strip(), dec_string.strip()))
                elif cell_str_lower.startswith('stare'):
                    value = command[len('stare'):].strip()
                    repeats_string, this_fov_name = tuple(value.split(maxsplit=1))
                    plan_actions.append(('stare', repeats_string.strip(), this_fov_name.strip()))
                else:
                    # Treat as a fov_name, with warning if not in fov list.
                    fov_name = cell_str_as_read.strip()
                    plan_actions.append(('fov', fov_name))
                print(plan_actions[-1:])

    parsed_list.append(plan_actions)  # Close out the last plan.
    return parsed_list, an


def make_lines_from_one_action(action, an, fov_dict, instrument):
    """
    :param action: one tuple from parsed_list, representing one action.
    :param an: Astronight object.
    :param fov_dict: FOV dictionary.
    :param instrument: Instrument object.
    :return: 3-tuple of summary_lines (list of strings, lines for summary file),
                        acp_plan_lines (list of strings, >=1 lines for acp_plan now in action),
                        duration (seconds, estimated time this action will require).
    """
    summary_lines = []  # default
    acp_plan_lines = []  # default
    raw_duration = 0  # in seconds, default
    action_type = action[0]
    if action_type == 'Plan':
        # Construct first line of each.
        if action[2] is None:
            text = action[1]
        else:
            text = action[1] + ' ; ' + action[2]
        summary_lines = ['', 55*'-', 'Begin Plan ' + text]
        acp_plan_lines.extend(an.acp_header_string().split('\n'))
        acp_plan_lines.extend([';'])
        raw_duration = 0
    elif action_type == 'chill':
        summary_lines = ['CHILL  ' + action[1]]
        acp_plan_lines = ['#CHILL  ' + action[1]]
        raw_duration = 0
    elif action_type == 'waituntil':
        if action[1] == 'hhmm':
            summary_lines = ['WAITUNTIL ' + action[2] + ' utc']
            acp_plan_lines = ['#WAITUNTIL 1, ' + action[2] + ' ; utc']
        elif action[1] == 'sun_degrees':
            summary_lines = ['WAITUNTIL sun reaches ' + action[2] + ' degrees alt']
            acp_plan_lines = ['#WAITUNTIL 1, ' + action[2] + ' ; deg sun alt']
        else:
            print("***** ERROR: WAITUNTIL action" + str(action) + ' is not understood.')
    elif action_type == 'quitat':
        dt = datetime_utc_from_hhmm(action[1], an)
        formatted_time = '{:%m/%d/%Y %H:%M}'.format(dt)
        summary_lines = ['QUITAT ' + action[1] + ' utc']
        acp_plan_lines = ['#QUITAT ' + formatted_time + ' ; utc']
        raw_duration = 0
    elif action_type == 'afinterval':
        summary_lines = ['AFINTERVAL ' + action[1]]
        acp_plan_lines = [';', '#AFINTERVAL  ' + action[1]]
        raw_duration = 0
    elif action_type == 'burn':
        summary_lines = ['BURN ' + action[1]]
        acp_plan_lines = [';', '#DITHER 0 ;', '#FILTER V,I ;', '#BINNING 1,1 ;', '#COUNT 1,1 ;',
                          '#INTERVAL 240,240 ;', ';----> BURN for new FOV file.',
                          action[1] + '\t' + action[2] + '\t' + action[3] + ' ;']
        raw_duration = 660
    elif action_type == 'autofocus':
        summary_lines = ['AUTOFOCUS']
        acp_plan_lines = [';', '#AUTOFOCUS']
        raw_duration = 160
    elif action_type == 'chain':
        summary_lines = ['Chain to \'' + action[1] + '\'']
        acp_plan_lines = [';', '#CHAIN ' + action[1]]
    elif action_type == 'flats':
        summary_lines = ['Flats']
        acp_plan_lines = [';', '#SCREENFLATS flats_VRI_16.txt ;']
    elif action_type == 'comment':
        summary_lines = acp_plan_lines = [';' + action[1]]
    elif action_type == 'stare':
        summary_lines = ['Stare ' + action[1] + ' repeats at ' + action[2]]
        acp_plan_lines, raw_duration = make_plan_obs_lines(fov_name=action[2], fov_dict=fov_dict,
                                                           an=an, instrument=instrument,
                                                           num_repeats=int(action[1]))
    elif action_type == 'fov':
        summary_lines = ['Observe ' + action[1]]
        acp_plan_lines, raw_duration = \
            make_plan_obs_lines(fov_name=action[1], fov_dict=fov_dict, an=an,
                                instrument=instrument, num_repeats=1)
    else:
        print("***** ERROR: parsed_list action " + str(action) + 'not understood.')
    return summary_lines, acp_plan_lines, raw_duration


def make_plan_obs_lines(fov_name, fov_dict, an, instrument, num_repeats=1):
    """
    Returns a list of strings representing one fov observation in an ACP plan, and duration/sec.
    :param fov_name: name of FOV to observe.
    :param fov_dict: fov_dictionary including this FOV.
    :param an: Astronight object for observing night.
    :param instrument: Instrument object for telescope on which to make object.
    :param num_repeats: number of repeats (for Stare observation style only).
    :return: 2-ple(list of strings for direct inclusion in ACP plan, entry raw duration in seconds).
    """
    this_fov = fov_dict[fov_name]
    obs_style = this_fov.observing_style
    filters = []
    binnings = []
    counts = []
    exp_times = []
    mags = dict()
    startup_duration = 60  # for slew + start guiding
    repeat_duration = 0
    for obs in this_fov.observing_list:
        filter, mag, count = obs
        filters.append(filter)
        binnings.append('1')
        counts.append(str(count))
        if obs_style.lower() in ['standard', 'monitor', 'stare']:
            exp_time = calc_exp_time(mag, filter, instrument, this_fov)
        elif obs_style.lower() == 'lpv':
            if len(mags) == 0:
                mags = this_fov.estimate_lpv_mags(an.local_middark_jd)  # dict (get on 1st obs only)
            exp_time = calc_exp_time(mags[filter], filter, instrument, this_fov)
        else:
            error_string = '****** WARNING: fov \'' + fov_name + \
                           '\' has unrecognized observing style \'' + obs_style + '\'.'
            return error_string, 0
        exp_time = round(exp_time)
        exp_times.append('{:.0f}'.format(exp_time))
        repeat_duration += 15 + count * (15 + exp_time)

    print("test " + fov_name)
    total_duration = startup_duration + num_repeats * repeat_duration
    if num_repeats > 1:
        duration_comment = str(round(repeat_duration / 60.0, 1)) + ' min/repeat --> ' +\
                           str(round(total_duration / 60.0, 1)) + ' min'
    else:
        duration_comment = ' --> ' + str(round(total_duration / 60.0, 1)) + ' min'

    acp_string = [';', '#DITHER 0 ;',
                  '#FILTER ' + ','.join(filters) + ' ;',
                  '#BINNING ' + ','.join(binnings) + ' ;',
                  '#COUNT ' + ','.join(counts) + ' ;',
                  '#INTERVAL ' + ','.join(exp_times) + ' ; ' + duration_comment,
                  ';----' + this_fov.acp_comments,
                  fov_name + '\t' + ra_as_hours(this_fov.ra) + '\t' + dec_as_hex(this_fov.dec)]
    if num_repeats > 1:
        acp_string.insert(1, '#REPEAT ' + str(num_repeats) + ' ;')
    return acp_string, total_duration


def calc_exp_time(mag, filter, instrument, fov):
    exp_time_from_mag = instrument.filters[filter]['reference_exposure_mag10'] *\
                        10.0 ** ((mag - 10.0) / 2.5)
    if fov.max_exposure is None:
        exp_time_asymptote = ABSOLUTE_MAX_EXPOSURE_TIME
    else:
        exp_time_asymptote = min(ABSOLUTE_MAX_EXPOSURE_TIME, fov.max_exposure)
    exp_time = 1.0 / (1.0 / exp_time_from_mag + 1.0 / exp_time_asymptote)
    return exp_time
