import os
import os.path
import numpy as np
import pandas as pd
from datetime import datetime, timedelta, timezone
from collections import defaultdict, namedtuple

from photrix.util import RaDec, datetime_utc_from_jd, time_hhmm
from photrix.fov import make_fov_dict, FovError
from photrix.user import Astronight
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
    '''
    Returns df_fov filtered to contain only specified observing styles.
    :param df_fov: input fov dataframe.
    :param obs_style_list: list of observing styles to include (or  string for one style).
       None simply returns the input df_fov.
    :return: filtered df_fov.
    '''
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
    print("\n\nan dark        >", an.ts_dark)
    print("an dark_no_moon>", an.ts_dark_no_moon)
    print("moon phase >", an.moon_phase)
    print("min alt > ", an.site.min_altitude)
    print("min_moon_degrees > ", min_moon_degrees)
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
    '''
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
    '''

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
        '''
        Updates cache's entry for one fov: this class's workhorse method.
        '''
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


def parse_an_sketch(excel_path='c:/24hrs/Scratch Plan.xlsx', site_name='BDO_Kansas', fov_dict=None,
                    output_directory=DEFAULT_PLAN_DIRECTORY):
    df = pd.read_excel(excel_path, sheetname='Sheet1').\
        dropna(axis=0, how='all').dropna(axis=1, how='all')
    nrow = len(df)
    ncol = len(df.columns)
    parsed_list = []  # nested list, one element per ACP plan.
    this_plan_id = ''
    plan_items = []

    an_date_string = str(df.iloc[0, 0]).strip()
    if 20170101 < int(an_date_string) < 20201231:  # max prob should be related to today's date.
        an = Astronight(an_date_string, site_name)
    else:
        print('>>>>> STOPPING: an_date_string '" + an_date_string + "' SEEMS UNREASONABLE.')
        return
    print('an_date_string: ' + an_date_string)  # TEST
    if fov_dict is None:
        fov_dict = make_fov_dict()

    for irow in range(1, nrow):
        for icol in range(ncol):
            cell = df.iloc[irow, icol]
            if isinstance(cell, str):
                do_this_cell = True
            else:
                do_this_cell = ~np.isnan(cell)
            if do_this_cell:
                txt_as_read = str(cell).strip()
                txt_lower = txt_as_read.lower()
                print(txt_as_read)
                if txt_lower.startswith('plan'):
                    if icol != 0:
                        print('>>>>>ERROR: Plan defined in a column other than the first one: "' +
                              txt_lower + '".')
                        return
                    else:
                        # Close previous plan if any, probably with chain to next plan:
                        parsed_list.append(plan_items)
                        plan_items = []
                        # Start next plan:
                        split_txt = txt_lower.split(';', maxsplit=1)
                        directive = split_txt[0]
                        this_plan_id = 'plan_' + an_date_string + '_' + \
                                       directive[len('plan'):].strip()
                        if len(split_txt) > 1:
                            comments = split_txt[1].strip()
                        else:
                            comments = None
                        plan_items.append(('plan', this_plan_id, comments))  # append tuple
                if txt_lower.startswith(';'):
                    comment_text = txt_lower.strip()[1:]
                    plan_items.append(('comment', comment_text))
                if txt_lower.startswith('afinterval'):
                    directive = txt_lower.split(';')[0].strip()
                    minutes = float(directive[len('afinterval'):].strip())
                    plan_items.append(('afinterval', minutes))
                if txt_lower.startswith('autofocus'):
                    plan_items.append(('autofocus',))
                if txt_lower.startswith('chill'):
                    directive = txt_lower.split(';')[0].strip()
                    degrees = float(directive[len('chill'):].strip())
                    plan_items.append(('chill', degrees))
                if txt_lower.startswith('burn'):
                    directive = txt_lower.split(';')[0].strip()
                    this_fov_name = directive[len('burn'):].strip()
                    plan_items.append(('burn', this_fov_name))
                if txt_lower.startswith('quitat'):
                    directive = txt_lower.split(';')[0].strip()
                    hhmm = directive[len('quitat'):].strip()
                    plan_items.append(('quitat', hhmm))
                if txt_lower.startswith('waituntil'):
                    directive = txt_lower.split(';')[0].strip()
                    value = directive[len('waituntil'):].strip()
                    if float(value) < 0:
                        plan_items.append(('waituntil', 'sun_degrees', value))
                    else:
                        plan_items.append(('waituntil', 'hhmm', value))
                if txt_lower.startswith('chain'):
                    directive = txt_lower.split(';')[0].strip()
                    next_plan_filename = directive[len('chain'):].strip()
                    if not next_plan_filename.endswith('.txt'):
                        next_plan_filename += '.txt'
                    plan_items.append(('chain', next_plan_filename))
                if txt_lower.startswith('flats'):
                    if this_plan_id[-2:].lower() != '_z':
                        print('>>>>> WARNING: flats directive encountered but plan_id is' +
                              this_plan_id + ', not the usual "Z".')
                    plan_items.append(('flats',))
                else:
                    # Treat as a fov_name, with warning if not in fov list.
                    fov_name = txt_as_read.strip()
                    plan_items.append(('FOV: ', fov_name))

    parsed_list.append(plan_items)  # Close out last plan.

    summary_lines = []
    # Add summary header here.
    for plan in parsed_list:
        acp_plan_lines = []
        plan_id = (plan[0])[1]
        for item in plan:
            new_summary_line, new_acp_plan_lines = make_lines_from_item(item, an, fov_dict)
            summary_lines.append(new_summary_line)
            acp_plan_lines.extend(new_acp_plan_lines)

        # Done with this plan, so write acp_plan_lines to ACP plan file.
        output_fullpath = os.path.join(output_directory, plan_id + '.txt')
        print('PRINT lines for plan ' + plan_id + ' to ', output_fullpath)
        # with open(output_fullpath) as this_file:
        #     for line in acp_plan_lines:
        #         this_file.write("{}\n".format(line))

    # Done with all parsing, so write summary_lines to Summary file.
    output_fullpath = os.path.join(output_directory, 'Summary_' + an.an_date_string + '.txt')
    print('PRINT all summary lines to ', output_fullpath)
    # with open(output_fullpath) as this_file:
    #     for line in summary_lines:
    #         this_file.write("{}\n".format(line))


def make_lines_from_item(item, an, fov_dict):
    summary_line = 'summary: ' + item[0]  # TEST DUMMY
    acp_plan_lines = list('acp_plan: ' + item[0])  # TEST DUMMY
    return summary_line, acp_plan_lines


def make_plan_header(an):
    return an.acp_header_string().split('\n')  # Grab an's header string and return it as lines.


def make_plan_fov_entry(an, fov_name, fov_dict):
    # Returns text plan lines for one fov observation.
    # For now, obs_style LPV needs estimated V mag to follow fov_name in its excel cell,
    #    e.g., "AU Aur / 12.5; then compute other mags and thus exp times from that.
    #    If this V mag is absent: (1) emit warning in summary, (2) assume mean V mag in summary,
    #    (3) write "XXX" for all exp times in its ACP plan entry.
    pass


def make_flat_lines():
    # Return list of plan text lines for panel flats.
    pass


def make_an_summary(an, plans, fov_dict):
    # This needs to translate  to make a summary.
    return ['PLACEHOLDER TEXT']
