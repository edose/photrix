# Python system imports:
import os
import os.path
from collections import OrderedDict
from datetime import datetime, timezone, timedelta
from math import floor, sqrt, ceil

# External library imports:
import ephem
import numpy as np
import pandas as pd

# Internal (photrix) imports:
from .fov import make_fov_dict, FovError, Fov
from .user import Astronight, Instrument, MOON_PHASE_NO_FACTOR
from .util import RaDec, datetime_utc_from_jd, hhmm_from_datetime_utc, \
    ra_as_hours, dec_as_hex, az_alt_at_datetime_utc, \
    degrees_as_hex, jd_from_datetime_utc, Timespan, event_utcs_in_timespan
from .web import get_aavso_webobs_raw_table

__author__ = "Eric Dose :: New Mexico Mira Project, Albuquerque"

# USAGE: *******************************************************************
# pl.make_an_roster('20170525', 'c:/Astro/ACP/AN20170525', user_update_tolerance_days=0.1,
#      exp_time_factor=0.75)
# pl.make_an_plan('c:/Astro/ACP/AN20170525/planning.xlsx', exp_time_factor=0.75)

# ROSTER Target Statement types:
# AZ Her  ; standard FOV target
# STARE 6 ST Tri  ;  standard stare FOV target (6 reps)
# BURN AA Aur 11:00:00 +34:00:00  ;  Burn target (240 sec in V and I)
# IMAGE target_name V=12 B=12.5(2) 12:00:00 +23:34:45  ;  arbitrary image, exp time from magnitude
# IMAGE target_name Clear=240sec(5) 12:00:00 +23:34:45  ;  arbitrary image, exp time requested directly


FOV_DIRECTORY = "C:/Dev/Photometry/FOV/"
STARE_EVENT_TYPES = {"eclipser": "minima", "exoplanet": "minima",
                     "delta scuti": "maxima", 'rr lyrae': 'maxima'}
MIN_AVAILABLE_SECONDS_DEFAULT = 900
MIN_AVAILABLE_SECONDS_STARE = 5400
MIN_MOON_DEGREES_DEFAULT = 45
MIN_MOON_DEGREES_STARE = 60
STARE_AN_PRIORITY_DIVIDER = 7.5  # >= this goes into the normal Roster list; < goes to low-pri list.
FITS_DIRECTORY = "C:/Astro/Images"
# DEFAULT_PLAN_DIRECTORY = 'C:/Astro/Plans'
DT_FMT = '%Y-%m-%d %H:%M:%S.%f%z'  # kludge around py inconsistency in python's datetime formats

PHOTRIX_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCAL_OBS_CACHE_FULLPATH = os.path.join(PHOTRIX_ROOT_DIRECTORY, "local_obs_cache.csv")

EARLIEST_AN_DATE = '20170101'
LATEST_AN_DATE = '20221231'  # Update this later, I suppose.
AN_START_REL_UTC_0000 = 19  # timedelta(UTC of actual AN start - nominal AN @ 0000 hours UTC)
#    (19 is good for North America)

# ********** Roster & cache parameters:
AAVSO_WEBOBS_ROWS_TO_GET = 100
MIN_ROWS_ONE_STARE = 10
MAX_DAYS_ONE_STARE = 0.5
DEFAULT_UPDATE_TOLERANCE_DAYS = 0.0416667  # 1 hour
FORCE_AUTOGUIDE_TOKEN = 'AG+'

# ********** ACP Timing:
CHILL_DURATION = 360  # seconds
PLAN_START_DURATION = 30  # seconds
AUTOFOCUS_DURATION = 180  # seconds, includes slew & filter wheel changes
CHAIN_DURATION = 3  # seconds; a guess
QUITAT_DURATION = 3  # seconds
SHUTDOWN_DURATION = 480  # seconds; a guess

# ********** Mount (L-500) Timing:
NEW_TARGET_DURATION = 34.3  # seconds; slew + settle + ACP processing (no guider start etc)

# ********** Camera & filter wheel (STXL-6303E) Timing:
MAX_AGGREGATE_EXPOSURE_NO_GUIDING = 119  # seconds;
GUIDE_STAR_ACQUISITION = 17.2  # seconds (if needed) (was 14.2)
GUIDER_CHECK_DURATION = 7  # seconds (if needed)  (was 4)
NEW_FILTER_DURATION = 5  # seconds; filter change and focuser change
NEW_EXPOSURE_DURATION_EX_GUIDER_CHECK = 19.3  # seconds; image download, plate solving (excl exposure)

# ********** EVD Preferences:
BURN_EXPOSURE = 240  # seconds per exposure
V_MAG_WARNING = 16.5  # a predicted V magnitude > this will trigger a warning line in Summary file.
ABSOLUTE_MAX_EXPOSURE_TIME = 600  # seconds
ABSOLUTE_MIN_EXPOSURE_TIME = 2.5  # seconds [20190318, was 3 seconds]
MIN_TOTAL_EXP_TIME_PER_FILTER = 9  # seconds, thus 4 [was 3] exposures max per filter for LPVs

# ********** MP parameters:
# defining MP color sequence as tuple of tuples: (filter, seconds exposure, repeats). Exps at V mag = 14.
MP_COLOR_V14 = (('I', 240, 1), ('R', 160, 1), ('V', 160, 2), ('R', 160, 1), ('I', 240, 1))


def make_df_fov(fov_directory=FOV_DIRECTORY, fov_names_selected=None):
    """
    Returns new, basic fov data frame, by reading FOV files (all or selected) in given directory_path.
    :param fov_directory: the directory_path from which to read FOV files.
    :param fov_names_selected: default = all FOV files within given directory_path.
    :return: basic data frame with columns: fov_name, main_target, fov_priority, obs_style,
         ra, dec. Index == fov_name.
    """
    fov_dict = make_fov_dict(fov_directory, fov_names_selected)
    fov_names = list(fov_dict.keys())
    df_fov = pd.DataFrame({'fov_name': fov_names})  # 1 column ('fov_name') only.
    # Add column of Fov objects, matched to column fov_names:
    df_fov['fov'] = [fov_dict[name] for name in fov_names]

    # Add other columns (directly from fov) for later convenience:
    df_fov['main_target'] = [fov_dict[name].main_target for name in fov_names]
    df_fov['fov_priority'] = [fov_dict[name].priority for name in fov_names]
    df_fov['obs_style'] = [fov_dict[name].observing_style for name in fov_names]
    df_fov['ra'] = [fov_dict[name].ra for name in fov_names]
    df_fov['dec'] = [fov_dict[name].dec for name in fov_names]
    df_fov['period'] = [fov_dict[name].period for name in fov_names]
    df_fov['target_type'] = [fov_dict[name].target_type for name in fov_names]
    df_fov['max_exposure'] = [fov_dict[name].max_exposure for name in fov_names]

    # Construct column 'radec' from 'ra' and 'dec':
    df_fov['radec'] = RaDec(0, 0)  # dummy value to be replaced (needed to set column object type).
    for ind in df_fov.index:
        ra = df_fov.loc[ind, 'ra']
        dec = df_fov.loc[ind, 'dec']
        df_fov.loc[ind, 'radec'] = RaDec(ra, dec)

    # Sort by fov name, set index to fov name.
    df_fov.sort_values(by='fov_name', inplace=True)
    df_fov.index = df_fov['fov_name']
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


def complete_df_fov_an(df_fov, user_update_tolerance_days=DEFAULT_UPDATE_TOLERANCE_DAYS,
                       an_string=None, site_name="DSW",
                       min_available_seconds=MIN_AVAILABLE_SECONDS_DEFAULT,
                       min_moon_degrees=MIN_MOON_DEGREES_DEFAULT,
                       remove_zero_an_priority=True,
                       remove_unobservables=True):
    if an_string is None or site_name == "":
        return df_fov
    an = Astronight(an_string, site_name)

    # Construct columns (specific to night and site) for available obs time, this astronight.
    df_fov = df_fov.assign(moon_deg=0.0) \
        .assign(start=an.local_middark_utc) \
        .assign(end=an.local_middark_utc) \
        .assign(mid=an.local_middark_utc) \
        .assign(seconds=0.0) \
        .assign(available=' - '.join(2 * [4 * ' '])) \
        .assign(an_priority=0.0) \
        .assign(an_priority_bars='')  # all dummy values to be overwritten later.

    # Fill in most columns.
    for ind in df_fov.index:
        ts_obs = an.ts_observable(df_fov.loc[ind, 'radec'], min_alt=an.site.min_altitude,
                                  min_moon_dist=min_moon_degrees)
        df_fov.loc[ind, 'moon_deg'] = df_fov.loc[ind, 'radec'].degrees_from(an.moon_radec)
        df_fov.loc[ind, 'start'] = ts_obs.start
        df_fov.loc[ind, 'end'] = ts_obs.end
        df_fov.loc[ind, 'mid'] = ts_obs.midpoint
        df_fov.loc[ind, 'seconds'] = ts_obs.seconds
        if ts_obs.seconds > 0:
            df_fov.loc[ind, 'available'] = ' - '.join([hhmm_from_datetime_utc(ts_obs.start),
                                                       hhmm_from_datetime_utc(ts_obs.end)])

    # Remove targets that can't be observed this astronight, *before* getting data from AAVSO:
    if remove_unobservables:
        enough_dark_time = df_fov['seconds'] >= min_available_seconds
        moon_dist_ok = df_fov['moon_deg'] >= min_moon_degrees
        is_observable = enough_dark_time & moon_dist_ok
        # print('Querying AAVSO for', str(sum(is_observable)), 'of', str(len(df_fov)), 'targets.')
        df_fov = df_fov[is_observable]

    # Update observations cache from AAVSO:
    loc = LocalObsCache()
    loc.update_fov_entries(df_fov, user_update_tolerance_days=user_update_tolerance_days)

    # Compute each target's priority for this astronight:
    for ind in df_fov.index:
        this_fov = df_fov.loc[ind, 'fov']
        df_fov.loc[ind, 'an_priority'] = loc.calc_an_priority(this_fov, an,
                                                              user_update_tolerance_days)
        max_bars = 16
        int_an_priority = int(round(df_fov.loc[ind, 'an_priority']))
        df_fov.loc[ind, 'an_priority_bars'] = \
            (8 * '.' + (max_bars - 8) * '#')[0: min(max_bars, int_an_priority)].ljust(max_bars)

    if remove_zero_an_priority:
        df_fov = df_fov[df_fov['an_priority'] > 0.0]

    return df_fov.sort_values(by=['mid', 'an_priority'], ascending=[True, False])


class LocalObsCache:
    """
    Holds a cache dataframe of most recent relevant observations for ~all FOVs.
    Can hold only one dataframe row per fov (however many filters constitute a previous obs).
    Will query AAVSO webobs site to refresh a database row if fov's main target looks too old.
    Cache dataframe columns are:
        fov_name [string]
        main_target [string]
        obs_style [string]
        cache_datetime: datetime this row was updated [datetime.datetime UTC]
        obs_datetime: datetime of most recent known observation [datetime.datetime UTC]
        obs_mag: magnitude of most recent observation [float]
        obs_mag_filter: filter in which obs_mag was measured [string]
    Typical usage: pl.make_an_plan('c:/Astro/ACP/AN20170525/planning.xlsx', exp_time_factor=0.75)
    """

    def __init__(self):
        # Read in local cache if it exists.
        if os.path.isfile(LOCAL_OBS_CACHE_FULLPATH):
            self.df_cache = self._read_cache_from_csv()
            need_to_create_empty_cache = self.df_cache is None
        else:
            need_to_create_empty_cache = True
        if need_to_create_empty_cache:
            #  Create *empty* dataframe with dtypes (incl. utc datetimes), write to cache file:
            self.df_cache = pd.DataFrame.from_dict(OrderedDict([
                ('fov_name', ['dummy']),
                ('main_target', ['dummy']),
                ('obs_style', ['dummy']),
                ('cache_datetime', [datetime.now(timezone.utc)]),
                ('obs_datetime', [datetime.now(timezone.utc)]),
                ('obs_mag', [0.0]),
                ('obs_mag_filter', ['dummy'])]))[:0]
            self.df_cache.index.name = 'row_index'
            csv_fullpath = self._write_cache_to_csv()  # empty cache to csv
            print('LocalObsCache: wrote new, empty cache file to ' + csv_fullpath)
        print('LocalObsCache opened; ' + str(len(self.df_cache)) + ' fovs.')

    def update_fov_entries(self, df_fov,
                           user_update_tolerance_days=DEFAULT_UPDATE_TOLERANCE_DAYS,
                           max_fovs_since_write=6):
        """
        For each fov available this night (in df_fov_list), update the cache.
        :param df_fov: df_fov (typically of fovs available this night) [pandas DataFrame].
        :param user_update_tolerance_days: pass-through parm [float].
        :param max_fovs_since_write: controls frequence of writes to cache.
        :return: number of fovs updated (fn effect is to update this class's cache dataframe.
        """
        fovs_since_write = 0
        for fov in df_fov['fov']:
            need_to_write_csv = (fovs_since_write >= max_fovs_since_write - 1)
            self.update_one_fov_entry(fov, user_update_tolerance_days, write_csv=need_to_write_csv)
            if need_to_write_csv:
                fovs_since_write = 0
            else:
                fovs_since_write += 1
        self._write_cache_to_csv()  # ensure cache written at end.

    def update_one_fov_entry(self, fov, user_update_tolerance_days=DEFAULT_UPDATE_TOLERANCE_DAYS,
                             write_csv=False):
        """
        This class's engine. Updates cache's entry for one fov, if entry is too aged.
        :param fov: fov to update in cache now [Fov object]
        :param user_update_tolerance_days: pass-through parm [float]
        :param write_csv:
        :return: if cache was updated, datetime (UTC) of new obs; else None.
        """
        # TODO: If query to AAVSO yields no latest obs, put some placeholder with cache_dt at least.
        if fov is None:
            raise FovError
        main_target = fov.main_target
        # self._curate_df_cache(fov_name, main_target)

        #  Determine whether update is needed, return if not.
        cache_row_pre_exists = fov.fov_name.lower() in list(self.df_cache['fov_name'].str.lower())
        if cache_row_pre_exists:
            now = datetime.now(timezone.utc)
            current_cache_datetime = self.df_cache.loc[fov.fov_name, 'cache_datetime']
            update_age = (now - current_cache_datetime).total_seconds() / (24 * 3600)
            if user_update_tolerance_days is None:
                update_tolerance_days = DEFAULT_UPDATE_TOLERANCE_DAYS
            else:
                update_tolerance_days = user_update_tolerance_days
            entry_fresh_enough = update_age <= update_tolerance_days
            if entry_fresh_enough:
                return self.df_cache.loc[fov.fov_name, 'obs_datetime']  # skip updating

        # Update fov's cache entry, from AAVSO webobs.
        obs_style = fov.observing_style
        obs_style_lower = obs_style.lower()
        target_type_lower = fov.target_type.lower()
        if target_type_lower == 'standard':
            return None
        if obs_style_lower == 'stare':
            num_obs = 200
        else:
            num_obs = 100
        print('AAVSO webobs query ' + fov.target_type +
              ' \'' + main_target + '\'...', end='', flush=True)
        recent_observations = AavsoWebobs(star_id=main_target, num_obs=num_obs)  # from AAVSO
        print('ok.', end='', flush=True)
        latest_obs_df = None  # default if no matches.
        if (obs_style_lower, target_type_lower) == ('lpv', 'mira'):
            latest_obs_df = self._latest_single_obs(fov, obs_style, recent_observations,
                                                    allow_filters=['V'])
        elif (obs_style_lower, target_type_lower) == ('lpv', 'lpv'):
            latest_obs_df = self._latest_single_obs(fov, obs_style, recent_observations,
                                                    allow_filters=['V', 'R'])
        elif obs_style_lower == 'monitor' and target_type_lower != 'astrometric':
            latest_obs_df = self._latest_single_obs(fov, obs_style, recent_observations,
                                                    allow_filters=['V', 'R'])
        elif obs_style_lower == 'stare':
            latest_obs_df = self._latest_stare_obs(fov, recent_observations,
                                                   allow_filters=['V', 'R'])
        else:
            print('\n*** WARNING: for fov \'' + fov.fov_name + '(obs_style, target_type) = (' +
                  obs_style + ', ' + fov.target_type + ') not understood.', end='', flush=True)
        if cache_row_pre_exists:
            self.df_cache = latest_obs_df.combine_first(self.df_cache)  # overwrites.
        else:
            #  This else-block is kludge for pandas' mis-handling of append to empty DataFrame.
            if len(self.df_cache) >= 1:
                self.df_cache = self.df_cache.append(latest_obs_df, sort=True)
            else:
                self.df_cache = latest_obs_df.copy()
        if write_csv:
            self._write_cache_to_csv()
            print('..csv written.', end='', flush=True)
        print('')
        if latest_obs_df is None:
            return None
        return latest_obs_df.iloc[0].loc['obs_datetime']  # obs datetime, to signal OK.

    def _latest_single_obs(self, fov, obs_style, recent_observations, allow_filters):
        """
        Takes a AavsoWebObs object and returns a pandas dataframe ready for inclusion
           in LocalCacheObs dataframe df_cache.
           Single-observation case (not stare).
        :param fov: fov to investigate for recent single observations [Fov object]
        :param obs_style: [string] ('Monitor' or 'LPV')
        :param recent_observations: recent observations for fov_name [AavsoWebObs object].
        :param allow_filters: list of filters [string] to include in finding latest observation.
        :return: 1-row dataframe of relevant data about latest stare observation for this fov_name;
                 return (with some None values) if no qualifying observation is found.
        """
        allow_filters_lower = [f.lower() for f in allow_filters]
        table_filters_lower = recent_observations.table['filter'].str.lower()
        rows_to_keep = [f.lower() in allow_filters_lower for f in table_filters_lower]
        if sum(rows_to_keep) <= 0:
            latest_obs = None
        else:
            latest_obs = recent_observations.table[rows_to_keep].nlargest(1, 'jd').iloc[0]
        if latest_obs is None:
            #  If no qualified observation found within webobs query,
            #  construct placeholder row in df_cache, to prevent repeating query needlessly.
            latest_obs_df = pd.DataFrame.from_dict(OrderedDict([
                ('fov_name', fov.fov_name),
                ('main_target', fov.main_target),
                ('obs_style', fov.observing_style),
                ('cache_datetime', [datetime.now(timezone.utc)]),
                ('obs_datetime', [None]),
                ('obs_mag', [None]),
                ('obs_mag_filter', [None])]))
            for column_name in ['cache_datetime']:
                latest_obs_df[column_name] = [x.to_pydatetime()
                                              for x in latest_obs_df[column_name]]
        else:
            latest_obs_df = pd.DataFrame.from_dict(OrderedDict([
                ('fov_name', fov.fov_name),
                ('main_target', fov.main_target),
                ('obs_style', fov.observing_style),
                ('cache_datetime', [datetime.now(timezone.utc)]),
                ('obs_datetime', [datetime_utc_from_jd(latest_obs.jd)]),
                ('obs_mag', [latest_obs.mag]),
                ('obs_mag_filter', [latest_obs.loc['filter']])]))
            for column_name in ['cache_datetime', 'obs_datetime']:
                latest_obs_df[column_name] = [x.to_pydatetime()
                                              for x in latest_obs_df[column_name]]
        latest_obs_df.index = latest_obs_df['fov_name'].copy()
        latest_obs_df.index.name = 'row_index'
        return latest_obs_df

    def _latest_stare_obs(self, fov, recent_observations, allow_filters):
        """
        Takes a AavsoWebObs object and returns a 1-row pandas dataframe ready for inclusion
            in LocalCacheObs dataframe df_cache.
            Stare case (multiple observations in one night), typically for eclipsers.
        :param fov: fov to investigate for recent stare observations [Fov object]
        :param recent_observations: recent observations for fov_name [AavsoWebObs object].
        :param allow_filters: list of filters [string] to include in finding latest observation.
        :return: dataframe of relevant data about latest stare observation for this fov_name;
                 return (with some None values) if no qualifying stare observation is found.
        """
        if len(recent_observations.table) <= MIN_ROWS_ONE_STARE:
            return None

        # Find latest qualifying stare in each filter, return latest observation of latest stare.
        latest_stare_obs_df = None
        for this_filter in allow_filters:
            stare_already_found_this_filter = False
            this_filter_lower = this_filter.lower()
            table_filters_lower = recent_observations.table['filter'].str.lower()
            rows_to_keep = [f.lower() == this_filter_lower for f in table_filters_lower]
            table_this_filter = recent_observations.table[rows_to_keep].sort_values(by='jd',
                                                                                    ascending=False)
            num_tests = len(table_this_filter) - MIN_ROWS_ONE_STARE + 1
            if num_tests >= 1:
                for first_test_irow in range(0, num_tests):
                    if not stare_already_found_this_filter:
                        test_latest_jd = table_this_filter['jd'] \
                            .iloc[first_test_irow]
                        test_earliest_jd = table_this_filter['jd'] \
                            .iloc[first_test_irow + MIN_ROWS_ONE_STARE - 1]

                        if test_latest_jd - test_earliest_jd <= MAX_DAYS_ONE_STARE:
                            stare_already_found_this_filter = True
                            if latest_stare_obs_df is None:
                                need_to_replace = True
                            else:
                                candidate_datetime = datetime_utc_from_jd(test_latest_jd)
                                existing_datetime = latest_stare_obs_df.iloc[0].loc['obs_datetime']
                                need_to_replace = (candidate_datetime > existing_datetime)

                            if need_to_replace:
                                latest_stare_obs = table_this_filter.iloc[first_test_irow]
                                latest_stare_obs_df = pd.DataFrame.from_dict(OrderedDict([
                                    ('fov_name', fov.fov_name),
                                    ('main_target', fov.main_target),
                                    ('obs_style', fov.observing_style),
                                    ('cache_datetime', [datetime.now(timezone.utc)]),
                                    ('obs_datetime', [datetime_utc_from_jd(latest_stare_obs.jd)]),
                                    ('obs_mag', [latest_stare_obs.mag]),
                                    ('obs_mag_filter', [latest_stare_obs.loc['filter']])]))
                                for column_name in ['cache_datetime', 'obs_datetime']:
                                    latest_stare_obs_df[column_name] = \
                                        [x.to_pydatetime()
                                         for x in latest_stare_obs_df[column_name]]
                                latest_stare_obs_df.index = latest_stare_obs_df['fov_name'].copy()
                                latest_stare_obs_df.index.name = 'row_index'

        if latest_stare_obs_df is None:
            #  If no qualified stare observation found within webobs query,
            #  construct placeholder row in df_cache, to prevent repeating query needlessly.
            latest_stare_obs_df = pd.DataFrame.from_dict(OrderedDict([
                ('fov_name', fov.fov_name),
                ('main_target', fov.main_target),
                ('obs_style', fov.observing_style),
                ('cache_datetime', [datetime.now(timezone.utc)]),
                ('obs_datetime', [None]),
                ('obs_mag', [None]),
                ('obs_mag_filter', [None])]))
            for column_name in ['cache_datetime']:
                latest_stare_obs_df[column_name] = [x.to_pydatetime()
                                                    for x in latest_stare_obs_df[column_name]]
        latest_stare_obs_df.index = latest_stare_obs_df['fov_name'].copy()
        latest_stare_obs_df.index.name = 'row_index'
        return latest_stare_obs_df

    @staticmethod
    def _read_cache_from_csv():
        cache = pd.read_csv(LOCAL_OBS_CACHE_FULLPATH, index_col=0)
        if len(cache) <= 0:
            return None
        for column_name in ['cache_datetime', 'obs_datetime']:
            if column_name not in cache.columns:
                return None
        # Parse cache_datetime column.
        cache['cache_datetime'] = [datetime.strptime(s, DT_FMT) for s in cache['cache_datetime']]

        # Parse obs_datetime column.
        for row_index in cache.index:
            if str(cache.loc[row_index, 'obs_datetime']).lower() != 'none':
                cache.loc[row_index, 'obs_datetime'] = \
                    datetime.strptime(cache.loc[row_index, 'obs_datetime'], DT_FMT)
                cache.loc[row_index, 'obs_mag'] = float(cache.loc[row_index, 'obs_mag'])
            else:
                cache.loc[row_index, 'obs_datetime'] = None
                cache.loc[row_index, 'obs_mag'] = None
                cache.loc[row_index, 'obs_mag_filter'] = None
        return cache

    def _write_cache_to_csv(self):
        # Very specifically writes datetimes in format: '2017-02-07 03:34:45.786374+0000'
        dt_format = '{:' + DT_FMT + '}'
        lines = [','.join(['row_index', 'fov_name', 'main_target', 'obs_style',
                           'cache_datetime', 'obs_datetime', 'obs_mag', 'obs_mag_filter']) + '\n']
        for row_index in self.df_cache.index:
            row = self.df_cache.loc[row_index]
            if row['obs_datetime'] is None or isinstance(row['obs_datetime'], type(pd.NaT)):
                line = ','.join([row_index, row['fov_name'], row['main_target'], row['obs_style'],
                                 dt_format.format(row['cache_datetime']),
                                 'None', 'None', 'None']) + '\n'
            else:
                line = ','.join([row_index, row['fov_name'], row['main_target'], row['obs_style'],
                                 dt_format.format(row['cache_datetime']),
                                 dt_format.format(row['obs_datetime']),
                                 '{:.4f}'.format(row['obs_mag']), row['obs_mag_filter']]) + '\n'
            lines.append(line)

        with open(LOCAL_OBS_CACHE_FULLPATH, 'w') as f:
            f.writelines(lines)
        # print("Cache written: " + str(len(self.df_cache)) + ' fovs.')
        return LOCAL_OBS_CACHE_FULLPATH

    def calc_an_priority(self, fov, an, user_update_tolerance_days=DEFAULT_UPDATE_TOLERANCE_DAYS):
        """
        Calculates astronight priority for one fov.
        :param fov:
        :param user_update_tolerance_days: pass-through parm [float].
        :param an: Astronight object for the night in question.
        :return: an_priority, from fov_priority and age of most recent obs [float].
        """
        if fov is None:
            print('LOC.calc_an_priority: fov \'' + fov + '\' not found in fov_dict.')
            return None
        if (fov.priority is None) or (fov.target_type.lower() == 'standard'):
            return 0
        if fov.priority <= 0:
            return 0

        # self.update_one_fov_entry(fov, user_update_tolerance_days, write_csv=True)
        if fov.fov_name not in self.df_cache.index:
            return 2 * fov.priority  # the maximum, since no latest obs was accessible.
        latest_obs = self.df_cache.loc[fov.fov_name]
        if latest_obs.obs_datetime is None:
            return 2 * fov.priority  # the maximum, since no latest obs was accessible.
        jd_latest_obs = jd_from_datetime_utc(latest_obs.obs_datetime)
        age_days = an.local_middark_jd - jd_latest_obs
        return fov.calc_priority_score(age_days)

    def _curate_df_cache(self, fov_name, main_target):
        """
        Cull damaged records from self.df_cache.
           *** Deactivated 2017-02-07 pending manual debugging. ***
        :param fov_name:
        :param main_target:
        :return: [nothing]
        """
        # Curation: If main_target is corrupt, delete all cache lines for that fov.
        if main_target is None:
            rows_to_delete = self.df_cache['fov_name'].str.lower() == fov_name.lower()
            self.df_cache = self.df_cache[rows_to_delete == False]
            return

        # Curation: If fov and target names don't match, delete all such fov and target lines.
        rows_with_wrong_target = \
            (self.df_cache['fov_name'].str.lower() == fov_name.lower()) & \
            (self.df_cache['main_target'].str.lower() != main_target.lower())
        rows_to_keep = [not row for row in rows_with_wrong_target]
        self.df_cache = self.df_cache[rows_to_keep]
        rows_with_wrong_fov = \
            (self.df_cache['main_target'].str.lower() == main_target.lower()) & \
            (self.df_cache['fov_name'].str.lower() != fov_name.lower())
        rows_to_keep = [not row for row in rows_with_wrong_fov]
        self.df_cache = self.df_cache[rows_to_keep]

    def __str__(self):
        return 'LocalObsCache object with ' + str(len(self.df_cache)) + \
               ' observations.'

    def __repr__(self):
        return 'planning.LocalObsCache()'


class AavsoWebobs:
    """
    Simple class: one object:one star. Holds dataframe for one star from AAVSO's webobs database.
    Also updates local cache file (to unneeded future calls to webobs).
    For Observation Styles: LPV, Monitor, and Stare; no need for Standard or Burn.
    Usage: table = AavsoWebobs("AU Aur") [for one obs/night], or
           table = AavsoWebobs("ST Tri", stare=True) [for at least 10 obs/night in filter].
    """

    def __init__(self, star_id=None, num_obs=AAVSO_WEBOBS_ROWS_TO_GET, dataframe=None):
        if dataframe is not None:
            self.table = dataframe  # typically for testing only.
            self.star_id = self.table['target_name'].iloc[0]
        else:
            self.table = get_aavso_webobs_raw_table(star_id, num_obs=num_obs)  # normal case
            self.star_id = star_id


# def get_local_aavso_reports(report_dir=None, earliest_an=None):
#     pass
    #     report_dict = {}
    #     for root, dirs, files in os.walk('C:/Astro/Images/Borea Photrix/'):
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
#         # TODO: finish writing get_local_obs_age_dict()
#         """
#         report_dir: directory_path in which all relevant AAVSO reports reside, as
#           "C:/Astro/2016/Photometry".
#         target_an: target astronight from which to count days, as "20151216".
#         limit_days: days into the past to look up old AAVSO reports.
#         Returns dict of (fov_name, days_since_last_local_obs).
#         """
#         pass
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


# ---------------------------------------------

def make_an_roster(an_date_string, output_directory, site_name='DSW', instrument_name='Borea',
                   user_update_tolerance_days=DEFAULT_UPDATE_TOLERANCE_DAYS,
                   exp_time_factor=1, min_an_priority=6):
    # TODO: recode download loop to only download those FOVs visible tonight & for which priority might
    #  be high enough (read from csv: some might already be known to have been recently observed).
    """
    Generates new .csv file containing info on each fov available this astronight.
       Typical usage: pl.make_an_roster("20170127", "C:/Astro/ACP/AN20170127/",
       user_update_tolerance_days=0.1, exp_time_factor=0.8)
    :param an_date_string: as '20170127'. Date of the evening to plan for [string]
    :param output_directory: directory_path in which to write Roster csv file [string]
    :param site_name: [string]
    :param instrument_name: [string]
    :param user_update_tolerance_days: esp for user to force update [float]
    :param exp_time_factor: multiply *raw* exp times by this; typically 0.6-0.9 [float]
    :param min_an_priority: hide Monitor and LPV targets with an_priority < this [float]
    :return: tuple of number of fovs, each obs style: (n_std, n_monitor_lpv, n_stare). [ints]
    """

    an = Astronight(an_date_string=an_date_string, site_name=site_name)
    df_fov = make_df_fov(fov_directory=FOV_DIRECTORY, fov_names_selected=None)
    print(str(len(df_fov)), 'FOVs read.')
    instrument = Instrument(instrument_name)
    an_year = int(an_date_string[0:4])
    an_month = int(an_date_string[4:6])
    an_day = int(an_date_string[6:8])
    day_of_week = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'] \
        [datetime(an_year, an_month, an_day).weekday()]
    lines_header = ['ROSTER file for     ' + an_date_string + '   ' + day_of_week,
                    '     as generated by photrix ' +
                    '{:%Y-%m-%d %H:%M  UTC}'.format(datetime.now(timezone.utc)),
                    '     using exposure time factor = ' + '{:5.3f}'.format(exp_time_factor),
                    an.acp_header_string().replace(',', ' '),
                    '; Site=' + site_name + '   Instrument=' + instrument_name +
                    '    min.alt = ' + '{:.1f}'.format(an.site.min_altitude) + u'\N{DEGREE SIGN}']

    # Handle obs_style = 'Standard':
    lines_std = ['\n\n\nSTANDARD roster for ' + an_date_string + ': ' + 50 * '-',
                 ',fov,fov, avail_utc,transit,minutes,   stars']
    df_fov_std = filter_df_fov_by_obs_styles(df_fov, obs_style_list=['Standard'])
    df_fov_std = complete_df_fov_an(df_fov_std, user_update_tolerance_days,
                                    an_string=an_date_string, site_name=site_name,
                                    min_available_seconds=MIN_AVAILABLE_SECONDS_DEFAULT,
                                    min_moon_degrees=MIN_MOON_DEGREES_DEFAULT,
                                    remove_zero_an_priority=False, remove_unobservables=True)
    for fov_index in df_fov_std.index:
        fov_name = df_fov_std.loc[fov_index, 'fov_name']
        available = df_fov_std.loc[fov_index, 'available']
        this_fov = Fov(fov_name)
        transit_hhmm = hhmm_from_datetime_utc(an.transit(RaDec(this_fov.ra, this_fov.dec)))
        exp_data = make_fov_exposure_data(fov_name, an, fov_dict=None, instrument=instrument,
                                          exp_time_factor=exp_time_factor,
                                          force_autoguide=False)  # default=autoguide iff exp-times warrant.
        if exp_data is None:
            return  # fail
        _, _, _, target_overhead, repeat_duration = exp_data
        minutes = (target_overhead + repeat_duration) / 60.0
        n_stars = len(this_fov.aavso_stars)
        this_fov_line = ',' + fov_name + ',' + fov_name + ', ' + available + ',' + \
                        "=\"" + transit_hhmm + "\"" + ',' + str(int(minutes)) + \
                        ',' + '{:3d}'.format(n_stars)  # formatting to placate Excel csv weirdness.
        lines_std.append(this_fov_line)

    # Handle obs_style = 'Stare':
    lines_stare_high_priority = \
        ['\n\n\nSTARE roster for ' + an_date_string + ': ' + 50 * '-',
         ',fov,fov, avail_utc,transit,min/rpt,   an_priority,,period,  events']
    lines_stare_low_priority = \
        ['\n\n\nSTARE roster (alternate; low-priority) for ' + an_date_string + ': ' + 50 * '-',
         ',fov,fov, avail_utc,transit,min/rpt,   an_priority,,period,  events']
    df_fov_stare = filter_df_fov_by_obs_styles(df_fov, obs_style_list=['Stare'])
    # Process each fov equally through most of the code,
    # then only in the last if block, write to one list or the other.
    df_fov_stare = complete_df_fov_an(df_fov_stare, user_update_tolerance_days,
                                      an_string=an_date_string, site_name=site_name,
                                      min_available_seconds=MIN_AVAILABLE_SECONDS_STARE,
                                      min_moon_degrees=MIN_MOON_DEGREES_STARE,
                                      remove_zero_an_priority=False, remove_unobservables=True)
    for fov_index in df_fov_stare.index:
        row = df_fov_stare.loc[fov_index]
        fov_name = row.loc['fov_name']
        available = row.loc['available']
        this_fov = Fov(fov_name)
        transit_hhmm = hhmm_from_datetime_utc(an.transit(RaDec(this_fov.ra, this_fov.dec)))
        exp_data = make_fov_exposure_data(fov_name, an, fov_dict=None, instrument=instrument,
                                          exp_time_factor=exp_time_factor,
                                          force_autoguide=False)  # default=autoguide iff exp-times warrant.
        if exp_data is None:
            return  # fail
        _, _, _, target_overhead, repeat_duration = exp_data
        minutes = (target_overhead + repeat_duration) / 60.0
        an_priority = row.loc['an_priority']
        an_priority_bars = row.loc['an_priority_bars']
        period = row.loc['period']
        row_ts = Timespan(row.loc['start'], row.loc['end'])  # Timespan object for this row.

        # For now, we will consider that each Stare FOV wants either minima or maxima but not both.
        event_type_string = STARE_EVENT_TYPES.get(this_fov.target_type.lower(), None)
        if event_type_string is None:
            print(this_fov.fov_name + ': probable bad target_type in FOV.')
        do_minima = event_type_string.lower().startswith("min")
        do_maxima = event_type_string.lower().startswith("max")

        #  Start with an *empty* dataframe of events, with correct dtypes:
        df_events = pd.DataFrame.from_dict(OrderedDict(
            [('event_type', 'dummy_type'),
             ('utc', [datetime.now(timezone.utc)])]))[:0]
        if do_minima:
            list_primary_mins = event_utcs_in_timespan(this_fov.JD_faint, this_fov.period, row_ts)
            if list_primary_mins is None:
                primaries_exist = False
            else:
                primaries_exist = len(list_primary_mins) >= 1
            if primaries_exist:
                df_primary_mins = pd.DataFrame.from_dict(dict([('utc', list_primary_mins,)]))
                df_primary_mins['event_type'] = "1'"
                df_events = df_events.append(df_primary_mins, sort=True)

            list_secondary_mins = event_utcs_in_timespan(this_fov.JD_second, this_fov.period,
                                                         row_ts)
            if list_secondary_mins is None:
                secondaries_exist = False
            else:
                secondaries_exist = len(list_secondary_mins) >= 1
            if secondaries_exist:
                df_secondary_mins = pd.DataFrame.from_dict(dict([('utc', list_secondary_mins,)]))
                df_secondary_mins['event_type'] = "2'"
                df_events = df_events.append(df_secondary_mins, sort=True)

        if do_maxima:
            list_maxima = event_utcs_in_timespan(this_fov.JD_bright, this_fov.period, row_ts)
            if list_maxima is None:
                maxima_exist = False
            else:
                maxima_exist = len(list_maxima) >= 1
            if maxima_exist:
                df_maxima = pd.DataFrame.from_dict(dict([('utc', list_maxima,)]))
                df_maxima['event_type'] = "max"
                df_events = df_events.append(df_maxima, sort=True)

        if len(df_events) >= 1:
            motive = this_fov.motive
            df_events.sort_values(by='utc', inplace=True)
            events_string = '   '
            for row in df_events.itertuples():
                events_string += str(row.event_type) + "=" + hhmm_from_datetime_utc(row.utc) + '  '
            this_fov_line = ',' + fov_name + ',' + fov_name + ',' + available + ',' + \
                            "=\"" + transit_hhmm + "\"" + ',' + str(int(minutes)) + ',' + \
                            str(int(round(an_priority))) + ' ,' + an_priority_bars + ',' + \
                            '{:7.3f}'.format(period) + ' ,' + events_string + ',' + \
                            "\"  " + motive + "\""  # formatting to placate Excel csv weirdness.
            if an_priority >= STARE_AN_PRIORITY_DIVIDER:
                lines_stare_high_priority.append(this_fov_line)
            else:
                lines_stare_low_priority.append(this_fov_line)

    # Handle obs_style = 'Monitor' or 'LPV':
    lines_mon_lpv = ['\n\n\nMONITOR / LPV roster for ' + an_date_string + ': ' + 50 * '-',
                     ',fov,fov,avail_utc,transit,minutes,   an_priority']
    df_fov_mon_lpv = filter_df_fov_by_obs_styles(df_fov, obs_style_list=['Monitor', 'LPV'])
    df_fov_mon_lpv = filter_df_fov_by_fov_priority(df_fov_mon_lpv,
                                                   min_fov_priority=0.5, include_std_fovs=False)
    df_fov_mon_lpv = complete_df_fov_an(df_fov_mon_lpv, user_update_tolerance_days,
                                        an_string=an_date_string,
                                        site_name=site_name,
                                        min_available_seconds=MIN_AVAILABLE_SECONDS_DEFAULT,
                                        min_moon_degrees=MIN_MOON_DEGREES_DEFAULT,
                                        remove_zero_an_priority=True, remove_unobservables=True)
    for fov_index in df_fov_mon_lpv.index:
        an_priority = df_fov_mon_lpv.loc[fov_index, 'an_priority']
        if an_priority >= min_an_priority:
            fov_name = df_fov_mon_lpv.loc[fov_index, 'fov_name']
            available = df_fov_mon_lpv.loc[fov_index, 'available']
            this_fov = Fov(fov_name)
            transit_hhmm = hhmm_from_datetime_utc(an.transit(RaDec(this_fov.ra, this_fov.dec)))
            exp_data = make_fov_exposure_data(fov_name, an, fov_dict=None, instrument=instrument,
                                              exp_time_factor=exp_time_factor,
                                              force_autoguide=False)  # so autoguide iff exp-times warrant.
            if exp_data is None:
                return  # fail
            _, _, _, target_overhead, repeat_duration = exp_data
            minutes = (target_overhead + repeat_duration) / 60.0
            an_priority_bars = df_fov_mon_lpv.loc[fov_index, 'an_priority_bars']
            motive = Fov(fov_name).motive
            this_fov_line = ',' + fov_name + ',' + fov_name + ', ' + available + ',' + \
                            "=\"" + transit_hhmm + "\"" + ',' + str(int(minutes)) + ',' + \
                            str(int(round(an_priority))) + ' ,' + an_priority_bars + ',' + \
                            "\"  " + motive + "\""  # formatting to placate Excel csv weirdness.
            lines_mon_lpv.append(this_fov_line)

    # Assemble all output lines:
    lines_all = lines_header + \
                lines_std + \
                lines_stare_high_priority + lines_stare_low_priority + \
                lines_mon_lpv

    # Write all lines to file:
    os.makedirs(output_directory, exist_ok=True)
    output_fullpath = os.path.join(output_directory, 'Roster_' + an.an_date_string + '.csv')
    csv_written = False
    while not csv_written:
        try:
            with open(output_fullpath, 'w') as this_file:
                this_file.write('\n'.join(lines_all))
            csv_written = True
        except PermissionError:
            input('***** CLOSE file \'' + output_fullpath + '\' and hit Enter.')
    print('Done.')


class Plan:
    """ Holds all data for one ACP plan.
    :param plan_id: name of this plan, minimalist, as 'C' [string]
    """

    def __init__(self, plan_id=None, plan_comment=None):
        self.plan_id = plan_id
        self.plan_comment = plan_comment.strip()
        self.directives = []  # as parsed from user input in parse_excel().

        # Values populated in make_events():
        self.utc_quitat = None  # forced stop time (at end of last event)
        self.afinterval = None  # in minutes; default=None if no afinterval requested for this plan.
        self.sets_requested = 1  # default
        self.chain_destination = None  # next plan filename, None if no chaining requested for this plan.
        self.events = []  # holds only ONE element per intended event, no matter how many SETS in a plan.

        # Values populated in make_timeline():
        self.utc_start = None  # actual start time (computed later)
        self.utc_end = None  # actual end time (computed later)
        self.sets_completed = 0  # actual number of set cycles completed (integer)
        self.afinterval_autofocus_count = 0  # count of autofocuses caused by AFINTERVAL, this set.

        # Lists of text lines to go before and after main body of lines (from make_events()):
        self.summary_pre_lines = []
        self.summary_post_lines = []
        self.acp_pre_lines = []
        self.end_warning_lines = []
        self.acp_post_lines = []

    def quitat_reached_at(self, utc):
        if self.utc_quitat is None:
            return False
        return utc >= self.utc_quitat

    def __str__(self):
        return 'Plan object: ' + self.plan_id


class Directive:
    """ Holds all initial data for one user-given directive (e.g., one cell in Excel spreadsheet).
    :param type: type of directive, from approved list, e.g., 'CHILL' or 'fov'  [string, case-insens.]
    :param spec: dictionary of specification data to hold, depends on directive type [directory_path].
    """

    def __init__(self, type, spec_dict):
        self.type = type.lower()
        self.spec = spec_dict
        self.comment = ''

    def __str__(self):
        return 'Directive object: ' + self.type


class Event:
    """ Holds all data for one event to be executed (per set).
        Each Event object will result in at least one line in summary doc and in ACP plan file.
    """

    def __init__(self, event_type, summary_text, acp_lines, duration_total=0, duration_dict=None,
                 target_name=None, ra=None, dec=None):
        self.type = event_type.lower()
        self.summary_text = summary_text  # starting text for summary document
        self.summary_lines = []  # final lines for summary document
        self.acp_lines = acp_lines  # list (always) of lines for ACP plan file [list of strings]
        self.duration_total = duration_total  # total for event; 0 for waituntil, quitat, etc.
        # .duration_dict exists only for exposure-event types: burn, stare, fov, and image,
        #    as: {'target_overhead': in sec, 'repeat_count': n, 'counts': [n], 'exp_times': [in sec]} :
        self.duration_dict = duration_dict  # dict describing durations of indiv exposures, incl overheads.
        self.utc_end = None  # for waituntil.
        self.target_name = target_name  # for any exposure event type (= fov name for most)
        self.ra = ra  # string, for exposure-event types
        self.dec = dec  # "

        # For this event's summary line. Values later populated by make_timeline():
        self.status = None  # any of: 'ok', 'chain', 'quitat', 'wait'
        self.utc_summary_display = None  # stored only for first SET.
        self.min_altitude = None  # in degrees, for ALL sets in this plan.

    def calc_actual_duration(self, utc_start, utc_quitat, afinterval, utc_most_recent_autofocus):
        if utc_quitat is None:
            return self.duration_total, 0, utc_most_recent_autofocus
        if utc_start >= utc_quitat:
            return 0, 0, utc_most_recent_autofocus
        if self.duration_dict is None:
            return None
        total_exp_time = self.duration_dict['repeat_count'] * sum([c * e for (c, e) in
                                                                   zip(self.duration_dict['counts'],
                                                                       self.duration_dict['exp_times'])])
        exposure_count = self.duration_dict['repeat_count'] * sum(self.duration_dict['counts'])
        overhead_per_exposure = (self.duration_total - total_exp_time -
                                 self.duration_dict['target_overhead']) / exposure_count
        utc_running = utc_start + timedelta(seconds=self.duration_dict['target_overhead'])
        event_autofocus_count = 0  # accumulator for this event.
        for i_repeat in range(self.duration_dict['repeat_count']):
            for c, e in zip(self.duration_dict['counts'], self.duration_dict['exp_times']):
                for i_exp in range(c):
                    # Update clock for any AFINTERVAL-triggered autofocus:
                    if afinterval is not None:
                        minutes_since_last_autofocus = \
                            (utc_running - utc_most_recent_autofocus).total_seconds() / 60.0
                        if minutes_since_last_autofocus > afinterval:
                            utc_running += timedelta(seconds=AUTOFOCUS_DURATION)
                            utc_most_recent_autofocus = utc_running
                            event_autofocus_count += 1
                    # Update clock for exposure itself:
                    utc_running += timedelta(seconds=overhead_per_exposure + e)
                    # Terminate event if QUITAT time has passed:
                    if utc_running >= utc_quitat:
                        return (utc_running - utc_start).total_seconds(), \
                               event_autofocus_count, utc_most_recent_autofocus
        return (utc_running - utc_start).total_seconds(), event_autofocus_count, utc_most_recent_autofocus

    def calc_lower_altitude(self, an, utc1, utc2):
        longitude, latitude = an.site.longitude, an.site.latitude
        longitude_hex, latitude_hex = degrees_as_hex(longitude), degrees_as_hex(latitude)
        target_radec = RaDec(self.ra, self.dec)
        _, alt_deg_utc1 = az_alt_at_datetime_utc(longitude_hex, latitude_hex, target_radec, utc1)
        _, alt_deg_utc2 = az_alt_at_datetime_utc(longitude_hex, latitude_hex, target_radec, utc2)
        return min(alt_deg_utc1, alt_deg_utc2)

    def __str__(self):
        return 'Event object: ' + self.summary_text


def make_an_plan(plan_excel_path='c:/24hrs/Planning.xlsx', site_name='DSW', instrument_name='Borea',
                 fov_dict=None, earliest_an_start_hhmm=None, exp_time_factor=1):
    """  Main user fn to take sketch Excel file and generate Summary and ACP Plan files.
    :param plan_excel_path: full path to Excel file holding all info for one night's observations.
    :param site_name: a Site object for location of observations.
    :param instrument_name: an Instrument object for scope to be used.
    :param fov_dict: fov_dict if available, default=None to generate new fov_dict (normal case).
    :param earliest_an_start_hhmm: 'hhmm' time to start plan, default=None for 'earliest possible' (normal case).
    :param exp_time_factor: multiply *raw* exp times by this; typically 0.6-0.9 [float]
    :return: Writes out Summary file with dateline, and one or more ACP plan files.
    Typical usage: pl.make_an_plan('c:/Astro/ACP/AN20170525/planning.xlsx', exp_time_factor=0.7)
    """

    # TODO: LocalObsCache updates only for fovs actually used, not including Burns.  ???meant for roster???
    plan_list, an = parse_excel(plan_excel_path, site_name)

    reorder_directives(plan_list)

    if fov_dict is None:
        fov_dict = make_fov_dict()
    instrument = Instrument(instrument_name)

    make_events(plan_list, instrument, fov_dict, an=an, exp_time_factor=exp_time_factor)

    output_directory = os.path.split(plan_excel_path)[0]

    make_timeline(plan_list, an=an, earliest_hhmm=earliest_an_start_hhmm)

    make_acp_plan_files(plan_list, an, output_directory, exp_time_factor)

    make_summary_file(plan_list, fov_dict, an, output_directory, exp_time_factor)


def parse_excel(excel_path, site_name='DSW'):
    """
    Parses sketch Excel file, returns a list of Plan objects containing all directives and the
        relevant Astronight object.
    :param excel_path: full path to Excel file holding all info for one night's observations [str].
    :param site_name: a Site object for location of observations [string]
    :return: list of Plan objects, astronight object (2-tuple)

    ----- Target types & their syntax:
    FOV_name  ::  for LPV, standards, and other once-per-night targets having FOV files,
       e.g., "FF Lyr" and "Std_SA32".
    STARE nnn FOV_name  ::  for stare targets; nnn=number of obs cycles,
       e.g., "STARE 100 ST Tri" (typically time-limited by QUITAT, not by number of obs cycles).
    BURN  FOV_name  RA  Dec  ::  for 240-second images in V and I only; no FOV file necessary,
       e.g., "BURN FF Lyr 12:00:00 +23:34:45".
    IMAGE  target_name  filter_mag_or_sec_string  RA  Dec  ::  arbitrary imaging,
       either in magnitudes for the current instrument, or in explicity exposure times (may be mixed on
       a single line):
       *** Magnitudes syntax: "IMAGE New target V=12 B=12.5(2) 12:00:00 +23:34:45" to image New target in
            V filter (once) at targeted mag 12, and B filter twice at targeted mag 12.5.
       *** Exposure syntax:  "IMAGE New target V=120s B=240s(2) 12:00:00 +23:34:45" to image New target in
           V filter (once) at 120 seconds, and B filter twice at 240 seconds (exposure times are NOT
           limited, so be careful!)
       All text between "IMAGE" and first word including a "=" character will make up the target name.
    COLOR  target_name  multiplier  RA  Dec  :: for MP color imaging (mp_color.py),
        e.g., "COLOR MP_1626 1.1x 21:55:08 +24:24:45. 'x' in multipler is optional but recommended.

    ----- Legal directives:
    PLAN  plan_id  ::  starts a plan section and names it.
    ;   comment_text  :: semicolon at beginning of cell makes cell a comment only.
    AFINTERVAL nnn  ::  autofocus interval in minutes
    SETS  nn  ::  number of times to repeat all targets, autofocuses, chills, etc
    AUTOFOCUS       :: force autofocus
    CHILL -nn  :: chill the cooler to -nn deg C
    QUITAT nn:nn  ::  quit plan at nn:nn UTC
    WAITUNTIL nn:nn  ::  wait to start plan (first Set) until nn:nn UTC
    SKIPFILTER filter_name  ::  skip filter for following targets; omit filter_name to restore all filters.
    SHUTDOWN  ::  perform ACP shutdown of camera and park scope
    CHAIN plan_id  ::  chain to next plan
    BURN target_id RA Dec  ::  shorthand for IMAGE target_id V=240sec(1) I=240sec(1) RA Dec
    IMAGE target_id exp_specs RA Dec  ::  take images of target at RA, Dec; exp_specs define the
       filters and exposures, e.g., V=12.8 R=120sec(2) I=11(3) where 12.8 is a magnitude, 120sec
       is 120 seconds, and (2) specifies 2 exposures (and resulting images).
    COLOR target_id multipler RA Dec :: take MP color sequence defined by MP_COLOR_V14.
    fov_name  ::  if line begins with none of the above, it's a FOV name and takes filters, exposures,
        RA, and Dec from its FOV file.

    """
    df = pd.read_excel(excel_path, header=None).dropna(axis=0, how='all').dropna(axis=1, how='all')
    nrow = len(df)
    ncol = len(df.columns)
    parsed_list = []  # nested list, one element per ACP plan.
    this_plan_id = ''
    an_date_string = str(df.iloc[0, 0]).strip()
    if int(EARLIEST_AN_DATE) < int(an_date_string) < int(LATEST_AN_DATE):
        an = Astronight(an_date_string, site_name)
    else:
        print('>>>>> STOPPING: an_date_string '" + an_date_string + "
              ""' SEEMS UNREASONABLE (update LATEST_AN_DATE?).')
        return

    plan_list = []
    this_plan = None
    macro_dict = dict()
    macro_field_keys = ['^' + str(i + 1) for i in range(9)]
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

                # Add MACRO directive that stores text in a dict for later use; then continue loop:
                if cell_str_lower.startswith('macro '):
                    _, macro_key, macro_command = tuple(cell_str_as_read.split(maxsplit=2))
                    macro_dict[macro_key] = macro_command
                    continue

                # If first word of command is in macro dict, substitute expanded macro for command.
                words = cell_str_as_read.split()
                macro_command = macro_dict.get(words[0], None)
                if macro_command is not None:
                    macro_misused = False
                    insert_strings = words[1:]
                    for i_key, macro_field_key in enumerate(macro_field_keys):
                        iloc = macro_command.find(macro_field_key)
                        if iloc >= 0:
                            if i_key < len(insert_strings):
                                insert_string = insert_strings[i_key]
                            else:
                                macro_misused = True
                                insert_string = '???'
                            macro_command = macro_command.replace(macro_field_key, insert_string)
                    if macro_misused:
                        print(' >>>>> ERROR: Macro misused in cell \'' + cell_str_as_read + '\'')
                    cell_str_as_read = macro_command
                    cell_str_lower = cell_str_as_read.lower()  # refresh this variable.

                # Handle any comment after first semi-colon:
                split_str = cell_str_as_read.split(';', maxsplit=1)
                command = split_str[0].strip()
                if len(split_str) > 1:
                    comment = split_str[1].rstrip()
                else:
                    comment = None

                # Determine action type and add action to directive_list:
                if cell_str_lower.startswith('plan'):
                    if this_plan is not None:
                        plan_list.append(this_plan)  # save previous plan, if any
                    this_plan_id = an_date_string + '_' + command[len('plan'):].strip()
                    this_plan = Plan(this_plan_id, comment)
                elif cell_str_lower.startswith('sets'):
                    set_count = command[len('sets'):].strip()
                    this_plan.directives.append(Directive('sets', {'count': int(set_count)}))
                elif cell_str_as_read.startswith(';'):
                    this_plan.directives.append(Directive('comment', {'text': comment}))
                elif cell_str_lower.startswith('afinterval'):
                    minutes = command[len('afinterval'):].strip()
                    this_plan.directives.append(Directive('afinterval', {'minutes': int(minutes)}))
                elif cell_str_lower.startswith('autofocus'):
                    this_plan.directives.append(Directive('autofocus', {}))
                elif cell_str_lower.startswith('chill'):
                    tempC = command[len('chill'):].strip()
                    this_plan.directives.append(Directive('chill', {'tempC': float(tempC)}))
                elif cell_str_lower.startswith('quitat'):
                    hhmm_utc = command[len('quitat'):].strip().replace(':', '')
                    this_plan.directives.append(Directive('quitat', {'utc': hhmm_utc}))
                elif cell_str_lower.startswith('waituntil'):
                    value = command[len('waituntil'):].strip().replace(':', '')
                    spec_dict = {'sun_degrees': None, 'utc': None}  # overwrite one of these, just below.
                    if float(value) < 0:
                        spec_dict['sun_degrees'] = float(value)
                    else:
                        spec_dict['utc'] = value
                    this_plan.directives.append(Directive('waituntil', spec_dict))
                elif cell_str_lower.startswith('skipfilter'):
                    value = command[len('skipfilter'):].strip()
                    if cell_str_lower.startswith('skipfilters'):
                        value = command[len('skipfilters'):].strip()  # deprecated SKIPFILTERS (plural) case
                    skipfilter_list = [item.strip() for item in value.split()]
                    this_plan.directives.append(Directive('skipfilter', {'filters': skipfilter_list}))
                elif cell_str_lower.startswith('shutdown'):
                    this_plan.directives.append(Directive('shutdown', {}))
                elif cell_str_lower.startswith('chain'):
                    next_plan_filename = 'plan_' + an_date_string + '_' + \
                                         command[len('chain'):].strip().upper()
                    if not next_plan_filename.endswith('.txt'):
                        next_plan_filename += '.txt'
                    this_plan.directives.append(Directive('chain', {'filename': next_plan_filename}))
                elif cell_str_lower.startswith('burn'):
                    value = command[len('burn'):].strip()
                    this_fov_name, ra_string, dec_string = extract_ra_dec(value)
                    # this_fov_name, ra_string, dec_string = tuple(value.rsplit(maxsplit=2))
                    this_plan.directives.append(Directive('burn', {'fov_name': this_fov_name.strip(),
                                                                   'ra': ra_string.strip(),
                                                                   'dec': dec_string.strip(),
                                                                   'force_autoguide': True}))
                elif cell_str_lower.startswith('stare'):
                    value = command[len('stare'):].strip()
                    repeats, this_fov_name = tuple(value.split(maxsplit=1))
                    this_plan.directives.append(Directive('stare', {'fov_name': this_fov_name.strip(),
                                                                    'repeat_count': int(repeats),
                                                                    'force_autoguide': True}))
                elif cell_str_lower.startswith('image'):
                    value = command[len('image'):].strip()
                    force_autoguide, value = get_and_remove_option(value, FORCE_AUTOGUIDE_TOKEN)
                    subvalue, ra_string, dec_string = extract_ra_dec(value)
                    # subvalue, ra_string, dec_string = tuple(value.rsplit(maxsplit=2))
                    filter_entries = []
                    target_name = "WARNING: NO TARGET NAME"
                    while True:
                        if len(subvalue) <= 0:
                            print(">>>>> WARNING: No target name for command '" +
                                  cell_str_as_read + "'.")
                            break
                        if len(subvalue.split()) == 1:
                            target_name = subvalue
                            break
                        subsubvalue, item = subvalue.rsplit(maxsplit=1)
                        is_filter_entry = '=' in item
                        if is_filter_entry:
                            filter_entries.append(item.strip())
                        else:
                            target_name = subvalue
                            break
                        subvalue = subsubvalue
                    filter_entries.reverse()
                    if len(filter_entries) >= 1:
                        this_plan.directives.append(Directive('image',
                                                              {'target_name': target_name,
                                                               'filter_entries': filter_entries,
                                                               'ra': ra_string,
                                                               'dec': dec_string,
                                                               'force_autoguide': force_autoguide}))
                elif cell_str_lower.startswith('color'):
                    value = command[len('image'):].strip()
                    force_autoguide, value = get_and_remove_option(value, FORCE_AUTOGUIDE_TOKEN)  # ignored.
                    subvalue, ra_string, dec_string = extract_ra_dec(value)
                    target_name, multiplier_string = tuple(subvalue.rsplit(maxsplit=1))
                    # target_name, multiplier_string, ra_string, dec_string =
                    # tuple(value.rsplit(maxsplit=3))
                    multiplier_string = multiplier_string.lower().split('x')[0]
                    multiplier = float(multiplier_string)
                    entries = tuple([(filt, multiplier * exp14, repeats)
                                     for (filt, exp14, repeats) in MP_COLOR_V14])
                    this_plan.directives.append(Directive('color',
                                                          {'target_name': target_name,
                                                           'entries': entries,
                                                           'multiplier_string': multiplier_string,
                                                           'ra': ra_string,
                                                           'dec': dec_string,
                                                           'force_autoguide': True}))  # always autoguide.
                else:
                    # Anything else we treat as a fov_name:
                    value = command  # use the whole command string (before comment); no directive string.
                    force_autoguide, value = get_and_remove_option(value, FORCE_AUTOGUIDE_TOKEN)
                    fov_name = value.strip()
                    if len(fov_name) >= 2:
                        this_plan.directives.append(Directive('fov', {'fov_name': fov_name,
                                                                      'force_autoguide': force_autoguide}))

    plan_list.append(this_plan)  # Ensure we save the last plan.
    return plan_list, an


def get_and_remove_option(string, option):
    """ From a value string (e.g., 'IMAGE MP_191 AG+ Clear=200sec(1) 12:34:45 -06:34:21') and an
    option string (e.g., 'AG+), determine whether option is in value, and return value string with
    all instances of option token (space-bounded) removed, for further processing.
    Used in parse_excel().
    :param string: value string of directive. [string]
    :param option: option string to locate. [string]
    :return: (flag value, value string with option removed. [2-tuple of (boolean, string)]
    """
    p = ' ' + string + ' '
    pu = p.upper()
    pu_option = ' ' + option.upper() + ' '
    flag = (pu.find(pu_option) >= 0)
    while True:
        i = pu.find(pu_option)
        if i == -1:
            break
        pu = pu[:i + 1] + pu[i + len(pu_option):]
        p = p[:i + 1] + p[i + len(pu_option):]
    return flag, p.strip()


def reorder_directives(plan_list):
    """
    Puts directives within each Plan object in the desired order, returns the updated plan list.
    :param plan_list: the plan list whose directives are to be reordered [list of Plan objects].
    :return: the plan list with reordered directives [list of Plan objects].
    """
    # Directives within each sublist retain user's given order.
    ideal_directive_ordering = [['quitat'],
                                ['afinterval'],
                                ['sets'],
                                ['waituntil', 'chill', 'stare', 'fov', 'burn',
                                 'image', 'color', 'autofocus', 'comment', 'skipfilter'],
                                ['shutdown'],
                                ['chain']]
    for plan in plan_list:
        reordered_directive_list = []
        for directive_order_sublist in ideal_directive_ordering:
            for i_directive in range(len(plan.directives)):
                this_directive = plan.directives[i_directive]
                if this_directive.type.lower() in directive_order_sublist:
                    reordered_directive_list.append(this_directive)
        plan.directives = reordered_directive_list
        num_omitted = len(plan.directives) - len(reordered_directive_list)
        if num_omitted > 0:
            print('>>>>> WARNING: ' + str(num_omitted) + ' actions in plan ' + plan.plan_id +
                  'were omitted during ordering.')
    # return plan_list


def make_events(plan_list, instrument, fov_dict, an, exp_time_factor):
    """ Translate user's directives into executable events (to be repeated if more than one set).
        For simplicity, handle all directives, even if not enough plan time to complete them all (common).
        Compute event durations here, but postpone creation of full plan timeline to later function.
    :param plan_list:
    :param instrument:
    :param fov_dict:
    :param an:
    :param exp_time_factor:
    :return: [nothing--it modifies plan_list in place].
    """

    for plan in plan_list:
        skipfilter_list = []  # default

        # For each directive: make event and add it to plan's event list:
        for directive in plan.directives:

            if directive.type == 'waituntil':  # NB there may be >1 waituntil, but only 1 active quitat.
                if directive.spec['sun_degrees'] is not None:
                    sun_degrees = directive.spec['sun_degrees']
                    site_obs = ephem.Observer()
                    site_obs.lat, site_obs.lon = str(an.site.latitude), str(an.site.longitude)
                    site_obs.elevation = an.site.elevation
                    sun = ephem.Sun(site_obs)
                    site_obs.horizon = str(sun_degrees)
                    utc_end = site_obs.previous_setting(sun, an.local_middark_utc) \
                        .datetime().replace(tzinfo=timezone.utc)
                    this_summary_text = 'WAITUNTIL sun reaches ' + \
                                        '{0:g}'.format(sun_degrees) + u'\N{DEGREE SIGN}' + ' alt'
                    this_acp_entry = ['#WAITUNTIL 1, ' + '{0:g}'.format(sun_degrees) + ' ; deg sun alt']
                else:
                    hhmm = ('0' + directive.spec['utc'])[-4:]
                    utc_end = an.datetime_utc_from_hhmm(hhmm)
                    formatted_time = '{:%m/%d/%Y %H:%M}'.format(utc_end)
                    this_summary_text = 'WAITUNTIL ' + hhmm + ' utc'
                    this_acp_entry = ['#WAITUNTIL 1, ' + formatted_time + ' ; utc']
                this_event = Event('waituntil', this_summary_text, this_acp_entry)
                this_event.utc_end = utc_end
                plan.events.append(this_event)

            elif directive.type == 'chill':
                this_summary_text = 'CHILL  ' + '{0:g}'.format(directive.spec['tempC'])
                this_acp_entry = ['#CHILL  ' + '{0:g}'.format(directive.spec['tempC'])]
                this_event = Event('chill', this_summary_text, this_acp_entry, CHILL_DURATION)
                plan.events.append(this_event)

            elif directive.type == 'stare':
                n_repeats = directive.spec['repeat_count']
                fov_name = directive.spec['fov_name']
                this_summary_text = 'Stare ' + str(n_repeats) + ' repeats at ' + fov_name
                if directive.spec['force_autoguide'] is True:
                    this_summary_text += ' AG+'
                exp_data = make_fov_exposure_data(fov_name, an, fov_dict, instrument,
                                                  exp_time_factor=exp_time_factor,
                                                  skipfilter_list=skipfilter_list,
                                                  force_autoguide=directive.spec['force_autoguide'])
                if exp_data is None:
                    return  # fail
                filters, counts, exp_times, target_overhead, repeat_duration = exp_data
                event_duration = target_overhead + n_repeats * repeat_duration
                duration_comment = str(round(repeat_duration / 60.0, 1)) + ' min/repeat --> ' + \
                                   str(round(event_duration / 60.0, 1)) + ' min (nominal)'
                this_fov = fov_dict[fov_name]
                this_acp_entry = [';', '#REPEAT ' + str(n_repeats) + ';',
                                  '#DITHER 0 ;',
                                  '#FILTER ' + ','.join(filters) + ' ;',
                                  '#BINNING ' + ','.join(len(filters) * ['1']) + ' ;',
                                  '#COUNT ' + ','.join([str(c) for c in counts]) + ' ;',
                                  '#INTERVAL ' + ','.join([str(e).split('.0')[0]
                                                           for e in exp_times]) +
                                  ' ; ' + duration_comment,
                                  ';----' + this_fov.acp_comments, fov_name + '\t' +
                                  ra_as_hours(this_fov.ra) + '\t' + dec_as_hex(this_fov.dec)]
                if directive.spec['force_autoguide'] is True:
                    this_acp_entry.insert(1, '#AUTOGUIDE  ;  Automatic for stare target')
                duration_dict = {'target_overhead': target_overhead,
                                 'repeat_count': n_repeats,
                                 'counts': counts,
                                 'exp_times': exp_times}
                this_event = Event('stare', this_summary_text, this_acp_entry,
                                   event_duration, duration_dict,
                                   ra=ra_as_hours(this_fov.ra), dec=dec_as_hex(this_fov.dec))
                this_event.target_name = fov_name
                plan.events.append(this_event)

            elif directive.type == 'fov':
                fov_name = directive.spec['fov_name']
                this_summary_text = fov_name
                if directive.spec['force_autoguide'] is True:
                    this_summary_text += ' AG+'
                exp_data = make_fov_exposure_data(fov_name, an, fov_dict, instrument,
                                                  exp_time_factor=exp_time_factor,
                                                  skipfilter_list=skipfilter_list,
                                                  force_autoguide=directive.spec['force_autoguide'])
                if exp_data is None:
                    return  # fail
                filters, counts, exp_times, target_overhead, repeat_duration = exp_data
                event_duration = target_overhead + 1 * repeat_duration
                duration_comment = ' --> ' + str(round(event_duration / 60.0, 1)) + ' min'
                this_fov = fov_dict[fov_name]
                this_acp_entry = [';', '#DITHER 0 ;',
                                  '#FILTER ' + ','.join(filters) + ' ;',
                                  '#BINNING ' + ','.join(len(filters) * ['1']) + ' ;',
                                  '#COUNT ' + ','.join([str(c) for c in counts]) + ' ;',
                                  '#INTERVAL ' + ','.join([str(e).split('.0')[0]
                                                           for e in exp_times]) +
                                  ' ; ' + duration_comment,
                                  ';----' + this_fov.acp_comments, fov_name + '\t' +
                                  ra_as_hours(this_fov.ra) + '\t' + dec_as_hex(this_fov.dec)]
                if directive.spec['force_autoguide'] is True:
                    this_acp_entry.insert(1, '#AUTOGUIDE  ;  Forced')
                duration_dict = {'target_overhead': target_overhead,
                                 'repeat_count': 1,
                                 'counts': counts,
                                 'exp_times': exp_times}
                this_event = Event('fov', this_summary_text, this_acp_entry,
                                   event_duration, duration_dict,
                                   ra=ra_as_hours(this_fov.ra), dec=dec_as_hex(this_fov.dec))
                this_event.target_name = fov_name
                plan.events.append(this_event)

            elif directive.type == 'burn':
                future_fov_name = directive.spec['fov_name']
                ra = directive.spec['ra']
                dec = directive.spec['dec']
                this_summary_text = 'BURN ' + future_fov_name + '  ' + ra + '  ' + dec
                if directive.spec['force_autoguide'] is True:
                    this_summary_text += ' AG+'
                this_acp_entry = [';', '#DITHER 0 ;', '#FILTER V,I ;', '#BINNING 1,1 ;',
                                  '#COUNT 1,1 ;', '#INTERVAL ' +
                                  str(BURN_EXPOSURE) + ',' + str(BURN_EXPOSURE) +
                                  ' ;----> BURN for new FOV file.',
                                  future_fov_name + '\t' + ra + '\t' + dec + ' ;']
                if directive.spec['force_autoguide'] is True:
                    this_acp_entry.insert(1, '#AUTOGUIDE  ;  Automatic for burn target')
                target_overhead, repeat_duration = tabulate_target_durations(
                    filters=['V', 'I'], counts=[1, 1],
                    exp_times=[BURN_EXPOSURE, BURN_EXPOSURE], force_autoguide=True)
                event_duration = target_overhead + 1 * repeat_duration
                duration_dict = {'target_overhead': event_duration - 2 * BURN_EXPOSURE,
                                 'repeat_count': 1,
                                 'counts': [1, 1],
                                 'exp_times': 2 * [BURN_EXPOSURE]}
                this_event = Event('burn', this_summary_text, this_acp_entry,
                                   event_duration, duration_dict,
                                   ra=ra, dec=dec)
                this_event.target_name = future_fov_name
                plan.events.append(this_event)

            elif directive.type == 'image':
                target_name = directive.spec['target_name']
                filter_entries = directive.spec['filter_entries']
                ra = directive.spec['ra']
                dec = directive.spec['dec']
                filters, counts, exp_times, target_overhead, repeat_duration = \
                    make_image_exposure_data(filter_entries, instrument, exp_time_factor=exp_time_factor,
                                             force_autoguide=directive.spec['force_autoguide'])
                event_duration = target_overhead + 1 * repeat_duration
                this_summary_text = 'Image ' + target_name +\
                                    ' ' + ''.join([f + '=' + '{0:g}'.format(e) + 's(' + str(c) + ')'
                                                   for (f, e, c) in zip(filters, exp_times, counts)]) +\
                                    ' ' + ra + ' ' + dec
                if directive.spec['force_autoguide'] is True:
                    this_summary_text += ' AG+'
                duration_comment = ' --> ' + str(round(event_duration / 60.0, 1)) + ' min'
                this_acp_entry = [';', '#DITHER 0 ;',
                                  '#FILTER ' + ','.join(filters) + ' ;',
                                  '#BINNING ' + ','.join(len(filters) * ['1']) + ' ;',
                                  '#COUNT ' + ','.join([str(c) for c in counts]) + ' ;',
                                  '#INTERVAL ' + ','.join([str(e).split('.0')[0]
                                                           for e in exp_times]) +
                                  ' ; ' + duration_comment,
                                  ';---- from IMAGE directive -----',
                                  target_name + '\t' + ra + '\t' + dec]
                if directive.spec['force_autoguide'] is True:
                    this_acp_entry.insert(1, '#AUTOGUIDE  ;  Forced')
                duration_dict = {'target_overhead': target_overhead,
                                 'repeat_count': 1,
                                 'counts': counts,
                                 'exp_times': exp_times}
                this_event = Event('image', this_summary_text, this_acp_entry,
                                   event_duration, duration_dict,
                                   ra=ra, dec=dec)
                this_event.target_name = target_name
                plan.events.append(this_event)

            elif directive.type == 'color':
                target_name = directive.spec['target_name']
                entries = directive.spec['entries']
                ra = directive.spec['ra']
                dec = directive.spec['dec']
                filters, counts, exp_times, target_overhead, repeat_duration = \
                    make_color_exposure_data(entries, force_autoguide=True)
                event_duration = target_overhead + 1 * repeat_duration
                this_summary_text = 'Color ' + target_name +\
                                    '  ' + '.'.join(filters) +\
                                    '  ' + directive.spec['multiplier_string'] +\
                                    'x  ' + '{0:.1f}'.format(event_duration / 60.0) + ' min.' +\
                                    '  ' + ra + '  ' + dec
                if directive.spec['force_autoguide'] is True:
                    this_summary_text += ' AG+'
                duration_comment = ' --> ' + str(round(event_duration / 60.0, 1)) + ' min'
                this_acp_entry = [';', '#DITHER 0 ;',
                                  '#FILTER ' + ','.join(filters) + ' ;',
                                  '#BINNING ' + ','.join(len(filters) * ['1']) + ' ;',
                                  '#COUNT ' + ','.join([str(c) for c in counts]) + ' ;',
                                  '#INTERVAL ' + ','.join([str(e).split('.0')[0]
                                                           for e in exp_times]) +
                                  ' ; ' + duration_comment,
                                  ';---- from COLOR directive -----',
                                  target_name + '\t' + ra + '\t' + dec]
                if directive.spec['force_autoguide'] is True:
                    this_acp_entry.insert(1, '#AUTOGUIDE  ;  Forced')
                duration_dict = {'target_overhead': target_overhead,
                                 'repeat_count': 1,
                                 'counts': counts,
                                 'exp_times': exp_times}
                this_event = Event('color', this_summary_text, this_acp_entry,
                                   event_duration, duration_dict,
                                   ra=ra, dec=dec)
                this_event.target_name = target_name
                plan.events.append(this_event)

            elif directive.type == 'autofocus':
                this_summary_text = 'AUTOFOCUS'
                this_acp_entry = [';', '#AUTOFOCUS']
                event_duration = AUTOFOCUS_DURATION
                this_event = Event('autofocus', this_summary_text, this_acp_entry, event_duration)
                plan.events.append(this_event)

            elif directive.type == 'comment':
                comment_text = directive.spec['text']
                this_summary_text = ';' + comment_text
                this_acp_entry = [';' + comment_text]
                event_duration = 0
                this_event = Event('comment', this_summary_text, this_acp_entry, event_duration)
                plan.events.append(this_event)

            elif directive.type == 'skipfilter':
                new_skipfilter_list = directive.spec['filters']
                if len(new_skipfilter_list) == 0:
                    skipfilter_list_text = 'none'
                else:
                    skipfilter_list_text = ' '.join(new_skipfilter_list)
                this_summary_text = 'SKIPFILTER ' + skipfilter_list_text
                this_acp_entry = [';', '; (skipfilter: ' + skipfilter_list_text + ')']
                skipfilter_list = new_skipfilter_list  # changing this state variable
                event_duration = 0
                this_event = Event('skipfilter', this_summary_text, this_acp_entry, event_duration)
                plan.events.append(this_event)

            elif directive.type == 'shutdown':
                this_summary_text = 'SHUTDOWN'
                this_acp_entry = [';', '#SHUTDOWN']
                event_duration = SHUTDOWN_DURATION
                this_event = Event('shutdown', this_summary_text, this_acp_entry, event_duration)
                plan.events.append(this_event)

            elif directive.type == 'quitat':
                plan.utc_quitat = an.datetime_utc_from_hhmm(directive.spec['utc'])

            elif directive.type == 'afinterval':
                plan.afinterval = float(directive.spec['minutes'])

            elif directive.type == 'sets':
                plan.sets_requested = int(directive.spec['count'])

            elif directive.type == 'chain':
                plan.chain_destination = directive.spec['filename']

            else:
                print(">>>>> ERROR: in plan", plan.plan_id,
                      ', directive', directive.type, 'not understood.')


def make_timeline(plan_list, an, earliest_hhmm):
    # TODO: SHUTDOWN needs repair, to make it function & stop (1) in mid-plan, (2) even with SETS.
    # For now, SHUTDOWN must go in it's own (last) plan.

    # Initialize times & intervals to state before first plan:
    utc_running = None
    if earliest_hhmm is not None:
        utc_running = an.datetime_utc_from_hhmm(earliest_hhmm)
    else:
        utc_running = an.datetime_utc_from_hhmm('0000') + timedelta(hours=AN_START_REL_UTC_0000)
        if utc_running > an.ts_dark.start:
            utc_running -= timedelta(hours=24)
    utc_most_recent_autofocus = utc_running - timedelta(days=1000)  # keep python happy with a prev value.
    shutdown_performed = False

    for plan in plan_list:
        plan.utc_start = utc_running
        no_plan_exposures_yet_encountered = True

        for i_set in range(1, plan.sets_requested + 1):  # i_set = 1 to sets_requested, inclusive.
            skipfilter_list = []  # reset at beginning of set execution.
            for event in plan.events:
                # First, do autofocus if AFINTERVAL since latest autofocus has passed or at plan startup:
                # TODO: rewrite (move?) this, so that long Stare & Image events can have > 1 autofocus.
                if plan.afinterval is not None:
                    minutes_since_last_autofocus = \
                        (utc_running - utc_most_recent_autofocus).total_seconds() / 60.0
                    if event.type in ['burn', 'stare', 'fov', 'image', 'color']:
                        if minutes_since_last_autofocus > plan.afinterval or \
                            no_plan_exposures_yet_encountered:
                            # Perform AFINTERVAL autofocus:
                            utc_running += timedelta(seconds=AUTOFOCUS_DURATION)
                            utc_most_recent_autofocus = utc_running
                            plan.afinterval_autofocus_count += 1
                            if plan.sets_requested == 1:
                                event.summary_text += ' (af)'
                    if plan.quitat_reached_at(utc_running):
                        break  # if quitat time reached during afinterval autofocus, do not run event.

                utc_start_event = utc_running

                # Store event's actual end time (incl quitat if active):
                if event.type == 'waituntil':
                    # WAITUNTIL only works in first set (set 1):
                    if i_set == 1:
                        if plan.utc_quitat is not None:
                            utc_end_event = min(event.utc_end, plan.utc_quitat)  # not later than QUITAT.
                        else:
                            utc_end_event = event.utc_end
                        # But definitely not before utc_running (time goes not backward):
                        utc_end_event_actual = max(utc_end_event, utc_running)
                elif event.type in ['comment', 'skipfilter']:
                    utc_end_event_actual = utc_start_event  # zero duration
                elif event.type in ['chill', 'shutdown', 'autofocus']:
                    utc_end_event_actual = utc_start_event + timedelta(seconds=event.duration_total)
                elif event.type in ['burn', 'stare', 'fov', 'image', 'color']:
                    actual_duration, event_autofocus_count, utc_most_recent_autofocus = \
                        event.calc_actual_duration(utc_start_event, plan.utc_quitat,
                                                   plan.afinterval, utc_most_recent_autofocus)
                    if event_autofocus_count >= 1:
                        event.summary_text += ' (' + str(event_autofocus_count) + ' af)'
                        plan.afinterval_autofocus_count += event_autofocus_count
                    utc_end_event_actual = utc_start_event + timedelta(seconds=actual_duration)
                    no_plan_exposures_yet_encountered = False
                else:
                    print('make_timeline() doesn\'t recognize event type"" ' + event.type)

                # Store event's summary display time (hhmm on summary line, usually for set 1):
                if i_set == 1:
                    event.utc_summary_display = utc_start_event

                # Update event's minimum altitude (all sets, target event types only):
                if event.type in ['burn', 'stare', 'fov', 'image', 'color']:
                    this_lower_alt = event.calc_lower_altitude(an, utc_start_event, utc_end_event_actual)
                    if event.min_altitude is None:
                        event.min_altitude = this_lower_alt
                    else:
                        event.min_altitude = min(event.min_altitude, this_lower_alt)

                # Store event's status:
                if event.type == 'chain':
                    event.status = 'CHAIN'
                elif plan.quitat_reached_at(utc_running):
                    event.status = 'QUITAT'
                elif utc_start_event < an.ts_dark.start or utc_end_event_actual > an.ts_dark.end:
                    event.status = 'LIGHT'
                elif event.type in ['burn', 'stare', 'fov', 'image', 'color', 'autofocus', 'chill']:
                    event.status = str(i_set)  # default
                elif event.type in ['shutdown', 'waituntil']:
                    event.status = 'ok'
                else:
                    event.status = ''

                # Finally, update master clock at end of this event:
                utc_running = utc_end_event_actual

                # For SHUTDOWN, signal end of entire run:
                if event.type == 'shutdown':
                    shutdown_performed = True
                    break  # out of event loop (to next set).

                # Stop events if shutdown run or quitat reached:
                if plan.quitat_reached_at(utc_running) or shutdown_performed:
                    break  # out of event loop (to next set, which will also stop)

            # Quit set if shutdown run or quitat reached:
            if plan.quitat_reached_at(utc_running) or shutdown_performed:
                break  # out of set loop (to next plan)

            plan.sets_completed = i_set

        # Finish any end-of-plan business (incl saving statistics):
        plan.utc_end = utc_running

        # Quit plan if shutdown run or quitat reached:
        if shutdown_performed:
            break  # out of plan loop to end of timeline.

    # Finish end-of-night business (or could go outside this function, instead):
    pass


def make_acp_plan_files(plan_list, an, output_directory, exp_time_factor):
    # First, delete old ACP plan files:
    filenames = os.listdir(output_directory)
    for filename in filenames:
        if filename.startswith("plan_") and filename.endswith(".txt"):
            fullpath = os.path.join(output_directory, filename)
            os.remove(fullpath)

    # Then, make an ACP-format plan file for each plan:
    for plan in plan_list:
        if plan.plan_comment is None:
            plan_comment = ''
        else:
            plan_comment = '; ' + plan.plan_comment
        # noinspection PyListCreation
        plan_acp_lines = ['; ACP PLAN ' + plan.plan_id + plan_comment,
                          ';     as generated by photrix at ' +
                          '{:%Y-%m-%d %H:%M  UTC}'.format(datetime.now(timezone.utc)),
                          ';     using exposure time factor = ' + '{:5.3f}'.format(exp_time_factor)]
        plan_acp_lines.append(an.acp_header_string())

        # Add SETS ACP directive if one exists:
        if plan.sets_requested > 1:
            plan_acp_lines.extend([';', '#SETS ' + str(int(plan.sets_requested))])

        # Add QUITAT ACP directive if one exists:
        if plan.utc_quitat is not None:
            formatted_time = '{:%m/%d/%Y %H:%M}'.format(plan.utc_quitat)
            plan_acp_lines.extend([';', '#QUITAT ' + formatted_time + ' ; utc'])

        # Add AFINTERVAL ACP directive if one exists:
        if plan.afinterval is not None:
            if plan.afinterval > 0:
                plan_acp_lines.append('#AFINTERVAL ' + '{0:g}'.format(plan.afinterval))

        if plan.utc_quitat is not None or plan.afinterval is not None:
            plan_acp_lines.append(';')

        # Add event lines:
        for event in plan.events:
            plan_acp_lines.extend(event.acp_lines)

        # Add CHAIN ACP directive if one exists:
        if plan.chain_destination is not None:
            plan_acp_lines.extend([';', '#CHAIN ' + plan.chain_destination])

        # Write this ACP plan file:
        filename = 'plan_' + plan.plan_id + '.txt'
        output_fullpath = os.path.join(output_directory, filename)
        print('PRINT plan ' + plan.plan_id)
        with open(output_fullpath, 'w') as this_file:
            this_file.write('\n'.join(plan_acp_lines))


def make_summary_file(plan_list, fov_dict, an, output_directory, exp_time_factor):
    # First, delete old summary files:
    filenames = os.listdir(output_directory)
    for filename in filenames:
        if filename.startswith("Summary_") and filename.endswith(".txt"):
            fullpath = os.path.join(output_directory, filename)
            os.remove(fullpath)

    # Unpack summary_lines:
    an_year = int(an.an_date_string[0:4])
    an_month = int(an.an_date_string[4:6])
    an_day = int(an.an_date_string[6:8])
    day_of_week = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'] \
        [datetime(an_year, an_month, an_day).weekday()]
    header_lines = ['SUMMARY for AN' + an.an_date_string + '   ' + day_of_week.upper(),
                    '     as generated by photrix at ' +
                    '{:%Y-%m-%d %H:%M  UTC}'.format(datetime.now(timezone.utc)),
                    '     using exposure time factor = ' + '{:5.3f}'.format(exp_time_factor) +
                    '     min.alt = ' + '{:.1f}'.format(an.site.min_altitude) + u'\N{DEGREE SIGN}',
                    an.acp_header_string(), '\n']
    moon_is_a_factor = an.moon_phase > MOON_PHASE_NO_FACTOR  # for this astronight
    radec_dict = dict()  # collector to find and warn against Image events of same target w/ diff RA,Dec.

    # Local function:
    def make_summary_line(status_text, hhmm_text, utc_day_indicator, min_altitude, summary_text):
        if status_text is None:
            status_text = ''
        if hhmm_text is None:
            hhmm_text = 4 * ' '
        if utc_day_indicator is None:
            utc_day_indicator = ' '
        if min_altitude is not None:
            altitude_text = str(int(round(min_altitude)))
        else:
            altitude_text = '  '
        if status_text == '1':
            status_text = 'ok'
        return ' '.join([status_text.rjust(8), hhmm_text + utc_day_indicator,
                         altitude_text, summary_text])

    # Construct summary_lines for every event:
    for i_plan, plan in enumerate(plan_list):
        # Add lines to top of plan summary:
        plan.summary_pre_lines.append(
            make_summary_line(None, None, None, None, 60 * '-'))
        hhmm_start = hhmm_from_datetime_utc(plan.utc_start)
        hhmm_end = hhmm_from_datetime_utc(plan.utc_end)
        if i_plan == 0:
            display_start = 'dusk to '
        else:
            display_start = hhmm_start + '-'
        plan.summary_pre_lines.append(
            make_summary_line(None, None, None, None,
                              'Begin Plan ' + plan.plan_id + ' :: ' + display_start + hhmm_end + ' utc'))
        if plan.plan_comment is not None:
            if len(plan.plan_comment.strip()) > 0:
                plan.summary_pre_lines.append(
                    make_summary_line(None, None, None, None, plan.plan_comment))
        if i_plan > 0:
            plan.summary_pre_lines.append(
                make_summary_line(None, hhmm_start, None, None, 'Plan entered.'))
        if plan.sets_requested > 1:
            plan.summary_pre_lines.append(
                make_summary_line(None, None, None, None,
                                  'SETS ' + '{0:g}'.format(plan.sets_requested)))
        if plan.utc_quitat is not None:
            plan.summary_pre_lines.append(
                make_summary_line(None, None, None, None,
                                  'QUITAT ' + hhmm_from_datetime_utc(plan.utc_quitat) + ' utc'))
        if plan.afinterval is not None:
            plan.summary_pre_lines.append(
                make_summary_line(None, None, None, None,
                                  'AFINTERVAL ' + '{0:g}'.format(plan.afinterval)))

        # Add lines to end of plan summary:
        if plan.chain_destination is not None:
            plan.summary_post_lines.append(
                make_summary_line('CHAIN', hhmm_end, ' ', None,
                                  'Chain to \'' + plan.chain_destination + '\''))
        if plan.afinterval is not None:
            plan.summary_post_lines.append(
                make_summary_line(None, None, None, None,
                                  str(plan.afinterval_autofocus_count) + ' AFINTERVAL autofocuses done.'))
        plan.summary_post_lines.append('\n')
        for event in plan.events:
            # Add warning line if time wasted by waiting to start this plan (prev plan ended early).
            if event.type == 'waituntil':
                if i_plan > 0:
                    gap_minutes = (event.utc_end - plan.utc_start).total_seconds() / 60.0
                    if gap_minutes > 1.0:
                        plan.summary_pre_lines.append(
                            '      >>>>>>>>>> WARNING: WAITUNTIL gap = ' +
                            str(int(gap_minutes)) + ' minutes.')
            # Construct main summary text line for this event, write into its Event object:
            if event.type in ['waituntil', 'comment', 'skipfilter']:
                hhmm_text, utc_day_indicator = '    ', ' '
            else:
                if event.utc_summary_display is not None:
                    hhmm_text = hhmm_from_datetime_utc(event.utc_summary_display)
                    if event.utc_summary_display < an.datetime_utc_from_hhmm('0000'):
                        utc_day_indicator = '-'
                    elif event.utc_summary_display > an.datetime_utc_from_hhmm('0000') + timedelta(days=1):
                        utc_day_indicator = '+'
                    else:
                        utc_day_indicator = ' '
                else:
                    utc_day_indicator = ' '
            if event.type in ['fov', 'stare', 'image', 'color', 'burn', 'autofocus', 'chill']:
                if event.status is None:
                    event.status = 'SKIPPED'
                    hhmm_text = None
                    utc_day_indicator = None
                    event.min_altitude = None
            summary_text_line = make_summary_line(event.status, hhmm_text, utc_day_indicator,
                                                  event.min_altitude, event.summary_text)
            event.summary_lines = [summary_text_line]

            # Add warning line if Image event whose RA,Dec differs from previous Image event of same target.
            if event.type == 'image':
                previous_list = radec_dict.get(event.target_name, None)
                this_dict_value = (event.ra, event.dec)
                if previous_list is None:
                    radec_dict[event.target_name] = [this_dict_value]  # start list & skip warning.
                else:
                    if any([v != this_dict_value for v in previous_list]):
                        event.summary_lines.append(
                            '      >>>>>>>>>> WARNING: ' +
                            ' Previous Image entry for ' + event.target_name + ' has different RA, Dec.')
                    radec_dict[event.target_name].append(this_dict_value)  # add value to list.

            # Add warning line if moon is too close to this object and moon is up:
            if moon_is_a_factor:
                if event.type in ['burn', 'image', 'color', 'fov', 'stare']:
                    moon_dist = an.moon_radec.degrees_from(RaDec(event.ra, event.dec))  # in degrees
                    if moon_dist < MIN_MOON_DEGREES_DEFAULT:
                        if event.utc_summary_display is not None:
                            if not (an.ts_dark_no_moon.start <= event.utc_summary_display
                                    <= an.ts_dark_no_moon.end):
                                event.summary_lines.append(
                                    '      >>>>>>>>>> WARNING: ' + event.target_name +
                                    ' MOON DISTANCE = ' +
                                    str(int(round(moon_dist))) +
                                    u'\N{DEGREE SIGN}' + ', vs. max ' +
                                    str(MIN_MOON_DEGREES_DEFAULT) +
                                    u'\N{DEGREE SIGN}')

            # Add warning line if fov target is estimated too faint in V:
            if event.type == 'fov':
                this_fov = fov_dict[event.target_name]
                if this_fov.observing_style.lower() == 'lpv':
                    mags = this_fov.estimate_lpv_mags(an.local_middark_jd)
                    v_mag = mags.get('V', None)
                    if v_mag is not None:
                        if v_mag >= V_MAG_WARNING:
                            event.summary_lines.append(
                                '      >>>>>>>>>> WARNING: above target estim. V Mag ~ ' +
                                '{:.2f}'.format(v_mag) +
                                ' very faint (>=' + '{0:g}'.format(V_MAG_WARNING) + ').')

            # Add warning line if autofocus and more than one sets requested (causing too many autofocuses):
            if event.type == 'autofocus' and plan.sets_requested > 1:
                event.summary_lines.append(
                    '      >>>>>>>>>> WARNING: autofocus not recommended when sets > 1.')

        if plan.chain_destination is not None:
            # Add plan warning line if plan chains to itself:
            if plan.chain_destination.lower() == 'plan_' + plan.plan_id.lower() + '.txt':
                plan.end_warning_lines.append(
                    '      >>>>>>>>>> ERROR: this plan attempts to chain to itself.')
            # Add plan warning line if chained-to plan does not exist:
            elif i_plan != len(plan_list) - 1:
                if plan.chain_destination.lower() != \
                    ('plan_' + plan_list[i_plan + 1].plan_id + '.txt').lower():
                    plan.end_warning_lines.append(
                        '      >>>>>>>>>> ERROR: this plan attempts to chain,'
                        ' but not to next plan.')

        # Add plan warning if no autofocus (or afinterval) given:
        if plan.afinterval is None and all([e.type != 'autofocus' for e in plan.events]):
            if any([e.type in ['burn', 'image', 'fov', 'stare', 'color'] for e in plan.events]):
                plan.end_warning_lines.append(
                    '      >>>>>>>>>> WARNING: this plan has no autofocus or afinterval.')
        # Add plan warning if autofocus and afinterval) both in same plan:
        if plan.afinterval is not None and any([e.type == 'autofocus' for e in plan.events]):
            if any([e.type in ['burn', 'image', 'fov', 'stare', 'color'] for e in plan.events]):
                plan.end_warning_lines.append(
                    '      >>>>>>>>>> WARNING: this plan has both autofocus and afinterval.')

    # Construct file contents by appending all required text lines:
    all_summary_lines = header_lines
    for plan in plan_list:
        all_summary_lines.extend(plan.summary_pre_lines)
        for event in plan.events:
            all_summary_lines.extend(event.summary_lines)
        all_summary_lines.extend(plan.end_warning_lines)
        all_summary_lines.extend(plan.summary_post_lines)

    # Write Summary file:
    output_fullpath = os.path.join(output_directory, 'Summary_' + an.an_date_string + '.txt')
    print('PRINT summary to ', output_fullpath)
    with open(output_fullpath, 'w') as this_file:
        this_file.write('\n'.join(all_summary_lines))


def make_fov_exposure_data(fov_name, an, fov_dict=None, instrument=None, exp_time_factor=1,
                           skipfilter_list=[], force_autoguide=None):
    """
    Calculates exposure data for ONE REPEAT of one given fov.
    :param fov_name: [string]
    :param an: [Astronight object]
    :param fov_dict:
    :param instrument: instrument data [Instrument object]
    :param exp_time_factor:
    :param skipfilter_list: list of filter names to omit from this FOV observation [list of strings].
    :return: tuple: (filters [list of str], counts [list of int], exp_times [list of float],
        target_overhead [float], repeat_duration [float])
    """
    if fov_dict is not None:
        # this_fov = fov_dict[fov_name]
        this_fov = fov_dict.get(fov_name, None)
        if this_fov is None:
            print(' >>>>> ERROR: FOV file not found for \'' + fov_name + '\'')
            return
    else:
        this_fov = Fov(fov_name)
    if not isinstance(instrument, Instrument):
        print(" >>>>> ERROR: make_fov_exposure_data() parm 'instrument' must be " +
              "a valid Instrument object")
        return None
    if force_autoguide is None:
        print(" >>>>> ERROR in make_fov_exposure_data(): force_autoguide is None but must be boolean.")
        return None
    obs_style = this_fov.observing_style
    filters = []
    counts = []
    exp_times = []
    mags = dict()
    omit_list = [f.lower() for f in skipfilter_list]
    for obs in this_fov.observing_list:
        filter, mag, count = obs
        if filter.lower().strip() not in omit_list:
            filters.append(filter)
            counts.append(count)
            if obs_style.lower() in ['standard', 'monitor', 'stare']:
                exp_time = calc_exp_time(mag, filter, instrument, this_fov.max_exposure,
                                         exp_time_factor=exp_time_factor)
            elif obs_style.lower() == 'lpv':
                if len(mags) == 0:
                    mags = this_fov.estimate_lpv_mags(an.local_middark_jd)  # dict (get on 1st obs only)
                exp_time = calc_exp_time(mags[filter], filter, instrument, this_fov.max_exposure,
                                         exp_time_factor=exp_time_factor)
            else:
                print('****** WARNING: fov \'' + fov_name +
                      '\' has unrecognized observing style \'' + obs_style + '\'.')
                return None
            exp_times.append(exp_time)
    if obs_style.lower() != 'stare':
        counts, exp_times = repeat_short_exp_times(counts, exp_times)
    target_overhead, repeat_duration = tabulate_target_durations(filters, counts, exp_times,
                                                                 force_autoguide=force_autoguide)
    # return types (3 lists, two floats): [str], [int], [float], float, float
    return filters, counts, exp_times, target_overhead, repeat_duration


def make_image_exposure_data(filter_entries, instrument, exp_time_factor=1, force_autoguide=None):
    """
    Calculates exposure data for given user-defined target ("IMAGE" directive).
    :param exp_time_factor: user-supplied multiplier of exp time from nominal, usually 0.5-1. [float]
    :param filter_entries: list of exposure-defining strings, as ['I=12','V=13(3)'], where I and V
        are filter names, 12 and 13 are target magnitudes, and (3) is an image count (=1 if absent).
    :param instrument: the instrument for which these exposures are wanted. [Instrument object]
    :param force_autoguide: True iff user wants to force autoguiding for this target. [boolean]
    :return: tuple of equal-length lists: (filters [str], counts [int], exp_times [float])
    """
    if force_autoguide is None:
        print(" >>>>> ERROR in make_image_exposure_data(): force_autoguide is None but must be boolean.")
        return None
    filters = []
    counts = []
    exp_times = []
    for entry in filter_entries:
        this_filter, this_mag, this_count = None, None, None
        raw_filter, mag_string = entry.split("=", maxsplit=1)
        this_filter = raw_filter.strip()
        bits = mag_string.split("(")
        if len(bits) == 1:  # case e.g. "V=13.2"
            this_count = 1
        elif len(bits) == 2:  # case e.g. "V=13.2(1)"
            try:
                this_count = int(bits[1].replace(")", ""))
            except ValueError:
                print(' >>>> PARSING ERROR:', entry)
                return
        # TODO: I'm not crazy about the next if-statement's condition.
        if 's' in bits[0].lower():
            this_exp_time = float(bits[0].lower().split('s')[0])
        else:
            this_mag = float(bits[0])
            this_exp_time = calc_exp_time(this_mag, this_filter, instrument, max_exp_time=None,
                                          exp_time_factor=exp_time_factor)
        filters.append(this_filter)
        counts.append(this_count)
        exp_times.append(this_exp_time)
    counts, exp_times = repeat_short_exp_times(counts, exp_times)
    target_overhead, repeat_duration = tabulate_target_durations(filters, counts, exp_times,
                                                                 force_autoguide=force_autoguide)
    # return types (3 lists, two floats): [str], [int], [float], float, float
    return filters, counts, exp_times, target_overhead, repeat_duration


def make_color_exposure_data(entries, force_autoguide=True):
    filters = [e[0] for e in entries]
    exp_times = [e[1] for e in entries]
    counts = [e[2] for e in entries]
    target_overhead, repeat_duration = tabulate_target_durations(filters, counts, exp_times,
                                                                 force_autoguide=force_autoguide)
    return filters, counts, exp_times, target_overhead, repeat_duration


def tabulate_target_durations(filters, counts, exp_times, force_autoguide):
    aggregate_exposure = sum([counts[i] * exp_times[i] for i in range(len(counts))])
    guiding_is_active = force_autoguide or aggregate_exposure > MAX_AGGREGATE_EXPOSURE_NO_GUIDING
    # TODO: get some of the next base values from instrument object.
    target_overhead = NEW_TARGET_DURATION + (GUIDE_STAR_ACQUISITION if guiding_is_active else 0)
    repeat_duration = aggregate_exposure + \
                      len(filters) * NEW_FILTER_DURATION + \
                      sum(counts) * NEW_EXPOSURE_DURATION_EX_GUIDER_CHECK + \
                      (sum(counts) * GUIDER_CHECK_DURATION if guiding_is_active else 0)
    return target_overhead, repeat_duration


def repeat_short_exp_times(counts, exp_times):
    for i in range(len(counts)):
        if counts[i] * exp_times[i] < MIN_TOTAL_EXP_TIME_PER_FILTER:
            counts[i] = ceil(MIN_TOTAL_EXP_TIME_PER_FILTER / exp_times[i])
    return counts, exp_times


def calc_exp_time(mag, filter, instrument, max_exp_time, exp_time_factor=1):
    # Raw exposure time from mag + properties of instrument (camera & filters).
    exp_time_from_mag = instrument.filter_data[filter]['reference_exposure_mag10'] * \
                        10.0 ** ((mag - 10.0) / 2.5)

    # Apply exposure time factor (from user, for this night) (before asymptotes and limits):
    exp_time = exp_time_factor * exp_time_from_mag

    # Apply absolute maximum as soft asymptote:
    exp_time = sqrt(1.0 / (1.0 / exp_time ** 2 + 1.0 / ABSOLUTE_MAX_EXPOSURE_TIME ** 2))

    # Apply absolute minimum as soft asymptote:
    # as of 20170406, absolute minimum is from this module, not necessarily from instrument object.
    # i.e., use more stringent of the two minima.
    effective_minimum = max(ABSOLUTE_MIN_EXPOSURE_TIME, instrument.camera['shortest_exposure'])
    exp_time = sqrt(exp_time ** 2 + effective_minimum ** 2)

    # Apply rounding (at least 2 significant digits):
    if exp_time >= 10.0:
        exp_time = round(exp_time, 0)  # round to nearest second
    else:
        exp_time = round(exp_time, 1)  # round to nearest 0.1 second

    # Apply fov's hard maximum:
    if max_exp_time is not None:
        exp_time = min(max_exp_time, exp_time)

    return exp_time


def extract_ra_dec(value_string):
    """ Split value string into subvalue, ra, dec, whether ra and dec are in standard hex format
        (e.g., 12:34:56 -11:33:42) or in TheSkyX format (e.g., 06h 49m 40.531s  +63° 00' 06.920").
    :param value_string: input string from parse_excel(), as above.
    :return: 3-tuple: subvalue (= everything but RA and Dec), ra, dec. [3-tuple of strings].
    """
    split_string = tuple(value_string.rsplit(maxsplit=2))
    if split_string[1].endswith('\'') and split_string[2].endswith('\"'):
        # RA and Dec are in TheSkyX format, e.g., 06h 49m 40.531s  +63° 00' 06.920":
        split_string = tuple(value_string.rsplit(maxsplit=6))
        if len(split_string) != 7:
            print(' >>>>> ERROR: cannot parse value string', value_string)
            return None
        subvalue = split_string[0]
        ra_items = [s.replace('h', '').replace('m', '').replace('s', '').strip() 
                    for s in split_string[1:4]]
        ra = ':'.join(ra_items)        
        dec_items = [s.replace('°', '').replace('\'', '').replace('"', '').strip() 
                     for s in split_string[4:7]]
        dec = ':'.join(dec_items)
    else:
        # RA and Dec are given directly in std hex:
        subvalue, ra, dec = split_string
    return subvalue, ra, dec