import os
import os.path
from collections import namedtuple
from datetime import datetime, timezone, timedelta
from math import floor

import ephem
import numpy as np
import pandas as pd

from photrix.fov import make_fov_dict, FovError, Fov
from photrix.user import Astronight, Instrument
from photrix.util import RaDec, datetime_utc_from_jd, hhmm_from_datetime_utc, \
    ra_as_hours, dec_as_hex, az_alt_at_datetime_utc, \
    degrees_as_hex, jd_from_datetime_utc, Timespan, find_minima_in_timespan
from photrix.web import get_aavso_webobs_raw_table

__author__ = "Eric Dose :: New Mexico Mira Project, Albuquerque"

FOV_DIRECTORY = "C:/Dev/Photometry/FOV/"
FOV_OBSERVING_STYLES = ["Standard", "Stare", "Monitor", "LPV", "Burn"]
MIN_AVAILABLE_SECONDS_DEFAULT = 900
MIN_AVAILABLE_SECONDS_STARE = 7200
MIN_MOON_DEGREES_DEFAULT = 45
MIN_MOON_DEGREES_STARE = 60
FITS_DIRECTORY = "J:/Astro/Images"
DEFAULT_PLAN_DIRECTORY = 'C:/Astro/Plans'
DT_FMT = '%Y-%m-%d %H:%M:%S.%f%z'  # kludge around py inconsistency in datetime formats


PHOTRIX_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCAL_OBS_CACHE_FULLPATH = os.path.join(PHOTRIX_ROOT_DIRECTORY, "local_obs_cache.csv")

AAVSO_WEBOBS_ROWS_TO_GET = 100
MIN_ROWS_ONE_STARE = 10
MAX_DAYS_ONE_STARE = 0.5
UPDATE_TOLERANCE_DAYS = 0.25

ABSOLUTE_MAX_EXPOSURE_TIME = 600  # seconds
AUTOFOCUS_DURATION = 180  # seconds, includes slew & filter wheel changes


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


def complete_df_fov_an(df_fov, an_string=None, site_name="DSW",
                       min_available_seconds=MIN_AVAILABLE_SECONDS_DEFAULT,
                       min_moon_degrees=MIN_MOON_DEGREES_DEFAULT,
                       remove_zero_an_priority=True,
                       remove_unobservables=True):
    if an_string is None or site_name == "":
        return df_fov
    an = Astronight(an_string, site_name)
    # print("\n\nan dark        >", an.ts_dark)
    # print("an dark_no_moon>", an.ts_dark_no_moon)
    # print("moon phase >", an.moon_phase)
    # print("min alt > ", an.site.min_altitude)
    # print("min_moon_degrees > ", min_moon_degrees)

    # Construct columns (specific to night and site) for available obs time, this astronight.
    df_fov = df_fov.assign(moon_deg=0.0) \
        .assign(start=an.local_middark_utc) \
        .assign(end=an.local_middark_utc) \
        .assign(mid=an.local_middark_utc) \
        .assign(seconds=0.0) \
        .assign(available=' - '.join(2*[4*' '])) \
        .assign(an_priority=0.0) \
        .assign(an_priority_bars='') # all dummy values to be overwritten later.

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
    # Fill in night's priorities.
    loc = LocalObsCache()
    loc.update_fov_entries(fov_list=df_fov['fov_name'])
    for ind in df_fov.index:
        df_fov.loc[ind, 'an_priority'] = loc.calc_an_priority(df_fov.loc[ind, 'fov_name'], an)
        max_bars = 16
        int_an_priority = int(round(df_fov.loc[ind, 'an_priority']))
        df_fov.loc[ind, 'an_priority_bars'] = \
            (8*'*' + (max_bars-8)*'#')[0: min(max_bars, int_an_priority)].ljust(max_bars)

    if remove_zero_an_priority:
        df_fov = df_fov[df_fov['an_priority'] > 0.0]

    if remove_unobservables:
        enough_dark_time = df_fov['seconds'] >= min_available_seconds
        moon_dist_ok = df_fov['moon_deg'] >= min_moon_degrees
        is_observable = enough_dark_time & moon_dist_ok
        df_fov = df_fov[is_observable]

    return df_fov.sort_values(by=['mid', 'an_priority'])


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
    """
    def __init__(self, fov_dict=None):
        # Read in local cache if it exists.
        if os.path.isfile(LOCAL_OBS_CACHE_FULLPATH):
            self.df_cache = self._read_cache_from_csv()
            need_to_create_empty_cache = self.df_cache is None
        else:
            need_to_create_empty_cache = True
        if need_to_create_empty_cache:
            #  Create *empty* dataframe with dtypes (incl. utc datetimes), write to cache file:
            # self.df_cache = pd.DataFrame({'fov_name': ['dummy'],
            #                               'main_target': ['dummy'],
            #                               'obs_style': ['dummy'],
            #                               'cache_datetime': [datetime.now(timezone.utc)],
            #                               'obs_datetime': [datetime.now(timezone.utc)],
            #                               'obs_mag': [0.0],
            #                               'obs_mag_filter': ['dummy']},
            #                              columns=['fov_name', 'main_target', 'obs_style',
            #                                       'cache_datetime', 'obs_datetime',
            #                                       'obs_mag', 'obs_mag_filter'])[:0]
            #
            self.df_cache = pd.DataFrame.from_items([('fov_name', ['dummy']),
                                            ('main_target', ['dummy']),
                                            ('obs_style', ['dummy']),
                                            ('cache_datetime', [datetime.now(timezone.utc)]),
                                            ('obs_datetime', [datetime.now(timezone.utc)]),
                                            ('obs_mag', [0.0]),
                                            ('obs_mag_filter', ['dummy'])])[:0]
            self.df_cache.index.name = 'row_index'
            csv_fullpath = self._write_cache_to_csv()  # empty cache to csv
            print('LocalObsCache: wrote new, empty cache file to ' + csv_fullpath)
        if fov_dict is None:
            self.fov_dict = make_fov_dict()
            print('LocalObsCache: made new fov_dict with ' +
                  str(len(self.fov_dict)) + ' fov entries.')
        else:
            self.fov_dict = fov_dict

    def update_one_fov_entry(self, fov_name, write_csv=False):
        """
        This class's engine. Updates cache's entry for one fov, if entry is too aged.
        :param write_csv:
        :param fov_name: [string]
        :return: if cache was updated, datetime (UTC) of new obs; else None.
        """
        # TODO: If query to AAVSO yields no latest obs, put some placeholder with cache_dt at least.
        this_fov = self.fov_dict.get(fov_name)  # Fov object for this fov_name.
        if this_fov is None:
            raise FovError
        main_target = this_fov.main_target
        # self._curate_df_cache(fov_name, main_target)

        #  Determine whether update is needed, return if not.
        cache_row_pre_exists = fov_name.lower() in list(self.df_cache['fov_name'].str.lower())
        if cache_row_pre_exists:
            now = datetime.now(timezone.utc)
            current_cache_datetime = self.df_cache.loc[fov_name, 'cache_datetime']
            update_age = (now - current_cache_datetime).total_seconds() / (24*3600)
            entry_fresh_enough = update_age <= UPDATE_TOLERANCE_DAYS
            if entry_fresh_enough:
                return self.df_cache.loc[fov_name, 'obs_datetime']  # skip updating

        # Update fov's cache entry, from AAVSO webobs.
        obs_style = this_fov.observing_style
        obs_style_lower = obs_style.lower()
        target_type_lower = this_fov.target_type.lower()
        if target_type_lower == 'standard':
            return None
        if obs_style_lower == 'stare':
            num_obs = 200
        else:
            num_obs = 100
        print('AAVSO webobs query \'' + main_target + '\'...', end='', flush=True)
        recent_observations = AavsoWebobs(star_id=main_target, num_obs=num_obs)  # from AAVSO
        print('ok.', end='', flush=True)
        latest_obs_df = None  # default if no matches.
        if (obs_style_lower, target_type_lower) == ('lpv', 'mira'):
            latest_obs_df = self._latest_single_obs(fov_name, obs_style, recent_observations,
                                                    allow_filters=['V'])
        elif (obs_style_lower, target_type_lower) == ('lpv', 'lpv'):
            latest_obs_df = self._latest_single_obs(fov_name, obs_style, recent_observations,
                                                    allow_filters=['V', 'R'])
        elif obs_style_lower == 'monitor' and target_type_lower != 'astrometric':
            latest_obs_df = self._latest_single_obs(fov_name, obs_style, recent_observations,
                                                    allow_filters=['V', 'R'])
        elif obs_style_lower == 'stare':
            latest_obs_df = self._latest_stare_obs(fov_name, obs_style, recent_observations,
                                                   allow_filters=['V', 'R'])
        else:
            print('\n*** WARNING: for fov \'' + fov_name + '(obs_style, target_type) = (' +
                  obs_style + ', ' + this_fov.target_type + ') not understood.', end='', flush=True)
        if cache_row_pre_exists:
            self.df_cache = latest_obs_df.combine_first(self.df_cache)  # overwrites.
        else:
            #  This else-block is kludge for pandas' mis-handling of append to empty DataFrame.
            if len(self.df_cache) >= 1:
                self.df_cache = self.df_cache.append(latest_obs_df)
            else:
                self.df_cache = latest_obs_df.copy()
        if write_csv:
            self._write_cache_to_csv()
            print('..csv written.', end='', flush=True)
        print('')
        if latest_obs_df is None:
            return None
        return latest_obs_df.iloc[0].loc['obs_datetime']  # obs datetime, to signal OK.

    def update_fov_entries(self, fov_list, max_fovs_since_write=6):
        """
        For each fov available this night (in df_fov_list), update the cache.
        :param fov_list: list of fovs (typically fovs available this night, as df_fov['fov_name']).
        :param max_fovs_since_write: controls frequence of writes to cache.
        :return: number of fovs updated (fn effect is to update this class's cache dataframe.
        """
        fovs_since_write = 0
        for fov_name in fov_list:
            need_to_write_csv = (fovs_since_write >= max_fovs_since_write - 1)
            self.update_one_fov_entry(fov_name, write_csv=need_to_write_csv)
            if need_to_write_csv:
                fovs_since_write = 0
            else:
                fovs_since_write += 1
        self._write_cache_to_csv()

    def _latest_single_obs(self, fov_name, obs_style, recent_observations, allow_filters):
        """
        Takes a AavsoWebObs object and returns a pandas dataframe ready for inclusion
           in LocalCacheObs dataframe df_cache.
           Single-observation case (not stare).
        :param fov_name: [string]
        :param obs_style: [string] ('Monitor' or 'LPV')
        :param recent_observations: recent observations for fov_name [AavsoWebObs object].
        :param allow_filters: list of filters [string] to include in finding latest observation.
        :return: 1-row dataframe of relevant data about latest stare observation for this fov_name;
                 return (with some None values) if no qualifying observation is found.
        """
        latest_obs_df = None
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
            latest_obs_df = pd.DataFrame.from_items([
                ('fov_name', [fov_name]),
                ('main_target', [self.fov_dict[fov_name].main_target]),
                ('obs_style', [obs_style]),
                ('cache_datetime', [datetime.now(timezone.utc)]),
                ('obs_datetime', [None]),
                ('obs_mag', [None]),
                ('obs_mag_filter', [None])])
            for column_name in ['cache_datetime']:
                latest_obs_df[column_name] = [x.to_pydatetime()
                                              for x in latest_obs_df[column_name]]
        else:
            latest_obs_df = pd.DataFrame.from_items([
                ('fov_name', [fov_name]),
                ('main_target', [self.fov_dict[fov_name].main_target]),
                ('obs_style', [obs_style]),
                ('cache_datetime', [datetime.now(timezone.utc)]),
                ('obs_datetime', [datetime_utc_from_jd(latest_obs.jd)]),
                ('obs_mag', [latest_obs.mag]),
                ('obs_mag_filter', [latest_obs.loc['filter']])])
            for column_name in ['cache_datetime', 'obs_datetime']:
                latest_obs_df[column_name] = [x.to_pydatetime()
                                              for x in latest_obs_df[column_name]]
        latest_obs_df.index = latest_obs_df['fov_name'].copy()
        latest_obs_df.index.name = 'row_index'
        return latest_obs_df

    def _latest_stare_obs(self, fov_name, obs_style, recent_observations, allow_filters):
        """
        Takes a AavsoWebObs object and returns a 1-row pandas dataframe ready for inclusion
            in LocalCacheObs dataframe df_cache.
            Stare case (multiple observations in one night), typically for eclipsers.
        :param fov_name: [string]
        :param obs_style: [string] (almost certainly == 'stare')
        :param recent_observations: recent observations for fov_name [AavsoWebObs object].
        :param allow_filters: list of filters [string] to include in finding latest observation.
        :return: dataframe of relevant data about latest stare observation for this fov_name;
                 return (with some None values) if no qualifying stare observation is found.
        """
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
                        test_latest_jd = table_this_filter['jd']\
                            .iloc[first_test_irow]
                        test_earliest_jd = table_this_filter['jd']\
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
                                latest_stare_obs_df = pd.DataFrame.from_items([
                                     ('fov_name', [fov_name]),
                                     ('main_target', [self.fov_dict[fov_name].main_target]),
                                     ('obs_style', [obs_style]),
                                     ('cache_datetime', [datetime.now(timezone.utc)]),
                                     ('obs_datetime', [datetime_utc_from_jd(latest_stare_obs.jd)]),
                                     ('obs_mag', [latest_stare_obs.mag]),
                                     ('obs_mag_filter', [latest_stare_obs.loc['filter']])])
                                for column_name in ['cache_datetime', 'obs_datetime']:
                                    latest_stare_obs_df[column_name] = \
                                        [x.to_pydatetime()
                                         for x in latest_stare_obs_df[column_name]]
                                latest_stare_obs_df.index = latest_stare_obs_df['fov_name'].copy()
                                latest_stare_obs_df.index.name = 'row_index'

        if latest_stare_obs_df is None:
            #  If no qualified stare observation found within webobs query,
            #  construct placeholder row in df_cache, to prevent repeating query needlessly.
            latest_stare_obs_df = pd.DataFrame.from_items([
                ('fov_name', [fov_name]),
                ('main_target', [self.fov_dict[fov_name].main_target]),
                ('obs_style', [obs_style]),
                ('cache_datetime', [datetime.now(timezone.utc)]),
                ('obs_datetime', [None]),
                ('obs_mag', [None]),
                ('obs_mag_filter', [None])])
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
            if row['obs_datetime'] is None or isinstance(row['obs_datetime'], pd.tslib.NaTType):
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

    def calc_an_priority(self, fov_name, an):
        """
        Calculates astronight priority for one fov.
        :param fov_name:
        :param an: Astronight object for the night in question.
        :return: an_priority, from fov_priority and age of most recent obs.
        """
        this_fov = self.fov_dict.get(fov_name, None)
        if this_fov is None:
            print('LOC.calc_an_priority: fov \'' + fov_name + '\' not found in fov_dict.')
            return None
        if (this_fov.priority is None) or (this_fov.target_type.lower() == 'standard'):
            return 0
        if this_fov.priority <= 0:
            return 0

        self.update_one_fov_entry(fov_name, write_csv=True)
        if fov_name not in self.df_cache.index:
            return 2 * this_fov.priority  # the maximum, since no latest obs was accessible.
        latest_obs = self.df_cache.loc[fov_name]
        if latest_obs.obs_datetime is None:
            return 2 * this_fov.priority  # the maximum, since no latest obs was accessible.
        jd_latest_obs = jd_from_datetime_utc(latest_obs.obs_datetime)
        age_days = an.local_middark_jd - jd_latest_obs
        return this_fov.calc_priority_score(age_days)

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


# def go(fov_name):
#     from photrix.user import Astronight
#     an = Astronight('20170127', 'DSW')
#     loc = LocalObsCache()
#     return loc.calc_an_priority(fov_name=fov_name, an=an)


# def stares():
#     an = Astronight('20170211', 'DSW')
#     loc = LocalObsCache()
#     stare_fov_names = [fov.fov_name for fov in loc.fov_dict.values()
#                        if fov.observing_style.lower() == 'stare']
#     stare_fov_names.sort(key=str.lower)  # in place; sort is mostly for debugging reproducibility.
#     fov_name_list = []
#     priority_list = []
#     bar_list = []
#     for fov_name in stare_fov_names:
#         fov_name_list.append(fov_name)
#         an_priority = loc.calc_an_priority(fov_name=fov_name, an=an)
#         int_an_priority = int(round(an_priority))
#         priority_list.append(int_an_priority)
#         bar_list.append(int_an_priority * '#')
#     df_stares = pd.DataFrame.from_items([('fov', fov_name_list),
#                                          ('an_priority', priority_list),
#                                          ('bars', bar_list)])
#     df_stares = df_stares.sort_values(by='fov')
#     # print(df_stares)
#     return df_stares


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


def get_local_obs_age_dict(fov_dict=None, report_dir=None, target_an=None, limit_days=366):
        # TODO: finish writing get_local_obs_age_dict()
        """
        report_dir: directory in which all relevant AAVSO reports reside, as
          "J:/Astro/2016/Photometry".
        target_an: target astronight from which to count days, as "20151216".
        limit_days: days into the past to look up old AAVSO reports.
        Returns dict of (fov_name, days_since_last_local_obs).
        """
        pass
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

def make_an_roster(an_date_string, output_directory, site_name='DSW', instrument_name='Borea'):
    """
    Generates new .csv file containing info on each fov available this astronight.
    :param an_date_string: as '20170127. Date of the evening to plan for [string]
    :param output_directory: directory in which to write Roster csv file [string]
    :param site_name: [string]
    :param instrument_name: [string]
    :return: tuple of number of fovs, each obs style: (n_std, n_monitor_lpv, n_stare). [ints]
    """

    an = Astronight(an_date_string=an_date_string, site_name=site_name)
    df_fov = make_df_fov(fov_directory=FOV_DIRECTORY, fov_names_selected=None)
    instrument = Instrument(instrument_name)
    lines_header = ['ROSTER file for     ' + an_date_string + '     generated by photrix ' +
                    '{:%Y-%m-%d %H:%M  UTC}'.format(datetime.now(timezone.utc)),
                    an.acp_header_string().replace(',',' ')]

    # Handle obs_style = 'Standard':
    lines_std = ['\n\n\nSTANDARD roster for ' + an_date_string + ':' + 50*'-',
                 ',fov,fov, avail_utc,transit,minutes,stars']
    df_fov_std = filter_df_fov_by_obs_styles(df_fov, obs_style_list=['Standard'])
    df_fov_std = complete_df_fov_an(df_fov_std, an_string=an_date_string, site_name=site_name,
                                    min_moon_degrees=MIN_MOON_DEGREES_DEFAULT,
                                    remove_zero_an_priority=False, remove_unobservables=True)
    for fov_index in df_fov_std.index:
        fov_name = df_fov_std.loc[fov_index, 'fov_name']
        available = df_fov_std.loc[fov_index, 'available']
        this_fov = Fov(fov_name)
        transit_hhmm = hhmm_from_datetime_utc(an.transit(RaDec(this_fov.ra, this_fov.dec)))
        _, _, _, target_overhead, repeat_duration = make_obs_exposure_data(fov_name,
                                                                           an, fov_dict=None,
                                                                           instrument=instrument)
        minutes = (target_overhead + repeat_duration) / 60.0
        n_stars = len(this_fov.aavso_stars)
        this_fov_line = ',' + fov_name + ',' + fov_name + ', ' + available + ',' +\
                        "=\"" + transit_hhmm + "\"" + ',' + '{:.1f}'.format(minutes) +\
                        ',' + '{:3d}'.format(n_stars)
        lines_std.append(this_fov_line)

    # Handle obs_style = 'Stare':
    lines_stare = ['\n\n\nSTARE roster for ' + an_date_string + ':' + 50*'-',
                   ',fov,fov, avail_utc,transit,min/rpt,an_priority,,period,minima']
    df_fov_stare = filter_df_fov_by_obs_styles(df_fov, obs_style_list=['Stare'])
    df_fov_stare = filter_df_fov_by_fov_priority(df_fov_stare,
                                                 min_fov_priority=0.5, include_std_fovs=False)
    df_fov_stare = complete_df_fov_an(df_fov_stare, an_string=an_date_string, site_name=site_name,
                                      min_moon_degrees=MIN_MOON_DEGREES_STARE,
                                      remove_zero_an_priority=True, remove_unobservables=True)
    for fov_index in df_fov_stare.index:
        row = df_fov_stare.loc[fov_index]
        fov_name = row.loc['fov_name']
        available = row.loc['available']
        this_fov = Fov(fov_name)
        transit_hhmm = hhmm_from_datetime_utc(an.transit(RaDec(this_fov.ra, this_fov.dec)))
        if str(fov_name).startswith('VSX'):
            iii = 4
        _, _, _, target_overhead, repeat_duration = make_obs_exposure_data(fov_name,
                                                                           an, fov_dict=None,
                                                                           instrument=instrument)
        minutes = (target_overhead + repeat_duration) / 60.0
        an_priority = row.loc['an_priority']
        an_priority_bars = row.loc['an_priority_bars']
        period = row.loc['period']
        row_ts = Timespan(row.loc['start'], row.loc['end'])  # Timespan object for this row.

        list_primary_mins = find_minima_in_timespan(this_fov.JD_faint, this_fov.period, row_ts)
        if list_primary_mins is None:
            primaries_exist = False
        else:
            df_primary_mins = pd.DataFrame.from_items([('utc', list_primary_mins, )])
            df_primary_mins['type'] = 1
            primaries_exist = len(df_primary_mins) >= 1

        list_secondary_mins = find_minima_in_timespan(this_fov.JD_second, this_fov.period, row_ts)
        if list_secondary_mins is None:
            secondaries_exist = False
        else:
            df_secondary_mins = pd.DataFrame.from_items([('utc', list_secondary_mins,)])
            df_secondary_mins['type'] = 2
            secondaries_exist = len(df_secondary_mins) >= 1

        if primaries_exist and secondaries_exist:
            df_mins = df_primary_mins.append(df_secondary_mins)
        elif primaries_exist and not secondaries_exist:
            df_mins = df_primary_mins
        elif not primaries_exist and secondaries_exist:
            df_mins = df_secondary_mins
        else:
            df_mins = None
        if df_mins is not None:
            motive = this_fov.motive
            df_mins.sort_values(by='utc', inplace=True)
            minima = 'minima: '
            for row in df_mins.itertuples():
                minima += str(row.type) + "=" + hhmm_from_datetime_utc(row.utc) + '  '
            this_fov_line = ',' + fov_name + ',' + fov_name + ',' + available + ',' + \
                            "=\"" + transit_hhmm + "\"" + ',' + '{:.1f}'.format(minutes) + ',' +\
                            str(int(round(an_priority))) + ' ,' + an_priority_bars + ',' +\
                            '{:7.3f}'.format(period) + ' ,' + minima + ',' + motive
            lines_stare.append(this_fov_line)

    # Handle obs_style = 'Monitor' or 'LPV':
    lines_mon_lpv = ['\n\n\nMONITOR / LPV roster for ' + an_date_string + ':' + 50*'-',
                     ',fov,fov,avail_utc,transit,minutes,an_priority']
    df_fov_mon_lpv = filter_df_fov_by_obs_styles(df_fov, obs_style_list=['Monitor', 'LPV'])
    df_fov_mon_lpv = filter_df_fov_by_fov_priority(df_fov_mon_lpv,
                                                   min_fov_priority=0.5, include_std_fovs=False)
    df_fov_mon_lpv = complete_df_fov_an(df_fov_mon_lpv, an_string=an_date_string,
                                        site_name=site_name,
                                        min_moon_degrees=MIN_MOON_DEGREES_DEFAULT,
                                        remove_zero_an_priority=True, remove_unobservables=True)
    for fov_index in df_fov_mon_lpv.index:
        fov_name = df_fov_mon_lpv.loc[fov_index, 'fov_name']
        available = df_fov_mon_lpv.loc[fov_index, 'available']
        this_fov = Fov(fov_name)
        transit_hhmm = hhmm_from_datetime_utc(an.transit(RaDec(this_fov.ra, this_fov.dec)))
        _, _, _, target_overhead, repeat_duration = make_obs_exposure_data(fov_name,
                                                                           an, fov_dict=None,
                                                                           instrument=instrument)
        minutes = (target_overhead + repeat_duration) / 60.0
        an_priority = df_fov_mon_lpv.loc[fov_index, 'an_priority']
        an_priority_bars = df_fov_mon_lpv.loc[fov_index, 'an_priority_bars']
        motive = Fov(fov_name).motive
        this_fov_line = ',' + fov_name + ',' + fov_name + ', ' + available + ',' + \
                        "=\"" + transit_hhmm + "\"" + ',' + '{:.1f}'.format(minutes) + ',' +\
                        str(int(round(an_priority))) + ' ,' + an_priority_bars + ',' + motive
        lines_mon_lpv.append(this_fov_line)

    # Assemble all lines, and write to file.
    lines_all = lines_header + lines_std + lines_stare + lines_mon_lpv
    output_fullpath = os.path.join(output_directory, 'Roster_' + an.an_date_string + '.csv')
    csv_written = False
    while not csv_written:
        try:
            with open(output_fullpath, 'w') as this_file:
                this_file.write('\n'.join(lines_all))
            csv_written = True
        except PermissionError:
            input('***** CLOSE file \'' + output_fullpath + '\' and hit Enter.')


def make_an_plan(plan_excel_path='c:/24hrs/Planning.xlsx', site_name='DSW', instrument_name='Borea',
                 fov_dict=None, an_start_hhmm=None):
    """
    Main user fn to take sketch Excel file and generate Summary and ACP Plan files.
    :param plan_excel_path: full path to Excel file holding all info for one night's observations.
    :param site_name: a Site object for location of observations.
    :param instrument_name: an Instrument object for scope to be used.
    :param fov_dict: fov_dict if available, default=None to generate new fov_dict (normal case).
    :param an_start_hhmm: 'hhmm' time to start plan, default=None for 'use excel file (normal case).
    :return: Writes out Summary file with dateline, and one or more ACP plan files.
    """

    # TODO: LocalObsCache updates only for fovs actually used, not including Burns.
    parsed_list, an = parse_excel(plan_excel_path, site_name)

    if fov_dict is None:
        fov_dict = make_fov_dict()
    instrument = Instrument(instrument_name)
    raw_plan_list = make_raw_plan_list(parsed_list, an)

    reordered_plan_list = reorder_actions(raw_plan_list)

    plan_list = add_raw_durations_and_lines(reordered_plan_list, an, fov_dict, instrument)

    plan_list = add_timetable(plan_list, an, an_start_hhmm)

    plan_list = add_altitudes(plan_list, an, fov_dict)

    output_directory = os.path.split(plan_excel_path)[0]  # output files go to same dir as excel input.
    write_acp_plans(plan_list, output_directory)
    write_summary(plan_list, an, output_directory)


def parse_excel(excel_path, site_name='DSW'):
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
    # print('an_date_string: ' + an_date_string)  # TEST

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
                # print(cell_str_as_read)

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
                    plan_actions.append(('autofocus', '#AUTOFOCUS'))  # i.e., directly from user
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
                elif cell_str_lower.startswith('shutdown'):
                        plan_actions.append(('shutdown',))
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
                # print(plan_actions[-1:])

    parsed_list.append(plan_actions)  # Close out the last plan.
    return parsed_list, an


def make_raw_plan_list(parsed_list, an):
    # Construct raw master plan_list (with as-yet incomplete actions):
    Plan = namedtuple('Plan', ['plan_id', 'action_list'])
    Action = namedtuple('Action', ['action_type', 'parsed_action', 'an_priority',
                                   'raw_duration', 'n_afinterval_autofocus',
                                   'status', 'start_utc', 'altitude_deg',
                                   'summary_lines', 'acp_plan_lines'])
    raw_plan_list = []  # will be the master container (a list of Plan namedtuples)
    for parsed_plan in parsed_list:
        # print('\n***' + str(parsed_plan[0]))
        # for plan_action in parsed_plan[1:]:
        #     print(str(plan_action))

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
                                 n_afinterval_autofocus=0,  # possibly overwrite later
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


def add_raw_durations_and_lines(plan_list, an, fov_dict, instrument):
    for i_plan in range(len(plan_list)):  # a Plan namedtuple
        this_plan = plan_list[i_plan]
        for i_action in range(len(this_plan.action_list)):  # an Action namedtuple
            action = this_plan.action_list[i_action]
            this_type = action.action_type.lower()
            parsed_action = action.parsed_action
            if this_type == 'plan':
                if parsed_action[2] is None:
                    text = parsed_action[1]
                else:
                    text = parsed_action[1] + ' ; ' + parsed_action[2]
                summary_lines = ['', 55 * '-', 'Begin Plan ' + text]
                acp_plan_lines = an.acp_header_string().split('\n') + [';']
                raw_duration = 60
            elif this_type == 'chill':
                summary_lines = ['CHILL  ' + parsed_action[1]]
                acp_plan_lines = ['#CHILL  ' + parsed_action[1]]
                raw_duration = 0
            elif this_type == 'waituntil':
                if parsed_action[1] == 'hhmm':
                    dt = an.datetime_utc_from_hhmm(parsed_action[2])
                    formatted_time = '{:%m/%d/%Y %H:%M}'.format(dt)
                    hhmm = hhmm_from_datetime_utc(dt)  # to let user verify correct parsing
                    summary_lines = ['WAITUNTIL ' + hhmm + ' utc']
                    acp_plan_lines = ['#WAITUNTIL 1, ' + formatted_time + ' ; utc']
                elif parsed_action[1] == 'sun_degrees':
                    summary_lines = ['WAITUNTIL sun reaches ' +
                                     parsed_action[2] + u'\N{DEGREE SIGN}' + ' alt']
                    acp_plan_lines = ['#WAITUNTIL 1, ' + parsed_action[2] + ' ; deg sun alt']
                else:
                    print("***** ERROR: WAITUNTIL action" + str(parsed_action) +
                          ' is not understood.')
                    summary_lines, acp_plan_lines = [], []
                raw_duration = 0  # special case: no duration but executes a delay
            elif this_type == 'quitat':
                dt = an.datetime_utc_from_hhmm(parsed_action[1])
                formatted_time = '{:%m/%d/%Y %H:%M}'.format(dt)
                hhmm = hhmm_from_datetime_utc(dt)  # to let user verify correct parsing
                summary_lines = ['QUITAT ' + hhmm + ' utc']
                acp_plan_lines = ['#QUITAT ' + formatted_time + ' ; utc']
                raw_duration = 0  # special case: no duration but may terminate the plan
            elif this_type == 'afinterval':
                summary_lines = ['AFINTERVAL ' + parsed_action[1]]
                acp_plan_lines = [';', '#AFINTERVAL  ' + parsed_action[1]]
                raw_duration = 0  # but prob inserts AUTOFOCUS lines which do have duration
            elif this_type == 'autofocus':
                summary_lines = ['AUTOFOCUS']
                acp_plan_lines = [';', '#AUTOFOCUS']
                raw_duration = AUTOFOCUS_DURATION
            elif this_type == 'comment':
                summary_lines = [';' + parsed_action[1]]
                acp_plan_lines = [';' + parsed_action[1]]
                raw_duration = 0
            elif this_type == 'burn':
                summary_lines = ['BURN ' + parsed_action[1] + '  ' +
                                 parsed_action[2] + '  ' + parsed_action[3]]
                acp_plan_lines = [';', '#DITHER 0 ;', '#FILTER V,I ;', '#BINNING 1,1 ;',
                                  '#COUNT 1,1 ;', '#INTERVAL 240,240 ;',
                                  ';----> BURN for new FOV file.',
                                  parsed_action[1] + '\t' +
                                  parsed_action[2] + '\t' + parsed_action[3] + ' ;']
                raw_duration = 11 * 60
            elif this_type == 'fov':
                fov_name = parsed_action[1]
                summary_lines = ['Observe ' + fov_name]
                filters, counts, exp_times, target_overhead, repeat_duration = \
                    make_obs_exposure_data(fov_name, an, fov_dict, instrument)
                raw_duration = target_overhead + 1 * repeat_duration
                duration_comment = ' --> ' + str(round(raw_duration / 60.0, 1)) + ' min'
                this_fov = fov_dict[fov_name]
                acp_plan_lines = [';', '#DITHER 0 ;',
                                  '#FILTER ' + ','.join(filters) + ' ;',
                                  '#BINNING ' + ','.join(len(filters)*['1']) + ' ;',
                                  '#COUNT ' + ','.join([str(c) for c in counts]) + ' ;',
                                  '#INTERVAL ' + ','.join([str(round(e)) for e in exp_times]) +
                                  ' ; ' + duration_comment,
                                  ';----' + this_fov.acp_comments, fov_name + '\t' +
                                  ra_as_hours(this_fov.ra) + '\t' + dec_as_hex(this_fov.dec)]
            elif this_type == 'stare':
                n_repeats, fov_name = int(parsed_action[1]), parsed_action[2]
                summary_lines = ['Stare ' + str(n_repeats) + ' repeats at ' + fov_name]
                filters, counts, exp_times, target_overhead, repeat_duration = \
                    make_obs_exposure_data(fov_name, an, fov_dict, instrument)
                raw_duration = target_overhead + n_repeats * repeat_duration
                duration_comment = str(round(repeat_duration / 60.0, 1)) + ' min/repeat --> ' + \
                                   str(round(raw_duration / 60.0, 1)) + ' min (raw)'
                this_fov = fov_dict[fov_name]
                acp_plan_lines = [';', '#REPEAT ' + str(n_repeats) + ';',
                                  '#DITHER 0 ;',
                                  '#FILTER ' + ','.join(filters) + ' ;',
                                  '#BINNING ' + ','.join(len(filters)*['1']) + ' ;',
                                  '#COUNT ' + ','.join([str(c) for c in counts]) + ' ;',
                                  '#INTERVAL ' + ','.join([str(round(e)) for e in exp_times]) +
                                  ' ; ' + duration_comment,
                                  ';----' + this_fov.acp_comments, fov_name + '\t' +
                                  ra_as_hours(this_fov.ra) + '\t' + dec_as_hex(this_fov.dec)]
            elif this_type == 'flats':
                flats_filename = 'flats_VRI_16.txt'
                summary_lines = ['Flats: ' + flats_filename]
                acp_plan_lines = [';', '#SCREENFLATS ' + flats_filename + ' ;']
                raw_duration = 777  # <--------- CODE THIS !!!
            elif this_type == 'darks':
                summary_lines = ['dummy']
                acp_plan_lines = ['dummy']
                raw_duration = 888  # <--------- CODE THIS !!!
            elif this_type == 'chain':
                summary_lines = ['Chain to \'' + parsed_action[1] + '\'']
                acp_plan_lines = [';', '#CHAIN ' + parsed_action[1]]
                raw_duration = 0  # though the plan chained to will have a duration
            elif this_type == 'shutdown':
                summary_lines = ['SHUTDOWN']
                acp_plan_lines = [';', '#SHUTDOWN']
                raw_duration = 480  # a good guess
            else:
                print('*** WARNING: insert_raw_duration() cannot understand action \'' +
                      action.action_type + 'in plan \'' + this_plan.plan_id + '\'')
                summary_lines, acp_plan_lines, raw_duration = [], [], None

            # Add warning if plan's last action is not chain (except last plan):
            this_plan_is_not_last = (i_plan < len(plan_list)-1)
            this_action_is_last = (i_action == len(this_plan.action_list)-1)
            this_action_should_be_chain = this_plan_is_not_last and this_action_is_last
            if this_action_should_be_chain and this_type != 'chain':
                summary_lines += \
                    ['*** WARNING: plan does not end on #CHAIN but probably should do.']
            new_action = action._replace(summary_lines=summary_lines,
                                         acp_plan_lines=acp_plan_lines,
                                         raw_duration=raw_duration)
            this_plan.action_list[i_action] = new_action
    return plan_list


def add_timetable(plan_list, an, an_start_hhmm):
    """
    Takes naive "raw durations" for Actions, and calculates and inserts realistic timetable data.
    Also inserts imputed Autofocus actions as required by #AFINTERVAL directives.
    :param plan_list: input list of Plan namedtuples
    :param an:
    :param an_start_hhmm: HHMM string denoting desired UTC start time, or None for dusk twilight.
    :return: list of Plan namedtuples with 'status', 'start_utc', and 'altitude_deg' filled in.
    """
    # Calculate timeline and action status, add to each action in each plan:
    # Develop projected timeline (starting time for each action in entire night):
    if an_start_hhmm is None:
        an_start_dt = an.ts_dark.start
    else:
        an_start_dt = an.datetime_utc_from_hhmm(an_start_hhmm)
    running_dt = an_start_dt
    for i_plan in range(len(plan_list)):
        this_plan = plan_list[i_plan]  # a Plan namedtuple

        # Scan this plan for #WAITUNTIL or #QUITAT, which alter timelines directly,
        #    and for #AFINTERVAL, which is likely to insert AUTOFOCUS actions:
        quitat_dt = None  # datetime utc
        waituntil_dt = None  # datetime utc
        afinterval_is_active = False
        previous_autofocus_dt = None
        for this_action in this_plan.action_list:  # an Action namedtuple
            if this_action.action_type == 'waituntil':
                if this_action.parsed_action[1] == 'hhmm':
                    waituntil_dt = an.datetime_utc_from_hhmm(this_action.parsed_action[2])
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
                quitat_dt = an.datetime_utc_from_hhmm(this_action.parsed_action[1])
            elif this_action.action_type == 'afinterval':
                afinterval_is_active = True
                afinterval_timedelta = timedelta(seconds=float(this_action.parsed_action[1]) * 60)

        # Set starting time for this plan:
        if waituntil_dt is not None:
            if waituntil_dt > running_dt:
                running_dt = waituntil_dt

        # Now construct starting time and completion status for each action in this plan:
        observation_is_first_in_plan = True
        for i_action in range(len(this_plan.action_list)):
            this_action = this_plan.action_list[i_action]  # an Action namedtuple
            action_start_dt = running_dt
            action_timedelta = timedelta(seconds=this_action.raw_duration)
            action_expected_end_dt = running_dt + action_timedelta
            n_afinterval_autofocus = 0  # default

            # If AFINTERVAL in play, then make an Afinterval_autofocus namedtuple & advance time.
            if afinterval_is_active:
                if this_action.action_type in ['fov', 'stare', 'burn']:
                    if observation_is_first_in_plan:
                        action_needs_pre_afinterval_autofocus = True
                        n_afinterval_autofocus_during_action = floor(action_timedelta /
                                                           (afinterval_timedelta -
                                                            timedelta(seconds=AUTOFOCUS_DURATION)))
                        n_afinterval_autofocus = 1 + n_afinterval_autofocus_during_action
                        action_timedelta += timedelta(seconds=
                                                      n_afinterval_autofocus *
                                                      AUTOFOCUS_DURATION)  # overwrite
                        action_expected_end_dt = running_dt + action_timedelta  # overwrite
                        previous_autofocus_dt = action_start_dt + \
                                                n_afinterval_autofocus_during_action *\
                                                timedelta(seconds=AUTOFOCUS_DURATION)
                        observation_is_first_in_plan = False
                    else:
                        n_afinterval_autofocus_during_action = floor((action_expected_end_dt -
                                                        previous_autofocus_dt) /
                                                        (afinterval_timedelta -
                                                        timedelta(seconds=AUTOFOCUS_DURATION)))
                        n_afinterval_autofocus = n_afinterval_autofocus_during_action
                        action_timedelta += n_afinterval_autofocus * \
                                            timedelta(seconds=AUTOFOCUS_DURATION)
                        action_expected_end_dt = running_dt + action_timedelta  # overwrite
                        previous_autofocus_dt += n_afinterval_autofocus *\
                                                 timedelta(seconds=AUTOFOCUS_DURATION)

            if quitat_dt is None:
                running_dt = action_expected_end_dt
                if action_expected_end_dt > an.ts_nosun.end:
                    action_new_status = 'SUN UP!'
                elif action_expected_end_dt > an.ts_dark.end:
                    action_new_status = 'twilight'
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
            new_action = this_action._replace(n_afinterval_autofocus=n_afinterval_autofocus,
                                              status=action_new_status,
                                              start_utc=action_start_dt)
            this_plan.action_list[i_action] = new_action
        new_plan = this_plan._replace(action_list=this_plan.action_list)
        plan_list[i_plan] = new_plan
    return plan_list


def add_altitudes(plan_list, an, fov_dict):
    for i_plan in range(len(plan_list)):
        this_plan  = plan_list[i_plan]
        for i_action in range(len(this_plan.action_list)):
            this_action = this_plan.action_list[i_action]
            if this_action.action_type.lower() in ['fov', 'stare', 'burn']:
                longitude, latitude = an.site.longitude, an.site.latitude
                if this_action.action_type.lower() == 'stare':
                    fov_name = this_action.parsed_action[2]
                else:
                    fov_name = this_action.parsed_action[1]
                longitude_hex, latitude_hex = degrees_as_hex(longitude), degrees_as_hex(latitude)
                if this_action.action_type.lower() == 'burn':
                    target_radec = RaDec(this_action.parsed_action[2],
                                         this_action.parsed_action[3])
                else:
                    this_fov = fov_dict[fov_name]
                    target_radec = RaDec(this_fov.ra, this_fov.dec)
                datetime_utc = this_action.start_utc
                az_deg, alt_deg = az_alt_at_datetime_utc(longitude_hex, latitude_hex,
                                                         target_radec, datetime_utc)
                new_action = this_action._replace(altitude_deg=alt_deg)
                this_plan.action_list[i_action] = new_action
        new_plan = this_plan._replace(action_list=this_plan.action_list)
        plan_list[i_plan] = new_plan
    return plan_list


def make_obs_exposure_data(fov_name, an, fov_dict=None, instrument=None):
    """
    Calculates exposure data for ONE REPEAT of given fov.
    :param fov_name:
    :param an:
    :param fov_dict:
    :param instrument: instrument data [Instrument object]
    :return: tuple of equal-length lists: (filters [str], counts [int], exp_times [float])
    """
    if fov_dict is not None:
        this_fov = fov_dict[fov_name]
    else:
        this_fov = Fov(fov_name)
    if not isinstance(instrument, Instrument):
        print("*** ERROR: make_obs_exposure_data() parm 'instrument' must be " +
              "a valid Instrument object")
        return None
    obs_style = this_fov.observing_style
    filters = []
    counts = []
    exp_times = []
    mags = dict()
    for obs in this_fov.observing_list:
        filter, mag, count = obs
        filters.append(filter)
        counts.append(count)
        if obs_style.lower() in ['standard', 'monitor', 'stare']:
            exp_time = calc_exp_time(mag, filter, instrument, this_fov)
        elif obs_style.lower() == 'lpv':
            if len(mags) == 0:
                mags = this_fov.estimate_lpv_mags(an.local_middark_jd)  # dict (get on 1st obs only)
            exp_time = calc_exp_time(mags[filter], filter, instrument, this_fov)
        else:
            print('****** WARNING: fov \'' + fov_name +
                  '\' has unrecognized observing style \'' + obs_style + '\'.')
            return None
        exp_times.append(exp_time)
    target_overhead = 75  # TODO: get this from instrument object, includes slew and guider start.
    filter_overhead = 15  # get this from instrument object, filter change time.
    exposure_overhead = 15  # get this from instrument object, ccd download & guider settle times.
    repeat_duration = len(filters) * filter_overhead + sum(counts) * exposure_overhead + \
                      sum([counts[i] * exp_times[i] for i in range(len(counts))])
    # types (3 lists, two floats): [str], [int], [float], float, float
    return filters, counts, exp_times, target_overhead, repeat_duration


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
    all_summary_lines = ['SUMMARY for AN ' + an.an_date_string,
                         'as generated by photrix at ' +
                         '{:%Y-%m-%d %H:%M  UTC}'.format(datetime.now(timezone.utc)),
                     an.acp_header_string()]
    for i_plan in range(len(plan_list)):
        this_plan = plan_list[i_plan]
        for i_action in range(len(this_plan.action_list)):
            this_action = this_plan.action_list[i_action]
            action_summary_lines = this_action.summary_lines.copy()  # modify copy, later replace.

            # Make summary line prefix parts:
            status_str = this_action.status.rjust(8)
            start_hhmm_str = hhmm_from_datetime_utc(this_action.start_utc)
            if this_action.status == 'SKIPPED':
                start_hhmm_str = len(start_hhmm_str) * ' '

            # Add '+' to hhmm if it actually refers to next day.
            next_day_string = ' '  # default
            if i_action >= 1:
                prev_start_utc = this_plan.action_list[i_action - 1].start_utc
            elif i_action == 0 and i_plan >= 1:
                prev_plan = plan_list[i_plan-1]
                prev_start_utc = prev_plan.action_list[len(prev_plan.action_list)-1].start_utc
            else:
                prev_start_utc = None
            if prev_start_utc is not None:
                if this_action.start_utc.date() > prev_start_utc.date():
                    next_day_string = '+'

            # Make degrees altitude prefix part.
            if this_action.action_type.lower() in ['fov', 'stare', 'burn']:
                alt_string = '{:2d}'.format(int(round(this_action.altitude_deg)))
            else:
                alt_string = '  '

            # Error if plan chains to itself, or if chained-to plan file does not exist:
            if this_action.action_type.lower() == 'chain':
                chain_target_filename = this_action.parsed_action[1].lower()
                if chain_target_filename == 'plan_' + this_plan.plan_id.lower() + '.txt':
                    action_summary_lines += ['*** ERROR: Plan attempts to chain to itself.']
                else:
                    target_plan_exists = False
                    for plan in plan_list:
                        this_plan_filename = 'plan_' + plan.plan_id.lower() + '.txt'
                        if this_plan_filename == chain_target_filename:
                            target_plan_exists = True
                    if not target_plan_exists:
                        action_summary_lines += ['*** ERROR: Chain-to plan does not exist.']

            # Construct and add prefix, add to summary lines:
            prefix = status_str + ' ' + start_hhmm_str + next_day_string + ' ' + alt_string + ' '
            lines_without_prefixes = len(action_summary_lines) - 1
            if this_action.action_type.lower() == 'plan':
                line_prefixes = lines_without_prefixes * [len(prefix) * ' '] + [prefix]  # last line
            else:
                line_prefixes = [prefix] + lines_without_prefixes * [len(prefix) * ' ']  # 1st line
            prefixed_summary_lines = [line_prefixes[i] + action_summary_lines[i]
                                      for i in range(len(action_summary_lines))]
            all_summary_lines.extend(prefixed_summary_lines)
    # Write Summary file:
    output_fullpath = os.path.join(output_directory, 'Summary_' + an.an_date_string + '.txt')
    print('PRINT all summary lines to ', output_fullpath)
    with open(output_fullpath, 'w') as this_file:
        this_file.write('\n'.join(all_summary_lines))


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
                  '#BINNING ' + ','.join(len(filters)*['1']) + ' ;',
                  '#COUNT ' + ','.join(counts) + ' ;',
                  '#INTERVAL ' + ','.join(exp_times) + ' ; ' + duration_comment,
                  ';----' + this_fov.acp_comments,
                  fov_name + '\t' + ra_as_hours(this_fov.ra) + '\t' + dec_as_hex(this_fov.dec)]
    if num_repeats > 1:
        acp_string.insert(1, '#REPEAT ' + str(num_repeats) + ' ;')
    return acp_string, total_duration


def calc_exp_time(mag, filter, instrument, fov):
    if filter == 'V11.7':
        iii = 4
    exp_time_from_mag = instrument.filters[filter]['reference_exposure_mag10'] *\
                        10.0 ** ((mag - 10.0) / 2.5)
    if fov.max_exposure is None:
        exp_time_asymptote = ABSOLUTE_MAX_EXPOSURE_TIME
    else:
        exp_time_asymptote = min(ABSOLUTE_MAX_EXPOSURE_TIME, fov.max_exposure)
    exp_time = 1.0 / (1.0 / exp_time_from_mag + 1.0 / exp_time_asymptote)
    return exp_time
