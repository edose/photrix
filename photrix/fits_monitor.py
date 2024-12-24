__author__ = "Eric Dose, Albuquerque"

""" This module: 
      
"""

# Python core:
import os
from time import sleep
import winsound
from datetime import datetime, timezone, timedelta
# from typing import Union, List, Tuple, Set
import sys


# External packages:

# Author's packages:
from util import hhmm_from_datetime_utc


THIS_PACKAGE_ROOT_DIRECTORY = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
INI_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, 'ini')

FITS_MASTER_PATH = 'C:\\Users\\Eric\\Dropbox\\Family Room\\Images'
SECONDS_BETWEEN_QUERIES = 60
MAX_DURATION_NO_NEW_FITS = 900  # Alert when no new file has appeared for this long.
UTC_OFFSET_AT_SITE = -7
UTC_OFFSET_OF_AN = 12 - UTC_OFFSET_AT_SITE  # UTC minus (0h AN of matching date).

SOUND_REPETITIONS_ON_NO_NEW_FITS = 20
# Closing: Win 10 Critical Stop; use Sonata/Windows Critical Stop.wav
SOUND_ON_NO_NEW_FITS = 'SystemHand'


class NoNewFitsError(Exception):
    pass


def play_sound_alias(sound: str, count: int):
    for i in range(count):
        winsound.PlaySound(sound, winsound.SND_ALIAS)


# noinspection DuplicatedCode
def monitor_fits_files():
    # Start and stop times for this night (from batch file calling this py file):
    n_args = len(sys.argv)
    if n_args != 3:
        raise ValueError('FITS monitor requires exactly 2 arguments but ' +
                         str(n_args) + ' given.')

    # COMPUTE Astronight date & FITS path:

    now_utc = datetime.now(timezone.utc)
    now_an_dt = now_utc - timedelta(hours=UTC_OFFSET_OF_AN)
    now_an = 10000 * now_an_dt.year + 100 * now_an_dt.month + now_an_dt.day
    now_an_str = str(now_an)
    path = os.path.join(FITS_MASTER_PATH, now_an_str)
    if not os.path.exists(path):
        os.makedirs(path)
    print('Reading', path)

    # COMPUTE monitoring start and stop times in UTC:
    start_hour_utc_arg = int(sys.argv[1])
    stop_hour_utc_arg = int(sys.argv[2])
    an_start_utc = datetime(now_an_dt.year, now_an_dt.month, now_an_dt.day,
                            0, 0, 0).replace(tzinfo=timezone.utc) + \
        timedelta(hours=UTC_OFFSET_OF_AN)
    hours_start_after_an_start = (start_hour_utc_arg - UTC_OFFSET_OF_AN) % 24
    start_utc = an_start_utc + timedelta(hours=hours_start_after_an_start)
    hhmm_start = hhmm_from_datetime_utc(start_utc)
    hours_stop_after_an_start = (stop_hour_utc_arg - UTC_OFFSET_OF_AN) % 24
    stop_utc = an_start_utc + timedelta(hours=hours_stop_after_an_start)

    # Prepare monitoring loop:
    ff_prev = FitsFilenames(path)
    print('Will start at', start_utc.strftime('%Y-%m-%d %H:%M:%S'))
    print('Will alert after', int(MAX_DURATION_NO_NEW_FITS),
          'seconds without new FITS.')
    print('Will stop at', stop_utc.strftime('%Y-%m-%d %H:%M:%S'))
    print()
    last_new_fits_utc = None

    # MONITORING LOOP:
    while True:
        now_utc = datetime.now(timezone.utc)
        hhmm_now = hhmm_from_datetime_utc(now_utc)

        # CASE: Start time has not yet been reached -> wait.
        if now_utc < start_utc:
            print(hhmm_now, 'waiting to start at', hhmm_start)
            seconds_before_start = (start_utc - now_utc).total_seconds()
            sleep(min(5 * SECONDS_BETWEEN_QUERIES, int(seconds_before_start)) + 1)
            continue

        # CASE: Stop time has been reached -> Stop and exit.
        if now_utc >= stop_utc:
            print('\n' + hhmm_now, 'STOPPED normally.')
            break

        if last_new_fits_utc is None:
            # If the first real iteration of monitoring loop:
            last_new_fits_utc = now_utc
            print(hhmm_now, 'Start monitoring')

        ff = FitsFilenames(path)
        if ff.new_fits_found_vs(ff_prev):
            # CASE: New FITS file found -> update.
            print('\n' + hhmm_from_datetime_utc(ff.dt_update),
                  ' new FITS file found', end='', flush=True)
            last_new_fits_utc = ff.dt_update
        else:
            # CASE: No new FITS file -> stop if system seems frozen, else skip.
            since_last_new_fits = (ff.dt_update - last_new_fits_utc).total_seconds()
            print('.', end='', flush=True)
            if since_last_new_fits > MAX_DURATION_NO_NEW_FITS:
                print('\n' + hhmm_now, ' ***** SYSTEM APPEARS FROZEN *****')
                print('for', since_last_new_fits, 'seconds (max=',
                      MAX_DURATION_NO_NEW_FITS, ')')
                play_sound_alias(SOUND_ON_NO_NEW_FITS, SOUND_REPETITIONS_ON_NO_NEW_FITS)
                raise NoNewFitsError()
        ff_prev = ff
        sleep(SECONDS_BETWEEN_QUERIES)


# def astronight_code(dt_now: datetime = None) -> str:
#     if dt_now is None:
#         dt_now = datetime.now(timezone.utc)
#     subtract_one_day = (dt_now.hour < HOUR_UTC_NEW_AN)
#     dt_an = dt_now + (timedelta(days=-1) if subtract_one_day else timedelta(days=0))
#     return f'{dt_an.year:04}' + f'{dt_an.month:02}' + f'{dt_an.day:02}'
#
#
# def make_fits_path() -> str:
#     """ Constructs path for directory holding new FITS files. """
#     return os.path.join(FITS_MASTER_PATH, astronight_code())


class FitsFilenames:
    def __init__(self, fits_path):
        self.fits_path = fits_path
        self.fits_filenames = None
        self.dt_update = None
        self.update()

    def update(self):
        """ Returns list of all FITS filenames in directory. """
        self.fits_filenames = [fname for fname in os.listdir(self.fits_path)
                               if os.path.isfile(os.path.join(self.fits_path, fname))
                               and (fname.endswith((".fts", ".fit", ".fts")))
                               # and (fname.startswith("MP_"))
                               ]
        self.dt_update = datetime.now(timezone.utc)

    def new_fits_found_vs(self, other) -> bool:
        """ Return True iff any FITS exists in other list that is not in this list. """
        this_set = set(self.fits_filenames)
        other_set = set(other.fits_filenames)
        new_fits = this_set - other_set
        return len(new_fits) >= 1


if __name__ == '__main__':
    # n = len(sys.argv)
    # for i in range(1, n):
    #     print(i, sys.argv[i])
    monitor_fits_files()
