__author__ = "Eric Dose, Albuquerque"

""" This module: 
      
"""

# Python core:
import os
from time import sleep
import winsound
from datetime import datetime, timezone, timedelta
from typing import Union
import sys

# Author's packages:
from util import hhmm_from_datetime_utc

THIS_PACKAGE_ROOT_DIRECTORY = \
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INI_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, 'ini')

HCRO_MONITOR_FULLPATH = 'C:\\Users\\Eric\\Dropbox\\Family Room/Roof\\roof.log'
SECONDS_BETWEEN_FILE_QUERIES = 30  # Set this short--queries are local, thus cheap.
MAX_DURATION_DISCONNECT = 720  # Alert when disconnected, etc. for this long.
MAX_BETWEEN_ROOF_STATUS = 180  # reads may normally go unchanged for this long.
# HOUR_UTC_TO_STOP = 12  # Typically, time when monitor no longer needed, not at dawn.

SOUND_REPETITIONS_ON_STATUS_CHANGE = 40
# Opening: Win 10 Asterisk; use Sonata/Windows Error.wav
SOUND_ON_OPENING = 'SystemAsterisk'
# Closing: Win 10 Critical Stop; use Sonata/Windows Critical Stop.wav
SOUND_ON_CLOSING = 'SystemHand'


def play_sound_alias(sound: str, count: int):
    for i in range(count):
        winsound.PlaySound(sound, winsound.SND_ALIAS)


# def play_sound_file(sound_file: str, count: int):
#     for i in range(count):
#         winsound.PlaySound(sound_file, winsound.SND_FILENAME|winsound.SND_NOWAIT)


class RoofStatusError(Exception):
    pass


class OddStatusWhileConnected(Exception):
    pass


class OddStatusWhileDisconnected(Exception):
    pass


class DisconnectedTooLong(Exception):
    pass


__________ASCOM_ACCESS________________________________________________________ = 0


def monitor_hcro_status_via_ascom():
    """ Make noise if HCRO Safety Monitor goes safe or unsafe, or is unknown status,
        according to HCRO Kronos safety monitor via ASCOM.
    """
    print('Playing OPENED sound twice, then CLOSED sound twice...')
    play_sound_alias(SOUND_ON_OPENING, 2)
    play_sound_alias(SOUND_ON_CLOSING, 2)
    print('HCRO roof status')
    dt_now = datetime.now(timezone.utc)

    # Start and stop times for this night (from batch file calling this py file):
    n_args = len(sys.argv)
    if n_args != 2:
        raise ValueError('FITS monitor requires exactly 1 argument but ' +
                         str(n_args) + ' given.')
    stop_hour_utc_arg = int(sys.argv[1])
    dt_to_stop = datetime(dt_now.year, dt_now.month, dt_now.day,
                          stop_hour_utc_arg, 0, 0).replace(tzinfo=timezone.utc) + \
        (timedelta(days=1) if dt_now.hour > stop_hour_utc_arg else timedelta(days=0))

    mon = RoofStatus(HCRO_MONITOR_FULLPATH)
    dt_previous_iteration = None
    is_first_iteration = True

    # Status Loop:
    while True:
        dt_now = datetime.now(timezone.utc)
        hhmm_now = hhmm_from_datetime_utc(dt_now)

        # CASE: Stop time has been reached -> Stop and exit.
        if dt_now >= dt_to_stop:
            print(hhmm_now, 'STOPPED normally.')
            break

        roof_status_change = mon.update()
        # print(mon.dump())
        # print()
        is_good_connection = (mon.status_latest_read_attempt in ('open', 'closed'))
        # is_unchanged_roof_status = not (roof_status_change in ('opened', 'closed'))
        hhmm_read = hhmm_from_datetime_utc(mon.dt_latest_read_attempt)
        if mon.dt_latest_valid_roof_status_reading is None:
            seconds_since_last_roof_status = None
        else:
            seconds_since_last_roof_status = \
                (dt_now - mon.dt_latest_valid_roof_status_reading).total_seconds()

        if is_good_connection:
            # MONITOR is CONNECTED on last attempt (normal case): ======================
            if is_first_iteration:
                print(hhmm_read + '   ' + str(mon.current_roof_status).upper())
            elif seconds_since_last_roof_status > MAX_BETWEEN_ROOF_STATUS:
                # WARNING (only): roof status is stale:
                print(' >>>>> roof presumed ' +
                      str(mon.current_roof_status).upper() + ' / ' +
                      '{0:.0f}'.format(seconds_since_last_roof_status) +
                      ' seconds since last valid update.')
            elif mon.dt_latest_read_attempt == dt_previous_iteration:
                # SKIP (no-op) CASE: Unchanged & not timed-out-> skip iteration:
                # print('# skip case at ', hhmm_now)
                pass
            elif roof_status_change == 'opened':
                print(32 * '*', '\n >>>>> OPENED at', hhmm_read)
                play_sound_alias(SOUND_ON_OPENING, SOUND_REPETITIONS_ON_STATUS_CHANGE)
            elif roof_status_change == 'closed':
                print(32 * '*', '\n >>>>> CLOSED at', hhmm_read)
                play_sound_alias(SOUND_ON_CLOSING, SOUND_REPETITIONS_ON_STATUS_CHANGE)
            elif roof_status_change is None:
                # NORMAL CASE: fresh reading, and roof status has not changed:
                print(hhmm_read, '  ' + str(mon.current_roof_status).upper(),
                      '     (since ' +
                      hhmm_from_datetime_utc(mon.dt_current_roof_status_began) + ')')
            else:
                play_sound_alias(SOUND_ON_CLOSING, 4)
                mon.dump()
                raise OddStatusWhileConnected()

        else:
            # MONITOR is DISCONNECTED (or file badly parsed, etc). =====================

            # EXIT CASE: Monitor has been disconnected longer than allowed:
            if mon.duration_current_disconnect > MAX_DURATION_DISCONNECT:
                play_sound_alias(SOUND_ON_CLOSING, 4)
                mon.dump()
                msg = 'Disconnected for ' + \
                      '{0:.0f}'.format(mon.duration_current_disconnect) + \
                      ' seconds vs max allowed = ' + str(MAX_DURATION_DISCONNECT) + '.'
                raise DisconnectedTooLong(msg)

            # NEVER YET CONNECTED CASE:
            elif mon.dt_current_roof_status_began is None:
                print(hhmm_now + '   no connection yet made...')

            # NORMAL DISCONNECTED CASE: roof was open or closed on latest disconnect,
            #   and presumed to remain in that roof state:
            elif mon.current_roof_status in ['open', 'closed']:
                print(hhmm_now + '   ' + 'Presumed',
                      str(mon.current_roof_status).upper() +
                      '  (since ' +
                      hhmm_from_datetime_utc(mon.dt_current_roof_status_began) + ')' +
                      ' but disconnected since ' +
                      hhmm_from_datetime_utc(mon.dt_current_disconnect_began))

            # ERROR CASE: should never get here.
            else:
                play_sound_alias(SOUND_ON_CLOSING, 4)
                mon.dump()
                raise OddStatusWhileDisconnected()

        dt_previous_iteration = mon.dt_latest_read_attempt
        is_first_iteration = False
        # print('# sleep', str(SECONDS_BETWEEN_FILE_QUERIES), 'sec.')
        sleep(SECONDS_BETWEEN_FILE_QUERIES)


class RoofStatus:
    """ Holds, updates, and reports on current roof status, as obtained from
        a one-line file updated elsewhere."""
    def __init__(self, fullpath: str):
        self.fullpath = fullpath
        # Most recently parsed dt and status:
        self.dt_latest_read_attempt = None
        self.status_latest_read_attempt = None
        # Current roof status:
        self.dt_latest_valid_roof_status_reading = None
        self.dt_current_roof_status_began = None
        self.current_roof_status = None
        # When current disconnect was first detected (or None):
        self.dt_current_disconnect_began = None

    def update(self) -> Union[str, None]:
        """  Interprets all state data with results of latest read attempt.
        :return: action ['opened', 'closed', None]
        """
        self.dt_latest_read_attempt, self.status_latest_read_attempt = \
            self._read_and_parse()
        # print('@', utc_string(datetime.now(timezone.utc)) + ', read:  ',
        #       utc_string(self.dt_latest_read_attempt), self.status_latest_read_attempt)

        if self.status_latest_read_attempt == 'open':
            self.dt_current_disconnect_began = None
            self.dt_latest_valid_roof_status_reading = self.dt_latest_read_attempt
            # CASE: roof remains open:
            if self.current_roof_status == 'open':
                return None
            # CASE: roof has just opened:
            if self.current_roof_status == 'closed':
                self.dt_current_roof_status_began = self.dt_latest_read_attempt
                self.current_roof_status = self.status_latest_read_attempt
                return 'opened'
            # CASE: this is the first known roof status:
            if self.current_roof_status is None:
                self.dt_current_roof_status_began = self.dt_latest_read_attempt
                self.current_roof_status = self.status_latest_read_attempt
                return None

        if self.status_latest_read_attempt == 'closed':
            self.dt_current_disconnect_began = None
            self.dt_latest_valid_roof_status_reading = self.dt_latest_read_attempt
            # CASE: roof remains closed:
            if self.current_roof_status == 'closed':
                return None
            # CASE: roof has just closed:
            if self.current_roof_status == 'open':
                self.dt_current_roof_status_began = self.dt_latest_read_attempt
                self.current_roof_status = self.status_latest_read_attempt
                return 'closed'
            # CASE: this is the first known roof status:
            if self.current_roof_status is None:
                self.dt_current_roof_status_began = self.dt_latest_read_attempt
                self.dt_latest_valid_roof_status_reading = self.dt_latest_read_attempt
                self.current_roof_status = self.status_latest_read_attempt
                return None

        if self.status_latest_read_attempt is None:
            # CASE: this is initial detection of disconnection or failed read:
            if self.dt_current_disconnect_began is None:
                self.dt_current_disconnect_began = self.dt_latest_read_attempt
            return None

        # ERROR CASE (sentinel): Should never reach here.
        print('RoofStatus.update() ERROR. Dump follows.')
        self.dump()
        raise RoofStatusError()

    def _read_and_parse(self) -> (datetime, Union[str, None]):
        """ Read and parse ONLY the one-line file as it exists at this moment.
            Does not attempt to interpret it, keeps no history whatever.
        Typical line to parse: "02:57:28 UTC     connected open".
        Perfect read returns: time_read, 'open' or 'closed'.
        Any failure to read or parse returns: time_attempt, None.
        :return: 2-tuple (dt: datetime, status: string | None),
            where dt (UTC) is from file if parsable, otherwise is Now; and
            where status must be in ['open', 'closed', None].
        """

        dt_attempt = datetime.now(timezone.utc)
        try:
            f = open(self.fullpath, 'r')
        # CASE: Cannot open file -> time_now, None.
        except OSError:
            return dt_attempt, None

        with f:
            lines = f.readlines()
        # CASE: Cannot parse (number of lines not 1) -> time_now, None.
        if len(lines) < 1:
            return dt_attempt, None

        line = lines[0]
        words = line.split('UTC')

        # CASE: Cannot parse (wrong 1-line format) -> time_now, None.
        if len(words) != 2:
            return dt_attempt, None

        try:
            dt_parsed = \
                datetime.strptime(words[0].strip(),
                                  '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)
        # CASE: Cannot parse (datetime in wrong format) -> time_now, None.
        except ValueError:
            return dt_attempt, None

        words = words[1].replace(',', ' ').split(maxsplit=1)
        connected_string = words[0].strip().lower()

        # CASE: Cannot parse for whatever reason -> time_now, None.
        if connected_string != 'connected':
            return dt_attempt, None
        status_parsed = words[1].strip().lower()
        if status_parsed not in ['open', 'closed']:
            return dt_attempt, None

        # NORMAL CASE: return parsed dt and status.
        return dt_parsed, status_parsed

    def dump(self) -> None:
        print('RoofStatus dump at:',
              utc_string(datetime.now(timezone.utc)))
        print('    .dt_latest_read_attempt:',
              utc_string(self.dt_latest_read_attempt))
        print('    .status_latest_read_attempt:',  self.status_latest_read_attempt)
        print('    .dt_current_roof_status_began:',
              utc_string(self.dt_current_roof_status_began))
        print('    .current_roof_status:', self.current_roof_status)
        print('    .dt_latest_valid_roof_status_reading: ',
              utc_string(self.dt_latest_valid_roof_status_reading))
        print('    .dt_current_disconnect_began:',
              utc_string(self.dt_current_disconnect_began))
        print()

    @property
    def duration_current_disconnect(self) -> float:
        if self.dt_current_disconnect_began is None:
            return 0
        return (self.dt_latest_read_attempt -
                self.dt_current_disconnect_began).total_seconds()


def utc_string(dt: datetime) -> str:
    return dt.strftime('%Y-%m-%d %H:%M:%S UTC') if dt is not None else 'None'


if __name__ == '__main__':
    monitor_hcro_status_via_ascom()
