__author__ = "Eric Dose :: New Mexico Mira Project, Albuquerque"

import winsound
from time import sleep
from datetime import datetime, timezone

import requests
# from bs4 import BeautifulSoup

from util import hhmm_from_datetime_utc

STATUS_URL = 'http://deepskyportal.com/weather/BetaRoofStatusFile.txt'
STATUS_PATH_FOR_TEST = 'C:/24hrs/RoofTest.txt'

HTTP_OK_CODE = 200  # "OK. The request has succeeded."
MAX_CONSECUTIVE_TIMEOUTS = 10
MAX_CONSECUTIVE_ERRORS = 10
SECONDS_BETWEEN_CONSECUTIVE_TIMEOUTS = 10  # retry cadence.
SECONDS_BETWEEN_STARTUP_QUERIES = 5        # initial polling cadence.
SECONDS_BETWEEN_ONGOING_QUERIES = 60       # normal polling cadence.
SECONDS_WHEN_CHANGE_SENSED = 5            # fast cadence when roof status change has been sensed.
CONSISTENT_CHANGE_QUERIES_REQUIRED = 5     # number of uniform query results to confirm status change.
TAG_TO_TEST = 'Roof Status:'

SOUND_HAS_OPENED = 'SystemAsterisk'  # = Win 10 Asterisk; use Sonata/Windows Error.wav
SOUND_HAS_CLOSED = 'SystemHand'      # = Win 10 Critical Stop; use Sonata/Windows Critical Stop.wav
SOUND_REPETITIONS = 50


def monitor_roof(url=STATUS_URL):
    """ Make noise if roof opens or closes, according to very small roof status file from Deep Sky West"""
    print('Playing OPENED sound twice, then CLOSED sound twice...')
    winsound.PlaySound(SOUND_HAS_OPENED, winsound.SND_ALIAS)
    winsound.PlaySound(SOUND_HAS_OPENED, winsound.SND_ALIAS)
    winsound.PlaySound(SOUND_HAS_CLOSED, winsound.SND_ALIAS)
    winsound.PlaySound(SOUND_HAS_CLOSED, winsound.SND_ALIAS)

    status_list = (CONSISTENT_CHANGE_QUERIES_REQUIRED + 1) * ['-']  # start without info (as dashes).
    last_event_string = ''
    while True:
        status_text = get_status_text(url)
        if not status_text.startswith('OK'):
            print(' >>>>> ERROR=', status_text, '--> STOPPING.')
            exit(0)
        immed_status = parse_immed_status(status_text)
        if immed_status == 'ERROR':
            print(' >>>>> ERROR=', immed_status, '--> STOPPING.')
            exit(0)
        status_list.append(immed_status)  # so that earliest status is first list item (time reads L->R).
        status_list = status_list[1:]     # pop earliest status off front of list.
        hhmm = hhmm_from_datetime_utc(datetime.now(timezone.utc))
        print(hhmm + ': is', status_list[-1] + last_event_string)

        if roof_has_opened(status_list):
            print(32 * '*', '\n >>>>> OPENED at', hhmm)
            last_event_string = ' (since ' + hhmm + ')'
            for i in range(SOUND_REPETITIONS):
                winsound.PlaySound(SOUND_HAS_OPENED, winsound.SND_ALIAS)
        elif roof_has_closed(status_list):
            print(32 * '*', '\n >>>>> CLOSED at', hhmm)
            last_event_string = ' (since ' + hhmm + ')'
            for i in range(SOUND_REPETITIONS):
                winsound.PlaySound(SOUND_HAS_CLOSED, winsound.SND_ALIAS)
        else:
            if change_is_sensed(status_list):
                sleep(SECONDS_WHEN_CHANGE_SENSED)
            elif any([status == '-' for status in status_list]):
                sleep(SECONDS_BETWEEN_STARTUP_QUERIES)
            else:
                sleep(SECONDS_BETWEEN_ONGOING_QUERIES)


def get_status_text(url):
    """ Get and return status text from roof status url.
    :param url: URL holding updated roof status in Deep Sky West format. [string]
    :return status_text = status text prepended by 'OK', 'ERROR', or 'TIMEOUT' (only). [string]
    """
    if url.lower() == 'test':
        with open(STATUS_PATH_FOR_TEST) as f:
            return 'OK' + f.readlines()[0]
    for n_timeouts in range(MAX_CONSECUTIVE_TIMEOUTS):
        try:
            r = requests.get(url)
        except requests.exceptions.Timeout:
            print(' >>>>> Warning:', str(n_timeouts), 'consecutive timeouts.')
            sleep(SECONDS_BETWEEN_CONSECUTIVE_TIMEOUTS)
            continue
        except requests.exceptions.RequestException as e:
            print(e)
            return 'ERROR'
        if r.status_code != HTTP_OK_CODE:
            print(' >>>>> Could not get', STATUS_URL)
            return 'ERROR'
        return 'OK' + r.text
    print(' >>>>> ERROR:', str(MAX_CONSECUTIVE_TIMEOUTS), 'consecutive timeouts.')
    return 'TIMEOUT'


def parse_immed_status(status_text):
    """ Parses immediate status from text and returns (only) 'open', 'closed', or 'error'. [string] """
    core_text = status_text.split(TAG_TO_TEST)[1].strip().upper()
    if core_text.startswith('OPEN'):
        return 'open'
    elif core_text.startswith('CLOSE'):
        return 'closed'
    print(' >>>>> ERROR: cannot parse status text >' + status_text + '<')
    return 'error'


def roof_has_opened(status_list):
    return all([s == 'closed' for s in status_list[:2]]) and all([s == 'open' for s in status_list[2:]])


def roof_has_closed(status_list):
    return all([s == 'open' for s in status_list[:2]]) and all([s == 'closed' for s in status_list[2:]])


def change_is_sensed(status_list):
    return ('open' in status_list) and ('closed' in status_list)


if __name__ == '__main__':
    monitor_roof()
