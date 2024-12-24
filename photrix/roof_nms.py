__author__ = "Eric Dose, Albuquerque"

""" This module: 
      
"""

# Python core:
import os
from time import sleep
import winsound
from collections import defaultdict
from datetime import datetime, timezone

# External packages:
import requests
from bs4 import BeautifulSoup

# Author's packages:
from util import hhmm_from_datetime_utc
from imap import get_most_recent_relevant_email


THIS_PACKAGE_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INI_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, 'ini')

NMS_STATUS_URL = 'https://nmskies.com/weather.php'
NMS_STATUS_DICT = defaultdict(lambda: 'unknown',  # default if image name is unknown.
                              {'daylight.jpg': 'closed',
                               'clouds.jpg': 'closed',
                               'snow.jpg': 'closed',
                               'fog.jpg': 'closed',
                               'wind.jpg': 'closed',
                               # 'gust.jpg': 'closed',
                               'rain.jpg': 'closed',
                               'open.jpg': 'open'})

# HTML CONSTANTS:
HTTP_OK_CODE = 200  # "OK. The request has succeeded."
SECONDS_BETWEEN_HTML_QUERIES = 120
MAX_CONSECUTIVE_HTML_TIMEOUTS = 3
SECONDS_BETWEEN_CONSECUTIVE_HTML_TIMEOUTS = 20  # retry cadence.
STR_LEFT_OF_IMAGE_NAME = 'images/'
STR_RIGHT_OF_IMAGE_NAME = '?image='

# E-MAIL CONSTANTS:
SECONDS_BETWEEN_EMAIL_QUERIES = 180
NMS_WEATHER_SUBJECT_START = '[NMSWX]'
NMS_TECHSUPPORT_FROM_STRINGS = ['tech', '@nmskies.com']
NMS_OPEN = ['OK to open']
NMS_CLOSED = ['opening delayed', 'closed due to']
NMS_EMAIL_SOUND = 'SystemExclamation'


SOUND_REPETITIONS_ON_STATUS_CHANGE = 40
SOUND_ON_OPENING = 'SystemAsterisk'  # = Win 10 Asterisk; use Sonata/Windows Error.wav
SOUND_ON_CLOSING = 'SystemHand'      # = Win 10 Critical Stop; use Sonata/Windows Critical Stop.wav
SOUND_REPETITIONS_ON_NEW_EMAIL_ONLY = 1


def play_sound_alias(sound, count):
    for i in range(count):
        winsound.PlaySound(sound, winsound.SND_ALIAS)


def play_sound_file(sound_file, count):
    for i in range(count):
        winsound.PlaySound(sound_file, winsound.SND_FILENAME|winsound.SND_NOWAIT)

__________HTML_ACCESS________________________________________________________ = 0


def get_nms_weather_html_request(url):
    """ Get and return HTML (request object) representing current NMS Weather webpage. """
    r = None  # default if r never set.
    for n_timeouts in range(MAX_CONSECUTIVE_HTML_TIMEOUTS):
        try:
            r = requests.get(url)
        except requests.exceptions.Timeout:
            print(' >>>>> Warning:', str(n_timeouts), 'consecutive timeouts.')
            sleep(SECONDS_BETWEEN_CONSECUTIVE_HTML_TIMEOUTS)
            continue
        except requests.exceptions.RequestException as e:
            print(e)
            return 'ERROR: RequestException (URL not found).'
        if r.status_code != HTTP_OK_CODE:
            print(' >>>>> Could not get', NMS_STATUS_URL)
        return r
    return r


def get_nms_status_image_name():
    """ Extract image name (e.g., 'daylight.jpg') from HTML of NMS Weather webpage. """
    image_name = ''
    r = get_nms_weather_html_request(NMS_STATUS_URL)
    if r.status_code == HTTP_OK_CODE:
        soup = BeautifulSoup(r.text, 'html.parser')
        img_tags = soup.find_all('div', class_='img')  # NB: "class_" not "class" (reserved).
        s = [x for x in img_tags][0].contents[1]['src']  # source code of page's very first image.
        i_left = s.find(STR_LEFT_OF_IMAGE_NAME, 0)
        i_right = s.find(STR_RIGHT_OF_IMAGE_NAME, len(STR_LEFT_OF_IMAGE_NAME))
        if i_left >= 0 and i_right >= len(STR_LEFT_OF_IMAGE_NAME):
            image_name = s[i_left + len(STR_RIGHT_OF_IMAGE_NAME):i_right].strip()
    return image_name


def monitor_nms_status_via_html():
    """ Make noise if NMS *observatory* (not just a roof) opens or closes,
        according to New Mexico Skies' weather web page.
        (logic here is simpler than for Deep Sky West polling, as NMS web page needs no fault tolerance.)
    """
    print('Playing OPENED sound twice, then CLOSED sound twice...')
    play_sound_alias(SOUND_ON_OPENING, 2)
    play_sound_alias(SOUND_ON_CLOSING, 2)

    # STATUS LOOP:
    previous_status = None
    last_event_string = ''
    while True:
        # status = get_nms_status_from_image_name()
        image_name = get_nms_status_image_name()
        status = NMS_STATUS_DICT[image_name]
        hhmm = hhmm_from_datetime_utc(datetime.now(timezone.utc))
        print(hhmm + ': NMS is', status.upper(), last_event_string)
        if previous_status is not None:
            if status == 'open' and previous_status == 'closed':
                print(32 * '*', '\n >>>>> OPENED at', hhmm)
                last_event_string = ' (since ' + hhmm + ')'
                play_sound_alias(SOUND_ON_OPENING, SOUND_REPETITIONS_ON_STATUS_CHANGE)
            elif status == 'closed' and previous_status == 'open':
                print(32 * '*', '\n >>>>> CLOSED at', hhmm)
                last_event_string = ' (since ' + hhmm + ')'
                play_sound_alias(SOUND_ON_CLOSING, SOUND_REPETITIONS_ON_STATUS_CHANGE)
            elif status not in ['closed', 'open']:  # only get here if there is a problem.
                play_sound_alias(SOUND_ON_CLOSING, 4)
                print('STOPPING ON ERROR: status is >' + status + '<' +
                      ' //  status_image_name is >' + image_name + '<')
                break
        previous_status = status
        sleep(SECONDS_BETWEEN_HTML_QUERIES)


__________E_MAIL_ACCESS________________________________________________________ = 0


def test_monitor():
    print()
    wx_dict = get_most_recent_relevant_email(subject_start=NMS_WEATHER_SUBJECT_START,
                                             from_strings=NMS_TECHSUPPORT_FROM_STRINGS)
    print('WX:', wx_dict)
    other_dict = get_most_recent_relevant_email(subject_start=None,
                                                from_strings=NMS_TECHSUPPORT_FROM_STRINGS)
    print('Other:', other_dict)


def monitor_nms_status_via_e_mail():
    """ Make specific noise whenever NMS techsupport sends an e-mail."""
    print('Playing OPENED sound twice, then CLOSED sound twice...')
    play_sound_alias(NMS_EMAIL_SOUND, 4)
    wx_dict = get_most_recent_relevant_email(subject_start=NMS_WEATHER_SUBJECT_START,
                                             from_strings=NMS_TECHSUPPORT_FROM_STRINGS)
    other_dict = get_most_recent_relevant_email(subject_start=None,
                                                from_strings=NMS_TECHSUPPORT_FROM_STRINGS)
    if other_dict == wx_dict:
        other_dict = None
    most_recent_wx_dict = wx_dict
    most_recent_other_dict = other_dict

    # Display most recent WX subject, if any WX message exists:
    if wx_dict is not None:
        print('\nMOST RECENT WEATHER E-MAIL at', wx_dict['date'].strip() +
              '\n' + wx_dict['subject'])
    else:
        print('\nNo most recent weather e-mail found.')

    # Display most recent Other subject, if any non-WX message exists:
    if other_dict != wx_dict:
        if other_dict is not None:
            print('\nMOST RECENT NON-WEATHER E-MAIL at', other_dict['date'].strip() +
                  '\n' + other_dict['subject'])

    # Monitoring loop:
    while True:
        sleep(SECONDS_BETWEEN_EMAIL_QUERIES)
        wx_dict = get_most_recent_relevant_email(subject_start=NMS_WEATHER_SUBJECT_START,
                                                 from_strings=NMS_TECHSUPPORT_FROM_STRINGS)
        other_dict = get_most_recent_relevant_email(subject_start=None,
                                                    from_strings=NMS_TECHSUPPORT_FROM_STRINGS)

        # If new weather e-mail, display subject and sound long audible alarm:
        if wx_dict is not None and wx_dict != most_recent_wx_dict:
            print('\nNEW WEATHER E-MAIL at', wx_dict['date'].strip() +
                  '\n' + wx_dict['subject'])
            play_sound_alias(NMS_EMAIL_SOUND, 20)
            most_recent_wx_dict = wx_dict.copy()

        # If new non-weather e-mail, display subject and sound short audible alarm:
        if other_dict is not None and other_dict != wx_dict and other_dict != most_recent_other_dict:

            print('\nNEW NON-WEATHER E-MAIL at', other_dict['date'].strip() +
                  '\n' + other_dict['subject'])
            play_sound_alias(NMS_EMAIL_SOUND, 5)
            most_recent_other_dict = other_dict.copy()


# def extract_nms_status_from_email(msg_dict):
#     """Return status string as indicated by subject line of NMS status e-mail; or None if no e-mail found.
#     :param msg_dict: dict of e-mail fields from probable NMS status e-mail. [py dict]
#     :return: nms_status_string, or None if no message was passed in. [string]
#     """
#     if msg_dict is None:
#         return None
#     for open_string in NMS_OPEN:
#         if open_string.lower() in msg_dict['subject'].lower():
#             return 'open', msg_dict['subject']
#     for closed_string in NMS_CLOSED:
#         if closed_string.lower() in msg_dict['subject'].lower():
#             return 'closed', msg_dict['subject']
#     raise StatusEmailParsingError('Subject line >' + msg_dict['subject'] +
#                                   '< cannot be parsed as open or closed.')


__________OTHER_FUNCTIONS___________________________________________________ = 0


def record_nms_status_image_name():
    """ Repeatedly get NMS observatory status image name from NMS weather page,
        so that we eventually learn all image names.
    """
    # slow_interval = 300  # seconds, for normal monitoring
    # fast_interval = 60  # seconds, for monitoring when roof expected to open or close
    # timer = RoofTimer(slow_interval, fast_interval)
    current_interval = 240  # seconds.
    fullpath = 'C:/Astro/NMS/record_nms_status_image_name.txt'
    with open(fullpath, 'a') as f:  # will append to file.
        start_string = '\nSTARTING:\n'
        print(start_string)
        f.write(start_string)
        while True:
            # timer.wait_interval()
            image_name = get_nms_status_image_name()
            utc_latest_capture = datetime.now(timezone.utc)
            utc_string = utc_latest_capture.strftime('%Y-%m-%d %H:%M:%S')
            write_string = utc_string + ': NMS image file >' + image_name + '<'
            print(write_string)
            f.write(write_string + '\n')
            f.flush()
            sleep(current_interval)


if __name__ == '__main__':
    monitor_nms_status_via_e_mail()
    # monitor_nms_status_via_html()
    # record_nms_status_image_name()
