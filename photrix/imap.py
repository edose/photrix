__author__ = "Eric Dose, Albuquerque"

""" This module:  web.py
    Various e-mail access routines.
"""

# Python core:
import os
import imaplib
import email
from email.header import decode_header


THIS_PACKAGE_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

HTTP_OK_CODE = 200  # "OK. The request has succeeded."
MAX_WEBOBS_LINES = 200  # for safety (& as dictated in any case by webobs API (as of Jan 2017)).


__________E_MAIL_ACCESS______________________________________________________ = 0


def get_field_string(msg, field_name):
    content, encoding = decode_header(msg[field_name])[0]
    if isinstance(content, bytes):
        return content.decode(encoding)
    return content


def get_imap_credentials():
    fullpath = os.path.join('C:/Astro/NMS/creds', 'creds.txt')
    with open(fullpath) as f:
        lines = f.readlines()
    creds_dict = {}
    for line in lines:
        items = [s.strip() for s in line.split(':')]
        if len(items) == 2:
            creds_dict[items[0]] = items[1]
    return creds_dict


def get_most_recent_relevant_email(subject_start, from_strings):
    if isinstance(from_strings, str):
        from_strings = [from_strings]
    creds = get_imap_credentials()
    imap = imaplib.IMAP4_SSL(creds['imap'])
    # authenticate
    imap.login(creds['login'], creds['fargelbnurrr'])
    status, messages = imap.select('INBOX')
    n_in_inbox = int(messages[0])
    n_to_retrieve = n_in_inbox
    for i in range(n_in_inbox, n_in_inbox - n_to_retrieve, -1):  # fetches latest to earliest.
        # TODO: try/except around imap.fetch() call, try again a few times before raising own exception.
        # TODO: **or** use imap.search() to avoid the problem.
        # res, msg = imap.fetch(str(i), "(RFC822)")
        res, msg = imap.fetch(str(i), '(BODY.PEEK[HEADER])')  # fetches headers only.
        for response in msg:
            if res == 'OK' and isinstance(response, tuple):
                this_msg = email.message_from_bytes(response[1])
                date = get_field_string(this_msg, 'Date')
                subject = get_field_string(this_msg, 'Subject')
                from_str = get_field_string(this_msg, 'From')
                subject_start_relevant = subject_start is None or \
                                         subject.lower().startswith(subject_start.lower())
                from_field_relevant = from_strings is None or \
                                      all([f.lower() in from_str.lower() for f in from_strings])
                if subject_start_relevant and from_field_relevant:
                    return {'i': i, 'date': date, 'subject': subject, 'from_str': from_str}
                # if subject.lower().startswith(subject_start.lower()) and\
                #     all([f.lower() in from_str.lower() for f in from_strings]):
                #     return {'i': i, 'date': date, 'subject': subject, 'from_str': from_str}
    return None
