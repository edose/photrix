import pandas as pd
import requests
from bs4 import BeautifulSoup

from photrix.user import Astronight
from photrix.util import float_or_none

__author__ = "Eric Dose :: New Mexico Mira Project, Albuquerque"
HTTP_OK_CODE = 200  # "OK. The request has succeeded."
MAX_WEBOBS_LINES = 200  # for safety (& as dictated in any case by webobs API (as of Jan 2017)).


def get_aavso_webobs_raw_table(star_id, num_obs=100, jd_start=None, jd_end=None):
    #  simply returns an pandas dataframe containing data,
    #   no parsing or cacheing. If star ID not in webobs, this returns a dataframe with no rows.
    """
    Downloads observations from AAVSO's webobs for ONE star (not fov), returns pandas dataframe.
       If star not in AAVSO's webobs site, return a dataframe with no rows.
       Columns: target_name, date_string, filter, observer, jd, mag, error.
    :param star_id: the STAR id (not the fov's name).
    :param num_obs: number of observations to get.
    :param jd_start: optional Julian date.
    :param jd_end: optional JD.
    :return: simple pandas dataframe containing data for 1 star, 1 row per observation downloaded.
    """

    star_safe_name = star_id.replace("+", "%2B").replace(" ", "+")
    num_obs = min(num_obs, MAX_WEBOBS_LINES)
    url = "https://www.aavso.org/apps/webobs/results/?star=" + star_safe_name + \
          "&num_results=" + str(num_obs) + \
          "&obs_types=vis+ccd"
    if jd_start is not None:
        url += "&start=" + str(jd_start)
    if jd_end is not None:
        url += "&end=" + str(jd_end)
    r = requests.get(url)
    obs_list = []
    if r.status_code == HTTP_OK_CODE:
        soup = BeautifulSoup(r.text, 'html.parser')
        obs_lines = soup.find_all('tr', class_='obs')  # NB: "class_" not "class" (reserved).
        for line in obs_lines:
            cells = line.find_all('td')
            cell_strings = [cell.text for cell in cells]
            obs_list.append(cell_strings)
    df = pd.DataFrame(obs_list, columns=['X0', 'target_name', 'jd_str', 'date_string',
                                         'mag_str', 'error_str',
                                         'filter', 'observer', 'X8'])
    df = df.assign(jd=[float_or_none(xx) for xx in df['jd_str']],
                   mag=[float_or_none(xx) for xx in df['mag_str']],
                   error=[float_or_none(xx) for xx in df['error_str']])
    # print(df.dtypes)
    df = df.drop(['X0', 'X8', 'jd_str', 'mag_str', 'error_str'], axis=1)
    # print(df.head())
    return df


# def most_recent_jd_mag(self, filter='V'):
#     df = self.table
#     df = df[df['filter'].str.lower() == filter.lower()]  # only rows with chosen filter.
#     df = df.sort_values('jd', ascending=False)  # latest row first
#     if len(df) <= 0:
#         return None
#     return df['jd'].iloc[0], df['mag'].iloc[0]  # return tuple of JD, mag
#
#
# def days_gap_jd(self, jd_target, filter='V'):
#     return jd_target - self.most_recent_jd_mag(filter)[0]
#
#
# def days_gap_an(self, an, filter='V'):
#     if not isinstance(an, Astronight):
#         return None
#     return self.days_gap_jd(jd_target=an.local_middark_jd, filter=filter)

