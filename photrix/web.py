import pandas as pd
import requests
from bs4 import BeautifulSoup

from photrix.user import Astronight

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


class WebObsAAVSO:
    def __init__(self, star_id):
        star_safe_name = star_id.replace("+", "%2B").replace(" ", "+")
        url = "https://www.aavso.org/apps/webobs/results/?star=" + \
              star_safe_name + "&num_results=50&obs_types=vis+ccd"
        r = requests.get(url)
        obs_list = []
        if r.status_code == 200:
            soup = BeautifulSoup(r.text, 'html.parser')
            obs_lines = soup.find_all('tr', class_='obs')  # NB: "class_" not "class" (reserved).
            for line in obs_lines:
                cells = line.find_all('td')
                cell_strings = [cell.text for cell in cells]
                obs_list.append(cell_strings)
        df = pd.DataFrame(obs_list, columns=['X0', 'target_name', 'jd_str', 'date_string',
                                             'mag_str', 'error_str',
                                             'filter', 'observer', 'X8'])
        df = df.assign(jd=[float(xx) for xx in df['jd_str']],
                       mag=[float(xx) for xx in df['mag_str']],
                       error=[float(xx) for xx in df['error_str']])
        # print(df.dtypes)
        df = df.drop(['X0', 'X8', 'jd_str', 'mag_str', 'error_str'], axis=1)
        # print(df.head())
        self.table = df

    def most_recent_jd_mag(self, filter='V'):
        df = self.table
        df = df[df['filter'].str.lower() == filter.lower()]  # only rows with chosen filter.
        df = df.sort_values('jd', ascending=False)  # latest row first
        if len(df) <= 0:
            return None
        return df['jd'].iloc[0], df['mag'].iloc[0]  # return tuple of JD, mag

    def days_gap_jd(self, jd_target, filter='V'):
        return jd_target - self.most_recent_jd_mag(filter)[0]

    def days_gap_an(self, an, filter='V'):
        if not isinstance(an, Astronight):
            return None
        return self.days_gap_jd(jd_target=an.local_middark_jd, filter=filter)

