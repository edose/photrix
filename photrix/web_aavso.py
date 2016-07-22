import requests
from bs4 import BeautifulSoup

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


class WebObsAAVSO:
    observations = []  # will be list of lists

    def __init__(self, star_id):
        star_safe_name = star_id.replace("+", "%2B").replace(" ", "+")
        url = "https://www.aavso.org/apps/webobs/results/?star=" + \
              star_safe_name + "&num_results=50&obs_types=vis+ccd"
        r = requests.get(url)
        if r.status_code == 200:
            soup = BeautifulSoup(r.text, 'html.parser')
            obs_lines = soup.find_all('tr', class_='obs')  # NB: "class_" not "class" (reserved).
            for line in obs_lines:
                cells = line.find_all('td')
                cell_text = [cell.text for cell in cells]
                self.observations.append(cell_text)

            # at this point, we have a list of lists.
            # --> Parse and save internally, probably as a pandas data frame.
            # --> Do we want to get the 2nd line/obs as well?

    def mag_range(self, filter="V"):
        pass
