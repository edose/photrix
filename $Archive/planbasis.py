import ephem
from user import Site, Instrument

from photrix.util import Timespan

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


class PlanBasis:
    def __init__(self, site_name, instrument_name, astronight_date_string):
        """
        Object: relevant info specific to one observing night, at one site, with one instrument.
        Usage: basis = PlanBasis("Bois dArc", "Borea", "20160102")
        """
        self.site_name = site_name
        self.site = Site(site_name)
        self.instrument_name = instrument_name
        self.instrument = Instrument(instrument_name)
        self.astronight_date_string = astronight_date_string

        # Construct ephem observer (at specified site site).
        obs = ephem.Observer()
        obs.lat = str(self.site.latitude)
        obs.lon = str(self.site.longitude)

        # Get local midnight at this site, for this Astronight.
        an_year = int(astronight_date_string[0:4])
        an_month = int(astronight_date_string[4:6])
        an_day = int(astronight_date_string[6:8])
        approx_hours_vs_ut = self.site.longitude / 15.0
        self.localMidnight_UT = ephem.Date((an_year, an_month, an_day, 24, 0, 0)) \
            - (approx_hours_vs_ut * ephem.hour)  # the 24 means midnight *after* evening begins.
        self.localMidnight_JD = ephem.julian_date(self.localMidnight_UT)

        # Prepare Timespan for fully dark.
        self.sun = ephem.Sun()
        obs.horizon = str(self.instrument.twilightSunAlt)
        self.fully_dark = Timespan(obs.previous_setting(self.sun, self.localMidnight_UT),
                                   obs.next_rising(self.sun, self.localMidnight_UT))

        # Prepare Timespan for moon up while fully dark.
        self.moon = ephem.Moon()
        obs.horizon = '0'
        if self.moon.alt > 0:
            self.moon_up_and_fully_dark = Timespan(
                obs.previous_rising(self.moon, start=self.localMidnight_UT),
                obs.next_setting(self.moon, start=self.localMidnight_UT))\
                .intersect(self.fully_dark)

        else:
            moon_up_west_and_fully_dark = Timespan(
                self.fully_dark.start,
                obs.previous_setting(self.moon, start=self.localMidnight_UT))
            moon_up_east_and_fully_dark = Timespan(
                obs.next_rising(self.moon, start=self.localMidnight_UT), self.time_dark.end)
            if moon_up_west_and_fully_dark.seconds > moon_up_east_and_fully_dark.seconds:
                self.moon_up_and_fully_dark = moon_up_west_and_fully_dark
            else:
                self.moon_up_and_fully_dark = moon_up_east_and_fully_dark

        # Prepare locations of moon and planets (simplify: one RA Dec for entire night)
        self.jupiter = ephem.Jupiter(self.localMidnight_UT)
        self.saturn = ephem.Saturn(self.localMidnight_UT)
        self.venus = ephem.Venus(self.localMidnight_UT)

    def __repr__(self):
        return "PlanBasis('" + self.astronight_date_string + "', '" + \
               self.site_name + "', '" + self.instrument_name + "')"

    def ACP_header_string(self):
        """ Returns an info string to include atop an ACP plan file, as:
            ; LST = cdt + 1128
            ; sunset 1946 cdt, -10 deg alt 2038-0507 cdt
            ; moon 15% @ (~22h04,-10) rise 0436 cdt
        Usage: s = an.ACP_header_string()
        """
        site = self.site
        site.date = self.localMidnight_UT

        # Time info line.
        time_from_UT = self.instrument.timezone
        LST_diff = self.localMidnightLST
        header_string = "; " + self.an_date_string + "  local time = " + str(time_from_UT) + \
                        "  //  LST = local + " + str(LST_diff) + "\n"

        # # Sun info line.
        # site.horizon = '0'
        # sunset_local = site.previous_setting(self.sun, self.localMidnightUT) + time_from_UT
        # sunrise_local = site.next_rising(self.sun, self.localMidnightUT) + time_from_UT
        # site.horizon = self.instrument.twilightSunAlt  # this needs to be divided by degrees/radian
        # twilight_angle = site.horizon
        # evening_twilight_local = site.previous_setting(self.sun,
        #                                                self.localMidnightUT) + time_from_UT
        # morning_twilight_local = site.next_rising(self.sun, self.localMidnightUT) + time_from_UT
        # site.horizon = '0'
        # header_string += "; sunset " + str(sunset_local) + " local,   " + str(twilight_angle) + \
        #                  " deg alt @ " + str(evening_twilight_local) + " - " + \
        #                  str(morning_twilight_local) + "   sunrise " + str(sunrise_local) + "\n"
        #
        # # Moon info line.
        # moon_pct = self.moon.phase
        # moon_up = self.moon_up_and_dark_timespan
        # # --> Here, do logic or whether moonrise or moonset is actually within nighttime
        # header_string += "; moon " + str(moon_pct) + "% @ (" + str(self.moon.ra) + ", " + \
        #                  str(self.moon.dec) + ")   rise/set=??? local" + "\n"

        return header_string

    def when_FOV_observable(self, FOV):
        """ Returns Timespan object containing when FOV is observable during this Astronight.
        Usage: ts = an.when_FOV_observable(FOV)
        """
        pass  # probably uses same logic as for moon, just above at end of AN constructor.









        # ---------------------------------------------------------------------


if __name__ == '__main__':
    an = Astronight("20160409", "Borea")
    print(repr(an))
    acp = an.ACP_header_string()
    print(acp)
