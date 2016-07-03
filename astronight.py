import ephem
from util import Timespan
from user import Instrument

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"

# THIS MUST WAIT UNTIL module site IS SORTED OUT.

class Astronight:
    def __init__(self, an_date_string, instrument_name):
        """
        Object: relevant info specific to one observing night.
        Usage: an = Astronight("20160102", "Borea")
        """
        self.an_date_string = an_date_string
        self.instrument_name = instrument_name

        instrument = Instrument(instrument_name)  # instrument is of pyphot class instrument.
        self.instrument = instrument

        site = ephem.Observer()                   # site is an ephem-package object.
        site.lat = str(instrument.latitude)
        site.lon = str(instrument.longitude)
        self.site = site

        # get local midnight for requested Astronight.
        an_year = int(an_date_string[0:4])
        an_month = int(an_date_string[4:6])
        an_day = int(an_date_string[6:8])
        self.localMidnightUT = ephem.Date((an_year, an_month, an_day, 24, 0, 0)) \
            - instrument.timezone * ephem.hour
        self.localMidnightJD = ephem.julian_date(self.localMidnightUT)
        site.date = self.localMidnightUT
        self.localMidnightLST = site.sidereal_time()

        # Prepare all solar system body objects.
        self.sun = ephem.Sun()
        self.moon = ephem.Moon()
        self.jupiter = ephem.Jupiter()
        self.saturn = ephem.Saturn()
        self.venus = ephem.Venus()

        # Prepare timespans for fully dark, this AN.
        site.horizon = str(instrument.twilightSunAlt)
        self.time_dark = Timespan(site.previous_setting(self.sun, self.localMidnightUT),
                                  site.next_rising(self.sun, self.localMidnightUT))

        # Prepare timespan for moon up, this AN.
        site.horizon = '0'
        site.date = self.localMidnightUT
        self.moon.compute(site)
        if self.moon.alt > 0:
            self.moon_up_and_dark_timespan = Timespan(
                site.previous_rising(self.moon, start=self.localMidnightUT),
                site.next_setting(self.moon, start=self.localMidnightUT)).intersect(self.time_dark)

        else:
            moon_up_west_and_dark_timespan = Timespan(
                self.time_dark.start, site.previous_setting(self.moon, start=self.localMidnightUT))
            moon_up_east_and_dark_timespan = Timespan(
                site.next_rising(self.moon, start=self.localMidnightUT), self.time_dark.end)
            if moon_up_west_and_dark_timespan.seconds > moon_up_east_and_dark_timespan.seconds:
                self.moon_up_and_dark_timespan = moon_up_west_and_dark_timespan
            else:
                self.moon_up_and_dark_timespan = moon_up_east_and_dark_timespan

    def __repr__(self):
        return "Astronight '" + self.an_date_string + "' at site '" + self.instrument.sitename + "'."

    def ACP_header_string(self):
        """ Returns an info string to include atop an ACP plan file, as:
            ; sunset 1946 cdt, -10 deg alt @ 2038 cdt  //  LST = cdt + 1128.
            ; moon 15% @ (~22h04,-10) rise 0436 cdt
        Usage: s = an.ACP_header_string()
        """
        site = self.site
        site.date = self.localMidnightUT

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
