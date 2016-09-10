from datetime import datetime, timedelta
import ephem
from photrix.user import Site
from photrix.util import Timespan, RaDec

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


class Astronight:
    def __init__(self, an_date_string, site_name, offset_from_UTC):
        """
        Object: relevant info specific to one observing night at one site.
        Usage: an = Astronight("20160102", "BDO_Kansas")
        """
        self.an_date_string = an_date_string
        self.site_name = site_name
        self.site = Site(site_name)

        observer = ephem.Observer()  # observer is an ephem-package object.
        observer.lat = str(self.site.latitude)
        observer.lon = str(self.site.longitude)
        observer.elevation = self.site.elevation
        self.observer = observer

        # get local midnight for requested Astronight.
        an_year = int(an_date_string[0:4])
        an_month = int(an_date_string[4:6])
        an_day = int(an_date_string[6:8])
        self.localMidnightUT = datetime(an_year, an_month, an_day, 12, 0, 0) + \
            timedelta(hours=12-offset_from_UTC)
        self.localMidnightJD = ephem.julian_date(self.localMidnightUT)
        observer.date = self.localMidnightUT
        self.localMidnightLST = observer.sidereal_time()  # in radians

        # Prepare all solar system body objects, exposing RaDecs as of local midnight.
        observer.date = self.localMidnightUT
        sun = ephem.Sun(observer)
        moon = ephem.Moon(observer)
        self.moon_radec = RaDec(str(moon.ra), str(moon.dec))
        self.moon_phase = moon.moon_phase
        jupiter = ephem.Jupiter(observer)
        self.jupiter_radec = RaDec(str(jupiter.ra), str(jupiter.dec))
        saturn = ephem.Saturn(observer)
        self.saturn_radec = RaDec(str(saturn.ra), str(saturn.dec))
        venus = ephem.Venus(observer)
        self.venus_radec = RaDec(str(venus.ra), str(venus.dec))

        # Prepare timespans for fully dark, this AN.
        observer.horizon = str(self.site.twilight_sun_alt)
        twilight_dusk = observer.previous_setting(sun, self.localMidnightUT).datetime()
        twilight_dawn = observer.next_rising(sun, self.localMidnightUT).datetime()
        self.dark = Timespan(twilight_dusk, twilight_dawn)

        # Prepare timespan for moon above horizon, this AN & site.
        observer.horizon = '0'
        observer.date = self.localMidnightUT
        moon.compute(observer)
        if moon.alt > 0:
            pr = observer.previous_rising(moon, start=self.localMidnightUT).datetime()
            ns = observer.next_setting(moon, start=self.localMidnightUT).datetime()
            moon_up = Timespan(pr, ns)
            self.dark_and_moon_up = moon_up.intersect(self.dark)
            #self.dark_and_moon_up = Timespan(
            #    observer.previous_rising(moon, start=self.localMidnightUT).datetime(),
            #    observer.next_setting(moon, start=self.localMidnightUT).datetime()) \
            #    .intersect(self.dark)
        else:
            dark_and_moon_up_west = Timespan(
                self.dark.start,
                observer.previous_setting(moon, start=self.localMidnightUT))
            dark_and_moon_up_east = Timespan(
                observer.next_rising(moon, start=self.localMidnightUT),
                self.dark.end)
            # choose the longer in unlikely case: moon up at both twilights but not at midnight.
            if dark_and_moon_up_west.seconds > dark_and_moon_up_east.seconds:
                self.moon_up_and_dark = dark_and_moon_up_west
            else:
                self.moon_up_and_dark = dark_and_moon_up_east

    def available(self, site_long, site_lat, min_alt, ra, dec):
        """
        Returns Timespan object defining which  may be observed during this astronight.
        Usage: ts = astronight.Astronight.available(site.long, site.lat, +28, fov.ra, fov.dec)
        """
        pass

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
