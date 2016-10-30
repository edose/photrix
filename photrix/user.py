from math import pi
from datetime import datetime, timedelta
import os
import json
import ephem
from copy import copy
from photrix.util import Timespan, RaDec, hex_degrees_as_degrees

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"

PHOTRIX_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SITE_DIRECTORY = os.path.join(PHOTRIX_ROOT_DIRECTORY, "site")
MOON_PHASE_NO_FACTOR = 0.05


class Site:
    """
    Object: holds site info (no instrument or AN info), read from a json site file.
    Usage: site = Site("BDO_Kansas")
    Attributes (string unless otherwise noted):
        .name : site's name
        .filename : json file name
        .description : site description
        .longitude, .latitude : of site, in degrees (float)
        .elevation : of site, in meters (float)
        .min_altitude : in degrees (float)
        .twilight_sun_alt : in degrees (float)
        .is_valid : True if attribute values appear valid (boolean)
    """
    def __init__(self, site_name, site_directory=SITE_DIRECTORY):
        site_fullpath = os.path.join(site_directory, site_name + ".json")
        with open(site_fullpath) as data_file:
            data = json.load(data_file)
        self.name = data.get("name", site_name)
        self.filename = site_name + ".json"
        self.description = data.get("description", "")
        self.longitude = data.get("longitude", None)  # West longitude is negative.
        if self.longitude is not None:
            self.longitude = hex_degrees_as_degrees(str(self.longitude))
        self.latitude = data.get("latitude", None)  # South latitude is negative.
        if self.latitude is not None:
            self.latitude = hex_degrees_as_degrees(str(self.latitude))
        self.elevation = data.get("elevation", 500)  # in meters
        self.min_altitude = data.get("min_altitude", 0)
        self.twilight_sun_alt = data.get("twilight_sun_alt", -10)

        is_valid = (self.longitude is not None) and (self.latitude is not None)  # default
        if is_valid:
            if not ((self.longitude >= -180) and (self.longitude <= +180)):
                is_valid = False
            if not ((self.latitude >= -90) and (self.longitude <= +90)):
                is_valid = False
        self.is_valid = is_valid

    def __repr__(self):
        return "site('" + self.filename + "')"

    def __str__(self):
        return self.__repr__() + " valid=" + str(self.is_valid)


class Instrument:
    """
    Object: holds instrument info (no site or AN info).
    Usage: inst = Instrument("Borea")
    Attributes (string unless otherwise noted):
        .name :
        .filename :
        .description :
        .min_distance_full_moon : in degrees (float)
        .mount["model"] :
        .mount["slew_rate_ra"] : in deg/sec (float)
        .mount["slew_rate_dec"] : in degrees (float)
        .mount["sec_to_speed_ra"] : in seconds (float)
        .mount["sec_to_speed_dec"] : in seconds (float)
        .ota["model"] :
        .ota["focal_length_mm"] : in millimeters (float)
        .camera["model"] :
        .camera["pixels_x"] : (int)
        .camera["pixels_y"] : (int)
        .camera["microns_per_pixel"] : (float)
        .camera["shortest_exposure"] : in seconds (float)
        .camera["saturation_adu"] : maximum ADUs while linear (float)
        .filters[filter(string)]["reference_exposure_mag10"] : possibly several (float)
        .filters[filter(string)]["transform"][color(string)] : possibly several (float)
        .is_valid : True if attribute values appear valid (boolean)
    """
    def __init__(self, instrument_name):
        photrix_root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        instrument_fullpath = os.path.join(photrix_root_dir, "instrument",
                                           (instrument_name + ".json"))
        with open(instrument_fullpath) as data_file:
            data = json.load(data_file)

        self.name = data.get("name", instrument_name)
        self.filename = instrument_name + ".json"
        self.description = data.get("description", "")
        self.min_distance_full_moon = data.get("min_distance_full_moon", 60)  # degrees
        mount = data.get("mount")
        mount["model"] = mount.get("model", "")
        mount["slew_rate_ra"] = mount.get("slew_rate_ra", 4)
        mount["slew_rate_dec"] = mount.get("slew_rate_dec", 4)
        mount["sec_to_speed_ra"] = mount.get("sec_to_speed_ra", 1)
        mount["sec_to_speed_dec"] = mount.get("secto_speed_dec", 1)
        self.mount = mount

        ota = data.get("ota")
        ota["model"] = ota.get("model", "")
        ota["focal_length_mm"] = ota.get("focal_length_mm", 0)
        self.ota = ota

        camera = data.get("camera")
        camera["model"] = camera.get("model", "")
        camera["pixels_x"] = camera.get("pixels_x", 0)
        camera["pixels_y"] = camera.get("pixels_y", 0)
        camera["microns_per_pixel"] = camera.get("microns_per_pixel", 0)
        camera["shortest_exposure"] = camera.get("shortest_exposure", 0)
        camera["saturation_adu"] = camera.get("saturation_adu", 64000)
        self.camera = camera

        self.filters = data.get("filters")

        is_valid = True  # default to be falsified if any error.

        self.is_valid = is_valid

    def filter_list(self):
        return list(self.filters.keys())

    def __repr__(self):
        return "Instrument('" + self.filename + "')"

    def __str__(self):
        return self.__repr__() + " valid=" + str(self.is_valid)


class Astronight:
    def __init__(self, an_date_string, site_name):
        """
        Object: relevant info specific to one observing night at one site.
        Usage: an = Astronight("20160102", "BDO_Kansas")
        Attributes (string unless otherwise noted):
            .an_date_string : as "20161011"
            .site_name : as "BDO_Kansas"
            .site : this site (Site object)
            ._observer : this site and this astronight (ephem.Observer object)
            .local_middark_ut = local mid-dark time, this astronight, in UTC time (datetime object)
            .local_middark_jd = local mid-dark time, this astronight, in Julian date (float)
            .local_middark_lst = local mid-dark time, this astronight, local sidercial time in degrees (float)
            .moon_radec = moon location at middark (RaDec object)
            .moon_phase = moon phase at middark (float)
            .dark = evening twilight and morning twilight (Timespan object)
        """
        self.an_date_string = an_date_string
        self.site_name = site_name
        self.site = Site(site_name)

        obs = ephem.Observer()  # for local use (within __init__()).
        obs.lat, obs.lon = str(self.site.latitude), str(self.site.longitude)
        obs.elevation = self.site.elevation
        self._observer = obs  # for access by other methods in this class.

        # get local middark times for requested Astronight.
        an_year = int(an_date_string[0:4])
        an_month = int(an_date_string[4:6])
        an_day = int(an_date_string[6:8])
        approx_midnight_utc = datetime(an_year, an_month, an_day, 24, 0, 0) + \
            timedelta(hours=-self.site.longitude / 15.0)
        sun = ephem.Sun()
        obs.horizon = str(self.site.twilight_sun_alt)
        sun.compute(obs)
        twilight_dusk = obs.previous_setting(sun, start=approx_midnight_utc)
        twilight_dawn = obs.next_rising(sun, start=approx_midnight_utc)
        self.dark = Timespan(twilight_dusk, twilight_dawn)
        self.local_middark_utc = self.dark.midpoint
        self.local_middark_jd = ephem.julian_date(self.local_middark_utc)
        obs.date = self.local_middark_utc
        self.local_middark_lst = obs.sidereal_time() * 180. / pi  # in degrees

        # Prep all other solar-system bodies, as RaDecs at local mid-dark.
        moon = ephem.Moon()
        obs.date = self.local_middark_utc
        moon.compute(obs)
        self.moon_radec = RaDec(str(moon.ra), str(moon.dec))
        self.moon_phase = moon.moon_phase
        # jupiter = ephem.Jupiter(obs)
        # self.jupiter_radec = RaDec(str(jupiter.ra), str(jupiter.dec))
        # saturn = ephem.Saturn(obs)
        # self.saturn_radec = RaDec(str(saturn.ra), str(saturn.dec))
        # venus = ephem.Venus(obs)
        # self.venus_radec = RaDec(str(venus.ra), str(venus.dec))

    def available(self, min_alt, ra, dec):
        """
        Returns Timespan object defining when this RA,Dec may be observed during this astronight.
        Site data are taken from stored ephem Observer object self._observer.
        Usage: ts = an.available(+28, fov.ra, fov.dec)
        """
        obs = copy(self._observer)
        obs.horizon = str(min_alt)
        b = ephem.FixedBody()
        b._ra = str(ra)
        b._dec = str(dec)
        b._epoch = '2000'
        b.compute(obs)
        previous_transit_utc = self._observer.previous_rising(b, start=self.local_middark_utc)
        next_transit_utc = self._observer.next_setting(b, start=self.local_middark_utc)
        if next_transit_utc - self.local_middark_utc < self.local_middark_utc - previous_transit_utc:
            nearest_transit_utc = next_transit_utc
        else:
            nearest_transit_utc = previous_transit_utc
        obs.date = nearest_transit_utc
        b.compute(obs)
        if b.alt > min_alt * pi / 180.0:
            start = obs.previous_rising(b, start=nearest_transit_utc)
            end = obs.next_setting(b, start=nearest_transit_utc)
            above_min_alt = Timespan(start, end)
        else:
            above_min_alt = Timespan(nearest_transit_utc, nearest_transit_utc)  # null Timespan
        return self.dark.intersect(above_min_alt)

    def __repr__(self):
        return "Astronight '" + self.an_date_string + "' at site '" + self.site_name + "'."

    def acp_header_string(self):
        """ Returns an info string to include atop an ACP plan file, as:
            ; sunset-rise 0033-1145 UT , -11 deg alt @ 0119-1055 UT  //  LST = cdt + 1128.
            ; moon 15% @ (~22h04,-10) rise 0936 UT
        Usage: s = an.acp_header_string()
        """
        obs = copy(self._observer)
        sun = ephem.Sun(obs)
        moon = ephem.Moon(obs)

        obs.horizon = '-0.25'  # top edge below horizon
        sunset_local = obs.previous_setting(sun, self.local_middark_utc)
        sunrise_local = obs.next_rising(sun, self.local_middark_utc)

        obs.horizon = str(self.site.min_altitude)
        evening_twilight_local = obs.previous_setting(sun, self.local_middark_utc)
        morning_twilight_local = obs.next_rising(sun, self.local_middark_utc)

        obs.horizon = '-0.25'  # top edge below horizon
        moon_pct = self.moon_phase * 100.0
        moon_ra = self.moon_radec.ra
        moon_dec = self.moon_radec.dec
        # need moon logic here  [ rise vs set; also add "no factor" if phase < threshold ]
        previous_transit_utc = obs.previous_rising(moon, start=self.local_middark_utc)
        next_transit_utc = obs.next_setting(moon, start=self.local_middark_utc)
        if next_transit_utc-self.local_middark_utc < self.local_middark_utc-previous_transit_utc:
            nearest_transit_utc = next_transit_utc
        else:
            nearest_transit_utc = previous_transit_utc
        obs.date = nearest_transit_utc
        moon.compute(obs)
        if moon.alt > -0.25 * pi / 180.0:
            start = obs.previous_rising(moon, start=nearest_transit_utc).datetime()
            end = obs.next_setting(moon, start=nearest_transit_utc).datetime()
            above_min_alt = Timespan(start, end)
        else:
            above_min_alt = Timespan(nearest_transit_utc.datetime(),
                                     nearest_transit_utc.datetime())  # null Timespan
        return self.dark.intersect(above_min_alt)



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


