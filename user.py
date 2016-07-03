import json
import os
from util import hex_degrees_as_degrees

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


class Site:
    """
    Object: holds site info (no instrument or AN info).
    Usage: site = Site("BDO_Kansas")
    """
    def __init__(self, site_name):
        site_fullpath = os.path.dirname(__file__) + "/site/" + site_name + ".json"
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
        self.altitude = data.get("altitude", 500)  # in meters
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
    """
    def __init__(self, instrument_name):
        instrument_fullpath = os.path.dirname(__file__) + "/instrument/" + instrument_name + ".json"
        # print (">>>" + instrument_fullpath + "<<<")
        with open(instrument_fullpath) as data_file:
            data = json.load(data_file)

        self.name = data.get("name", instrument_name)
        self.filename = instrument_name + ".json"
        self.description = data.get("description", "")
        self.min_altitude = data.get("min_altitude", 0)
        self.twilight_sun_alt = data.get("twilight_sun_alt", -10)
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
