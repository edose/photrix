import json
import os
from photrix.util import hex_degrees_as_degrees

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"

PHOTRIX_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SITE_DIRECTORY = os.path.join(PHOTRIX_ROOT_DIRECTORY, "site")

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
        instrument_fullpath = os.path.join(photrix_root_dir, "instrument", \
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
