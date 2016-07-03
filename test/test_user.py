from user import *
from util import hex_degrees_as_degrees

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


def test_Site():
    site_name = "Site_test1"
    s = Site(site_name)
    assert s.is_valid
    assert s.name == site_name
    assert s.filename == s.name + ".json"
    assert s.description.startswith("Bois d'Arc Obs")

    assert s.longitude == hex_degrees_as_degrees("-95:53:18")
    assert s.latitude == hex_degrees_as_degrees("+38:55:29")
    assert s.altitude == 350


def test_Instrument():
    instrument_name = "Instrument_test1"
    i = Instrument(instrument_name)
    assert i.is_valid
    assert i.name == instrument_name
    assert i.filename == i.name + ".json"
    assert i.min_altitude == 25
    assert i.twilight_sun_alt == -9
    assert i.min_distance_full_moon == 50
    assert i.mount["model"].startswith("Paramount MX")
    assert i.mount["slew_rate_ra"] == 4
    assert i.mount["slew_rate_dec"] == 4
    assert i.mount["sec_to_speed_ra"] == 1
    assert i.mount["sec_to_speed_dec"] == 1
    assert i.ota["model"].startswith("Celestron C14 ")
    assert i.ota["focal_length_mm"] == 2710
    assert i.camera["model"].startswith("SBIG STXL 6303")
    assert i.camera["pixels_x"] == 3072
    assert i.camera["pixels_y"] == 2048
    assert i.camera["microns_per_pixel"] == 9
    assert i.camera["shortest_exposure"] == 0.6
    assert i.camera["saturation_adu"] == 54000
    assert i.filters["V"]["reference_exposure_mag10"] == 22
    assert i.filters["R"]["reference_exposure_mag10"] == 30
    assert i.filters["V"]["transform"]["V-I"] == 0.02
    assert i.filters["V"]["transform"]["B-V"] == 0
    assert i.filters["V"]["transform"]["V-R"] == 0
    assert i.filters["I"]["transform"]["V-I"] == 0.025

    instrument_name = "Instrument_test2"
    i = Instrument(instrument_name)
    assert i.is_valid
    assert i.name == "XXX"
    assert i.filename == instrument_name + ".json"
    assert i.min_altitude == 0  # absent -> default
    assert i.twilight_sun_alt == -10  # absent -> default
    assert i.min_distance_full_moon == 60  # absent -> default
    assert i.mount["model"] == ""
    assert i.mount["slew_rate_ra"] == 7
    assert i.mount["slew_rate_dec"] == 4
    assert i.mount["sec_to_speed_ra"] == 1
    assert i.mount["sec_to_speed_dec"] == 1
    assert i.ota["model"] == ""
    assert i.ota["focal_length_mm"] == 0
    assert i.camera["model"] == ""
    assert i.camera["pixels_x"] == 0
    assert i.camera["pixels_y"] == 0
    assert i.camera["microns_per_pixel"] == 0
    assert i.camera["shortest_exposure"] == 0
    assert i.camera["saturation_adu"] == 64000
    assert i.filters["V"]["reference_exposure_mag10"] == 22
    assert i.filters["V"]["transform"]["V-I"] == 0.02
