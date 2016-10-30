from datetime import datetime, timezone
from photrix import user                 # call: user.fn() & user.Class()
from photrix.util import hex_degrees_as_degrees

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


def test_Site():
    site_name = "Site_test1"
    s = user.Site(site_name)
    assert s.is_valid
    assert s.name == site_name
    assert s.filename == s.name + ".json"
    assert s.description.startswith("Bois d'Arc Obs")
    assert s.longitude == hex_degrees_as_degrees("-95:53:18")
    assert s.latitude == hex_degrees_as_degrees("+38:55:29")
    assert s.elevation == 350
    assert s.min_altitude == 25
    assert s.twilight_sun_alt == -9


def test_Instrument():
    instrument_name = "Instrument_test1"
    i = user.Instrument(instrument_name)
    assert i.is_valid
    assert i.name == instrument_name
    assert i.filename == i.name + ".json"
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
    assert set(i.filter_list()) == {"V", "R", "I"}  # set

    instrument_name = "Instrument_test2"
    i = user.Instrument(instrument_name)
    assert i.is_valid
    assert i.name == "XXX"
    assert i.filename == instrument_name + ".json"
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


def test_Astronight():
    # Case: moon up at midnight.
    an_date_string = "20160910"
    site_string = "BDO_Kansas"
    an = user.Astronight(an_date_string, site_string, -5)
    assert an.an_date_string == an_date_string
    assert an.site_name == site_string  # the name
    assert an.site.name == site_string  # the object (need to expose though Astronight()?)
    assert an.localMidnightUT == datetime(2016, 9, 11, 5, 0, 0, tzinfo=timezone.utc)
    assert abs(an.localMidnightJD - 2457642.708338) < 1 / 24 / 3600  # one sec tolerance
    target_lst_seconds = 21*3600 + 59*60 + 3
    an_lst_seconds = an.localMidnightLST * 240.0
    assert abs(target_lst_seconds - an_lst_seconds) < 1  # one sec tolerance
    assert abs(an.moon_radec.ra - 278.16825) < 1/3600
    assert abs(an.moon_radec.dec - -19.1090833) < 1/3600
    assert an.dark.start == datetime(2016, 9, 11, 1, 22, 54, 830573, tzinfo=timezone.utc)
    assert an.dark.end == datetime(2016, 9, 11, 11, 17, 48, 383649, tzinfo=timezone.utc)

    # Case: full moon, mid-winter.
    an_date_string = "20161213"
    site_string = "BDO_Kansas"
    an = user.Astronight(an_date_string, site_string, -6)  # not daylight savings time
    assert an.an_date_string == an_date_string
    assert an.site_name == site_string  # the name
    assert an.site.name == site_string  # the object (need to expose though Astronight()?)
    assert an.localMidnightUT == datetime(2016, 12, 14, 6, 0, 0, tzinfo=timezone.utc)
    assert abs(an.localMidnightJD - 2457736.75) < 1 / 24 / 3600  # one sec tolerance
    target_lst_seconds = 5 * 3600 + 9 * 60 + 49
    an_lst_seconds = an.localMidnightLST * 240.0
    assert abs(target_lst_seconds - an_lst_seconds) < 1  # one sec tolerance
    assert abs(an.moon_radec.ra - 86.071) < 1 / 3600
    assert abs(an.moon_radec.dec - 18.297444) < 1 / 3600
    assert an.dark.start == datetime(2016, 12, 13, 23, 50, 19, 965103, tzinfo=timezone.utc)
    assert an.dark.end == datetime(2016, 12, 14, 12, 46, 22, 225678, tzinfo=timezone.utc)

    # Case: full moon, mid-summer
    an_date_string = "20160619"
    site_string = "BDO_Kansas"
    an = user.Astronight(an_date_string, site_string, -5)
    assert an.an_date_string == an_date_string
    assert an.site_name == site_string  # the name
    assert an.site.name == site_string  # the object (need to expose though Astronight()?)
    assert an.localMidnightUT == datetime(2016, 6, 20, 5, 0, tzinfo=timezone.utc)
    assert abs(an.localMidnightJD - 2457559.7083333) < 1 / 24 / 3600  # one sec tolerance
    target_lst_seconds = 16 * 3600 + 31 * 60 + 49
    an_lst_seconds = an.localMidnightLST * 240.0
    assert abs(target_lst_seconds - an_lst_seconds) < 1  # one sec tolerance
    assert abs(an.moon_radec.ra - 266.443917) < 1 / 3600
    assert abs(an.moon_radec.dec - -19.225917) < 1 / 3600
    assert an.dark.start == datetime(2016, 6, 20, 2, 45, 38, 994307, tzinfo=timezone.utc)
    assert an.dark.end == datetime(2016, 6, 20, 10, 4, 39, 440813, tzinfo=timezone.utc)

    # Case: new moon.
    an_date_string = "20160930"
    site_string = "BDO_Kansas"
    an = user.Astronight(an_date_string, site_string, -5)
    assert an.an_date_string == an_date_string
    assert an.site_name == site_string  # the name
    assert an.site.name == site_string  # the object (need to expose though Astronight()?)
    assert an.localMidnightUT == datetime(2016, 10, 1, 5, 0, 0, tzinfo=timezone.utc)
    assert abs(an.localMidnightJD - 2457662.7083333) < 1 / 24 / 3600  # one sec tolerance
    target_lst_seconds = 23 * 3600 + 17 * 60 + 54
    an_lst_seconds = an.localMidnightLST * 240.0
    assert abs(target_lst_seconds - an_lst_seconds) < 1  # one sec tolerance
    assert abs(an.moon_radec.ra - 190.5197916) < 1 / 3600
    assert abs(an.moon_radec.dec - -2.498417) < 1 / 3600
    assert an.dark.start == datetime(2016, 10, 1, 0, 50, 14, 943153, tzinfo=timezone.utc)
    assert an.dark.end == datetime(2016, 10, 1, 11, 36, 34, 354740, tzinfo=timezone.utc)
