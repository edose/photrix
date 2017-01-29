from datetime import datetime, timezone, timedelta
from pytest import approx
from photrix import user                 # call: user.fn() & user.Class()
from photrix.util import hex_degrees_as_degrees, ra_as_hours, dec_as_hex, RaDec
from photrix.fov import Fov

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
    print()
    # Test constructor, case = moon up at midnight. ----------------------------------------------
    an_date_string = "20160910"
    site_string = "BDO_Kansas"
    an = user.Astronight(an_date_string, site_string)
    assert an.an_date_string == an_date_string
    assert an.site_name == site_string  # the name

    assert abs((an.ts_dark.start -
                datetime(2016, 9, 11, 1, 33, 30, 944563, tzinfo=timezone.utc)).total_seconds()) < 1
    assert abs((an.ts_dark.end -
                datetime(2016, 9, 11, 11, 7, 11, 687656, tzinfo=timezone.utc)).total_seconds()) < 1

    target_local_middark_utc = datetime(2016, 9, 11, 6, 20, 21, tzinfo=timezone.utc)
    assert abs((an.local_middark_utc - target_local_middark_utc).total_seconds()) <= 1
    assert an.local_middark_jd == approx(2457642.764131, 1/(24*3600))  # one sec tolerance
    target_lst_seconds = 23*3600 + 19*60 + 37
    an_lst_seconds = an.local_middark_lst * 240.0
    assert target_lst_seconds == approx(an_lst_seconds, abs=1)

    target_moon_ra = (15/3600) * (18*3600+35*60+7.2)  # degrees
    target_moon_dec = -(19+1/60+24/3600)  # degrees
    assert an.moon_radec.ra == approx(target_moon_ra, abs=1/3600)  # tol = 1 arcsecond
    assert an.moon_radec.dec == approx(target_moon_dec, abs=1/3600)  # tol = 1 arcsecond
    assert an.moon_phase == approx(0.6722, abs=0.005)
    assert abs((an.ts_dark_no_moon.start -
                datetime(2016, 9, 11, 6, 36, 39, 829350, tzinfo=timezone.utc)).total_seconds()) < 1
    assert abs((an.ts_dark_no_moon.end -
                datetime(2016, 9, 11, 11, 7, 11, 687656, tzinfo=timezone.utc)).total_seconds()) < 1

    # Test constructor, case = full moon, mid-winter. --------------------------------------------
    an_date_string = "20161213"
    site_string = "BDO_Kansas"
    an = user.Astronight(an_date_string, site_string)
    assert an.an_date_string == an_date_string
    assert an.site_name == site_string  # the name

    assert abs((an.ts_dark.start -
                datetime(2016, 12, 14, 0, 1, 23, 877599,
                tzinfo=timezone.utc)).total_seconds()) < 1
    assert abs((an.ts_dark.end -
                datetime(2016, 12, 14, 12, 35, 17, 958638,
                tzinfo=timezone.utc)).total_seconds()) < 1

    target_local_middark_utc = datetime(2016, 12, 14, 6, 18, 20, 877599, tzinfo=timezone.utc)
    assert abs((an.local_middark_utc - target_local_middark_utc).total_seconds()) <= 1
    assert an.local_middark_jd == approx(2457736.762742, 1/(24*3600))  # one sec tolerance
    target_lst_seconds = 5*3600 + 28*60 + 13
    an_lst_seconds = an.local_middark_lst * 240.0
    assert target_lst_seconds == approx(an_lst_seconds, abs=1)

    target_moon_ra = (15/3600) * (5*3600+44*60+50)  # degrees
    target_moon_dec = +(18+18/60+39/3600)  # degrees
    assert an.moon_radec.ra == approx(target_moon_ra, abs=1/3600)  # tol = 1 arcsecond
    assert an.moon_radec.dec == approx(target_moon_dec, abs=1/3600)  # tol = 1 arcsecond
    assert an.moon_phase == approx(0.9973, abs=0.005)
    assert an.ts_dark_no_moon.seconds == 0

    # Test constructor, case = full moon, mid-summer ---------------------------------------------
    an_date_string = "20160619"
    site_string = "BDO_Kansas"
    an = user.Astronight(an_date_string, site_string)
    assert an.an_date_string == an_date_string
    assert an.site_name == site_string  # the name

    assert abs((an.ts_dark.start -
                datetime(2016, 6, 20, 2, 59, 26, 172446,
                tzinfo=timezone.utc)).total_seconds()) < 1
    assert abs((an.ts_dark.end -
                datetime(2016, 6, 20, 9, 50, 52, 262724,
                tzinfo=timezone.utc)).total_seconds()) < 1

    target_local_middark_utc = datetime(2016, 6, 20, 6, 25, 9, 172446, tzinfo=timezone.utc)
    assert abs((an.local_middark_utc - target_local_middark_utc).total_seconds()) <= 1
    assert an.local_middark_jd == approx(2457559.767467, 1/(24*3600))  # one sec tolerance
    target_lst_seconds = 17*3600 + 57*60 + 12
    an_lst_seconds = an.local_middark_lst * 240.0
    assert target_lst_seconds == approx(an_lst_seconds, abs=1)

    target_moon_ra = (15/3600) * (17*3600+47*60+48)  # degrees
    target_moon_dec = -(19+16/60+5.50/3600)  # degrees
    assert an.moon_radec.ra == approx(target_moon_ra, abs=1/3600)  # tol = 1 arcsecond
    assert an.moon_radec.dec == approx(target_moon_dec, abs=1/3600)  # tol = 1 arcsecond
    assert an.moon_phase == approx(0.9978, abs=0.005)
    assert an.ts_dark_no_moon.seconds == 0

    # Test constructor, case = new moon. ---------------------------------------------------------
    an_date_string = "20160930"
    site_string = "BDO_Kansas"
    an = user.Astronight(an_date_string, site_string)
    assert an.an_date_string == an_date_string
    assert an.site_name == site_string  # the name

    assert abs((an.ts_dark.start -
                datetime(2016, 10, 1, 1, 0, 32, 860858,
                tzinfo=timezone.utc)).total_seconds()) < 1
    assert abs((an.ts_dark.end -
                datetime(2016, 10, 1, 11, 26, 15, 525660,
                tzinfo=timezone.utc)).total_seconds()) < 1

    target_local_middark_utc = datetime(2016, 10, 1, 6, 13, 23, 860858, tzinfo=timezone.utc)
    assert abs((an.local_middark_utc - target_local_middark_utc).total_seconds()) <= 1
    assert an.local_middark_jd == approx(2457662.759304, 1/(24*3600))  # one sec tolerance
    target_lst_seconds = 0*3600 + 31*60 + 30
    an_lst_seconds = an.local_middark_lst * 240.0
    assert target_lst_seconds == approx(an_lst_seconds, abs=1)

    target_moon_ra = (15/3600) * (12*3600+45*60+16.08)  # degrees
    target_moon_dec = -(2+41/60+18.0/3600)  # degrees
    assert an.moon_radec.ra == approx(target_moon_ra, abs=1/3600)  # tol = 1 arcsecond
    assert an.moon_radec.dec == approx(target_moon_dec, abs=1/3600)  # tol = 1 arcsecond
    assert an.moon_phase == approx(0.0011, abs=0.005)
    assert an.ts_dark_no_moon == an.ts_dark

    # Test ts_observable(), set up.  =============================================================
    an_date_string = "20160919"
    site_string = "BDO_Kansas"
    an = user.Astronight(an_date_string, site_string)
    hip_116928 = RaDec('23:42:02.662', '+01:46:45.557')

    # Test ts_observable(), case = object farther than min dist from moon (common case).  --------
    ts_obs = an.ts_observable(hip_116928, min_moon_dist=45)
    assert abs((ts_obs.start - datetime(2016, 9, 20, 2, 13, 12, 660671,
                tzinfo=timezone.utc)).total_seconds()) <= 60
    assert abs((ts_obs.end - datetime(2016, 9, 20, 10, 3, 0, 540779,
                tzinfo=timezone.utc)).total_seconds()) <= 60

    # Test ts_observable(), case = object closer than min dist from moon. ------------------------
    ts_obs = an.ts_observable(hip_116928, min_moon_dist=90)
    # print("b", ts_obs, "\n\n\n")
    assert abs((ts_obs.start - datetime(2016, 9, 20, 2, 13, 12, 660671,
                tzinfo=timezone.utc)).total_seconds()) <= 60
    assert abs((ts_obs.end - datetime(2016, 9, 20, 2, 38, 2, 264284,
                tzinfo=timezone.utc)).total_seconds()) <= 60

    # Test ts_observable(), case = ignore moon altogether (set min_moon_dist to 0). --------------
    ts_obs = an.ts_observable(hip_116928, min_moon_dist=0)
    # print("c", ts_obs, "\n\n\n")
    assert abs((ts_obs.start - datetime(2016, 9, 20, 2, 13, 12, 660671,
                tzinfo=timezone.utc)).total_seconds()) <= 60
    assert abs((ts_obs.end - datetime(2016, 9, 20, 10, 3, 11, 540779,
                tzinfo=timezone.utc)).total_seconds()) <= 60

    # Test ts_observable(), case = disable observing any time moon is up at all. -----------------
    ts_obs = an.ts_observable(hip_116928, min_moon_dist=200)
    # print("d", ts_obs, "\n\n\n")
    assert abs((ts_obs.start - datetime(2016, 9, 20, 2, 13, 12, 660671,
                tzinfo=timezone.utc)).total_seconds()) <= 60
    assert abs((ts_obs.end - datetime(2016, 9, 20, 2, 38, 2, 264284,
                tzinfo=timezone.utc)).total_seconds()) <= 60

    # Test ts_observable(), case = object farther than min dist from moon, higher min_alt. ------
    ts_obs = an.ts_observable(hip_116928, min_moon_dist=45, min_alt=35)
    # print("e", ts_obs, "\n\n\n")
    assert abs((ts_obs.start - datetime(2016, 9, 20, 3, 9, 53, 907299,
                tzinfo=timezone.utc)).total_seconds()) <= 60
    assert abs((ts_obs.end - datetime(2016, 9, 20, 9, 6, 30, 294148,
                tzinfo=timezone.utc)).total_seconds()) <= 60

    # Test ts_observable(), case = object closer than min dist from moon, higher min_alt. --------
    ts_obs = an.ts_observable(hip_116928, min_moon_dist=90, min_alt=35)
    # print("f", ts_obs, "\n\n\n")
    assert ts_obs.seconds == 0  # start and end times are unimportant (indeed they are undefined).

    # For remaining tests, assume Astronight object's returned values are ok (as were tested above).
    # Wide range of sky positions, to exercise all functions and exception handling.
    altais = RaDec('19:12:33.405', '+67:39:43.092')  # in NW sky (from North America)
    hip_22783 = RaDec('04:54:03.012', '+66:20:33.763')  # in NE sky
    mira = RaDec('02:19:20.804', '-02:58:43.518')  # in SE sky
    algedi = RaDec('20:18:03.324', '-12:32:41.419')  # in SW sky
    ankaa = RaDec('00:26:17.310', '-42:18:27.446')  # too far south to observe
    polaris = RaDec('02:31:49.133', '+89:15:50.598')  # circumpolar north

    # All targets, allow any moon.
    ts_obs = an.ts_observable(altais, min_moon_dist=0, min_alt=25)
    # print("altais / any moon >> ", ts_obs, "\n\n\n")
    assert abs((ts_obs.start - an.ts_dark.start).total_seconds()) <= 60
    assert abs((ts_obs.end - datetime(2016, 9, 20, 9, 47, 48, 586041,
                tzinfo=timezone.utc)).total_seconds()) <= 60
    ts_obs = an.ts_observable(hip_22783, min_moon_dist=0, min_alt=25)
    # print("hip_22783 / any moon >> ", ts_obs, "\n\n\n")
    assert abs((ts_obs.start - datetime(2016, 9, 20, 3, 23, 32, 838379,
                tzinfo=timezone.utc)).total_seconds()) <= 60
    assert abs((ts_obs.end - an.ts_dark.end).total_seconds()) <= 60
    ts_obs = an.ts_observable(mira, min_moon_dist=0, min_alt=25)
    # print("mira / any moon >> ", ts_obs, "\n\n\n")
    assert abs((ts_obs.start - datetime(2016, 9, 20, 5, 8, 37, 879958,
                tzinfo=timezone.utc)).total_seconds()) <= 60
    assert abs((ts_obs.end - an.ts_dark.end).total_seconds()) <= 60
    ts_obs = an.ts_observable(algedi, min_moon_dist=0, min_alt=25)
    # print("algedi / any moon >> ", ts_obs, "\n\n\n")
    assert abs((ts_obs.start - an.ts_dark.start).total_seconds()) <= 60
    assert abs((ts_obs.end - datetime(2016, 9, 20, 5, 35, 16, 643651,
                tzinfo=timezone.utc)).total_seconds()) <= 60
    ts_obs = an.ts_observable(ankaa, min_moon_dist=0, min_alt=25)
    # print("ankaa / any moon >> ", ts_obs, "\n\n\n")
    assert ts_obs.seconds == 0
    ts_obs = an.ts_observable(polaris, min_moon_dist=0, min_alt=25)
    # print("polaris / any moon >> ", ts_obs, "\n\n\n")
    assert ts_obs == an.ts_dark

    # All targets, allow NO moon.
    ts_obs = an.ts_observable(altais, min_moon_dist=220, min_alt=25)
    # print("altais / NO moon >> ", ts_obs, "\n\n\n")
    assert ts_obs == an.ts_dark_no_moon
    ts_obs = an.ts_observable(hip_22783, min_moon_dist=220, min_alt=25)
    # print("hip_22783 / NO moon >> ", ts_obs, "\n\n\n")
    assert ts_obs.seconds == 0
    ts_obs = an.ts_observable(mira, min_moon_dist=220, min_alt=25)
    # print("mira / NO moon >> ", ts_obs, "\n\n\n")
    assert ts_obs.seconds == 0
    ts_obs = an.ts_observable(algedi, min_moon_dist=220, min_alt=25)
    # print("algedi / NO moon >> ", ts_obs, "\n\n\n")
    assert ts_obs == an.ts_dark_no_moon
    ts_obs = an.ts_observable(ankaa, min_moon_dist=220, min_alt=25)
    # print("ankaa / NO moon >> ", ts_obs, "\n\n\n")
    assert ts_obs.seconds == 0
    ts_obs = an.ts_observable(polaris, min_moon_dist=220, min_alt=25)
    # print("polaris / NO moon >> ", ts_obs, "\n\n\n")
    assert ts_obs == an.ts_dark_no_moon

    # Continue testing .ts_observable():  --------------------------------------------------------
    # Same targets, new astronight with ~ opposite moon phase to previous astronight.
    an_date_string = "20161008"
    site_string = "BDO_Kansas"
    an = user.Astronight(an_date_string, site_string)

    # Case: allow any moon.
    ts_obs = an.ts_observable(altais, min_moon_dist=0, min_alt=25)
    # print("altais / any moon >> ", ts_obs, "\n\n\n")
    assert abs((ts_obs.start - an.ts_dark.start).total_seconds()) <= 60
    assert abs((ts_obs.end - datetime(2016, 10, 9, 8, 33, 5, 419491,
               tzinfo=timezone.utc)).total_seconds()) <= 60
    ts_obs = an.ts_observable(hip_22783, min_moon_dist=0, min_alt=25)
    # print("hip_22783 / any moon >> ", ts_obs, "\n\n\n")
    assert abs((ts_obs.start - datetime(2016, 10, 9, 2, 8, 51, 471332,
               tzinfo=timezone.utc)).total_seconds()) <= 60
    assert abs((ts_obs.end - an.ts_dark.end).total_seconds()) <= 60
    ts_obs = an.ts_observable(mira, min_moon_dist=0, min_alt=25)
    # print("mira / any moon >> ", ts_obs, "\n\n\n")
    assert abs((ts_obs.start - datetime(2016, 10, 9, 3, 53, 55, 971837,
               tzinfo=timezone.utc)).total_seconds()) <= 60
    assert abs((ts_obs.end - datetime(2016, 10, 9, 11, 6, 47, 395358,
               tzinfo=timezone.utc)).total_seconds()) <= 60
    ts_obs = an.ts_observable(algedi, min_moon_dist=0, min_alt=25)
    # print("algedi / any moon >> ", ts_obs, "\n\n\n")
    assert abs((ts_obs.start - an.ts_dark.start).total_seconds()) <= 60
    assert abs((ts_obs.end - datetime(2016, 10, 9, 4, 20, 34, 71309,
               tzinfo=timezone.utc)).total_seconds()) <= 60
    ts_obs = an.ts_observable(ankaa, min_moon_dist=0, min_alt=25)
    # print("ankaa / any moon >> ", ts_obs, "\n\n\n")
    assert ts_obs.seconds == 0
    ts_obs = an.ts_observable(polaris, min_moon_dist=0, min_alt=25)
    # print("polaris / any moon >> ", ts_obs, "\n\n\n")
    assert ts_obs == an.ts_dark

    # Case: allow NO moon:
    ts_obs = an.ts_observable(altais, min_moon_dist=220, min_alt=25)
    # print("altais / NO moon >> ", ts_obs, "\n\n\n")
    assert abs((ts_obs.start - an.ts_dark_no_moon.start).total_seconds()) <= 60
    assert abs((ts_obs.end - datetime(2016, 10, 9, 8, 33, 5, 419491,
               tzinfo=timezone.utc)).total_seconds()) <= 60
    ts_obs = an.ts_observable(hip_22783, min_moon_dist=220, min_alt=25)
    # print("hip_22783 / NO moon >> ", ts_obs, "\n\n\n")
    assert ts_obs == an.ts_dark_no_moon
    ts_obs = an.ts_observable(mira, min_moon_dist=220, min_alt=25)
    # print("mira / NO moon >> ", ts_obs, "\n\n\n")
    assert abs((ts_obs.start - an.ts_dark_no_moon.start).total_seconds()) <= 60
    assert abs((ts_obs.end - datetime(2016, 10, 9, 11, 6, 47, 395358,
               tzinfo=timezone.utc)).total_seconds()) <= 60
    ts_obs = an.ts_observable(algedi, min_moon_dist=220, min_alt=25)
    # print("algedi / NO moon >> ", ts_obs, "\n\n\n")
    assert ts_obs.seconds == 0
    ts_obs = an.ts_observable(ankaa, min_moon_dist=220, min_alt=25)
    # print("ankaa / NO moon >> ", ts_obs, "\n\n\n")
    assert ts_obs.seconds == 0
    ts_obs = an.ts_observable(polaris, min_moon_dist=220, min_alt=25)
    # print("polaris / NO moon >> ", ts_obs, "\n\n\n")
    assert ts_obs == an.ts_dark_no_moon

    # Case: object near sun (wholly unobservable) and during new moon (moon matters little):
    an_date_string = "20160930"
    site_string = "BDO_Kansas"
    an = user.Astronight(an_date_string, site_string)
    porrima = RaDec('12:41:38.954', '-01:26:56.733')
    ts_obs = an.ts_observable(porrima, min_moon_dist=0, min_alt=25)
    assert ts_obs.seconds == 0
    ts_obs = an.ts_observable(porrima, min_moon_dist=220, min_alt=25)
    assert ts_obs.seconds == 0
    ts_obs = an.ts_observable(porrima, min_moon_dist=0, min_alt=2)
    assert ts_obs.seconds == 0
    ts_obs = an.ts_observable(porrima, min_moon_dist=220, min_alt=2)
    assert ts_obs.seconds == 0

    # Case: object circumpolar and during new moon:
    polaris = RaDec('02:31:49.133', '+89:15:50.598')  # circumpolar north
    ts_obs = an.ts_observable(polaris, min_moon_dist=0, min_alt=25)
    assert ts_obs == an.ts_dark
    ts_obs = an.ts_observable(polaris, min_moon_dist=220, min_alt=25)
    assert ts_obs == an.ts_dark_no_moon
    ts_obs = an.ts_observable(polaris, min_moon_dist=0, min_alt=2)
    assert ts_obs == an.ts_dark
    ts_obs = an.ts_observable(polaris, min_moon_dist=220, min_alt=2)
    assert ts_obs == an.ts_dark_no_moon
    ts_obs = an.ts_observable(polaris, min_moon_dist=0, min_alt=60)  # it's never this high.
    assert ts_obs.seconds == 0
    ts_obs = an.ts_observable(polaris, min_moon_dist=220, min_alt=60)
    assert ts_obs.seconds == 0

    # Test ts_fov_observable (based on ts_observable() having been tested above):
    an_date_string = "20160930"
    site_string = "BDO_Kansas"
    an = user.Astronight(an_date_string, site_string)
    fov = Fov('ST Tri')
    ts_from_radec = an.ts_observable(RaDec(fov.ra, fov.dec), min_moon_dist=50, min_alt=40)
    ts_from_fov = an.ts_fov_observable(fov, min_moon_dist=50, min_alt=40)
    assert ts_from_fov == ts_from_radec

    # Test .datetime_utc_from_hhmm():
    an_date_string = "20160930"
    site_string = "BDO_Kansas"
    an = user.Astronight(an_date_string, site_string)
    assert datetime_utc_from_hhmm_OK('2323', an)
    assert datetime_utc_from_hhmm_OK('2020', an)
    assert datetime_utc_from_hhmm_OK('2000', an)
    assert datetime_utc_from_hhmm_OK('0000', an)
    assert datetime_utc_from_hhmm_OK('0600', an)
    assert datetime_utc_from_hhmm_OK('0900', an)

    an_date_string = "20170101"
    site_string = "BDO_Kansas"
    an = user.Astronight(an_date_string, site_string)
    assert datetime_utc_from_hhmm_OK('2323', an)
    assert datetime_utc_from_hhmm_OK('2020', an)
    assert datetime_utc_from_hhmm_OK('2000', an)
    assert datetime_utc_from_hhmm_OK('0000', an)
    assert datetime_utc_from_hhmm_OK('0600', an)
    assert datetime_utc_from_hhmm_OK('0900', an)


# ----- HELPER FUNCTIONS. --------------------------------------------------------------------
def datetime_utc_from_hhmm_OK(hhmm_string, an):
    dt = an.datetime_utc_from_hhmm(hhmm_string)
    hh = int(hhmm_string[0:2])
    mm = int(hhmm_string[2:4])
    days_diff = abs(dt-an.local_middark_utc).total_seconds() / (24*3600)
    return dt.hour == hh and dt.minute == mm and days_diff <= 0.5
