import os
from photrix import fov
from photrix import util

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"

PHOTRIX_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_FOV_DIRECTORY = os.path.join(PHOTRIX_ROOT_DIRECTORY, "test", "$fovs_for_test")

# For Schema 1.3 of Sept 2016 (as defined initially in R code "photometry").


def test_FOV_init_eclipser():
    """ Exhaustive test, #STARS section coverage (which is needed only this once)."""
    fov_name = "ST Tri"
    f = fov.Fov("ST Tri", TEST_FOV_DIRECTORY)
    assert f.fov_name == fov_name
    assert f.format_version == "1.3"
    assert f.ra == util.ra_as_degrees("02:42:00.0")
    assert f.dec == util.dec_as_degrees("+35:43:31")
    assert f.chart == "X15646NP"
    assert f.fov_date == "12/21/2015"
    assert f.main_target == "ST Tri"
    assert f.target_type == "Eclipser"
    assert f.period == float("0.47905145")
    assert f.JD_bright == float("2457360.45")
    assert f.JD_faint == float("2457360.69")
    assert f.mag_V_bright == 14
    assert f.mag_V_faint == 14.77
    assert f.color_VI_bright == 0.5
    assert f.color_VI_faint == 0.55
    assert f.observing_style == "Stare"
    assert f.alert is None
    assert len(f.observing_list) == 2
    assert f.observing_list == [("I", 12, 1), ("V", 13, 6)]
    assert f.stare_reference == "MIN"
    assert (f.stare_start, f.stare_stop) == (-2, 2)
    assert f.max_exposure == 240
    assert f.priority == 8
    assert f.gap_score_days == [60, 90, 150]
    assert f.acp_directives == ['#REPEAT 50', '#DITHER 0', '#FILTER I,V',
                                '#BINNING 1,1', "#COUNT 1,5", "#INTERVAL 120,240"]
    assert f.acp_comments == 'Eclipser EB  V=14-14.77  NO R COMPS'
    assert f.punches == [("143", -7.04, 11.94)]
    assert len(f.aavso_stars) == 11
    assert all([star.is_valid for star in f.aavso_stars])
    s = f.aavso_stars[0]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("156", 40.60386, 35.67683, "comp", 16.31, 15.608, None, None, None)
    s = f.aavso_stars[1]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("127", 40.55951, 35.68621, "comp", 13.556, 12.673, None, 11.697, None)
    s = f.aavso_stars[2]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("137", 40.55655, 35.71275, "check", 14.255, 13.676, None, 13.113, None)
    s = f.aavso_stars[3]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("147", 40.61729, 35.67576, "comp", 15.466, 14.673, None, None, None)
    s = f.aavso_stars[4]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("101", 40.46313, 35.76667, "comp", 10.578, 10.075, None, 9.503, None)
    s = f.aavso_stars[5]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("143", 40.56636, 35.69035, "comp", 15.441, 14.288, None, None, None)
    s = f.aavso_stars[6]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("131", 40.64438, 35.76555, "comp", 13.597, 13.087, None, 12.463, None)
    s = f.aavso_stars[7]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("139", 40.64442, 35.69147, "comp", 14.533, 13.9, None, None, None)
    s = f.aavso_stars[8]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("151", 40.62392, 35.63252, "comp", 16.098, 15.073, None, None, None)
    s = f.aavso_stars[9]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("ST Tri", 40.38631, 35.72553, "target", None, None, None, None, None)
    s = f.aavso_stars[10]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("V0680 Per", 40.42055, 35.71548, "target", None, None, None, None, None)


def test_FOV_init_mira():
    """ NB: #STARS section covered not here, but in test_FOV_eclipser()."""
    fov_name = "AU Aur"
    f = fov.Fov("AU Aur", TEST_FOV_DIRECTORY)
    assert f.fov_name == fov_name
    assert f.format_version == "1.3"
    assert f.ra == util.ra_as_degrees("04:54:15.0")
    assert f.dec == util.dec_as_degrees("+49:54:00")
    assert f.chart == "X15865BK"
    assert f.fov_date == "1/30/2016"
    assert f.main_target == "AU Aur"
    assert f.target_type == "Mira"
    assert f.period == float("400")
    assert f.JD_bright == float("2456520")
    assert f.JD_faint == float("2456720")
    assert f.mag_V_bright == 10.2
    assert f.mag_V_faint == 12.7
    assert f.color_VI_bright == 3.8
    assert f.color_VI_faint == 5
    assert f.observing_style == "LPV"
    assert len(f.observing_list) == 2
    assert f.alert is None
    assert f.observing_list == [("V", None, 1), ("I", None, 1)]
    assert f.stare_reference is None
    assert (f.stare_start, f.stare_stop) == (None, None)
    assert f.max_exposure is None
    assert f.priority == 8
    assert f.gap_score_days == [4, 8, 20]
    assert f.acp_directives == ['#DITHER 0', '#FILTER V,I', '#BINNING 1,1',
                                "#COUNT 1,1", "#INTERVAL 240,30"]
    assert f.acp_comments == 'Mira CARBON C6-7(N) V~10.2-12.6? very uncertain color, NO R COMPS'
    assert f.punches == [("AU Aur", 9.33, 2.71)]
    assert len(f.aavso_stars) == 11
    assert all([star.is_valid for star in f.aavso_stars])
    assert [star.star_id for star in f.aavso_stars] == \
           ['117', 'AU Aur', '141', '154', '132', '146', '124', '118', '155', '107', '112']


def test_FOV_mira_init_multi_exposure():
    """ NB: #STARS section covered not here, but in test_FOV_eclipser()."""
    fov_name = "AU Aur Modified"
    f = fov.Fov("AU Aur Modified", TEST_FOV_DIRECTORY)
    assert f.fov_name == fov_name
    assert f.format_version == "1.3"
    assert f.observing_style == "LPV"
    assert len(f.observing_list) == 2
    assert f.alert is None
    assert f.observing_list == [("V", None, 5), ("I", None, 2)]


def test_FOV_init_delta_scuti():
    """ NB: #STARS section covered not here, but in test_FOV_eclipser()."""
    fov_name = "ASAS J162540_19123"
    f = fov.Fov("ASAS J162540_19123", TEST_FOV_DIRECTORY)
    assert f.fov_name == fov_name
    assert f.format_version == "1.3"
    assert f.ra == util.ra_as_degrees("16:25:39.0")
    assert f.dec == util.dec_as_degrees("+19:03:18")
    assert f.chart == "X16073AJD"
    assert f.fov_date == "5/20/2016"
    assert f.main_target == "ASAS J162540+1912.3"
    assert f.target_type == "Delta Scuti"
    assert f.period == float("0.0571721")
    assert f.JD_bright == float("2457614.69717")
    assert f.JD_faint == float("2457614.67845")
    assert f.mag_V_bright == 12.38
    assert f.mag_V_faint == 12.74
    assert f.color_VI_bright == 0.25
    assert f.color_VI_faint == 0.45
    assert f.observing_style == "Stare"
    assert len(f.observing_list) == 2
    assert f.observing_list == [("I", 11, 1), ("V", 12, 16)]
    assert f.alert is None
    assert f.stare_reference == "ANY"
    assert (f.stare_start, f.stare_stop) == (0, 3)
    assert f.max_exposure == 80
    assert f.priority == 6
    assert f.gap_score_days == [60, 90, 150]
    assert f.acp_directives == ['#REPEAT 50', '#DITHER 0', '#FILTER I,V', '#BINNING 1,1',
                                "#COUNT 1,16", "#INTERVAL 60,60"]
    assert f.acp_comments == 'Eclipser DSCT  V=12.38-12.74  P=82 min  NO R COMPS'
    assert f.punches == []
    assert len(f.aavso_stars) == 8
    assert all([star.is_valid for star in f.aavso_stars])
    assert [star.star_id for star in f.aavso_stars] == \
           ['94', '110_1', '110', 'ASAS J162540+1912.3', 'U Her', '70', '114', '132']


def test_FOV_init_Z_Cam():
    """ NB: #STARS section covered not here, but in test_FOV_eclipser()."""
    fov_name = "NSV 14581"
    f = fov.Fov("NSV 14581", TEST_FOV_DIRECTORY)
    assert f.fov_name == fov_name
    assert f.format_version == "1.3"
    assert f.ra == util.ra_as_degrees("23:26:50.3")
    assert f.dec == util.dec_as_degrees("+82:22:11")
    assert f.chart == "X16224AYI"
    assert f.fov_date == "9/2/2016"
    assert f.main_target == "NSV 14581"
    assert f.target_type == "Z Cam"
    assert f.period == float("0.194334535")
    assert f.JD_bright == float("2452888.328")
    assert f.JD_faint == float("2452888.4255")
    assert f.mag_V_bright == 14
    assert f.mag_V_faint == 17.2
    assert f.color_VI_bright == 0.7
    assert f.color_VI_faint == 0.7
    assert f.observing_style == "Monitor"
    assert f.alert == 2
    assert len(f.observing_list) == 2
    assert f.observing_list == [("I", 13, 1), ("V", 14, 1)]
    assert f.stare_reference is None
    assert (f.stare_start, f.stare_stop) == (None, None)
    assert f.max_exposure is None
    assert f.priority == 3
    assert f.gap_score_days == [3, 5, 13]
    assert f.acp_directives == ['#DITHER 0', '#FILTER I,V', '#BINNING 1,1',
                                "#COUNT 1,1", "#INTERVAL 240,240"]
    assert f.acp_comments == 'Z Cam: gap-fill=5 days  V~15??'
    assert f.punches == []
    assert len(f.aavso_stars) == 13
    assert all([star.is_valid for star in f.aavso_stars])
    assert [star.star_id for star in f.aavso_stars] == \
           ['157', 'NSV 14581', '136', '148', '145', '151', '110',
            '142', '155', '134', '140', '153', '112']

def test_FOV_init_exoplanet():
    """ NB: #STARS section covered not here, but in test_FOV_eclipser()."""
    fov_name = "HAT_P_3"
    f = fov.Fov("HAT_P_3", TEST_FOV_DIRECTORY)
    assert f.fov_name == fov_name
    assert f.format_version == "1.3"
    assert f.main_target == "HAT-P-3"
    assert f.target_type == "Exoplanet"
    assert f.observing_style == "Stare"
    assert f.alert is None
    assert len(f.observing_list) == 4
    assert f.observing_list == [("R", 11.1, 8), ("V", 11.55, 1), ("R", 11.1, 8), ("I", 10.5, 1)]
    assert f.stare_reference == 'MIN'
    assert (f.stare_start, f.stare_stop) == (-2, 2)
    assert f.max_exposure == 80
    assert f.priority == 0
    assert f.gap_score_days == [60, 90, 150]
    assert f.acp_directives == ['#REPEAT 100', '#DITHER 0', '#FILTER R,V,R,I', '#BINNING 1,1,1,1',
                                "#COUNT 8,1,8,1", "#INTERVAL 40,60,40,60"]
    assert f.acp_comments == 'Exoplanet  V=11.86-11.876  1h24m duration // OBSERVE in R filter'
    assert f.punches == [('133', -9.58, 4.50), ('HAT-P-3', -4.58, 8.53)]
    assert len(f.aavso_stars) == 10
    assert all([star.is_valid for star in f.aavso_stars])
    assert [star.star_id for star in f.aavso_stars] == \
           ['125', '127', 'HAT-P-3', '131', '128', '140', '123', '133', '112', '107']


def test_FOV_init_standard():
    """ NB: #STARS section covered not here, but in test_FOV_eclipser()."""
    fov_name = "Std_SA100"
    f = fov.Fov("Std_SA100", TEST_FOV_DIRECTORY)
    assert f.fov_name == fov_name
    assert f.format_version == "1.3"
    assert f.ra == util.ra_as_degrees("08:53:14.3")
    assert f.dec == util.dec_as_degrees("-00:37:56")
    assert f.chart == "X15687X"
    assert f.fov_date == "12/20/2015"
    assert f.period is None
    assert f.JD_bright is None
    assert f.JD_faint is None
    assert f.mag_V_bright is None
    assert f.mag_V_faint is None
    assert f.color_VI_bright is None
    assert f.color_VI_faint is None
    assert f.main_target == "Standard"
    assert f.target_type == "Standard"
    assert f.observing_style == "Standard"
    assert f.alert is None
    assert len(f.observing_list) == 3
    assert f.observing_list == [("V", 11.5, 1), ("R", 11.2, 1), ("I", 10.2, 1)]
    assert f.stare_reference is None
    assert (f.stare_start, f.stare_stop) == (None, None)
    assert f.max_exposure is None
    assert f.priority is None
    assert f.gap_score_days is None
    assert f.acp_directives == ['#DITHER 0', '#FILTER V,R,I', '#BINNING 1,1,1',
                                "#COUNT 1,1,1", "#INTERVAL 60,60,60"]
    assert f.acp_comments == 'Standard field'
    assert f.punches == []
    assert len(f.aavso_stars) == 6
    assert all([star.is_valid for star in f.aavso_stars])
    assert [star.star_id for star in f.aavso_stars] == \
           ['118', '130', '101', '114', '124', '92']


def test_calc_gap_score():
    fov_name = "ST Tri"
    f = fov.Fov("ST Tri", TEST_FOV_DIRECTORY)
    assert f.fov_name == fov_name
    assert f.format_version == "1.3"
    gsd = f.gap_score_days
    print(str(gsd))
    assert f.calc_gap_score(0) == 0
    assert f.calc_gap_score(gsd[0] / 2.0) == 0
    assert f.calc_gap_score(gsd[0]) == 0
    assert f.calc_gap_score((gsd[0] + gsd[1]) / 2.0) == 0.5
    assert f.calc_gap_score(gsd[1]) == 1
    assert f.calc_gap_score((gsd[1] + gsd[2]) / 2.0) == 1.5
    assert f.calc_gap_score(gsd[2]) == 2
    assert f.calc_gap_score(1.5 * gsd[2]) == 2
    assert f.calc_gap_score(10 * gsd[2]) == 2


def test_estimate_mira_mags():
    fov_name = "AU Aur"
    f = fov.Fov("AU Aur", TEST_FOV_DIRECTORY)
    assert f.fov_name == fov_name
    assert f.format_version == "1.3"
    mag = f.estimate_mira_mags(f.JD_bright)
    assert mag['V'] == f.mag_V_bright
    assert mag['I'] == f.mag_V_bright + f.color_VI_bright
    mag = f.estimate_mira_mags(f.JD_faint)
    assert mag['V'] == f.mag_V_faint
    assert mag['I'] == f.mag_V_faint + f.color_VI_faint
    #  TODO fill out this test with (1) 'R' mags, and (2) more JD examples.


def test_all_fov_names():
    all_names = fov.all_fov_names(fov_directory=TEST_FOV_DIRECTORY)
    assert isinstance(all_names, list)
    assert isinstance(all_names[0], str)
    assert len(all_names) == 7


def test_FovDict_init():
    fov_dict = fov.FovDict(fov_directory=TEST_FOV_DIRECTORY)
    inner_dict = fov_dict.fov_dict
    assert isinstance(fov_dict, fov.FovDict)
    assert fov_dict.count() == 7
    assert isinstance(inner_dict, dict)
    assert isinstance(list(inner_dict.keys())[0], str)
    assert isinstance(list(inner_dict.values())[0], fov.Fov)
    assert len(inner_dict) == 7
    non_zero = fov_dict.no_zero_priority()
    assert non_zero.count() == 6
