import os
import pytest
from photrix import fov
from photrix import util

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"

PHOTRIX_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_FOV_DIRECTORY = os.path.join(PHOTRIX_ROOT_DIRECTORY, "test", "$fovs_for_test")
CURRENT_SCHEMA_VERSION = "1.4"  # For Schema 1.4 of Jan. 2017.


def test_fov():
    """
    This is for observing style = "Stare".
    This test function includes all general-case tests
        (which are not repeated in test functions for other observing styles).
    """
    fovname1 = "ST Tri"
    fov1 = fov.Fov(fovname1, fov_directory=TEST_FOV_DIRECTORY)
    # Test fields.
    assert fov1.fov_name == fovname1
    assert fov1.format_version == CURRENT_SCHEMA_VERSION
    assert fov1.ra == 2*15 + 42*15/60
    assert fov1.dec == 35 + 43/60 + 31/3600
    assert fov1.chart == "X15646NP"
    assert fov1.fov_date == "12/21/2015"
    assert fov1.main_target == "ST Tri"
    assert fov1.target_type == "Eclipser"
    assert fov1.period == 0.47905145
    assert fov1.motive == 'Test only.'
    assert fov1.JD_bright == 2457360.57
    assert fov1.JD_faint == 2457360.69
    assert fov1.JD_second == 2457360.45
    assert fov1.mag_V_bright == 14
    assert fov1.mag_V_faint == 14.77
    assert fov1.color_VI_bright == 0.5
    assert fov1.color_VI_faint == 0.55
    assert fov1.observing_style == "Stare"
    assert len(fov1.observing_list) == 2
    assert fov1.observing_list == [("I", 12, 1), ("V", 13, 6)]
    assert fov1.alert is None
    assert fov1.stare_reference == "MIN"
    assert fov1.stare_start == -2
    assert fov1.stare_stop == 2
    assert fov1.max_exposure == 240
    assert fov1.priority == 8
    assert fov1.gap_score_days == [60, 90, 150]
    assert fov1.acp_directives == ["#REPEAT 50", "#DITHER 0", "#FILTER I,V", "#BINNING 1,1",
                                   "#COUNT 1,5", "#INTERVAL 120,240"]
    assert fov1.acp_comments == "Eclipser EB  V=14-14.77  NO R COMPS"
    assert fov1.punches == [("143", -7.04, 11.94)]
    assert str(fov1) == "FOV 'ST Tri' with 11 sequence stars."
    assert len(fov1.aavso_stars) == 11

    assert len(fov1.aavso_stars) == 11
    assert all([star.is_valid for star in fov1.aavso_stars])
    s = fov1.aavso_stars[0]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("156", 40.60386, 35.67683, "comp", 16.31, 15.608, None, None, None)
    s = fov1.aavso_stars[1]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("127", 40.55951, 35.68621, "comp", 13.556, 12.673, None, 11.697, None)
    s = fov1.aavso_stars[2]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("137", 40.55655, 35.71275, "check", 14.255, 13.676, None, 13.113, None)
    s = fov1.aavso_stars[3]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("147", 40.61729, 35.67576, "comp", 15.466, 14.673, None, None, None)
    s = fov1.aavso_stars[4]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("101", 40.46313, 35.76667, "comp", 10.578, 10.075, None, 9.503, None)
    s = fov1.aavso_stars[5]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("143", 40.56636, 35.69035, "comp", 15.441, 14.288, None, None, None)
    s = fov1.aavso_stars[6]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("131", 40.64438, 35.76555, "comp", 13.597, 13.087, None, 12.463, None)
    s = fov1.aavso_stars[7]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("139", 40.64442, 35.69147, "comp", 14.533, 13.9, None, None, None)
    s = fov1.aavso_stars[8]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("151", 40.62392, 35.63252, "comp", 16.098, 15.073, None, None, None)
    s = fov1.aavso_stars[9]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("ST Tri", 40.38631, 35.72553, "target", None, None, None, None, None)
    s = fov1.aavso_stars[10]
    assert (s.star_id, s.ra, s.dec, s.star_type, s.magB, s.magV, s.magR, s.magI, s.magU) == \
           ("V0680 Per", 40.42055, 35.71548, "target", None, None, None, None, None)

    # Test calc_gap_score().
    assert fov1.calc_gap_score(0) == 0
    assert fov1.calc_gap_score(60) == 0
    assert fov1.calc_gap_score(75) == 0.5
    assert fov1.calc_gap_score(90) == 1
    assert fov1.calc_gap_score(120) == 1.5
    assert fov1.calc_gap_score(150) == 2
    assert fov1.calc_gap_score(365) == 2
    assert fov1.calc_gap_score(-1) == 0
    # Test calc_priority_score().
    ps = fov1.priority
    assert fov1.calc_priority_score(0) == 0
    assert fov1.calc_priority_score(60) == 0
    assert fov1.calc_priority_score(75) == 0.5 * ps
    assert fov1.calc_priority_score(90) == 1 * ps
    assert fov1.calc_priority_score(120) == 1.5 * ps
    assert fov1.calc_priority_score(150) == 2 * ps
    assert fov1.calc_priority_score(365) == 2 * ps
    assert fov1.calc_priority_score(-1) == 0


def test_fov_mira():
    """ NB: For tests on #STARS section, see test_fov() above."""
    fovname2 = "AU Aur"
    fov2 = fov.Fov(fovname2, fov_directory=TEST_FOV_DIRECTORY)
    # Test fields.
    assert fov2.fov_name == fovname2
    assert fov2.format_version == CURRENT_SCHEMA_VERSION
    assert fov2.main_target == "AU Aur"
    assert fov2.target_type == "Mira"
    assert fov2.period == 400
    assert fov2.motive == 'Test only: AU Aur.'
    assert fov2.JD_bright == 2456520
    assert fov2.JD_faint == 2456720
    assert fov2.observing_style == "LPV"
    assert len(fov2.observing_list) == 2
    assert fov2.observing_list[0] == ("V", None, 1)
    assert fov2.observing_list[1] == ("I", None, 1)
    assert fov2.alert is None
    assert fov2.stare_reference is None
    assert fov2.stare_start is None
    assert fov2.stare_stop is None
    assert fov2.max_exposure is None  # this needs to be given in real FOV file
    assert fov2.priority == 8
    assert fov2.gap_score_days == [4, 8, 20]
    assert fov2.acp_directives == ["#DITHER 0", "#FILTER V,I", "#BINNING 1,1",
                                   "#COUNT 1,1", "#INTERVAL 240,30"]
    assert fov2.acp_comments == "Mira CARBON C6-7(N) V~10.2-12.6? very uncertain color, NO R COMPS"
    assert fov2.punches == [("AU Aur", 9.33, 2.71)]
    assert str(fov2) == "FOV 'AU Aur' with 11 sequence stars."
    assert len(fov2.aavso_stars) == 11
    assert all([star.is_valid for star in fov2.aavso_stars])
    assert [star.star_id for star in fov2.aavso_stars] == \
           ['117', 'AU Aur', '141', '154', '132', '146', '124', '118', '155', '107', '112']

    # Test .estimate_lpv_mags().
    jd = fov2.JD_bright
    mags_bright = fov2.estimate_lpv_mags(jd)
    assert mags_bright['V'] == fov2.mag_V_bright
    assert mags_bright['R'] == fov2.mag_V_bright - 0.5*fov2.color_VI_bright
    assert mags_bright['I'] == fov2.mag_V_bright - fov2.color_VI_bright
    jd = fov2.JD_faint
    mags_faint = fov2.estimate_lpv_mags(jd)
    assert mags_faint['V'] == fov2.mag_V_faint
    assert mags_faint['R'] == fov2.mag_V_faint - 0.5*fov2.color_VI_faint
    assert mags_faint['I'] == fov2.mag_V_faint - fov2.color_VI_faint
    jd = fov2.JD_bright + 0.5 * (fov2.JD_faint - fov2.JD_bright)  # mid-dimming
    mags_jd = fov2.estimate_lpv_mags(jd)
    assert mags_jd['V'] == mags_bright['V'] + 0.5*(mags_faint['V']-mags_bright['V'])
    assert mags_jd['R'] == mags_bright['R'] + 0.5*(mags_faint['R']-mags_bright['R'])
    assert mags_jd['I'] == mags_bright['I'] + 0.5*(mags_faint['I']-mags_bright['I'])
    jd = fov2.JD_faint + 0.5 * (fov2.JD_bright+fov2.period - fov2.JD_faint)  # mid-brightening
    mags_jd = fov2.estimate_lpv_mags(jd)
    assert mags_jd['V'] == mags_faint['V'] + 0.5*(mags_bright['V']-mags_faint['V'])
    assert mags_jd['R'] == mags_faint['R'] + 0.5*(mags_bright['R']-mags_faint['R'])
    assert mags_jd['I'] == mags_faint['I'] + 0.5*(mags_bright['I']-mags_faint['I'])

    # Test much later dates than JD_bright etc, but same phase.
    jd = fov2.JD_bright + 11*fov2.period
    mags_bright = fov2.estimate_lpv_mags(jd)
    assert mags_bright['V'] == fov2.mag_V_bright
    assert mags_bright['R'] == fov2.mag_V_bright - 0.5*fov2.color_VI_bright
    assert mags_bright['I'] == fov2.mag_V_bright - fov2.color_VI_bright
    jd = fov2.JD_faint + 11*fov2.period
    mags_faint = fov2.estimate_lpv_mags(jd)
    assert mags_faint['V'] == fov2.mag_V_faint
    assert mags_faint['R'] == fov2.mag_V_faint - 0.5*fov2.color_VI_faint
    assert mags_faint['I'] == fov2.mag_V_faint - fov2.color_VI_faint
    jd = fov2.JD_bright + 0.5 * (fov2.JD_faint - fov2.JD_bright) + 11*fov2.period  # mid-dimming
    mags_jd = fov2.estimate_lpv_mags(jd)
    assert mags_jd['V'] == mags_bright['V'] + 0.5*(mags_faint['V']-mags_bright['V'])
    assert mags_jd['R'] == mags_bright['R'] + 0.5*(mags_faint['R']-mags_bright['R'])
    assert mags_jd['I'] == mags_bright['I'] + 0.5*(mags_faint['I']-mags_bright['I'])
    jd = fov2.JD_faint + 0.5 * (fov2.JD_bright+fov2.period - fov2.JD_faint) + \
        11*fov2.period  # mid-brightening
    mags_jd = fov2.estimate_lpv_mags(jd)
    assert mags_jd['V'] == mags_faint['V'] + 0.5*(mags_bright['V']-mags_faint['V'])
    assert mags_jd['R'] == mags_faint['R'] + 0.5*(mags_bright['R']-mags_faint['R'])
    assert mags_jd['I'] == mags_faint['I'] + 0.5*(mags_bright['I']-mags_faint['I'])

    # Test much earlier dates than JD_bright etc, but same phase.
    jd = fov2.JD_bright - 23*fov2.period
    mags_bright = fov2.estimate_lpv_mags(jd)
    assert mags_bright['V'] == fov2.mag_V_bright
    assert mags_bright['R'] == fov2.mag_V_bright - 0.5*fov2.color_VI_bright
    assert mags_bright['I'] == fov2.mag_V_bright - fov2.color_VI_bright
    jd = fov2.JD_faint - 23*fov2.period
    mags_faint = fov2.estimate_lpv_mags(jd)
    assert mags_faint['V'] == fov2.mag_V_faint
    assert mags_faint['R'] == fov2.mag_V_faint - 0.5*fov2.color_VI_faint
    assert mags_faint['I'] == fov2.mag_V_faint - fov2.color_VI_faint
    jd = fov2.JD_bright + 0.5 * (fov2.JD_faint - fov2.JD_bright) - 23*fov2.period  # mid-dimming
    mags_jd = fov2.estimate_lpv_mags(jd)
    assert mags_jd['V'] == mags_bright['V'] + 0.5*(mags_faint['V']-mags_bright['V'])
    assert mags_jd['R'] == mags_bright['R'] + 0.5*(mags_faint['R']-mags_bright['R'])
    assert mags_jd['I'] == mags_bright['I'] + 0.5*(mags_faint['I']-mags_bright['I'])
    jd = fov2.JD_faint + 0.5 * (fov2.JD_bright+fov2.period - fov2.JD_faint) - \
        23*fov2.period  # mid-brightening
    mags_jd = fov2.estimate_lpv_mags(jd)
    assert mags_jd['V'] == mags_faint['V'] + 0.5*(mags_bright['V']-mags_faint['V'])
    assert mags_jd['R'] == mags_faint['R'] + 0.5*(mags_bright['R']-mags_faint['R'])
    assert mags_jd['I'] == mags_faint['I'] + 0.5*(mags_bright['I']-mags_faint['I'])

    # Test more phases.
    from math import sin, pi
    fract = 0.15
    jd = fov2.JD_bright + fract * (fov2.JD_faint - fov2.JD_bright)
    mags_jd = fov2.estimate_lpv_mags(jd)
    linear_part = fract
    sine_part = (1 + sin((fract-0.5) * pi)) / 2
    mag_fract = (0.5 * linear_part) + (0.5 * sine_part)
    assert mags_jd['V'] == pytest.approx(mags_bright['V'] + mag_fract*(mags_faint['V']-mags_bright['V']))
    assert mags_jd['R'] == pytest.approx(mags_bright['R'] + mag_fract*(mags_faint['R']-mags_bright['R']))
    assert mags_jd['I'] == pytest.approx(mags_bright['I'] + mag_fract*(mags_faint['I']-mags_bright['I']))

    fract = 0.66
    jd = fov2.JD_faint + fract * (fov2.JD_bright - fov2.JD_faint)
    mags_jd = fov2.estimate_lpv_mags(jd)
    linear_part = fract
    sine_part = (1 + sin((fract-0.5) * pi)) / 2
    mag_fract = (0.5 * linear_part) + (0.5 * sine_part)
    assert mags_jd['V'] == pytest.approx(mags_faint['V'] + mag_fract*(mags_bright['V']-mags_faint['V']))
    assert mags_jd['R'] == pytest.approx(mags_faint['R'] + mag_fract*(mags_bright['R']-mags_faint['R']))
    assert mags_jd['I'] == pytest.approx(mags_faint['I'] + mag_fract*(mags_bright['I']-mags_faint['I']))


def test_fov_mira_init_multi_exposure():
    """ NB: For tests on #STARS section, see test_fov() above."""
    fov_name = "AU Aur Modified"
    f = fov.Fov("AU Aur Modified", TEST_FOV_DIRECTORY)
    assert f.fov_name == fov_name
    assert f.format_version == CURRENT_SCHEMA_VERSION
    assert f.observing_style == "LPV"
    assert len(f.observing_list) == 2
    assert f.alert is None
    assert f.observing_list == [("V", None, 5), ("I", None, 2)]


def test_fov_init_delta_scuti():
    """ NB: For tests on #STARS section, see test_fov() above."""
    fov_name = "ASAS J162540_19123"
    f = fov.Fov("ASAS J162540_19123", TEST_FOV_DIRECTORY)
    assert f.fov_name == fov_name
    assert f.format_version == CURRENT_SCHEMA_VERSION
    assert f.ra == util.ra_as_degrees("16:25:39.0")
    assert f.dec == util.dec_as_degrees("+19:03:18")
    assert f.chart == "X16073AJD"
    assert f.fov_date == "5/20/2016"
    assert f.main_target == "ASAS J162540+1912.3"
    assert f.target_type == "Delta Scuti"
    assert f.period == float("0.0571721")
    assert f.motive == 'Test only.'
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


def test_fov_init_z_cam():
    """ NB: For tests on #STARS section, see test_fov() above."""
    fov_name = "NSV 14581"
    f = fov.Fov("NSV 14581", TEST_FOV_DIRECTORY)
    assert f.fov_name == fov_name
    assert f.format_version == CURRENT_SCHEMA_VERSION
    assert f.ra == util.ra_as_degrees("23:26:50.3")
    assert f.dec == util.dec_as_degrees("+82:22:11")
    assert f.chart == "X16224AYI"
    assert f.fov_date == "9/2/2016"
    assert f.main_target == "NSV 14581"
    assert f.target_type == "Z Cam"
    assert f.period == float("0.194334535")
    assert f.motive == 'Test only.'
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


def test_fov_init_exoplanet():
    """ NB: For tests on #STARS section, see test_fov() above."""
    fov_name = "HAT_P_3"
    f = fov.Fov("HAT_P_3", TEST_FOV_DIRECTORY)
    assert f.fov_name == fov_name
    assert f.format_version == CURRENT_SCHEMA_VERSION
    assert f.main_target == "HAT-P-3"
    assert f.target_type == "Exoplanet"
    assert f.motive == 'Test only.'
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


def test_fov_init_standard():
    """ NB: For tests on #STARS section, see test_fov() above."""
    fov_name = "Std_SA100"
    f = fov.Fov("Std_SA100", TEST_FOV_DIRECTORY)
    assert f.fov_name == fov_name
    assert f.format_version == CURRENT_SCHEMA_VERSION
    assert f.ra == util.ra_as_degrees("08:53:14.3")
    assert f.dec == util.dec_as_degrees("-00:37:56")
    assert f.chart == "X15687X"
    assert f.fov_date == "12/20/2015"
    assert f.period is None
    assert f.motive == 'Test only, Std_SA100.'
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


def test_all_fov_names():
    all_names = fov.all_fov_names(fov_directory=TEST_FOV_DIRECTORY)
    assert isinstance(all_names, list)
    assert isinstance(all_names[0], str)
    assert len(all_names) == 8
    assert "ST Tri" in all_names


def test_make_fov_dict(fov_directory=TEST_FOV_DIRECTORY):
    # Full dictionary.
    fov_dict = fov.make_fov_dict(fov_directory)
    assert isinstance(fov_dict, dict)
    assert all([isinstance(f, fov.Fov) for f in list(fov_dict.values())])
    assert len(fov_dict) == 8
    assert "ST Tri" in fov_dict.keys()
    # Partial dictionary.
    fov_dict = fov.make_fov_dict(fov_directory,
                                 fov_names_selected=["ST Tri", "AU Aur", "ST Tri"])
    assert isinstance(fov_dict, dict)
    assert all([isinstance(f, fov.Fov) for f in list(fov_dict.values())])
    assert len(fov_dict) == 2
    assert "ST Tri" in fov_dict.keys()
