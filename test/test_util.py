from util import *

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


def test_ra_as_degrees():
    assert ra_as_degrees("180") == 180.0
    assert ra_as_degrees("0") == 0.0
    assert ra_as_degrees("360") == 360.0
    assert ra_as_degrees("-0.1") is None
    assert ra_as_degrees("360.1") is None

    assert ra_as_degrees("12:00") == 180.0
    assert ra_as_degrees("12:00:00") == 180.0
    assert ra_as_degrees("0:00:00") == 0.0
    assert ra_as_degrees("11:16:30") == 169.125
    assert ra_as_degrees("24:00:01") is None


def test_hex_degrees_as_degrees():
    assert hex_degrees_as_degrees("12") == 12.0
    assert hex_degrees_as_degrees("-12") == -12.0
    assert hex_degrees_as_degrees("0") == 0.0
    assert hex_degrees_as_degrees("-0") == 0.0
    assert hex_degrees_as_degrees("90") == 90
    assert hex_degrees_as_degrees("-90") == -90
    assert hex_degrees_as_degrees("90.125") == 90.125
    assert hex_degrees_as_degrees("-90.125") == -90.125

    assert hex_degrees_as_degrees("88:45") == 88.75
    assert hex_degrees_as_degrees("-88:45") == -88.75
    assert hex_degrees_as_degrees("12:34:30") == 12.575
    assert hex_degrees_as_degrees("-12:34:30") == -12.575
    assert hex_degrees_as_degrees("91:34:30") == 91.575
    assert hex_degrees_as_degrees("-91:34:30") == -91.575
    assert hex_degrees_as_degrees("91:45") == 91.75
    assert hex_degrees_as_degrees("-91:45") == -91.75


def test_dec_as_degrees():
    assert dec_as_degrees("12") == 12.0
    assert dec_as_degrees("-12") == -12.0
    assert dec_as_degrees("0") == 0.0
    assert dec_as_degrees("-0") == 0.0
    assert dec_as_degrees("90") == 90
    assert dec_as_degrees("-90") == -90
    assert dec_as_degrees("90.1") is None
    assert dec_as_degrees("-90.1") is None

    assert dec_as_degrees("88:45") == 88.75
    assert dec_as_degrees("-88:45") == -88.75
    assert dec_as_degrees("12:34:30") == 12.575
    assert dec_as_degrees("-12:34:30") == -12.575
    assert dec_as_degrees("91:34:30") is None
    assert dec_as_degrees("-91:34:30") is None
    assert dec_as_degrees("91:45") is None
    assert dec_as_degrees("-91:45") is None


def test_ra_as_hours():
    assert ra_as_hours(0.0) == "00:00:00.000"
    assert ra_as_hours(20.0) == "01:20:00.000"
    assert ra_as_hours(169.125) == "11:16:30.000"
    assert ra_as_hours(180.0) == "12:00:00.000"
    assert ra_as_hours(360.0) == "00:00:00.000"
    assert ra_as_hours(359.99) == "23:59:57.600"
    assert ra_as_hours(359.999) == "23:59:59.760"
    assert ra_as_hours(359.9999) == "23:59:59.976"
    assert ra_as_hours(359.99999) == "23:59:59.998"
    assert ra_as_hours(359.999999) == "00:00:00.000"
    assert ra_as_hours(359.9999999) == "00:00:00.000"
    assert ra_as_hours(359.99999999) == "00:00:00.000"
    assert ra_as_hours(-0.01) is None
    assert ra_as_hours(360.01) is None
    assert ra_as_hours(-44) is None
    assert ra_as_hours(654) is None


def test_dec_as_hex():
    assert dec_as_hex(0.0) == "+00:00:00.00"
    assert dec_as_hex(+90.0) == "+90:00:00.00"
    assert dec_as_hex(-90.0) == "-90:00:00.00"
    assert dec_as_hex(0.001) == "+00:00:03.60"
    assert dec_as_hex(-0.001) == "-00:00:03.60"
    assert dec_as_hex(-69.125) == "-69:07:30.00"
    assert dec_as_hex(69.125) == "+69:07:30.00"
    assert dec_as_hex(90.001) is None
    assert dec_as_hex(-90.001) is None
    assert dec_as_hex(255) is None
    assert dec_as_hex(-255) is None

