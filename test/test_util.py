import pytest

##### Your choices for importing (at least in this test module):
from photrix import util                  # call: util.ra_as_degrees()
# from photrix.util import *              # call: ra_as_degrees()
# from photrix.util import ra_as_degrees  # call: ra_as_degrees()
# import photrix.util as u                # call: u.ra_as_degrees()
# from .context import photrix            # XXX I can't make this work with a context.py file.

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


def test_ra_as_degrees():
    assert util.ra_as_degrees("180") == 180.0
    assert util.ra_as_degrees("0") == 0.0
    assert util.ra_as_degrees("360") == 360.0
    assert util.ra_as_degrees("-0.1") is None
    assert util.ra_as_degrees("360.1") is None

    assert util.ra_as_degrees("12:00") == 180.0
    assert util.ra_as_degrees("12:00:00") == 180.0
    assert util.ra_as_degrees("0:00:00") == 0.0
    assert util.ra_as_degrees("11:16:30") == 169.125
    assert util.ra_as_degrees("24:00:01") is None


def test_hex_degrees_as_degrees():
    assert util.hex_degrees_as_degrees("12") == 12.0
    assert util.hex_degrees_as_degrees("-12") == -12.0
    assert util.hex_degrees_as_degrees("0") == 0.0
    assert util.hex_degrees_as_degrees("-0") == 0.0
    assert util.hex_degrees_as_degrees("90") == 90
    assert util.hex_degrees_as_degrees("-90") == -90
    assert util.hex_degrees_as_degrees("90.125") == 90.125
    assert util.hex_degrees_as_degrees("-90.125") == -90.125

    assert util.hex_degrees_as_degrees("88:45") == 88.75
    assert util.hex_degrees_as_degrees("-88:45") == -88.75
    assert util.hex_degrees_as_degrees("12:34:30") == 12.575
    assert util.hex_degrees_as_degrees("-12:34:30") == -12.575
    assert util.hex_degrees_as_degrees("91:34:30") == 91.575
    assert util.hex_degrees_as_degrees("-91:34:30") == -91.575
    assert util.hex_degrees_as_degrees("91:45") == 91.75
    assert util.hex_degrees_as_degrees("-91:45") == -91.75


def test_dec_as_degrees():
    assert util.dec_as_degrees("12") == 12.0
    assert util.dec_as_degrees("-12") == -12.0
    assert util.dec_as_degrees("0") == 0.0
    assert util.dec_as_degrees("-0") == 0.0
    assert util.dec_as_degrees("90") == 90
    assert util.dec_as_degrees("-90") == -90
    assert util.dec_as_degrees("90.1") is None
    assert util.dec_as_degrees("-90.1") is None

    assert util.dec_as_degrees("88:45") == 88.75
    assert util.dec_as_degrees("-88:45") == -88.75
    assert util.dec_as_degrees("12:34:30") == 12.575
    assert util.dec_as_degrees("-12:34:30") == -12.575
    assert util.dec_as_degrees("91:34:30") is None
    assert util.dec_as_degrees("-91:34:30") is None
    assert util.dec_as_degrees("91:45") is None
    assert util.dec_as_degrees("-91:45") is None


def test_ra_as_hours():
    assert util.ra_as_hours(0.0) == "00:00:00.000"
    assert util.ra_as_hours(20.0) == "01:20:00.000"
    assert util.ra_as_hours(169.125) == "11:16:30.000"
    assert util.ra_as_hours(180.0) == "12:00:00.000"
    assert util.ra_as_hours(360.0) == "00:00:00.000"
    assert util.ra_as_hours(359.99) == "23:59:57.600"
    assert util.ra_as_hours(359.999) == "23:59:59.760"
    assert util.ra_as_hours(359.9999) == "23:59:59.976"
    assert util.ra_as_hours(359.99999) == "23:59:59.998"
    assert util.ra_as_hours(359.999999) == "00:00:00.000"
    assert util.ra_as_hours(359.9999999) == "00:00:00.000"
    assert util.ra_as_hours(359.99999999) == "00:00:00.000"
    assert util.ra_as_hours(-0.01) is None
    assert util.ra_as_hours(360.01) is None
    assert util.ra_as_hours(-44) is None
    assert util.ra_as_hours(654) is None


def test_dec_as_hex():
    assert util.dec_as_hex(0.0) == "+00:00:00.00"
    assert util.dec_as_hex(+90.0) == "+90:00:00.00"
    assert util.dec_as_hex(-90.0) == "-90:00:00.00"
    assert util.dec_as_hex(0.001) == "+00:00:03.60"
    assert util.dec_as_hex(-0.001) == "-00:00:03.60"
    assert util.dec_as_hex(-69.125) == "-69:07:30.00"
    assert util.dec_as_hex(69.125) == "+69:07:30.00"
    assert util.dec_as_hex(90.001) is None
    assert util.dec_as_hex(-90.001) is None
    assert util.dec_as_hex(255) is None
    assert util.dec_as_hex(-255) is None


def test_weighted_mean():
    assert util.weighted_mean([3], [7]) == 3
    assert util.weighted_mean([1, 3, 8], [0, 3, 9]) == 81/12
    with pytest.raises(ValueError) as e:
        util.weighted_mean([], [])                 # zero-length
    assert 'lengths of values & weights must be equal & non-zero' in str(e)
    print('>>>' + str(e) + '<<<')
    with pytest.raises(ValueError) as e:
        util.weighted_mean([2, 3], [4, 5, 3])      # unequal lengths
    assert 'lengths of values & weights must be equal & non-zero' in str(e)
    print('>>>' + str(e) + '<<<')
    with pytest.raises(ValueError) as e:
        util.weighted_mean([2, 3, 4], [1, 4, -5])  # sum(weights)==0
    assert 'sum of weights must be positive' in str(e)
    print('>>>' + str(e) + '<<<')
    assert util.weighted_mean([3], [7], True) == (3, 0)
    assert util.weighted_mean([1, 3, 8], [0, 3, 9], True) == (81/12, 2.9296875)
    assert util.weighted_mean([1, 3, 8], [0, 3, 9], True) == util.weighted_mean([3, 8], [3, 9], True)

