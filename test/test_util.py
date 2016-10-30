from datetime import datetime, timezone, timedelta
import pytest

# Your choices for importing (at least in this test module):
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
    assert util.weighted_mean([1, 3, 8], [0, 3, 9], True) == \
           util.weighted_mean([3, 8], [3, 9], True)


def test_ladder_round():
    assert util.ladder_round(0) == 0
    assert util.ladder_round(0.12) == 0.125
    assert util.ladder_round(45) == 50
    assert util.ladder_round(10) == 10
    assert util.ladder_round(100) == 100
    assert util.ladder_round(64) == 64
    assert util.ladder_round(99.7) == 100
    assert util.ladder_round(-99.7) == -100
    assert util.ladder_round(-45) == -50
    assert util.ladder_round(-64) == -64


def test_Timespan():
    # Normal case.
    dt1 = datetime(2016, 9, 10, 0, 0, 0, tzinfo=timezone.utc)
    dt2 = dt1 + timedelta(hours=1.5)
    ts1 = util.Timespan(dt1, dt2)
    assert ts1.start == dt1
    assert ts1.end == dt2
    assert ts1.seconds == 1.5 * 3600
    assert ts1.midpoint == dt1 + (dt2-dt1) / 2
    # Test .contains_time():
    assert ts1.contains_time(dt1)
    assert ts1.contains_time(dt2)
    assert ts1.contains_time(dt1 + timedelta(hours=0.5))
    assert not ts1.contains_time(dt1 - timedelta(hours=0.5))
    assert not ts1.contains_time(dt2 + timedelta(hours=0.5))
    # Test .contains_timespan():
    ts_contained = util.Timespan(dt1+timedelta(hours=0.1), dt2+timedelta(hours=-0.1))
    assert ts1.contains_timespan(ts_contained)
    assert not ts_contained.contains_timespan(ts1)
    ts_identical = util.Timespan(dt1+timedelta(hours=0), dt2+timedelta(hours=0))
    assert ts1.contains_timespan(ts_identical)
    assert ts_identical.contains_timespan(ts1)
    ts_shift_later = util.Timespan(dt1+timedelta(hours=0.5), dt2+timedelta(hours=0.5))
    assert not ts1.contains_timespan(ts_shift_later)
    assert not ts_shift_later.contains_timespan(ts1)
    # Test .intersect():
    ts1a = util.Timespan(dt1+timedelta(hours=0.5), dt2+timedelta(hours=3))
    ts1a_int1 = ts1.intersect(ts1a)
    assert ts1a_int1.start == ts1a.start
    assert ts1a_int1.end == ts1.end
    ts1a_int2 = ts1a.intersect(ts1)
    assert ts1a_int1.start == ts1a_int2.start
    assert ts1a_int1.end == ts1a_int2.end

    # Case: input times equal
    dt3 = dt1
    ts2 = util.Timespan(dt1, dt3)
    assert ts2.start == dt1
    assert ts2.end == dt3
    assert ts2.seconds == 0
    assert ts2.midpoint == dt1 == dt3
    assert ts2.contains_time(dt1)
    assert ts2.contains_time(dt3)
    assert not ts2.contains_time(dt1 + timedelta(hours=0.5))
    assert not ts2.contains_time(dt1 - timedelta(hours=0.5))
    assert not ts2.contains_time(dt3 + timedelta(hours=0.5))

    # Case: input times inverted
    ts3 = util.Timespan(dt2, dt1)
    assert ts3.start == dt2
    assert ts3.end == dt2
    assert ts3.seconds == 0
    assert ts3.midpoint == dt2
    assert ts3.contains_time(dt2)  # but not necessarily dt1 which is ~ invalid
    assert not ts3.contains_time(dt2 + timedelta(hours=0.5))
    assert not ts3.contains_time(dt2 - timedelta(hours=0.5))
    assert not ts3.contains_time(dt1 + timedelta(hours=0.5))

# TODO write test for class RaDec
