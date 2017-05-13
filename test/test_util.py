from datetime import datetime, timezone, timedelta
from math import isnan
import pytest

from photrix import util                  # call: util.ra_as_degrees()

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


def test_degrees_as_hex():
    assert util.degrees_as_hex(0.0) == "+00:00:00.00"
    assert util.degrees_as_hex(0.0) == util.degrees_as_hex(0.0, 2)  # default
    assert util.degrees_as_hex(+90.0) == "+90:00:00.00"
    assert util.degrees_as_hex(-90.0) == "-90:00:00.00"
    assert util.degrees_as_hex(0.001) == "+00:00:03.60"
    assert util.degrees_as_hex(-0.001) == "-00:00:03.60"
    assert util.degrees_as_hex(-69.125) == "-69:07:30.00"
    assert util.degrees_as_hex(69.125) == "+69:07:30.00"
    assert util.degrees_as_hex(90.001) == "+90:00:03.60"
    assert util.degrees_as_hex(-90.001, 4) == "-90:00:03.6000"
    assert util.degrees_as_hex(255) == "+255:00:00.00"
    assert util.degrees_as_hex(-255, 6) == "-255:00:00.000000"


def test_weighted_mean():
    assert util.weighted_mean([3], [7]) == 3
    assert util.weighted_mean([1, 3, 8], [0, 3, 9]) == 81/12
    with pytest.raises(ValueError) as e:
        util.weighted_mean([], [])                 # zero-length
    assert 'lengths of values & weights must be equal & non-zero' in str(e)
    # print('>>>' + str(e) + '<<<')
    with pytest.raises(ValueError) as e:
        util.weighted_mean([2, 3], [4, 5, 3])      # unequal lengths
    assert 'lengths of values & weights must be equal & non-zero' in str(e)
    # print('>>>' + str(e) + '<<<')
    with pytest.raises(ValueError) as e:
        util.weighted_mean([2, 3, 4], [1, 4, -5])  # sum(weights)==0
    assert 'sum of weights must be positive' in str(e)
    # print('>>>' + str(e) + '<<<')
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
    # Set up tests.
    dt1 = datetime(2016, 9, 10, 0, 0, 0, tzinfo=timezone.utc)
    dt2 = dt1 + timedelta(hours=1.5)
    ts1 = util.Timespan(dt1, dt2)
    # Test fields:
    assert ts1.start == dt1
    assert ts1.end == dt2
    assert ts1.seconds == 1.5 * 3600
    assert ts1.midpoint == dt1 + (dt2-dt1) / 2
    # Test copy():
    tscopy = ts1.copy()
    assert tscopy.start == ts1.start
    assert tscopy.end == ts1.end
    # Test equality:
    ts1eq = util.Timespan(dt1,dt2)
    assert ts1eq == ts1
    ts1neq1 = util.Timespan(dt1, dt2+timedelta(hours=1))
    assert ts1neq1 != ts1
    ts1neq2 = util.Timespan(dt1+timedelta(hours=1), dt2)
    assert ts1neq2 != ts1
    ts1neq3 = util.Timespan(dt1+timedelta(hours=1), dt2+timedelta(hours=1))
    assert ts1neq3 != ts1
    # Test .delay_seconds():
    ts1delay = ts1.delay_seconds(120)
    assert ts1delay.start == ts1.start + timedelta(seconds=120)
    assert ts1delay.end == ts1.end + timedelta(seconds=120)
    # Test .intersect():
    ts1a = util.Timespan(dt1 + timedelta(hours=0.5), dt2 + timedelta(hours=3))
    ts1a_int1 = ts1.intersect(ts1a)
    assert ts1a_int1.start == ts1a.start
    assert ts1a_int1.end == ts1.end
    ts1a_int2 = ts1a.intersect(ts1)
    assert ts1a_int1.start == ts1a_int2.start
    assert ts1a_int1.end == ts1a_int2.end
    # Test .subtract():
    ts1_sub = ts1.subtract(ts1.delay_seconds(+10000))  # no overlap
    assert ts1_sub == ts1
    ts1_sub = ts1.subtract(ts1.delay_seconds(-10000))  # also no overlap
    assert ts1_sub == ts1
    ts1_sub = ts1.subtract(ts1.delay_seconds(+1800))  # partial overlap of 1 hour
    assert ts1_sub.start == ts1.start
    assert ts1_sub.end == ts1.start + timedelta(seconds=+1800)
    ts1_sub = ts1.subtract(ts1.delay_seconds(-1800))  # partial overlap of 1 hour
    assert ts1_sub.start == ts1.delay_seconds(-1800).end
    assert ts1_sub.end == ts1.end
    ts1_other = util.Timespan(ts1.start+timedelta(seconds=100),
                              ts1.end+timedelta(seconds=-150))  # ts1_other within ts1 (early side).
    ts1_sub = ts1.subtract(ts1_other)
    assert ts1_sub.start == ts1_other.end
    assert ts1_sub.end == ts1.end
    ts1_other = util.Timespan(ts1.start+timedelta(seconds=160),
                              ts1.end+timedelta(seconds=-70))  # ts1_other within ts1 (late side).
    ts1_sub = ts1.subtract(ts1_other)
    assert ts1_sub.start == ts1.start
    assert ts1_sub.end == ts1_other.start
    ts1_other = util.Timespan(ts1.start+timedelta(seconds=-160),
                              ts1.end+timedelta(seconds=+7000))  # ts1 contained in ts1_other.
    ts1_sub = ts1.subtract(ts1_other)
    assert ts1_sub.start == ts1.start
    assert ts1_sub.end == ts1.start
    assert ts1_sub.seconds == 0
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
    # Test str():
    s = str(ts1)
    assert s.startswith("Timespan ") and s.endswith(" = 5400 seconds.")

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

    # Test .longer():
    dt1 = datetime(2016, 9, 10, 0, 0, 0, tzinfo=timezone.utc)
    dt2 = dt1 + timedelta(hours=1.5)
    dt3 = dt2 + timedelta(hours=1)
    tsa = util.Timespan(dt1, dt2)
    tsb = util.Timespan(dt1, dt3)
    # simple case: equal start, unequal end: returns longer Timespan.
    assert util.Timespan.longer(tsa, tsb) == tsb
    assert util.Timespan.longer(tsb, tsa) == tsb
    # edge case: identical inputs: returns first input.
    ts_copy = tsa.copy()
    assert util.Timespan.longer(tsa, ts_copy) == tsa
    assert util.Timespan.longer(ts_copy, tsa) == ts_copy
    # edge case: one zero-length input: returns the non-zero-length Timespan.
    ts_zero = util.Timespan(dt1, dt1)
    assert util.Timespan.longer(tsa, ts_zero) == tsa
    # edge case: both zero-length inputs: returns earlier Timespan.
    ts_zero_2 = util.Timespan(dt2, dt2)
    assert util.Timespan.longer(ts_zero, ts_zero_2) == ts_zero
    assert util.Timespan.longer(ts_zero_2, ts_zero) == ts_zero
    # case if_tie=="earlier": returns earlier Timespan.
    tsc = tsa.delay_seconds(7200)
    assert util.Timespan.longer(tsa, tsc, on_tie="earlier") == tsa
    assert util.Timespan.longer(tsc, tsa, on_tie="earlier") == tsa
    # case if_tie=="first": returns first input.
    assert util.Timespan.longer(tsa, tsc, on_tie="first") == tsa
    assert util.Timespan.longer(tsc, tsa, on_tie="first") == tsc
    # case if_tie==some other string: returns first input.
    assert util.Timespan.longer(tsa, tsc, on_tie="whatever") == tsa
    assert util.Timespan.longer(tsc, tsa, on_tie="whatever") == tsc


def test_RaDec():
    # Set up tests.
    rd1 = util.RaDec(60, 70)
    rd2 = util.RaDec(70, +80)
    rd3 = util.RaDec("10:30:00", "-45:15:00")
    rd4 = util.RaDec("10:30:00", -45.25)
    rd5 = util.RaDec(157.5, "-45:15:00")
    # Test fields and __eq__().
    assert rd1.as_degrees == (60, +70)
    assert rd1.as_degrees != (60, -70)
    assert rd1.as_hex == ("04:00:00.000", "+70:00:00.00")
    assert rd3.as_hex == ("10:30:00.000", "-45:15:00.00")
    assert rd3.ra == 10.5*15
    assert rd3.dec == -45.25
    assert rd1 != rd2
    assert rd3 == rd4 == rd5
    # Test .degrees_from() and .farther_from().
    assert rd1.degrees_from(rd2) == pytest.approx(10.293451406994343)
    assert util.RaDec(0, 60).degrees_from(util.RaDec(180, 60)) == pytest.approx(60.0)
    assert util.RaDec(0, 0).degrees_from(util.RaDec(180, 0)) == pytest.approx(180.0)
    assert rd1.farther_from(rd2, 10)
    assert not rd1.farther_from(rd2, 11)
    assert rd1.farther_from(rd3, 130)
    assert not rd1.farther_from(rd3, 140)
    # Test __str__() and __repr__().
    assert str(rd1) == "RaDec object:  04:00:00.000  +70:00:00.00"
    assert str(rd5) == "RaDec object:  10:30:00.000  -45:15:00.00"
    assert repr(rd2) == "RaDec('04:40:00.000', '+80:00:00.00')"
    assert repr(rd5) == "RaDec('10:30:00.000', '-45:15:00.00')"
    assert eval('util.' + repr(rd1)) == rd1
    assert eval('util.' + repr(rd5)) == rd5
    assert eval('util.' + repr(rd1)) != rd5


def test_get_phase():
    assert util.get_phase(5, 2, 10) == pytest.approx(0.3)
    assert util.get_phase(-5, 2, 10) == pytest.approx(0.3)
    assert util.get_phase(-335, 2, 10) == pytest.approx(0.3)
    assert util.get_phase(335, 2, 10) == pytest.approx(0.3)
    assert util.get_phase(3352, 2, 10) == 0
    assert util.get_phase(2, 2, 10) == 0
    assert util.get_phase(1.99, 2, 10) == pytest.approx(0.999)
    assert util.get_phase(12, 2, 10) == 0
    assert util.get_phase(7, 2, 10) == 0.5


def test_jd_from_datetime_utc():
    one_second = 1.0 / (24.0 * 3600.0)  # tolerance of 1 second, in days (for JD)
    datetime_j2000 = datetime(2000, 1, 1, 0, 0, 0).replace(tzinfo=timezone.utc)
    assert util.jd_from_datetime_utc(datetime_j2000) == pytest.approx(2451544.5, abs=one_second)

    assert util.jd_from_datetime_utc(None) is None
    assert util.jd_from_datetime_utc() is None

    datetime_1 = datetime(2017, 1, 9, 15, 23, 53).replace(tzinfo=timezone.utc)
    assert util.jd_from_datetime_utc(datetime_1) == pytest.approx(2457763.14158398, abs=one_second)
    datetime_2 = datetime(2020, 7, 9, 6, 23, 53).replace(tzinfo=timezone.utc)
    assert util.jd_from_datetime_utc(datetime_2) == pytest.approx(2459039.76658403, abs=one_second)
    datetime_3 = datetime(1986, 10, 11, 3, 12, 7).replace(tzinfo=timezone.utc)
    assert util.jd_from_datetime_utc(datetime_3) == pytest.approx(2446714.63341273, abs=one_second)


def test_datetime_utc_from_jd():
    jd_j2000 = 2451544.5
    datetime_j2000 = datetime(2000, 1, 1, 0, 0, 0).replace(tzinfo=timezone.utc)
    assert (util.datetime_utc_from_jd(jd_j2000) - datetime_j2000).total_seconds() == \
           pytest.approx(0, abs=1.0)

    datetime_now = datetime.now(timezone.utc)
    jd_now = util.jd_from_datetime_utc(datetime_now)  # tested just above.
    assert (util.datetime_utc_from_jd(jd_now) - datetime_now).total_seconds() == \
           pytest.approx(0, abs=1.0)

    datetime_1 = datetime(2017, 1, 9, 15, 23, 53).replace(tzinfo=timezone.utc)
    jd_1 = 2457763.14158398
    assert (util.datetime_utc_from_jd(jd_1) - datetime_1).total_seconds() == \
           pytest.approx(0, abs=1.0)
    datetime_2 = datetime(2020, 7, 9, 6, 23, 53).replace(tzinfo=timezone.utc)
    jd_2 = 2459039.76658403
    assert (util.datetime_utc_from_jd(jd_2) - datetime_2).total_seconds() == \
           pytest.approx(0, abs=1.0)
    datetime_3 = datetime(1986, 10, 11, 3, 12, 7).replace(tzinfo=timezone.utc)
    jd_3 = 2446714.63341273
    assert (util.datetime_utc_from_jd(jd_3) - datetime_3).total_seconds() == \
           pytest.approx(0, abs=1.0)


def test_time_hhmm():
    dt = datetime(2016, 1, 1, 23, 34, 45, 454545).replace(tzinfo=timezone.utc)
    assert util.hhmm_from_datetime_utc(dt) == '2335'
    dt = datetime(2016, 1, 1, 23, 34, 29, 999999).replace(tzinfo=timezone.utc)
    assert util.hhmm_from_datetime_utc(dt) == '2334'
    dt = datetime(2016, 1, 1, 23, 59, 31, 454545).replace(tzinfo=timezone.utc)
    assert util.hhmm_from_datetime_utc(dt) == '0000'
    dt = datetime(2016, 1, 31, 0, 0, 0, 0).replace(tzinfo=timezone.utc)
    assert util.hhmm_from_datetime_utc(dt) == '0000'
    dt = datetime(2016, 1, 31, 0, 0, 30, 0).replace(tzinfo=timezone.utc)  # banker's rounding.
    assert util.hhmm_from_datetime_utc(dt) == '0000'
    dt = datetime(2016, 1, 31, 0, 1, 30, 0).replace(tzinfo=timezone.utc)  # banker's rounding.
    assert util.hhmm_from_datetime_utc(dt) == '0002'
    dt = datetime(2016, 1, 31, 0, 0, 30, 1).replace(tzinfo=timezone.utc)
    assert util.hhmm_from_datetime_utc(dt) == '0001'

    # Test .az_alt_at_datetime_utc():
    aldebaran = util.RaDec('4:35:55.31', '+16:30:30.249')
    mira = util.RaDec('02:19:20.804', '-2:58:43.606')
    hip_113116 = util.RaDec('22:54:25.073', '84:20:46.620')

    long, lat = '-95:53:18', '+38:55:29'  # Bois d'Arc Obs, Kansas
    dt = datetime(2016, 1, 31, 0, 9, 30, 780000).replace(tzinfo=timezone.utc)
    az, alt = util.az_alt_at_datetime_utc(long, lat, aldebaran, dt)  # in degrees
    assert az == pytest.approx(util.hex_degrees_as_degrees('118:27:06'), abs=0.5)
    assert alt== pytest.approx(util.hex_degrees_as_degrees('+53:30:20'), abs=0.5)
    az, alt = util.az_alt_at_datetime_utc(long, lat, mira, dt)  # in degrees
    assert az == pytest.approx(util.hex_degrees_as_degrees('181:40:22'), abs=0.5)
    assert alt == pytest.approx(util.hex_degrees_as_degrees('+48:09:18'), abs=0.5)
    az, alt = util.az_alt_at_datetime_utc(long, lat, hip_113116, dt)  # in degrees
    assert az == pytest.approx(util.hex_degrees_as_degrees('354:02:03'), abs=0.5)
    assert alt == pytest.approx(util.hex_degrees_as_degrees('+42:09:15'), abs=0.5)

    dt = datetime(2016, 6, 13, 3, 17, 53, 780000).replace(tzinfo=timezone.utc)
    az, alt = util.az_alt_at_datetime_utc(long, lat, aldebaran, dt)  # in degrees
    assert az == pytest.approx(util.hex_degrees_as_degrees('323:35:35'), abs=0.5)
    assert alt== pytest.approx(util.hex_degrees_as_degrees('-26:10:56'), abs=0.5)
    az, alt = util.az_alt_at_datetime_utc(long, lat, mira, dt)  # in degrees
    assert az == pytest.approx(util.hex_degrees_as_degrees('0:42:33'), abs=0.5)
    assert alt == pytest.approx(util.hex_degrees_as_degrees('-53:58:40'), abs=0.5)
    az, alt = util.az_alt_at_datetime_utc(long, lat, hip_113116, dt)  # in degrees
    assert az == pytest.approx(util.hex_degrees_as_degrees('5:22:23'), abs=0.5)
    assert alt == pytest.approx(util.hex_degrees_as_degrees('+35:21:32'), abs=0.5)


def test_isfloat():
    assert util.isfloat('-25') is True
    assert util.isfloat('+34.4') is True
    assert util.isfloat('  34.4  ') is True
    assert util.isfloat('-.6') is True
    assert util.isfloat('-123.E6') is True
    assert util.isfloat('-inf') is True
    assert util.isfloat('NaN') is True
    assert util.isfloat(True) is True
    assert util.isfloat(False) is True
    assert util.isfloat('') is False
    assert util.isfloat('12.2.2') is False
    assert util.isfloat('-') is False
    assert util.isfloat('+') is False
    assert util.isfloat('(3)') is False


def test_float_or_none():
    assert util.float_or_none('-25') == -25
    assert util.float_or_none('+34.4') == pytest.approx(34.4)
    assert util.float_or_none('  34.4  ') == pytest.approx(34.4)
    assert util.float_or_none('-.6') == pytest.approx(-0.6)
    assert util.float_or_none('-123.E6') == pytest.approx(-123000000.0)
    assert util.float_or_none('-inf') == float('-inf')
    assert isnan(util.float_or_none('NaN'))
    assert util.float_or_none(True) == 1.0
    assert util.float_or_none(False) == 0.0
    assert util.float_or_none('') is None
    assert util.float_or_none('12.2.2') is None
    assert util.float_or_none('-') is None
    assert util.float_or_none('+') is None
    assert util.float_or_none('(3)') is None


def test_event_utcs_in_timespan():
    # Test normal case:
    ts1 = util.Timespan(datetime(2017, 2, 10, 1, 30, 0).replace(tzinfo=timezone.utc),
                        datetime(2017, 2, 10, 10, 30, 0).replace(tzinfo=timezone.utc))
    min1_utc_list = util.event_utcs_in_timespan(2455336.44, 0.3641393156, ts1)  # DE CVn 1'min
    assert len(min1_utc_list) == 1
    assert util.jd_from_datetime_utc(min1_utc_list[0]) == pytest.approx(2457794.74452, abs=0.0001)
    min2_utc_list = util.event_utcs_in_timespan(2455336.26, 0.3641393156, ts1)  # DE CVn 2'min
    assert len(min2_utc_list) == 2
    assert util.jd_from_datetime_utc(min2_utc_list[0]) == pytest.approx(2457794.56452, abs=0.0001)
    assert util.jd_from_datetime_utc(min2_utc_list[1]) == pytest.approx(2457794.92866, abs=0.0001)

    # Case: no mins in timespan:
    ts1 = util.Timespan(datetime(2017, 2, 10, 7, 0, 0).replace(tzinfo=timezone.utc),
                        datetime(2017, 2, 10, 8, 30, 0).replace(tzinfo=timezone.utc))
    min1_utc_list = util.event_utcs_in_timespan(2455336.44, 0.3641393156, ts1)  # DE CVn 1'min
    assert len(min1_utc_list) == 0
    min2_utc_list = util.event_utcs_in_timespan(2455336.26, 0.3641393156, ts1)  # DE CVn 2'min
    assert len(min2_utc_list) == 0

    # Case: missing data:
    assert util.event_utcs_in_timespan(None, 0.3641393156, ts1) is None
    assert util.event_utcs_in_timespan(2455336.44, None, ts1) is None


def test_mixed_model_fit_class():
    import numpy as np
    import pandas as pd
    # First, construct test data frame:
    points = 80
    np.random.seed(1234)
    d = {'A': np.random.randn(points),
         'B': np.random.randn(points),
         'C': np.random.randn(points),
         'Ran': np.random.randint(0, 3, points),
         'Dep': 0}
    df = pd.DataFrame(d)
    df['Dep'] = 17 + 1*df.A + 2*df.B + 4*df.C + 5*(df.Ran-1) + 1*np.random.randn(len(df))
    categories = ['X', 'Y', 'Z']
    df['Ran'] = [categories[r] for r in df['Ran']]
    df.index = df.index + 200
    # Split test data into model and test blocks:
    df_model = df[0:int(3*points/4)]
    df_test = df[len(df_model):]

    # Construct fit object:
    fit = util.MixedModelFit(df_model, dep_var='Dep', fixed_vars=['A', 'B', 'C'], group_var='Ran')

    # Test object and scalar attributes:
    assert isinstance(fit, util.MixedModelFit)
    assert fit.converged is True
    assert fit.nobs == len(df_model)
    assert fit.likelihood == pytest.approx(-95.7673)
    assert fit.dep_var == 'Dep'
    assert fit.fixed_vars == ['A', 'B', 'C']
    assert fit.group_var == 'Ran'
    assert fit.sigma == pytest.approx(1.030723)

    # Test fixed-effect results dataframe:
    assert list(fit.df_fixed_effects.index) == ['Intercept', 'A', 'B', 'C']
    assert list(fit.df_fixed_effects['Name']) == list(fit.df_fixed_effects.index)
    assert list(fit.df_fixed_effects['Value']) == pytest.approx([16.648186, 0.946692,
                                                                 1.959923, 4.069383], abs=0.00001)
    assert list(fit.df_fixed_effects['Stdev']) == pytest.approx([2.844632, 0.142185,
                                                                 0.134386, 0.145358], abs=0.00001)
    assert list(fit.df_fixed_effects['Tvalue']) == pytest.approx([5.8525, 6.65818,
                                                                  14.58429, 27.99568], abs=0.0001)
    assert list(fit.df_fixed_effects['Pvalue'] * 10**9) == pytest.approx([4.8426, 0.027724,
                                                                          0, 0], abs=0.0001)

    # Test random-effect (group) results dataframe:
    assert list(fit.df_random_effects.index) == ['X', 'Y', 'Z']
    assert list(fit.df_random_effects['GroupName']) == list(fit.df_random_effects.index)
    assert list(fit.df_random_effects['GroupValue']) == pytest.approx([-5.164649, 0.543793,
                                                                       4.620857], abs=0.00001)

    # Test observation results dataframe:
    assert list(fit.df_observations['FittedValue'])[0:4] == pytest.approx([24.809899, 10.130408,
                                                                           19.834543, 7.758331],
                                                                          abs=0.00001)
    assert list(fit.df_observations['Residual'])[0:4] == pytest.approx([0.490179, -0.786949,
                                                                        0.58315, -1.23926],
                                                                       abs=0.00001)

    # Verify predictions on model data:
    assert list(fit.predict(df_model)) == pytest.approx(fit.df_observations['FittedValue'])
    # Verify predictions on test data:
    assert list(fit.predict(df_test[0:4])) == pytest.approx([16.706808844306536, 18.50420441664172,
                                                            20.533844119080726, 18.952844883150998])






