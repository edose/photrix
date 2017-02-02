from datetime import datetime, timedelta, timezone
import math
import ephem

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


class Timespan:
    """ Holds one (start, end) span of time. Immutable.
        Input: 2 python datetimes (in UTC), defining start and end of timespan.
        methods:
        ts2 = ts.copy()
        ts2 == ts  # only if both start and end are equal
        ts2 = ts.delay_seconds(120)  # returns new Timespan, offset in both start and end
        ts.intersect(other)  # returns True iff any overlap at all
        ts2 = ts.subtract(other)  # returns new Timespan; longer of 2 possible spans if ambiguous.
        ts.contains_time(t)  # returns True iff ts.start <= t <= ts.end
        ts.contains_timespan(other)  # returns True iff ts wholly contains other
        str(ts)  # returns string describing Timespan's start, end, and duration in seconds.
    """
    def __init__(self, start_utc, end_utc):
        self.start = start_utc
        self.end = max(start_utc, end_utc)
        self.seconds = (self.end-self.start).seconds
        self.midpoint = self.start + timedelta(seconds=self.seconds / 2)

    def copy(self):
        return Timespan(self.start, self.end)

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end

    def delay_seconds(self, seconds):
        delay = timedelta(seconds=seconds)
        return Timespan(self.start+delay, self.end+delay)

    def intersect(self, other):
        new_start = max(self.start, other.start)
        new_end = min(self.end, other.end)
        return Timespan(new_start, new_end)

    def subtract(self, other):
        if self.intersect(other).seconds == 0:  # case: no overlap/intersection.
            return self
        if other.contains_timespan(self):  # case: self entirely subtracted away.
            return Timespan(self.start, self.start)
        if self.contains_timespan(other):  # case: 2 timespans -> take the longer.
            diff_early = Timespan(self.start, other.start)
            diff_late = Timespan(other.end, self.end)
            if diff_early.seconds >= diff_late.seconds:
                return diff_early
            else:
                return diff_late
        if self.start < other.start:  # remaining case: partial overlap.
            return Timespan(self.start, other.start)
        else:
            return Timespan(other.end, self.end)

    def contains_time(self, time_utc):
        return self.start <= time_utc <= self.end

    def contains_timespan(self, other):
        return (self.start <= other.start) & (self.end >= other.end)

    @staticmethod
    def longer(ts1, ts2, on_tie="earlier"):
        """
        Returns Timespan with longer duration (larger .seconds).
        If equal duration:
            if_tie=="earlier", return earlier.
            if_tie=="first", return ts1.
            [TODO: add "random" option later to return randomly chosen ts1 or ts2.]
        :param ts1: input Timespan object.
        :param ts2: input Timespan object.
        :param on_tie: "earlier" or "first". Any other string behaves as "first".
        :return: the Timespan object with longer duration.
        """
        if ts1.seconds > ts2.seconds:
            return ts1
        if ts2.seconds > ts1.seconds:
            return ts2
        # here: equal length cases. First, try to break duration tie with earlier midpoint.
        if on_tie.lower() == "earlier" and ts1.midpoint != ts2.midpoint:
            if ts1.midpoint < ts2.midpoint:
                return ts1
            return ts2
        # here, tie-breaking has failed. So simply return first of 2 input Timespans.
        return ts1

    def __str__(self):
        return "Timespan '" + str(self.start) + "' to '" + str(self.end) + "' = " + \
               str(self.seconds) + " seconds."


class RaDec:
    """
    Holds one Right Ascension, Declination sky position (internally as degrees).
    Parameters:
        ra : (hours hex string, or degrees)
        dec : (degrees hex string, or degrees)
    Methods:

    """
    def __init__(self, ra, dec):
        if isinstance(ra, str):
            self.ra = ra_as_degrees(ra)
        else:
            self.ra = ra
        if isinstance(dec, str):
            self.dec = dec_as_degrees(dec)
        else:
            self.dec = dec
        self.as_degrees = self.ra, self.dec  # stored internally as degrees
        self.as_hex = ra_as_hours(self.ra), dec_as_hex(self.dec)

    def degrees_from(self, other):
        deg_per_radian = 180.0 / math.pi
        diff_ra = abs(self.ra - other.ra) / deg_per_radian
        cos_dec_1 = math.cos(self.dec / deg_per_radian)
        cos_dec_2 = math.cos(other.dec / deg_per_radian)
        diff_dec = abs(self.dec - other.dec) / deg_per_radian
        arg = math.sqrt(math.sin(diff_dec/2.0)**2 + cos_dec_1*cos_dec_2*math.sin(diff_ra/2.0)**2)
        if arg > 0.001:
            return deg_per_radian * (2.0 * math.asin(arg))  # haversine formula
        else:
            # spherical law of cosines
            sin_dec_1 = math.sin(self.dec / deg_per_radian)
            sin_dec_2 = math.sin(other.dec / deg_per_radian)
            return deg_per_radian * \
                math.acos(sin_dec_1*sin_dec_2 + cos_dec_1*cos_dec_2*math.cos(diff_ra))

    def farther_from(self, other_ra_dec, degrees_limit):
        return self.degrees_from(other_ra_dec) > degrees_limit

    def __eq__(self, other):
        return (self.ra == other.ra) and (self.dec == other.dec)

    def __str__(self):
        ra_hex, dec_hex = self.as_hex
        return "RaDec object:  " + ra_hex + "  " + dec_hex

    def __repr__(self):
        ra_hex, dec_hex = self.as_hex
        return "RaDec('" + ra_hex + "', '" + dec_hex + "')"


def ra_as_degrees(ra_string):
    """ Input: string in either full hex ("12:34:56.7777") or degrees ("234.55")
        Returns: float of Right Ascension in degrees between 0 and 360.
    """
    ra_list = ra_string.split(":")
    if len(ra_list) == 1:
        ra_degrees = float(ra_list[0])  # input assumed to be in degrees.
    elif len(ra_list) == 2:
        ra_degrees = 15 * (float(ra_list[0]) + float(ra_list[1])/60.0)  # input assumed in hex.
    else:
        ra_degrees = 15 * (float(ra_list[0]) + float(ra_list[1]) / 60.0 +
                           float(ra_list[2])/3600.0)  # input assumed in hex.
    if (ra_degrees < 0) | (ra_degrees > 360):
        ra_degrees = None
    return ra_degrees


def hex_degrees_as_degrees(hex_degrees_string):
    """ Input: string in either full hex ("-12:34:56.7777") or degrees ("-24.55")
        Returns: float of degrees (not limited)
    """
    dec_list = hex_degrees_string.split(":")
    dec_list = [dec.strip() for dec in dec_list]
    if dec_list[0].startswith("-"):
        sign = -1
    else:
        sign = 1
    if len(dec_list) == 1:
        dec_degrees = float(dec_list[0])  # input assumed to be in degrees.
    elif len(dec_list) == 2:
        dec_degrees = sign * (abs(float(dec_list[0])) + float(dec_list[1])/60.0)  # input is hex.
    else:
        dec_degrees = sign * (abs(float(dec_list[0])) + float(dec_list[1]) / 60.0 +
                              float(dec_list[2])/3600.0)  # input is hex.
    return dec_degrees


def dec_as_degrees(dec_string):
    """ Input: string in either full hex ("-12:34:56.7777") or degrees ("-24.55")
        Returns: float of Declination in degrees, required to be -90 to +90, inclusive.
    """
    dec_degrees = hex_degrees_as_degrees(dec_string)
    if (dec_degrees < -90) | (dec_degrees > +90):
        dec_degrees = None
    return dec_degrees


def ra_as_hours(ra_degrees):
    """ Input: float of Right Ascension in degrees.
        Returns: string of RA as hours, in hex, to the nearest 0.001 RA seconds.
    """
    if (ra_degrees < 0) | (ra_degrees > 360):
        return None
    n_ra_milliseconds = round((ra_degrees * 3600 * 1000) / 15)
    ra_hours, remainder = divmod(n_ra_milliseconds, 3600 * 1000)
    ra_minutes, remainder = divmod(remainder, 60 * 1000)
    ra_seconds = round(remainder / 1000, 3)
    format_string = "{0:02d}:{1:02d}:{2:06.3f}"
    ra_str = format_string.format(int(ra_hours), int(ra_minutes), ra_seconds)
    if ra_str[:3] == "24:":
        ra_str = format_string.format(0, 0, 0)
    return ra_str


def dec_as_hex(dec_degrees):
    """ Input: float of Declination in degrees.
        Returns: string of Declination in hex, to the nearest 0.01 arcsecond.
    """
    if (dec_degrees < -90) | (dec_degrees > +90):
        return None
    if dec_degrees < 0:
        sign_str = "-"
    else:
        sign_str = "+"
    abs_degrees = abs(dec_degrees)
    milliseconds = round(abs_degrees * 3600 * 1000)
    degrees, remainder = divmod(milliseconds, 3600 * 1000)
    minutes, remainder = divmod(remainder, 60 * 1000)
    seconds = round(remainder / 1000, 2)
    format_string = "{0}{1:02d}:{2:02d}:{3:05.2f}"
    dec_str = format_string.format(sign_str, int(degrees), int(minutes), seconds)
    return dec_str


def weighted_mean(values, weights, return_stdev=False):
    """
    Returns weighted mean, and optionally the weighted std deviation of the *mean*.
    :param values: list (or other iterable) of values to be averaged
    :param weights: list (or other iterable) of weights; length must = length of values
    :param return_stdev: True-> return tuple (w_mean, w_stdev); False-> return w_mean only (float).
    :return: return_stdev==True-> return tuple (w_mean, w_stdev); False-> return w_mean only.
    """
    if (len(values) != len(weights)) or len(values) == 0 or len(weights) == 0:
        raise ValueError('lengths of values & weights must be equal & non-zero.')
    if sum(weights) <= 0:
        raise ValueError('sum of weights must be positive.')
    w_range = range(len(weights))
    norm_weights = [weights[i]/sum(weights) for i in w_range]
    w_mean = sum([norm_weights[i]*values[i] for i in w_range])
    if not return_stdev:
        return w_mean
    if len(values) == 1:
        w_stdev = 0
    else:
        v2 = sum([norm_weights[i]**2 for i in w_range])
        w_stdev = v2 * sum([norm_weights[i]*(values[i]-w_mean)**2 for i in w_range])
    return w_mean, w_stdev


DEFAULT_LADDER = [1.0, 1.25, 1.6, 2.0, 2.5, 3.2, 4.0, 5.0, 6.4, 8.0, 10.0]


def ladder_round(raw_value, ladder=DEFAULT_LADDER, direction="nearest"):
    """
    Rounds to a near-log scale value. May be useful for familiar exposure times.
    Can handle negative numbers, too. Zero returns zero.
    :param raw_value: the value we want to round
    :param ladder: ascending list of values from 1 to 10 to which to round.
    :param direction: "nearest" or "down" or "up"
    :return: raw_valued rounded to nearest ladder value, not counting powers of 10,
    e.g., 32.5 -> 32, 111 -> 100, 6321 -> 6400, -126 -> -125
    """
    if raw_value == 0:
        return 0
    base = math.copysign(10**(math.floor(math.log10(math.fabs(raw_value)))), raw_value)
    target = math.fabs(raw_value / base)
    if target in ladder:
        return raw_value
    for i, val in enumerate(ladder[1:]):
        if target < val:
            ratio_below = target / ladder[i]
            ratio_above = ladder[i+1] / target
            if direction == "down":
                return base * ladder[i]
            if direction == "up":
                return base * ladder[i+1]
            if ratio_below <= ratio_above:  # default case "nearest"
                return base * ladder[i]  # round downward
            else:
                return base * ladder[i+1]  # round upward


def get_phase(jd, jd_epoch, period):
    phase = math.modf((jd - jd_epoch) / period)[0]
    if phase < 0:
        phase += 1
    return phase


def jd_from_datetime_utc(datetime_utc=None):
    if datetime_utc is None:
        datetime_utc = datetime.now(timezone.utc)
    datetime_j2000 = datetime(2000, 1, 1, 0, 0, 0).replace(tzinfo=timezone.utc)
    jd_j2000 = 2451544.5
    seconds_since_j2000 = (datetime_utc - datetime_j2000).total_seconds()
    return jd_j2000 + seconds_since_j2000 / (24*3600)


def datetime_utc_from_jd(jd=None):
    if jd is None:
        return datetime.now(timezone.utc)
    datetime_j2000 = datetime(2000, 1, 1, 0, 0, 0).replace(tzinfo=timezone.utc)
    jd_j2000 = 2451544.5
    seconds_since_j2000 = 24 * 3600 * (jd - jd_j2000)
    return datetime_j2000 + timedelta(seconds=seconds_since_j2000)


def time_hhmm(datetime_utc):
    minutes = round(datetime_utc.hour*60  # NB: banker's rounding (nearest even)
                    + datetime_utc.minute
                    + datetime_utc.second/60
                    + datetime_utc.microsecond/(60*1000000)) % 1440
    hh = minutes // 60
    mm = minutes % 60
    return '{0:0>4d}'.format(100 * hh + mm)


def az_alt_at_datetime_utc(longitude, latitude, target_radec, datetime_utc):
    obs = ephem.Observer()  # for local use.
    if isinstance(longitude, str):
        obs.lon = longitude
    else:
        obs.lon = str(longitude * math.pi / 180)
    if isinstance(latitude, str):
        obs.lat = latitude
    else:
        obs.lat = str(latitude * math.pi / 180)
    obs.date = datetime_utc
    target_ephem = ephem.FixedBody()  # so named to suggest restricting its use to ephem.
    target_ephem._epoch = '2000'
    target_ephem._ra, target_ephem._dec = target_radec.as_hex  # text: RA in hours, Dec in deg
    target_ephem.compute(obs)
    return target_ephem.az * 180 / math.pi, target_ephem.alt * 180 / math.pi


def isfloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def float_or_none(string):
    try:
        return float(string)
    except ValueError:
        return None

