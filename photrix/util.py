from datetime import datetime, timezone, timedelta
import math

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


class Timespan:
    """ Holds one (start, end) span of time. Immutable.
        Input: 2 python datetimes (in UTC), defining start and end of timespan.
    """
    def __init__(self, start_utc, end_utc):
        self.start = start_utc
        self.end = max(start_utc, end_utc)
        self.seconds = (self.end - self.start).seconds
        self.midpoint = self.start + timedelta(seconds=self.seconds / 2)

    def intersect(self, timespan2):
        new_start = max(self.start, timespan2.start)
        new_end = min(self.end, timespan2.end)
        return Timespan(new_start, new_end)

    def contains_time(self, time_utc):
        return self.start <= time_utc <= self.end

    def __str__(self):
        return "Timespan '" + str(self.start) + "' to '" + str(self.end) + "' = " + \
               str(self.seconds) + " seconds."


class RaDec:
    """
    Holds one Right Ascension, Declination sky position (internally as degrees).
    ra : (hours hex string, or degrees)
    dec : (degrees hex string, or degrees)
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

    def as_degrees(self):
        return ra_as_degrees(self.ra), dec_as_degrees(self.dec)

    def as_hex(self):
        return ra_as_hours(self.ra), dec_as_hex(self.dec)

    def alt_az(self, latitude, longitude, time_utc):
        pass

    def degrees_from(self, other_ra_dec):
        pass

    def farther_from(self, other_ra_dec, degrees_limit):
        return self.degrees_from(other_ra_dec) > degrees_limit

    def __str__(self):
        ra_hex, dec_hex = self.as_hex()
        return "RaDec object:  " + ra_hex + "  " + dec_hex




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
    ra_str = format_string.format(ra_hours, ra_minutes, ra_seconds)
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
    dec_str = format_string.format(sign_str, degrees, minutes, seconds)
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

