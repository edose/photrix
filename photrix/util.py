from datetime import datetime, timedelta, timezone
import math
from math import floor, sqrt

import numpy as np
import pandas as pd
import statsmodels.regression.mixed_linear_model as sm  # statsmodels version >= 0.8 !
import ephem

__author__ = "Eric Dose :: New Mexico Mira Project, Albuquerque"


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
    dec_string = degrees_as_hex(dec_degrees, seconds_decimal_places=2)
    return dec_string


def degrees_as_hex(angle_degrees, seconds_decimal_places=2):
    """
    :param angle_degrees: any angle as degrees
    :return: same angle in hex notation, unbounded.
    """
    if angle_degrees < 0:
        sign = "-"
    else:
        sign = "+"
    abs_degrees = abs(angle_degrees)
    milliseconds = round(abs_degrees * 3600 * 1000)
    degrees, remainder = divmod(milliseconds, 3600 * 1000)
    minutes, remainder = divmod(remainder, 60 * 1000)
    seconds = round(remainder / 1000, 2)
    format_string = '{0}{1:02d}:{2:02d}:{3:0' + str(int(seconds_decimal_places)+3) + \
                    '.0' + str(int(seconds_decimal_places)) + 'f}'
    hex_string = format_string.format(sign, int(degrees), int(minutes), seconds)
    return hex_string


def weighted_mean(values, weights):
    """
    Returns weighted mean and weighted std deviation of the mean (not of observations).
    :param values: list (or other iterable) of values to be averaged
    :param weights: list (or other iterable) of weights; length must = length of values
    :return: 3-tuple (weighted mean, weighted std dev (population), weighted std dev of mean)
    """
    if (len(values) != len(weights)) or (len(values) == 0) or (len(weights) == 0):
        raise ValueError('lengths of values & weights must be equal & non-zero.')
    if sum(weights) <= 0:
        raise ValueError('sum of weights must be positive.')
    value_list = list(values)    # py list comprehension often misunderstands pandas Series indices.
    weight_list = list(weights)  # "
    norm_weights = [wt/sum(weights) for wt in weight_list]
    w_mean = sum([nwt * val for (nwt, val) in zip(norm_weights, value_list)])
    n_nonzero_weights = sum([w != 0 for w in weight_list])

    if n_nonzero_weights == 1:
        w_stdev_pop = 0
        w_stdev_w_mean = 0
    else:
        resid2 = [(val-w_mean)**2 for val in value_list]
        nwt2 = sum([nwt**2 for nwt in norm_weights])
        rel_factor = 1.0 / (1.0 - nwt2)  # reliability factor (better than N'/(N'-1))
        w_stdev_pop = sqrt(rel_factor * sum([nwt * r2 for (nwt, r2) in zip(norm_weights, resid2)]))
        w_stdev_w_mean = sqrt(nwt2) * w_stdev_pop
    return w_mean, w_stdev_pop, w_stdev_w_mean


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
        return None
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


def hhmm_from_datetime_utc(datetime_utc):
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


def event_utcs_in_timespan(jd_reference, period, timespan):
    """
    Returns a list of UTC times of period events within a given Timespan.
       A generalization of (and replacing) fn find_minima_in_timespan()
    :param jd_reference: Julian Date of any occurence of the period event (e.g., Mira max) [float]
    :param period: in days [float]
    :param timespan: target timespan (start and end datetimes) [Timespan object]
    :return: list of up to 10 UTCs of periodic events within the target timespan [list of datetimes]
       Return None if jd_reference or period are invalid. Return empty list if no such events.
    """
    if jd_reference is None or period is None:
        return None
    jd_ts_start = jd_from_datetime_utc(timespan.start)
    jd_ts_end = jd_from_datetime_utc(timespan.end)
    n_prior = floor((jd_ts_start - jd_reference) / period)
    jd_prior = jd_reference + n_prior * period
    utc_list = []
    for i in range(10):
        jd_test = jd_prior + i * period
        if jd_test > jd_ts_end:
            return utc_list
        if jd_test >= jd_ts_start:
            utc_list.append(datetime_utc_from_jd(jd_test))
    return utc_list


class MixedModelFit:
    """
    Object: holds info for one mixed-model (py::statsmodel) fit. 
    Generic in nature--NOT tied to astronomical usage.
    Uses formula form, i.e., statsmodel::sm.MixedLM.from_formula()
    """
    def __init__(self, data, dep_var=None, fixed_vars=None, group_var=None):
        """
        Executes mixed-model fit & makes data available.
        :param data: input data, one variable per column, one point per row [pandas Dataframe]
        :param dep_var: one column name as dependent 'Y' variable [string] 
        :param fixed_vars: one or more column names as independent 'X' variable [string or
                  list of strings]
        :param group_var: one column name as group (category; random-effect) variable [string]
        Usage: fit = MixedModel(df_input, 'Y', ['X1', 'X2'], 'a_group_type']
               fit = MixedModel(df_input, 'Y', 'X1', 'a_group_type'] (OK if only one indep var)
        """
        if not isinstance(data, pd.DataFrame):
            print('Parameter \'data\' must be a pandas Dataframe of input data.')
            return
        if dep_var is None or fixed_vars is None or group_var is None:
            print('Provide all parameters: dep_var, fixed_vars, and group_var.')
            return
        if not isinstance(dep_var, str) or not isinstance(group_var, str):
            print('Parameters \'dep_var\' and \'group_var\' must both be strings.')
            return
        fixed_vars_valid = False  # default if not validated
        if isinstance(fixed_vars, str):
            fixed_vars = list(fixed_vars)
            fixed_vars_valid = True
        if isinstance(fixed_vars, list):
            if len(fixed_vars) >= 1:
                if all([isinstance(var, str) for var in fixed_vars]):
                    fixed_vars_valid = True
        if not fixed_vars_valid:
            print('Parameter \'fixed_vars\' must be a string or a list of strings.')
            return
        formula = dep_var + ' ~ ' + ' + '.join(fixed_vars)

        model = sm.MixedLM.from_formula(formula, groups=data[group_var], data=data)
        fit = model.fit()

        self.statsmodels_object = fit  # instance of class MixedLMResults (py pkg statsmodels)

        # Scalar and naming attributes:
        self.converged = fit.converged  # bool
        self.nobs = fit.nobs  # number of observations used in fit
        self.likelihood = fit.llf
        self.dep_var = dep_var
        self.fixed_vars = fixed_vars
        self.group_var = group_var
        self.sigma = sqrt(sum(fit.resid**2)/(fit.nobs-len(fixed_vars)-2))

        # Fixed-effects dataframe (joins so we don't count on consistent input ordering):
        df = pd.DataFrame({'Value': fit.fe_params})
        df = df.join(pd.DataFrame({'Stdev': fit.bse_fe}))     # join on index (enforce consistency)
        df = df.join(pd.DataFrame({'Tvalue': fit.tvalues}))   # " & any random effect discarded
        df = df.join(pd.DataFrame({'Pvalue': fit.pvalues}))   # " & "
        df['Name'] = df.index
        self.df_fixed_effects = df.copy()

        # Random-effect dataframe, index=GroupName, cols=GroupName, GroupValue:
        df = pd.DataFrame(fit.random_effects).transpose()  # DataFrame, 1 row/group
        df = df.rename(columns={'groups': 'GroupValue'})
        df['GroupName'] = df.index
        self.df_random_effects = df.copy()

        # Observation dataframe (safe to count on consistent input ordering -> easier construction):
        df = pd.DataFrame({'FittedValue': fit.fittedvalues})
        df['Residual'] = fit.resid
        self.df_observations = df.copy()

        pass

    def predict(self, df_predict_input, include_random_effect=True):
        """
        Takes new_data and renders predicted dependent-variable values.
        Optionally includes effect of groups (random effects), unlike py::statsmodels.
        :param: new_data: new input data used to render predictions. 
           Extra (unused) columns OK; model selects only needed columns. [pandas DataFrame] 
        :param: include_random_effect: True to include them, False to omit/ignore [bool]
        :return: predictions of dependent-variable values matching rows of new data (pandas Series)
        """

        # Get predicted values on fixed effects only (per statsmodels' weird def. of 'predicted'):
        fixed_effect_inputs = df_predict_input[self.fixed_vars]  # 1 col per fixed effect variable
        predicted_on_fixed_only = self.statsmodels_object.predict(exog=fixed_effect_inputs)

        # If requested, add RE contibs (that were not included in MixedModels object 'fit'):
        if include_random_effect:
            df_random_effect_inputs = pd.DataFrame(df_predict_input[self.group_var])
            df_random_effect_values = self.df_random_effects[['GroupValue']]
            predicted_on_random_only = pd.merge(df_random_effect_inputs, df_random_effect_values,
                                                left_on=self.group_var,
                                                right_index=True, how='left',
                                                sort=False)['GroupValue']  # Series (left-join)
            # Random effect is ***SUBTRACTED***, because original fit was
            #    InstMag ~ CatMag + Random effect + offsets + other fixed effects
            #    and now its surrogate CatMag we're trying to estimate
            total_prediction = predicted_on_fixed_only - predicted_on_random_only
        else:
            total_prediction = predicted_on_fixed_only

        return total_prediction
