from datetime import datetime

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


class Timespan:
    """ Holds one (start, end) span of time. Immutable.
        Input: 2 UTC datetimes, defining start and end of timespan.
    """
    def __init__(self, start_utc, end_utc):
        if (not isinstance(start_utc, datetime)) or (not isinstance(end_utc, datetime)):
            pass  # TODO: freak out here.
        self.start = start_utc
        self.end = max(start_utc, end_utc)
        self.days = max(0, self.end - self.start)
        self.seconds = 24 * 3600 * self.days

    def intersect(self, timespan2):
        new_start = max(self.start, timespan2.start)
        new_end = min(self.end, timespan2.end)
        return Timespan(new_start, new_end)

    def contains_time(self, time_ut):
        return self.start < time_ut < self.end

    def start(self):
        return self.start_utc

    def end(self):
        return self.end_utc

    def midpoint(self):
        return (self.start_utc + self.end_utc) / 2

    def __str__(self):
        return "Timespan '" + str(self.start) + "' to '" + str(self.end) + "' = " + \
               str(self.seconds) + " seconds."


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
    # dec_list = dec_string.split(":")
    # dec_list = [dec.strip() for dec in dec_list]
    # if dec_list[0].startswith("-"):
    #     sign = -1
    # else:
    #     sign = 1
    # if len(dec_list) == 1:
    #     dec_degrees = float(dec_list[0])  # input assumed to be in degrees.
    # elif len(dec_list) == 2:
    #     dec_degrees = sign * (abs(float(dec_list[0])) + float(dec_list[1])/60.0)  # input is hex.
    # else:
    #     dec_degrees = sign * (abs(float(dec_list[0])) + float(dec_list[1]) / 60.0 +
    #                           float(dec_list[2])/3600.0)  # input is hex.
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


# ---------------------------------------------------------------------

if __name__ == '__main__':
    print(dec_as_hex(0.0))
    print(dec_as_hex(20.0))
    print(dec_as_hex(-69.125))

