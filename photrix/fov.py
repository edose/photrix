
import os
import pandas as pd
import json
from collections import Counter
from datetime import datetime, timezone

from .util import *
from .web import get_aavso_vsp_chart

__author__ = "Eric Dose :: New Mexico Mira Project, Albuquerque"

FOV_DIRECTORY = "C:/Dev/Photometry/FOV/"
CHART_DIRECTORY = 'C:/Dev/Photometry/Chart'
VALID_FOV_OBSERVING_STYLES = ["Standard", "Stare", "Monitor", "LPV"]
CURRENT_SCHEMA_VERSION = "1.5"
LPV_MAG_SINE_FRACTION = 0.5
VR_FRACTION_OF_VI = 0.5


class Fov:
    """
    Object: holds info for one Field Of View (generally maps to one AAVSO sequence or chart).
    For Schema 1.5 of April 2017 (first version defined initially in python/photrix).
    Observation types must be one of: "Stare", "Monitor", "LPV", or "Standard"
    Usage: fov = FOV("ST Tri") or fov = FOV("ST Tri", "C:/Dev/Photometry/FOV1_5/")
    """
    def __init__(self, fov_name, fov_directory=FOV_DIRECTORY):
        fov_fullpath = os.path.join(fov_directory, fov_name + ".txt")
        with open(fov_fullpath) as fov_file:
            lines = fov_file.readlines()
        lines = [line.split(";")[0] for line in lines]  # remove all comments
        lines = [line.strip() for line in lines]

        # ---------- Header section (all directives are required).
        self.fov_name = Fov._directive_value(lines, "#FOV_NAME")
        self.format_version = Fov._directive_value(lines, "#FORMAT_VERSION")
        if self.format_version != CURRENT_SCHEMA_VERSION:
            print(fov_name + ': Fov Version Error')
            raise FovVersionError
        ra_str, dec_str = Fov._directive_words(lines, "#CENTER")[:2]
        self.ra = ra_as_degrees(ra_str)
        self.dec = dec_as_degrees(dec_str)
        self.chart = Fov._directive_words(lines, "#CHART")[0]
        self.fov_date = Fov._directive_words(lines, "#DATE")[0]

        # ---------- Main-target section.
        self.main_target = Fov._directive_value(lines, "#MAIN_TARGET")
        self.target_type = Fov._directive_value(lines, "#TARGET_TYPE")
        self.motive = Fov._directive_value(lines, '#MOTIVE', default_value='')
        self.acp_comments = Fov._directive_value(lines, "#ACP_COMMENTS", default_value='')
        words = Fov._directive_words(lines, "#PERIOD")
        if words is not None:
            self.period = float(words[0])
        else:
            self.period = None
        # As of schema v 1.4, require 2 or 3 values for JD, Mag_V, Color_VI (lengths must match).
        self.JD_bright, self.JD_faint, self.JD_second = (None, None, None)  # default
        words = Fov._directive_words(lines, "#JD")
        if words is not None:
            self.JD_bright = float(words[0])
            self.JD_faint = float(words[1])
            if len(words) >= 3:
                self.JD_second = float(words[2])
        self.mag_V_bright, self.mag_V_faint, self.mag_V_second = (None, None, None)  # default
        words = Fov._directive_words(lines, "#MAG_V")
        if words is not None:
            self.mag_V_bright = float(words[0])
            self.mag_V_faint = float(words[1])
            if len(words) >= 3:
                self.mag_V_second = float(words[2])
        self.color_VI_bright, self.color_VI_faint, \
            self.color_VI_second = (None, None, None)  # default
        words = Fov._directive_words(lines, "#COLOR_VI")
        if words is not None:
            self.color_VI_bright = float(words[0])
            self.color_VI_faint = float(words[1])
            if len(words) >= 3:
                self.color_VI_second = float(words[2])

        # ---------- Observing section.
        obs_style_words = Fov._directive_words(lines, "#OBSERVING_STYLE")
        obs_style = obs_style_words[0]
        obs_values = obs_style_words[1:]
        if obs_style not in VALID_FOV_OBSERVING_STYLES:
            raise ValueError("In '" + self.fov_name + "', '" +
                             obs_style + "' is not a valid Observing Style.")
        self.observing_style = obs_style
        self.alert = None  # default
        self.observing_list = []
        for val in obs_values:
            items = val.split("=")
            tag = items[0]
            # Handle non-filter entries on #OBSERVING_STYLE line.
            if tag == "ALERT" and len(items) >= 2:
                self.alert = float(items[1])
                continue
            # (here, insert any additional future non-filter values like "ALERT=n")

            # Handle filter entries on #OBSERVING_STYLE line.
            this_filter, this_mag, this_count = None, None, None
            if len(items) == 1:  # cases without specified magnitude
                bits = items[0].split("(")
                if len(bits) == 1:  # case "V" for LPVs, one exposure
                    this_filter = bits[0]
                    this_mag = None
                    this_count = 1
                elif len(bits) == 2:  # case "V(2) for LPVs, 2 exposures, e.g.
                    this_filter = bits[0]
                    this_mag = None
                    this_count = int(bits[1].replace(")", ""))
            elif len(items) == 2:  # cases with specified magnitude
                bits = items[1].split("(")
                this_filter = items[0]
                if len(bits) == 1:  # case "V=13.2" for monitors etc, one exposure
                    this_mag = float(bits[0])
                    this_count = 1
                elif len(bits) == 2:  # case "V=13.2(2)" for stares
                    this_mag = float(bits[0])
                    this_count = int(bits[1].replace(")", ""))
            self.observing_list.append((this_filter, this_mag, this_count))
            if this_filter is None:
                raise FovError

        max_exp_value = Fov._directive_value(lines, "#MAX_EXPOSURE")
        if max_exp_value is not None:
            self.max_exposure = float(max_exp_value)
            if self.max_exposure <= 0:
                self.max_exposure = None
        else:
            self.max_exposure = None

        value = Fov._directive_value(lines, "#PRIORITY")
        if value is not None:
            self.priority = float(value)
        else:
            self.priority = None

        value = Fov._directive_value(lines, "#GAP_SCORE_DAYS")
        if value is not None:
            gap_score_words = value.split()
            if len(gap_score_words) >= 3:
                self.gap_score_days = [float(word) for word in gap_score_words[:3]]
            else:
                self.gap_score_days = [self.period * fraction for fraction in [0.01, 0.02, 0.05]]
        else:
            self.gap_score_days = None

    # ---------- AAVSO Sequence section.
        self.punches = Fov._get_punch_values(lines)
        self.aavso_stars = Fov._get_aavso_stars(lines)

    # ---------- Diagnostics and messages before finalizing object.
        if not self.observing_style.lower() == 'standard':
            star_is_check = [star.star_type == "check" for star in self.aavso_stars]
            if not any(star_is_check):
                print(">>>>> WARNING: FOV file ", self.fov_name,
                      " seems not to be a Standard FOV, but has NO CHECK STAR.")

    @staticmethod
    def _directive_value(lines, directive_string, default_value=None):
        for line in lines:
            if line.upper().startswith(directive_string):
                return line[len(directive_string):].strip()
        return default_value  # if directive absent.

    @staticmethod
    def _directive_words(lines, directive_string):
        value = Fov._directive_value(lines, directive_string, default_value=None)
        if value is None:
            return None
        return value.split()

    @staticmethod
    def _get_punch_values(lines):
        punch_values = []
        for line in lines:
            if line.upper().startswith("#PUNCH"):
                value_string = line[len("#PUNCH"):].strip()
                star_id = (value_string.split(":")[0])
                terms = (value_string.split(":")[1]).split()
                if len(terms) == 2:
                    punch_item = (star_id.strip(), float(terms[0]), float(terms[1]))  # tuple
                    punch_values.append(punch_item)
        return punch_values  # list of tuples

    @staticmethod
    def _get_aavso_stars(lines):
        aavso_stars = []
        stars_line_found = False
        for line in lines:
            if line.upper().startswith("#STARS"):
                stars_line_found = True
            else:
                if stars_line_found:
                    aavso_star = AavsoSequenceStar_WithMagErrors(line)
                    # aavso_star = AavsoSequenceStar_MagsOnly(line)  # legacy, from FOV Schema 1.4
                    if aavso_star.is_valid:
                        aavso_stars.append(aavso_star)
        return aavso_stars  # = list of AavsoSequenceStar_WithMagError objects, one per star.

    def calc_gap_score(self, days):
        if days < self.gap_score_days[0]:
            return 0
        elif days < self.gap_score_days[1]:
            return (days - self.gap_score_days[0]) / \
                   (self.gap_score_days[1] - self.gap_score_days[0])
        elif days < self.gap_score_days[2]:
            return 1 + (days - self.gap_score_days[1]) / \
                   (self.gap_score_days[2] - self.gap_score_days[1])
        else:
            return 2

    def calc_priority_score(self, days):
        return self.priority * self.calc_gap_score(days)

    def estimate_lpv_mags(self, jd):
        if self.period <= 0 or self.mag_V_bright >= self.mag_V_faint:
            return None
        jd_bright, jd_faint = self.JD_bright, self.JD_faint
        v_bright, v_faint = self.mag_V_bright, self.mag_V_faint
        color_bright, color_faint = self.color_VI_bright, self.color_VI_faint
        period = self.period
        # First, determine whether brightness is increasing or decreasing at jd.
        # phases must be in [0,1], where phase = 0 means max brightness.
        phase_jd = get_phase(jd, jd_bright, period)
        phase_faint = get_phase(jd_faint, jd_bright, period)
        if 0 <= phase_jd <= phase_faint:  # brightness decreasing at jd
            time_fract = phase_jd / phase_faint
            v_start = v_bright
            color_start = color_bright
            v_change = v_faint - v_bright
            color_change = color_faint - color_bright
        elif phase_faint <= phase_jd <= 1:  # brightness increasing at jd
            time_fract = (phase_jd - phase_faint) / (1 - phase_faint)
            v_start = v_faint
            color_start = color_faint
            v_change = v_bright - v_faint
            color_change = color_bright - color_faint
        else:
            return None  # phase_jd must be outside [0,1]
        #  Now, calculate linear and sine components and blend them (each filter).
        linear_mag_fract = time_fract
        sine_mag_fract = (1 + math.sin((time_fract-0.5)*math.pi)) / 2
        mag_fract = LPV_MAG_SINE_FRACTION * sine_mag_fract + \
                    (1 - LPV_MAG_SINE_FRACTION) * linear_mag_fract
        #  Render mag in each filter.
        mags = dict()
        mags['V'] = v_start + mag_fract * v_change
        mags['I'] = (v_start - color_start) + mag_fract * (v_change - color_change)
        mags['R'] = mags['V'] + VR_FRACTION_OF_VI * (mags['I'] - mags['V'])
        return mags

    def __str__(self):
        return "FOV '" + self.fov_name + "' with " + str(len(self.aavso_stars)) + " sequence stars."


class AavsoSequenceStar_MagsOnly:
    """
    Object: holds parsed mag info (NO mag error info) for one star of an AAVSO sequence.
    FOV 1.4 and earlier only.
    For FOV 1.5+, this is still needed to read (1.4-compatible) pre-FOVs to be merged with chart
       data (by insert_chart_data()).
    """
    def __init__(self, star_string):
        words = star_string.strip().split("\t")   # AAVSO sequence star lines have tab delimiters.
        words = [word.strip() for word in words]  # strip whitespace from all items.
        self.star_id, self.star_type = ("", "")   # invalid defaults.

        if len(words) >= 7:
            # Get raw data strings from line.
            self.star_id = words[0]
            self.ra = float(words[1])   # already in degrees RA
            self.dec = float(words[2])  # already in degrees Dec

            # Extract star type.
            star_type_char = words[4]
            star_type_dict = {"C": "comp", "H": "check", "T": "target"}
            if star_type_char in star_type_dict.keys():
                self.star_type = star_type_dict[star_type_char]
            else:
                self.star_type = None

            # Extract magnitudes. (This code probably sub-optimal.)
            mag_string = words[3]
            self.magB, self.magV, self.magR, self.magI, self.magU = (None, None, None, None, None)
            if self.star_type != "target":
                mag_words = mag_string.split("|")
                mag_words = [mag_word.strip() for mag_word in mag_words]  # strip whitespace
                for mag_word in mag_words:
                    mag_split = mag_word.split("_")
                    this_mag = float(mag_split[1])
                    if this_mag != 0:
                        if mag_split[0] == "1":
                            self.magB = this_mag
                        if mag_split[0] == "2":
                            self.magV = this_mag
                        if mag_split[0] == "4":
                            self.magR = this_mag
                        if mag_split[0] == "8":
                            self.magI = this_mag
                        if mag_split[0] == "1024":
                            self.magU = this_mag

        # Finally, validate this object, or not.
        self.is_valid = (len(self.star_id) >= 2) & (self.star_type is not None)

    def __str__(self):
        return "AAVSO Sequence Star '" + self.star_id + "', is_valid=" + str(self.is_valid)


class AavsoSequenceStar_WithMagErrors:
    """
    Object: holds parsed mag AND mag_error info for ONE star of an AAVSO sequence.
    For FOV 1.5+.
    """

    def __init__(self, star_string):
        words = star_string.strip().split("\t")  # AAVSO sequence star lines have tab delimiters.
        words = [word.strip() for word in words]  # strip whitespace from all items.
        self.star_id, self.star_type = "", ""  # invalid defaults.

        self.mags = {}  # dict of filter:(mag,error); empty will suffice for target stars.
        if len(words) >= 4:
            # Get raw data strings from line:
            self.star_id = words[0]
            self.ra = float(words[1])  # already in degrees RA
            self.dec = float(words[2])  # already in degrees Dec
            if words[3].lower() in ['check', 'target', 'comp']:
                self.star_type = words[3].lower()

            # Get magnitudes and errors for check and comp stars:
            if self.star_type in ['check', 'comp']:
                mag_words = words[4].split(" ")
                mag_words = [mag_word.strip() for mag_word in mag_words]  # strip whitespace
                mag_words = [mag_word for mag_word in mag_words if mag_word != '']  # delete empties
                for mag_word in mag_words:
                    mag_split = mag_word.split("=")
                    this_filter = mag_split[0]
                    mag_error = mag_split[1].split("(")
                    this_mag = float(mag_error[0])
                    this_error = float(mag_error[1].replace(")", "")) / 1000.0  # mMag to Mag
                    if this_mag > 0.0:  # mag missing from mags means no data available.
                        self.mags[this_filter] = (this_mag, this_error)
        # Finally, validate this object, or not. NB: mag_error[0] is the magnitude.
        num_positive_mags = sum([mag_error[0] > 0.0 for mag_error in self.mags.values()])
        self.is_valid = (len(self.star_id) >= 2) and \
                        (self.star_type is not None) and \
                        ((num_positive_mags >= 1) or (self.star_type == 'target'))

    def __str__(self):
        return "AAVSO Sequence Star '" + self.star_id + "', is_valid=" + str(self.is_valid)


class FovError(Exception):
    pass


class FovVersionError(Exception):
    pass


def all_fov_names(fov_directory=FOV_DIRECTORY):
    """ Returns list of all FOV names (from filenames in FOV_DIRECTORY). """
    fov_names = [fname[:-4] for fname in os.listdir(fov_directory)
                 if (fname.endswith(".txt")) and (not fname.startswith("$"))]
    return fov_names


def make_fov_dict(fov_directory=FOV_DIRECTORY, fov_names_selected=None):
    """
    Returns dict of FOVs, as:  FOV_name:FOV object.
    Usage: d = make_fov_dict()  --> returns dict of *all* FOVs.
    Usage: d = make_fov_dict(fov_names_selected=name_list) --> returns selected FOVs.
    """
    fov_all_names = all_fov_names(fov_directory)
    if fov_names_selected is None:
        fov_names = fov_all_names
    else:
        fov_names = list(set(fov_all_names) & set(fov_names_selected))
    fov_dict = {fov_name: Fov(fov_name, fov_directory) for fov_name in fov_names}
    return fov_dict


def print_fov_one_directive_line(fov_directory, all_names, directive, one_line_only=True):
    print("\n\nDirective = '" + directive + '":')
    for fov_name in all_names:
        fov_fullpath = os.path.join(fov_directory, fov_name + ".txt")
        with open(fov_fullpath) as fov_file:
            lines = fov_file.readlines()
        lines = [line.split(";")[0] for line in lines]  # remove all comments
        lines = [line.strip() for line in lines]

        directive_lines = [line for line in lines if line.startswith(directive)]
        error_prefix = "ERROR >>>>>"
        spaces = len(error_prefix) * ' '
        if one_line_only is True:
            if len(directive_lines) == 0:
                print(error_prefix, fov_name, "has NO directive lines for", directive)
            if len(directive_lines) > 1:
                for line in directive_lines:
                    print(error_prefix, fov_name, "MULTIPLE", line)
            if len(directive_lines) == 1:
                print(spaces, fov_name, directive_lines[0])
        else:
            for line in directive_lines:
                print(spaces, fov_name, line)


def fov_diag(fov_directory=FOV_DIRECTORY):
    """
    Comprehensive diagnostics for FOV files.
    This code is for FOV format 1.4.
    :param fov_directory: all .txt files in this directory will be checked [string]
    :return: Number of errors found [integer]
    """
    # Collect all FOVs into a dictionary:
    fd = make_fov_dict(fov_directory=fov_directory)
    fov_names = list(fd.keys())
    fov_names.sort()
    print(len(fov_names), ' FOV files found to test in directory \'', fov_directory, "\'.")

    # Make empty error dictionary:
    error_dict = []  # key=fov_name, value=list of error messages
    for name in fov_names:
        error_dict[name] = []

    # First, verify format versions.
    print("FOV format version required to be : \'" + CURRENT_SCHEMA_VERSION + '\'')
    for name in fov_names:
        fov = fd[name]
        if fov.format_version != CURRENT_SCHEMA_VERSION:
            error_dict[name].append('Format version \'' + fov.format_version + '\'')

    # Verify reasonable JD values:
    jd_min = 2451544.5  # January 1 2000
    jd_max = jd_from_datetime_utc(datetime.now(timezone.utc))  # time at this check
    print("JD limits applied: " + '{0:.3f}'.format(jd_min) + ' to ' + '{0:.3f}'.format(jd_max))
    for name in fov_names:
        fov = fd[name]
        if fov.target_type.lower() != "standard":
            if not jd_min <= fov.JD_bright <= jd_max:
                error_dict[name].append("JD_bright '" + '{0:.3f}'.format(fov.JD_bright) +
                                        "' unreasonable.")
            if not jd_min <= fov.JD_faint <= jd_max:
                error_dict[name].append(": JD_faint '" + '{0:.3f}'.format(fov.JD_faint) +
                                        "' unreasonable.")
            if fov.JD_second is not None:
                if not jd_min <= fov.JD_second <= jd_max:
                    error_dict[name].append(": JD_second '" + '{0:.3f}'.format(fov.JD_second) +
                                            "' unreasonable.")

    # Verify reasonable mag and color values:
    mag_bright = 5.0
    mag_faint = 18.0
    color_min = -0.2
    color_max = +7.5
    print("Mag limits: " + '{0:.3f}'.format(mag_bright) + ' to ' + '{0:.3f}'.format(mag_faint))
    print("Color limits: " + '{0:.3f}'.format(color_min) + ' to ' + '{0:.3f}'.format(color_max))
    for name in fov_names:
        fov = fd[name]
        if fov.target_type.lower() != "standard":
            if not mag_bright <= fov.mag_V_bright <= mag_faint:
                error_dict[name].append("mag_V_bright '" + '{0:.3f}'.format(fov.mag_V_bright) +
                                        "' unreasonable.")
            if not mag_bright <= fov.mag_V_faint <= mag_faint:
                error_dict[name].append("mag_V_faint '" + '{0:.3f}'.format(fov.mag_V_faint) +
                                        "' unreasonable.")
            if fov.mag_V_second is not None:
                if not mag_bright <= fov.mag_V_second <= mag_faint:
                    error_dict[name].append("mag_V_second '" + '{0:.3f}'.format(fov.mag_V_second) +
                                            "' seems unreasonable.")
            if not color_min <= fov.color_VI_bright <= color_max:
                error_dict[name].append("color_VI_bright '" +
                                        '{0:.3f}'.format(fov.color_VI_bright) +
                                        "' seems unreasonable.")
            if not color_min <= fov.color_VI_faint <= color_max:
                error_dict[name].append("color_VI_faint '" + '{0:.3f}'.format(fov.color_VI_faint) +
                                        "' seems unreasonable.")
            if f.color_VI_second is not None:
                if not color_min <= fov.color_VI_second <= color_max:
                    error_dict[name].append("color_VI_second '" +
                                            '{0:.3f}'.format(f.color_VI_second) +
                                            "' seems unreasonable.")

    # Ensure main target is in star list, as a target (skip standard FOVs):
    print(" Ensure main_target in star list, as a target [skip standard FOVs]")
    for name in fov_names:
        fov = fd[name]
        if fov.target_type.lower() != "standard":
            main_target_star_type = [star.star_type for star in fov.aavso_stars
                                     if star.star_id.lower() == fov.main_target.lower()]
            if len(main_target_star_type) <= 0:
                error_dict[name].append("main_target '" + fov.main_target +
                                        "' absent from star list.")
            if len(main_target_star_type) > 1:
                error_dict[name].append("main_target '" + fov.main_target +
                                        "' in star list more than once.")
            if len(main_target_star_type) == 1:
                if main_target_star_type[0] != "target":
                    error_dict[name].append("main_target '" + f.main_target +
                                            "' is in star list once but not as type 'target'.")

    # Ensure Observing styles are valid:
    print("Ensure Observing styles are valid.")
    valid_obs_styles_lower = \
        [valid_obs_style.lower() for valid_obs_style in VALID_FOV_OBSERVING_STYLES]
    for name in fov_names:
        fov = fd[name]
        if fov.observing_style.lower() not in valid_obs_styles_lower:
            error_dict[name].append("invalid obs_style \'" + fov.observing_style + "\'")

    # Ensure JD, mag, color are consistent (skip standard FOVs):
    print("\n", 10*"=", " Ensure mutual consistency of: JD, mag, color (skip standard FOVs)")
    for name in fov_names:
        fov = fd[name]
        if fov.target_type.lower() != "standard":
            # Ensure all present with at least 2 values:
            if None in [fov.JD_bright, fov.JD_faint]:
                error_dict[name].append('missing JD')
            if None in [fov.mag_V_bright, fov.mag_V_faint]:
                error_dict[name].append('missing mag.')
            if None in [fov.color_VI_bright, fov.color_VI_faint]:
                error_dict[name].append('missing color.')

            # Ensure secondary min values are either all present or all absent:
            all_present = None not in [fov.JD_second, fov.mag_V_second, fov.color_VI_second]
            all_absent = fov.JD_second is None and \
                         fov.mag_V_second is None and \
                         fov.color_VI_second is None
            if not (all_present or all_absent):
                error_dict[name].append('mismatched JD, mag, color (secondary min?).')

    # Alert on out-of-spec phases, or non-positive periods (Eclipser-like only):
    print("Alert on out-of-spec phases & non-positive periods (Eclipser-like only)")
    for name in fov_names:
        fov = fd[name]
        if fov.target_type.lower() in ['eclipser', 'exoplanet']:
            if fov.period <= 0:
                error_dict[name].append("PERIOD=" + '{0:8.3f}'.format(fov.period))
            else:
                # Verify that max JD is reasonable.
                phase_max = ((fov.JD_bright - fov.JD_faint) / fov.period) % 1.0
                if abs(phase_max-0.25) > 0.05 or abs(phase_max-0.75) > 0.05:
                    error_dict[name].append('Max phase of ' +
                                            '{0:.3f}'.format(phase_max) + ' unreasonable.')
                if fov.JD_second is not None:
                    phase_second = ((fov.JD_second - fov.JD_faint) / fov.period) % 1.0
                    if abs(phase_second-0.5) > 0.02:
                        error_dict[name].append('Secondary phase of ' +
                                                '{0:.3f}'.format(phase_second) + ' unreasonable.')

    # Finally, write out all errors, by fov name:
    num_errors = 0
    for name in fov_names:
        num_errors += len(error_dict[name])
    print(str(num_errors) + ' errors found.')

    for name in fov_names:
        fov_errors = error_dict[name]
        if len(fov_errors) >= 1:
            print('\n' + name + ':')
            for error in fov_errors:
                print(4*'' + error)


def fovs_by_ra(fov_directory=FOV_DIRECTORY):
    fov_dict = make_fov_dict(fov_directory=fov_directory)
    fov_names = list(fov_dict.keys())
    df_fov = pd.DataFrame({'fov_name': fov_names})  # 1 column ('fov_name') only.
    df_fov['obs_style'] = [fov_dict[name].observing_style for name in fov_names]
    df_fov['ra_hours'] = [fov_dict[name].ra/15.0 for name in fov_names]
    dict_category = {'Stare':['Stare'], 'LPV/Monitor':['LPV', 'Monitor']}
    print('\nFOV Counts as of ' + '{:%Y-%m-%d %H:%M  UTC}'.format(datetime.now(timezone.utc)))
    for category, styles in dict_category.items():
        styles_lower = [style.lower() for style in styles]
        matches_style = [s.lower() in styles_lower for s in df_fov['obs_style']]
        df_style = df_fov[matches_style]
        ra_list = df_style['ra_hours']
        hist = 24*[0]
        for ra in ra_list:
            int_ra = int(min(23.0, floor(ra)))
            hist[int_ra] += 1
        print('\n\nObs category \'' + category + '\' (' + str(len(ra_list)) + ' FOVs):')
        for i_ra in range(24):
            print('   ' + '{:2d}'.format(i_ra) + ': ' + '{:4d}'.format(hist[i_ra]))


def delete_directive(fov_directory=FOV_DIRECTORY, out_fov_directory=None,
                     directive_to_remove=None):
    """
    Not Tested.
    :param fov_directory: 
    :param out_fov_directory: 
    :param directive_to_remove: 
    :return: [None]
    """
    if out_fov_directory is None or directive_to_remove is None:
        print('\n\nPlease give a new FOV directory and a directive to remove.\n\n')
        return
    names = all_fov_names(fov_directory)
    print(str(len(names)) + ' FOVs to adjust.')
    os.makedirs(out_fov_directory, exist_ok=True)  # make output directory if doesn't exist.
    for name in names:
        in_fullpath = os.path.join(fov_directory, name + ".txt")
        with open(in_fullpath) as fov_file:
            lines = fov_file.readlines()
        out_lines = []
        for line in lines:
            if not line.startswith(directive_to_remove):
                out_lines.append(line)
        out_fullpath = os.path.join(out_fov_directory, name + ".txt")
        with open(out_fullpath, 'w') as out_file:
            out_file.writelines(out_lines)
    print('Done.')


def move_directive(fov_directory=FOV_DIRECTORY, out_fov_directory=None,
                   directive_to_move=None, directive_before_new_position=None):
    """
    :param fov_directory: 
    :param out_fov_directory: 
    :param directive_to_move: 
    :param directive_before_new_position: 
    :return: [None]
    """
    if out_fov_directory is None or directive_to_move is None or \
            directive_before_new_position is None:
        print('\n\nPlease give a new FOV directory. directive to move, and '
              'directive to insert after.\n\n')
        return
    names = all_fov_names(fov_directory)
    print(str(len(names)) + ' FOVs to adjust.')
    os.makedirs(out_fov_directory, exist_ok=True)  # make output directory if doesn't exist.
    for name in names:
        in_fullpath = os.path.join(fov_directory, name + ".txt")
        with open(in_fullpath) as fov_file:
            lines = fov_file.readlines()

        # Extract line to move:
        for line in lines:
            if line.startswith(directive_to_move):
                line_to_move = line
                break
        lines.remove(line_to_move)

        # Find new location and do insertion.
        new_lines = []
        for line in lines:
            new_lines.append(line)
            if line.startswith(directive_before_new_position):
                new_lines.append(line_to_move)

        out_fullpath = os.path.join(out_fov_directory, name + ".txt")
        with open(out_fullpath, 'w') as out_file:
            out_file.writelines(new_lines)
    print('Done.')


def change_directive_value(fov_directory=FOV_DIRECTORY, out_fov_directory=None,
                           directive_to_change=None, new_value=None, new_comment=None):
    """
    May be used to update FORMAT_VERSION, especially.
    :param fov_directory: source directory [string]
    :param out_fov_directory: output directory, may be same as source directory [string]
    :param directive_to_change: e.g., '#FORMAT_VERSION'
    :param new_value: e.g., '1.5' [string]
    :param new_comment:  [string, or None]
    :return: [nothing]
    """
    if fov_directory is None or directive_to_change is None or \
            new_value is None:
        print('\n\nPlease give a new FOV directory. directive to change, and '
              'new directive value.\n\n')
        return
    names = all_fov_names(fov_directory)
    os.makedirs(out_fov_directory, exist_ok=True)  # make output directory if doesn't exist.
    for name in names:
        new_lines = []
        in_fullpath = os.path.join(fov_directory, name + ".txt")
        with open(in_fullpath) as fov_file:
            lines = fov_file.readlines()
        for line in lines:
            if line.startswith(directive_to_change):
                if ';' in line:
                    halves = line.split(';', maxsplit=1)
                    value = (halves[0])[len(directive_to_change):]
                    comment = halves[1]
                else:
                    value = line[len(directive_to_change):]
                    comment = None
                num_leading_spaces = len(value) - len(value.lstrip())
                if new_comment is not None:
                    comment = new_comment
                new_line = directive_to_change + (num_leading_spaces * ' ') +\
                           new_value + ' ;' + comment + '\n'
                new_lines.append(new_line)
            else:
                new_lines.append(line)

        out_fullpath = os.path.join(out_fov_directory, name + ".txt")
        with open(out_fullpath, 'w') as out_file:
            out_file.writelines(new_lines)
        # print(name)
    print('Done.')


def insert_chart_data(fov_name, fov_directory=FOV_DIRECTORY):
    """
    Takes one pre-FOV file already having VPhot sequence mags, adds VSP chart data (mag errors)
        to render a valid FOV 1.5+ file which overwrites the pre-FOV file.
    Gets chart ID from FOV #CHART directive. Looks in CHART_DIRECTORY for chart JSON file, or 
        failing that downloads it from AAVSO and saves it.
    Handles only filters U,B,V,R,I as of 7/1/2017.
    :param fov_name: one FOV name [string]
    :param fov_directory: input and output directory (same) [string]
    :return: list of warning lines [list of strings]
    """
    # global mag_error_string
    warning_lines = []

    # Load FOV text file with VPhot data only, divide into top lines & #STARS lists:
    fov_fullpath = os.path.join(fov_directory, fov_name + ".txt")
    with open(fov_fullpath) as fov_file:
        lines = fov_file.readlines()
    top_lines = []
    stars_lines = []
    above_stars_directive = True
    for line in lines:
        if above_stars_directive:
            top_lines.append(line)  # top_lines includes the #STARS line.
            if line.startswith("#STARS"):
                above_stars_directive = False
        else:
            stars_lines.append(line)

    # Read JSON chart file, or download from VSP (JSON format) and cache it:
    os.makedirs(CHART_DIRECTORY, exist_ok=True)  # create directory if doesn't already exist.
    chart_id = None
    for line in top_lines:
        if line.startswith("#CHART"):
            chart_id = line.split(";")[0][6:].strip()
            break
    if chart_id is None:
        error_line = "No #CHART directive in fov " + fov_name
        print(error_line)
        warning_lines.append(error_line)
        return warning_lines
    chart_fullpath = os.path.join(CHART_DIRECTORY, chart_id + ".txt")
    if os.path.exists(chart_fullpath):
        with open(chart_fullpath, 'r') as chart_file:
            chart_json_text = chart_file.read()
    else:
        chart_json_text = get_aavso_vsp_chart(chart_id)
        print('Downloading chart \'' + chart_id + '\' for FOV \'' + fov_name + '\'')
        if chart_json_text == '':
            error_line = '>>>>> No chart \'' + chart_id +\
                         '\' in AAVSO VSP (or chart directory. No change made to fov \'' + \
                         fov_name + '\'.'
            print(error_line)
            warning_lines.append(error_line)
            return warning_lines
        with open(chart_fullpath, 'w') as fov_file:
            fov_file.write(chart_json_text)  # cache json
    json_obj = json.loads(chart_json_text)
    chart_stars = json_obj['photometry']

    # Match a chart star to each FOV star, make a new text line with both mags and errors:
    new_star_lines = []

    for line in stars_lines:
        if line.split(";")[0].strip() == '':
            new_star_lines.append(line)  # a comment line: copy as is.
            continue

        fov_star = AavsoSequenceStar_MagsOnly(line.strip())
        if fov_star.star_type.lower() in ['comp', 'check']:
            star_prefix = int(fov_star.star_id.split("_")[0])  # int from before any "_".
            star_radec = RaDec(fov_star.ra, fov_star.dec)
            filter_dict = dict()
            filter_dict['U'] = (fov_star.magU, 0.0)  # mag, default mag_error
            filter_dict['B'] = (fov_star.magB, 0.0)
            filter_dict['V'] = (fov_star.magV, 0.0)
            filter_dict['R'] = (fov_star.magR, 0.0)
            filter_dict['I'] = (fov_star.magI, 0.0)
            # Redefine mag=None to mag=0.0 (FOV 1.5+ convention).
            for filter_name in filter_dict.keys():
                mag, error = filter_dict[filter_name]
                filter_dict[filter_name] = (mag if mag is not None else 0.0, error)

            # Find chart star entry (if any) that matches FOV star line.
            chart_star_found = False
            for chart_star in chart_stars:
                if chart_star['label'] == star_prefix:
                    chart_star_radec = RaDec(chart_star['ra'], chart_star['dec'])
                    if star_radec.degrees_from(chart_star_radec) < 20.0 / 3600.0:
                        # Here, we've found the chart star matching the FOV star.
                        chart_star_found = True
                        # Take union of FOV filters & chart bands for this star, parse mags, errors:
                        # Then update with chart mag errors, adding a dict entry if absent:
                        for band in chart_star['bands']:  # which is a list of dicts
                            band_translate = {'Ic': 'I', 'Rc': 'R'}  # chart -> FOV band names
                            fov_filter_name = band_translate.get(band['band'], band['band'])
                            fov_filter_data = filter_dict.get(fov_filter_name, None)
                            if fov_filter_data is not None:
                                new_mag = filter_dict[fov_filter_name][0]
                            else:
                                new_mag = band['mag']
                                # Chart may have either null or 0 to signify 'no error data':
                            new_error = band['error'] if band['error'] is not None else 0.0
                            filter_dict[fov_filter_name] = (new_mag, new_error)
                        break  # go to next chart star.
            if not chart_star_found:
                warning_lines.append('FOV ' + fov_name +
                                     ', star ' + fov_star.star_id +
                                     ': no matching star found in chart ' + chart_id +
                                     ' w/ maglimit ' + str(json_obj['maglimit']))

            # Mags & errors for this star are ready; add text to this output line:
            mag_error_list = []
            for filter_name in filter_dict.keys():
                mag, error = filter_dict[filter_name]
                if error is None:
                    # TODO: add error handling
                    dummy = 0
                if mag == 0 and error == 0:
                    mag_error_list.append(filter_name + '=0(0)')  # abbreviated
                else:
                    mag_error_list.append(filter_name + '={0:.3f}'.format(mag) +
                                          '({0:d})'.format(round(1000 * error)))
            mag_error_string = ' '.join(mag_error_list)
        else:
            mag_error_string = ''  # for a target star.
        new_star_text_list = [str(fov_star.star_id),
                              '{0:9.5f}'.format(fov_star.ra),
                              '{0:9.5f}'.format(fov_star.dec),
                              fov_star.star_type,
                              mag_error_string]
        new_star_line = '\t'.join(new_star_text_list) + '\n'
        new_star_lines.append(new_star_line)

    # Build list of all output lines, then write to new FOV file:
    new_lines = top_lines.copy()
    new_lines.extend(new_star_lines)
    with open(fov_fullpath, 'w') as fov_file:
        fov_file.writelines(new_lines)
    return warning_lines


def get_chart_error(bands, chart_filter_name):
    """
    Helper fn for insert_chart_data().
    :param bands: list of band dicts, as extracted from chart json.
    :param chart_filter_name: chart-style filter name, as "B" or "Ic" (not "I").
    :return: error if found, zero otherwise [float]
    """
    chart_error = 0  # default if entry for filter not found.
    for i in range(len(bands)):
        if bands[i]['band'] == chart_filter_name:
            chart_error = bands[i]['error'] if bands[i]['error'] is not None else 0
            break
    return chart_error


# def all_fovs_1_4_to_1_5(fov_1_4_directory=FOV_DIRECTORY, fov_1_5_directory=None):
#     """
#     Convenience function to upgrade all FOVs in fov_1_4_directory from 1.4 to 1.5
#     :param fov_1_4_directory: directory from which to read (all) FOV 1.4 files  [string]
#     :param fov_1_5_directory: director to which to write FOV 1.5 files [string]
#     :return: summary of operations [string]
#     """
#     # Create new empty FOV 1.5 directory if needed:
#     if fov_1_5_directory is None:
#         if fov_1_4_directory.endswith('/FOV/'):
#             fov_1_5_directory = fov_1_4_directory.split('/FOV/')[0] + '/FOV1_5/'
#         else:
#             print('You must provide a valid fov_1_5_directory. Stopping.')
#             return
#     os.makedirs(fov_1_5_directory, exist_ok=True)  # create directory if doesn't already exist.
#
#     # Delete any pre-existing files in 1.5 directory:
#     filenames = os.listdir(fov_1_5_directory)
#     for filename in filenames:
#         fullpath = os.path.join(fov_1_5_directory, filename)
#         os.remove(fullpath)
#
#     # Populate new FOV 1.5 directory with all FOV 1.4 files, to begin:
#     import shutil
#     objectnames = os.listdir(fov_1_4_directory)
#     filenames = [name for name in objectnames if name.endswith(".txt")]
#     for filename in filenames:
#         source_path = os.path.join(fov_1_4_directory, filename)
#         shutil.copy2(source_path, fov_1_5_directory)
#     print('All ' + str(len(filenames)) + ' FOVs copied.')
#
#     # Perform all one-line 1.4 --> 1.5 changes (in new directory):
#     delete_directive(fov_1_5_directory, fov_1_5_directory, '#STARE_HOURS')
#     delete_directive(fov_1_5_directory, fov_1_5_directory, '#ACP_DIRECTIVES')
#     move_directive(fov_1_5_directory, fov_1_5_directory,
#                    directive_to_move='#ACP_COMMENTS', directive_before_new_position='#MOTIVE')
#     change_directive_value(fov_1_5_directory, fov_1_5_directory, '#FORMAT_VERSION', '1.5',
#                            ' FOV format version defined April 24 2017')
#     print('All one-line changes made.')
#
#     # For each FOV, update #STARS section to include mag errors from AAVSO VSP chart:
#     names = all_fov_names(fov_1_5_directory)
#     # names = [name for name in names if name.startswith('A')]  # <<< Remove this in production!
#     all_warning_lines = []
#     for name in names:
#         warning_lines = insert_chart_data(name, fov_1_5_directory)
#         all_warning_lines.extend(warning_lines)
#     print('\n'.join(all_warning_lines) + '\nDone.')
