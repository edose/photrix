from photrix.util import *
import os

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"

FOV_DIRECTORY = "C:/Dev/Photometry/FOV/"
CURRENT_SCHEMA_VERSION = "1.3"


class FOV_list:
    def __init__(self, fov_name_list):
        pass


class FOV:
    """
    Object: holds info for one Field Of View (generally maps to one AAVSO sequence or chart).
    For Schema 1.3 of Sept 2016 (as defined initially in R code "photometry").
    Observation types must be one of: "Stare", "Monitor", "LPV", "Standard"
    Usage: fov = FOV("ST Tri")
    """
    def __init__(self, fov_name, fov_directory=FOV_DIRECTORY):
        fov_fullpath = os.path.join(fov_directory, fov_name + ".txt")
        with open(fov_fullpath) as fov_file:
            lines = fov_file.readlines()
        lines = [line.split(";")[0] for line in lines]  # remove comments

        # ---------- Header section.
        self.fov_name = FOV._directive_value(lines, "#FOV_NAME")
        self.format_version = FOV._directive_value(lines, "#FORMAT_VERSION")
        if self.format_version != CURRENT_SCHEMA_VERSION:
            raise FOV_Error
        ra_str, dec_str = FOV._directive_words(lines, "#CENTER")[:2]
        self.ra = ra_as_degrees(ra_str)
        self.dec = dec_as_degrees(dec_str)
        self.chart = FOV._directive_words(lines, "#CHART")[0]
        self.fov_date = FOV._directive_words(lines, "#DATE")[0]

        # ---------- Main-target section.
        self.main_target = FOV._directive_value(lines, "#MAIN_TARGET")
        self.target_type = FOV._directive_value(lines, "#TARGET_TYPE")
        words = FOV._directive_words(lines, "#PERIOD")
        if words is not None:
            self.period = float(words[0])
        else:
            self.period = None
        # As of schema v 1.3, require both values for JD, Mag_V, Color_VI.
        self.JD_bright, self.JD_faint = (None, None)  # default
        words = FOV._directive_words(lines, "#JD")
        if words is not None:
            self.JD_bright = float(words[0])
            self.JD_faint = float(words[1])
        self.mag_V_bright, self.mag_V_faint = (None, None)  # default
        words = FOV._directive_words(lines, "#MAG_V")
        if words is not None:
            self.mag_V_bright = float(words[0])
            self.mag_V_faint = float(words[1])
        self.color_VI_bright, self.color_VI_faint = (None, None)  # default
        words = FOV._directive_words(lines, "#COLOR_VI")
        if words is not None:
            self.color_VI_bright = float(words[0])
            self.color_VI_faint = float(words[1])

        # ---------- Observing section.
        obs_style_words = FOV._directive_words(lines, "#OBSERVING_STYLE")
        obs_style = obs_style_words[0]
        obs_values = obs_style_words[1:]
        if obs_style not in ["Stare", "Monitor", "LPV", "Standard"]:
            raise ValueError("'" + obs_style + "' is not a valid Observing Style.")

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
            this_filter = None
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
                raise FOV_Error

        stare_hours_words = FOV._directive_words(lines, "#STARE_HOURS")
        if stare_hours_words is not None:
            if len(stare_hours_words) >= 1:
                if stare_hours_words[0] == "ANY":
                    if len(stare_hours_words) >= 2:
                        self.stare_reference = stare_hours_words[0]
                        self.stare_start = 0
                        self.stare_stop = float(stare_hours_words[1])
                elif stare_hours_words[0] in ["MIN", "MAX"]:
                    if len(stare_hours_words) >= 3:
                        self.stare_reference = stare_hours_words[0]
                        self.stare_start = float(stare_hours_words[1])
                        self.stare_stop = float(stare_hours_words[2])
        else:
            self.stare_reference = None  # defaults
            self.stare_start = None
            self.stare_stop = None

        max_exp_value = FOV._directive_value(lines, "#MAX_EXPOSURE")
        if max_exp_value is not None:
            self.max_exposure = float(max_exp_value)
            if self.max_exposure <= 0:
                self.max_exposure = None
        else:
            self.max_exposure = None

        value = FOV._directive_value(lines, "#PRIORITY")
        if value is not None:
            self.priority = float(value)
        else:
            self.priority = None
        value = FOV._directive_value(lines, "#GAP_SCORE_DAYS")
        if value is not None:
            gap_score_words = value.split()
            if len(gap_score_words) >= 3:
                self.gap_score_days = [float(word) for word in gap_score_words[:3]]
            else:
                self.gap_score_days = [self.period * fraction for fraction in [0.01, 0.02, 0.05]]
        else:
            self.gap_score_days = None
        self.acp_directives = FOV._directive_value(lines, "#ACP_DIRECTIVES").split("|")
        self.acp_comments = FOV._directive_value(lines, "#ACP_COMMENTS")

    # ---------- AAVSO Sequence section.
        self.punches = FOV._get_punch_values(lines)
        self.aavso_stars = FOV._get_aavso_stars(lines)

    # ---------- Diagnostics and messages before finalizing object.
        if not self.fov_name.startswith("Std_"):
            star_is_check = [star.star_type == "check" for star in self.aavso_stars]
            if not any(star_is_check):
                print(">>>>> WARNING: FOV file ", self.fov_name,
                      " seems not to be a Standard FOV, but has NO CHECK STAR.")

    @staticmethod
    def _directive_value(lines, directive_string):
        for line in lines:
            if line.upper().startswith(directive_string):
                return line[len(directive_string):].strip()
        return None  # if directive absent.

    @staticmethod
    def _directive_words(lines, directive_string):
        value = FOV._directive_value(lines, directive_string)
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
                    aavso_star = AavsoSequenceStar(line)
                    if aavso_star.is_valid:
                        aavso_stars.append(aavso_star)
        return aavso_stars

    def __str__(self):
        return "FOV '" + self.fov_name + "' with " + str(len(self.aavso_stars)) + " sequence stars."


class FOV_Error:
    pass


class AavsoSequenceStar:
    """
    Object: holds parsed info for one star of an AAVSO sequence.
    Internal use only (by class FOV).
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
        return "AAVSO Sequence Star '" + self.star_id + "', is_valid=" + self.is_valid





# ---------------------------------------------------------------------

if __name__ == '__main__':
    f = FOV("ST Cnc")
    print(f.name, "\n", repr(f))
