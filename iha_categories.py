# iha_categories.py
# -------------------
# Single shared source of truth for the five-category directional
# hydrologic alteration (HA) classification scheme (manuscript Methods
# 3.4), used by chess_scape_iha_analysis.py and, via figures/shared_utils.py,
# by all figure-generation scripts in figures/.
#
# Every category boundary (0.33, 0.67, -0.33, -0.67) is resolved
# consistently by assigning it to the lower-magnitude-of-alteration
# category, so bins are contiguous with no gaps or overlaps -- verified
# against all boundary values from -1.0 to +1.0 (the symmetric HA range
# that results from EXPECTED_FREQ = 0.5; see chess_scape_iha_analysis.py).

CAT_COLOURS = {
    "low":                "#cccccc",   # light grey    -- low alteration
    "mod_departure":      "#6baed6",   # medium blue   -- moderate departure
    "high_departure":     "#08519c",   # dark blue     -- high departure
    "mod_concentration":  "#fd8d3c",   # medium orange -- moderate concentration
    "high_concentration": "#a50f15",   # dark red      -- high concentration
}

CAT_ORDER = [
    "high_departure", "mod_departure", "low",
    "mod_concentration", "high_concentration",
]

CAT_LABELS = {
    "low":                "Low alteration  |HA| < 0.33",
    "mod_departure":      "Moderate departure  -0.67 <= HA < -0.33",
    "high_departure":     "High departure  HA < -0.67",
    "mod_concentration":  "Moderate concentration  0.33 < HA <= 0.67",
    "high_concentration": "High concentration  HA > 0.67",
}

# Concise labels for CSV cell values (verbose CAT_LABELS above are for plot
# legends, shown once per figure -- repeating the full threshold text in
# every one of 621 catchments x 32 parameters would be needlessly verbose).
CAT_LABELS_SIMPLE = {
    "low":                "low alteration",
    "mod_departure":      "moderate departure",
    "high_departure":     "high departure",
    "mod_concentration":  "moderate concentration",
    "high_concentration": "high concentration",
}


def classify_ha(ha_value: float) -> str:
    """
    Assign a five-category short key to a single HA value.

    Boundary convention: every threshold (0.33, 0.67, -0.33, -0.67)
    belongs to the lower-magnitude-of-alteration category, so this
    partitions the real line with no gaps or overlaps:

        ha < -0.67            -> "high_departure"
        -0.67 <= ha < -0.33   -> "mod_departure"
        -0.33 <= ha <= 0.33   -> "low"
        0.33 < ha <= 0.67     -> "mod_concentration"
        ha > 0.67             -> "high_concentration"
    """
    if ha_value < -0.67:
        return "high_departure"
    elif ha_value < -0.33:
        return "mod_departure"
    elif ha_value <= 0.33:
        return "low"
    elif ha_value <= 0.67:
        return "mod_concentration"
    else:
        return "high_concentration"


def is_departure(ha_value: float) -> bool:
    """True if ha_value falls in moderate or high departure (HA < -0.33)."""
    return ha_value < -0.33
