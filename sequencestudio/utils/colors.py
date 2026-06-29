"""Residue coloring schemes for sequence and alignment display."""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Nucleotide — ClustalX-style
# ---------------------------------------------------------------------------
NUCLEOTIDE_COLORS: dict[str, str] = {
    "A": "#64F73F",   # green
    "C": "#FFB340",   # orange
    "G": "#EB413C",   # red
    "T": "#3C88EE",   # blue
    "U": "#3C88EE",   # blue (RNA)
    "N": "#CCCCCC",   # grey
    "-": "#FFFFFF",
    ".": "#FFFFFF",
}

# ---------------------------------------------------------------------------
# Amino acid — ClustalX
# ---------------------------------------------------------------------------
CLUSTALX_AA_COLORS: dict[str, str] = {
    # Hydrophobic
    "A": "#80A0F0", "V": "#80A0F0", "I": "#80A0F0", "L": "#80A0F0",
    "M": "#80A0F0", "F": "#80A0F0", "W": "#80A0F0", "C": "#80A0F0",
    # Positive
    "K": "#F01505", "R": "#F01505", "H": "#F01505",
    # Negative
    "D": "#C048C0", "E": "#C048C0",
    # Polar
    "S": "#15C015", "T": "#15C015", "N": "#15C015", "Q": "#15C015",
    # Aromatic
    "Y": "#15A4A4",
    # Glycine
    "G": "#F09048",
    # Proline
    "P": "#C0C000",
    "-": "#FFFFFF", "*": "#FFFFFF",
}

# ---------------------------------------------------------------------------
# Amino acid — Zappo
# ---------------------------------------------------------------------------
ZAPPO_AA_COLORS: dict[str, str] = {
    "I": "#FFAAAA", "L": "#FFAAAA", "V": "#FFAAAA", "A": "#FFAAAA", "M": "#FFAAAA",
    "F": "#FF00AA", "W": "#FF00AA", "Y": "#FF00AA",
    "K": "#6600FF", "R": "#6600FF", "H": "#6600FF",
    "D": "#FF0000", "E": "#FF0000",
    "S": "#00FF00", "T": "#00FF00", "N": "#00FF00", "Q": "#00FF00",
    "P": "#FFFF00", "G": "#FFFF00", "C": "#FFFF00",
    "-": "#FFFFFF", "*": "#FFFFFF",
}

# ---------------------------------------------------------------------------
# Amino acid — Taylor
# ---------------------------------------------------------------------------
TAYLOR_AA_COLORS: dict[str, str] = {
    "A": "#CCFF00", "R": "#0000FF", "N": "#CC00FF", "D": "#FF0000",
    "C": "#FFFF00", "Q": "#FF00CC", "E": "#FF0066", "G": "#FF9900",
    "H": "#0066FF", "I": "#66FF00", "L": "#33FF00", "K": "#6600FF",
    "M": "#00FF00", "F": "#00FF66", "P": "#FFAA00", "S": "#FF3300",
    "T": "#FF6600", "W": "#00CCFF", "Y": "#00FFCC", "V": "#99FF00",
    "-": "#FFFFFF", "*": "#FFFFFF",
}

# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------
COLOR_SCHEMES: dict[str, dict[str, str]] = {
    "Nucleotide (ClustalX)": NUCLEOTIDE_COLORS,
    "Protein (ClustalX)":    CLUSTALX_AA_COLORS,
    "Protein (Zappo)":       ZAPPO_AA_COLORS,
    "Protein (Taylor)":      TAYLOR_AA_COLORS,
}


def get_text_color(hex_bg: str) -> str:
    """Return '#000000' or '#FFFFFF' for legible text on *hex_bg*."""
    h = hex_bg.lstrip("#")
    if len(h) != 6:
        return "#000000"
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    lum = (0.299 * r + 0.587 * g + 0.114 * b) / 255
    return "#000000" if lum > 0.45 else "#FFFFFF"
