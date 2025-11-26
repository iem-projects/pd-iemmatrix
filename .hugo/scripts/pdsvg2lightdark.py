#!/usr/bin/env python3

# convert an SVG (as created by Pd) to another SVG
# that can be used for both light (default) and dark modes
# this should kind of work if now Pd color-schemes were used

# © 2025, IOhannes m zmölnig, Institute of Electronic Music and Acoustics, KUG

import re

colpattern = re.compile(r"#[0-9a-fA-F]{6}\b")
hexcolpattern = re.compile(r"[0-9a-fA-F]{2}")

import logging

logging.basicConfig()
log = logging.getLogger(__name__)


def parseArgs():
    import argparse

    parser = argparse.ArgumentParser(
        description="""replace hardcoded colorvalues in SVGs (as exported by Pd)
        with variables that do something useful in light & dark modes"""
    )

    parser.add_argument(
        "--dark-background",
        default="#111111",
        help="reference colour for background in dark-mode (DEFAULT: %(default)r)",
    )
    parser.add_argument(
        "--dark-foreground",
        default="#BCC0C3",
        help="reference colour for foreground in dark-mode (DEFAULT: %(default)r)",
    )

    parser.add_argument(
        "infile",
        help="input SVG file",
    )
    parser.add_argument(
        "outfile",
        nargs="?",
        help="output SVG file (leave empty to output to stdout)",
    )

    args = parser.parse_args()
    return args


def svg2soup(filename):
    from bs4 import BeautifulSoup

    with open(filename) as f:
        return BeautifulSoup(f.read(), features="xml")


def getcolours(soup):
    colors = set()
    for element in soup.find("svg"):
        try:
            style = element["style"]
        except TypeError:
            continue
        colors.update(colpattern.findall(style))
    return colors


def hexcol2tuple(hexcol="#000000"):
    """convert a hexadecimal colour code to an (R, G, B) tuple"""
    return [int(_, 16) for _ in hexcolpattern.findall(hexcol[1:])]


def tuple2hexcol(col=(0x00, 0x00, 0x00)):
    hexcol = [f"{int(c):02X}" for c in col]
    return "#" + "".join(hexcol)


def colorchanmap(chan, src=(0x00, 0xFF), dst=(0xFF, 0x00)):
    """map a single color-channel <chan> in the <src> range to the <dst> range"""

    c1 = (chan - src[0]) / (src[1] - src[0])
    r = int(dst[0] + c1 * (dst[1] - dst[0]))
    return r


def colormap(
    col,
    src=((0xFF, 0xFF, 0xFF), (0x00, 0x00, 0x00)),
    dst=((0x00, 0x00, 0x00), (0xFF, 0xFF, 0xFF)),
):
    """convert the <col> RGB-tuple from light to black"""
    log.error(f"mapping {col} from {src} to {dst}")
    return [
        colorchanmap(c, (s0, s1), (d0, d1))
        for c, s0, s1, d0, d1 in zip(col, src[0], src[1], dst[0], dst[1])
    ]


def fix_styles_in_soup(soup, colormap: dict):
    """given a name:colour mapping in <colormap>, this replaces all <colour> occurences with their <name> in all 'style' attributes in soup's SVG"""

    #    rep = {"condition1": "", "condition2": "text"} # define desired replacements here
    #    # use these three lines to do the replacement
    #    rep = dict((re.escape(k), v) for k, v in rep.items())
    #    pattern = re.compile("|".join(rep.keys()))
    #    text = pattern.sub(lambda m: rep[re.escape(m.group(0))], text)

    replacedict = dict((re.escape(k), f"var({v})") for v, k in colormap.items())
    pattern = re.compile("|".join(replacedict.keys()))

    for element in soup.find("svg"):
        try:
            style = element["style"]
        except TypeError:
            continue

        element["style"] = pattern.sub(
            lambda m: replacedict[re.escape(m.group(0))], style
        )

    return soup


def make_style(lightcolours, darkcolours):
    """creates an SVG <style> element from the light and dark colours dicts"""

    # <style>:root { --pd-foreground: #000000; --pd-color1: #404040;}@media (prefers-color-scheme: dark) { :root { --pd-foreground: #BCC0C3; --pd-color1: #BFBFBF; } }</style>
    def dict2css(d):
        return " ".join((f"{k}: {v};" for k, v in d.items()))

    style = ":root { %s } @media (prefers-color-scheme: dark) { :root { %s } }" % (
        dict2css(lightcolours),
        dict2css(darkcolours),
    )
    return style


def insert_style(soup, style: str):
    style_tag = soup.new_tag("style", string=style)
    for element in soup.find("svg"):
        element.insert_before(style_tag)
        break
    return soup


def _main():
    args = parseArgs()
    lightcolours = {"foreground": "#000000"}

    light = [
        (0x00, 0x00, 0x00),
        (0xFF, 0xFF, 0xFF),
    ]
    dark = [
        hexcol2tuple(args.dark_foreground),
        hexcol2tuple(args.dark_background),
    ]

    soup = svg2soup(args.infile)
    colours = getcolours(soup)
    log.error(colours)
    for c in lightcolours.values():
        colours.discard(c)
    for idx, c in enumerate(sorted(colours)):
        lightcolours[f"col{idx}"] = c

    lcols = {k: hexcol2tuple(v) for k, v in lightcolours.items()}
    dcols = {k: colormap(v, light, dark) for k, v in lcols.items()}

    if "col9" in lcols:
        x = lcols["col9"]
        print(f"{x} -> {colormap(x, light, dark)}")

    log.error(lcols)
    log.error(dcols)

    darkcolours = {k: tuple2hexcol(v) for k, v in dcols.items()}

    lightcolours = {f"--pd-{k}": v for k, v in lightcolours.items()}
    darkcolours = {f"--pd-{k}": v for k, v in darkcolours.items()}

    style = make_style(lightcolours, darkcolours)

    soup = fix_styles_in_soup(soup, lightcolours)
    soup = insert_style(soup, style)

    if args.outfile is None:
        print(soup.prettify())
    else:
        with open(args.outfile, "w") as f:
            f.write(soup.prettify())


if __name__ == "__main__":
    _main()
