"""Rasterize the gallery's simple SVG with an existing Pillow installation.

Usage: python rasterize_contact_sheet.py path/to/contact-sheet.svg
The PNG is placed alongside its SVG source; no Maya installation is required.
"""

import argparse
from pathlib import Path
import xml.etree.ElementTree as ET

from PIL import Image, ImageDraw, ImageFont


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("svg", type=Path)
    args = parser.parse_args()
    root = ET.parse(args.svg).getroot()
    image = Image.new("RGB", (int(root.attrib["width"]), int(root.attrib["height"])))
    draw = ImageDraw.Draw(image)
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", 11)
    except OSError:
        font = ImageFont.load_default()
    for element in root:
        tag = element.tag.rsplit("}", 1)[-1]
        attrs = element.attrib
        if tag == "rect":
            draw.rectangle((0, 0, int(attrs["width"]), int(attrs["height"])), fill=attrs["fill"])
        elif tag == "text":
            draw.text(
                (float(attrs["x"]), float(attrs["y"]) - 11),
                element.text,
                fill=attrs["fill"],
                font=font,
            )
        elif tag == "polyline":
            points = [tuple(map(float, point.split(","))) for point in attrs["points"].split()]
            draw.line(points, fill=attrs["stroke"], width=int(attrs["stroke-width"]))
        else:
            raise ValueError("Unsupported SVG element: " + tag)
    image.save(args.svg.with_suffix(".png"))


if __name__ == "__main__":
    main()
