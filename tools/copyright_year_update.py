#!/usr/bin/env python3
"""
Replace the copyright year in the copyright string on all matching files
"""


import sys
import re
from datetime import datetime as dt
import pathlib


def replace_copyright_year(path, old_copyright_year, new_copyright_year):
    with open(path, "r") as fp:
        content = fp.read()

    if re.search(old_copyright_year, content):
        print("Copyright line found")
        corrected_content = re.sub(old_copyright_year, new_copyright_year, content)

        with open(path, "w") as fp:
            fp.write(corrected_content)


def main():
    old_copyright = re.escape("British Crown Copyright ") + ".... Met Office"

    new_year = dt.now().year
    new_copyright = "British Crown Copyright {new_year:d} Met Office".format(
        new_year=new_year
    )
    print("Old Copyright regex: ", old_copyright)
    print("New Copyright will be: ", new_copyright)

    # walk all files starting from current directory
    root = pathlib.Path()

    files = (
        [fl for fl in root.glob("**/*.txt") if fl.is_file()]
        + [fl for fl in root.glob("**/*.md") if fl.is_file()]
        + [fl for fl in root.glob("**/*.yaml") if fl.is_file()]
        + [fl for fl in root.glob("**/*.cc") if fl.is_file()]
        + [fl for fl in root.glob("**/*.h") if fl.is_file()]
    )

    for file_path in files:
        print(f"Replacing copyright year in: {file_path}")
        replace_copyright_year(file_path, old_copyright, new_copyright)


if __name__ == "__main__":
    main()
