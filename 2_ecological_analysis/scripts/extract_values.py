#!/usr/bin/env python3
"""
extract_values.py
Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
dynamics in an island plant radiation

Extracts the values of all .asc raster layers at each occurrence point
(Methods S3, "Ecological differentiation analysis"). Outputs a CSV with
Species, x, y and one column per environmental layer. Points and rasters must
use the same projected coordinate system (metres).

Usage:
    python3 scripts/extract_values.py <occurrence_csv> <asc_folder> <output_csv>

Used in the paper on the spatially thinned occurrences (001_spatial_thinning.R):
    python3 scripts/extract_values.py \\
        species_occurrence/spatial_thinning_10m/occurrences_thinned.csv \\
        data/ascii_layers \\
        species_occurrence/spatial_thinning_10m/extracted_data.csv
"""

import sys
import os
import glob
import csv


def read_asc(filepath):
    """Read an ASC raster file and return header info and data array."""
    with open(filepath, 'r') as f:
        header = {}
        for _ in range(6):
            key, val = f.readline().split()[:2]
            try:
                header[key.lower()] = int(val)
            except ValueError:
                header[key.lower()] = float(val)
        data = [[float(x) for x in line.split()] for line in f]
    return header, data


def sample_point(header, data, x, y):
    """Raster value of the cell containing (x, y); None if outside or NODATA."""
    col = int((x - header['xllcorner']) / header['cellsize'])
    row = int((header['yllcorner'] + header['nrows'] * header['cellsize'] - y) / header['cellsize'])
    if row < 0 or row >= header['nrows'] or col < 0 or col >= header['ncols']:
        return None
    if row >= len(data) or col >= len(data[row]):
        return None
    val = data[row][col]
    if abs(val - header.get('nodata_value', -9999)) < 0.001:
        return None
    return val


def main():
    if len(sys.argv) != 4:
        print("Usage: python3 extract_values.py <occurrence_csv> <asc_folder> <output_csv>")
        sys.exit(1)
    species_file, asc_folder, output_file = sys.argv[1:4]

    # Raster layers, named by file name
    asc_files = sorted(glob.glob(os.path.join(asc_folder, '*.asc')))
    if not asc_files:
        sys.exit(f"No .asc files found in {asc_folder}")
    layer_names = [os.path.splitext(os.path.basename(f))[0].replace('_aligned', '') for f in asc_files]
    rasters = [read_asc(f) for f in asc_files]

    # Occurrence points (comma- or tab-delimited; species and coordinate columns detected by name)
    with open(species_file, 'r') as f:
        delimiter = '\t' if '\t' in f.readline() else ','
        f.seek(0)
        reader = csv.DictReader(f, delimiter=delimiter)
        cols = {c.strip().lower(): c for c in reader.fieldnames}
        x_col = next(cols[c] for c in ['x', 'ddlongitude', 'longitude', 'lon', 'long'] if c in cols)
        y_col = next(cols[c] for c in ['y', 'ddlatitude', 'latitude', 'lat'] if c in cols)
        sp_col = next(cols[c] for c in ['species', 'species_name'] if c in cols)
        points = []
        for row in reader:
            try:
                points.append((row[sp_col].strip(), float(row[x_col]), float(row[y_col])))
            except (ValueError, KeyError):
                continue

    # Extract values and write output
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Species', 'x', 'y'] + layer_names)
        for species, x, y in points:
            vals = [sample_point(h, d, x, y) for h, d in rasters]
            writer.writerow([species, x, y] + ['' if v is None else v for v in vals])

    print(f"{len(points)} points x {len(layer_names)} layers written to {output_file}")


if __name__ == "__main__":
    main()
