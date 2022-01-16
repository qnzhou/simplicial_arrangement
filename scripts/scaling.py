#!/usr/bin/env python3

import csv
import os
import pathlib
import re
import subprocess

def run_scaling_benchmark(name):
    cwd = os.getcwd()
    exe = pathlib.Path(cwd) / "../build/simplicial_arrangement_tests"
    assert(exe.exists())
    r = subprocess.check_output([exe, name]).decode('ascii')
    return r

def parse_timing(name, text):
    header_pattern = re.compile("\s*(\d*)\s*planes")
    time_pattern = re.compile("\s*(\d*.\d*)\s*([um]?s)")

    plane_counts = []
    means = []
    std_devs = []

    chunks = text.split("3D {}".format(name))
    for c in chunks:
        lines = c.splitlines()
        r0 = header_pattern.match(lines[0])
        r1 = time_pattern.match(lines[1])
        r2 = time_pattern.match(lines[2])
        if r0 is None or r1 is None or r2 is None: continue

        num_planes = r0.group(1)
        mean = float(r1.group(1))
        std_dev = float(r2.group(1))
        if (r1.group(2) == "ms"): mean *= 1e-3
        elif (r1.group(2) == "us"): mean *= 1e-6
        elif (r1.group(2) == "ns"): mean *= 1e-9
        if (r2.group(2) == "ms"): std_dev *= 1e-3
        elif (r2.group(2) == "us"): std_dev *= 1e-6
        elif (r2.group(2) == "ns"): std_dev *= 1e-9

        print("{} planes: mean={}s  std_dev={}s".format(num_planes, mean,
            std_dev))
        plane_counts.append(num_planes)
        means.append(mean)
        std_devs.append(std_dev)

    return plane_counts, means, std_devs

def export_csv(filename, data):
    num_planes, means, std_devs = data
    with open(filename, 'w') as fout:
        writer = csv.writer(fout)
        writer.writerow(["num_planes", "mean", "std_dev"])
        num_entries = len(num_planes)
        for i in range(num_entries):
            writer.writerow([num_planes[i], means[i], std_devs[i]])

if __name__ == "__main__":
    r = run_scaling_benchmark("arrangement_scaling")
    data = parse_timing("ar", r)
    export_csv("arrangement_scaling.csv", data)

    r = run_scaling_benchmark("material_interface_scaling")
    data = parse_timing("mi", r)
    export_csv("material_interface_scaling.csv", data)

