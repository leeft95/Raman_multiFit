# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 00:15:54 2019

@author: leeva

"""
import os
from pathlib import Path

import numpy as np

import fittingV_ver2

"""
Dependency: lmfit,numpy,scipy,matplotlib
"""


class Data:
    def __init__(self, shift, signal, outfile, name, filename):
        self.shift = shift
        self.signal = signal
        self.outfile = outfile
        self.name = name
        self.filename = filename


files = []
name = []
rama = []

n = 10000
cond = 0

while True:
    try:
        data_path = Path(input("Input data path: \n"))
        if not data_path.is_dir():
            raise ValueError
        else:
            break
    except ValueError:
        print("Please input valid path")

while True:
    try:
        frmt = input('Input file format (form ".txt" or ".dat"): \n')
        if frmt not in [".dat", ".txt"]:
            raise ValueError
        else:
            break
    except ValueError:
        print("please input a valid file path")


while True:
    try:
        fit_t = int(
            input(
                "Input the number of peaks to fit \n 1 Peak == 1 \n 2 Peaks == 2 \n 3 Peaks == 3:\n"
            )
        )
        if fit_t > 3 or fit_t <= 0:
            raise ValueError
        else:
            break
    except ValueError:
        print("please select a valid option")

while True:
    try:
        fit_type = int(
            input(
                "Input the fit type:\n "
                "1 = Lorentzian + LinearBackground "
                "\n 2 = Voigt + LinearBackground\n "
            )
        )
        if fit_type > 2 or fit_type <= 0:
            raise ValueError
        else:
            break
    except ValueError:
        print("please select a valid option")

# data_path = Path("/mnt/c/Users/leeva/Desktop/Raman_Program/test")
# frmt = ".txt"
# fit_t = 3
# fit_type = 1
widths = 5

for fle in data_path.iterdir():
    if fle.suffix.endswith(frmt):
        files.append(fle)
        name.append(fle.name)

for i in range(len(files)):
    infile = open(files[i], "r")
    outfile = "new_" + os.path.splitext(name[i])[0] + ".png"
    line = infile.readline()
    data = np.loadtxt(infile)
    if frmt == ".dat":
        shift = data[:, 0]
        signal = data[:, 2]
    elif frmt == ".txt":
        shift = data[:, 0]
        signal = data[:, 1]
    rama.append(fittingV_ver2.Multi_fit(shift, signal, fit_type, widths))

"""
Change the sped == 'x' to where x is your own keyword for your file ie
its central measurement wavenumber
"""

if fit_t == 1:
    for to_fit in rama:
        summary_data, fit_report = to_fit.peaks_1()
        print(summary_data)
        print(fit_report)
elif fit_t == 2:
    for to_fit in rama:
        summary_data, fit_report = to_fit.peaks_2()
        print(summary_data)
        print(fit_report)
elif fit_t == 3:
    for to_fit in rama:
        summary_data, fit_report = to_fit.peaks_3()
        print(summary_data)
        print(fit_report)
