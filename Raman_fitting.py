# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 00:15:54 2019

@author: leeva

"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fittingV_ver2

"""
Dependency: lmfit, numpy, scipy, matplotlib, pandas
"""


def main():
    while True:
        try:
            data_path = Path(input("Input full path to data directory eg C:\\Path\\to\\Data: \n"))
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
            print("please input a valid extension")

    # # parse files from direcotry into a dictionary
    files = {fle.stem: fle for fle in data_path.iterdir() if fle.suffix == frmt}

    while True:
        try:
            fit_t = int(
                input(
                    "Input the number of peaks to fit \n 1 Peak == 1"
                    "\n 2 Peaks == 2 \n 3 Peaks == 3:\n"
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

    while True:
        try:
            widths = int(
                input(
                    "peak width to use in peak search algo, usually 4-6 works best\n"
                )
            )
            if not isinstance(widths, int):
                raise ValueError
            else:
                break
        except ValueError:
            print("please input a valid width")

    while True:
        try:
            generate_plots = str(input("Generate Plots y/n\n"))
            if generate_plots not in ["y", "n", "yes", "no", "Y", "N"]:
                raise ValueError
            else:
                break
        except ValueError:
            print("please select a valid option")

    def parse_file(file_path):
        infile = open(file_path, "r")
        data = np.loadtxt(infile)
        if frmt == ".dat":
            shift = data[:, 0]
            signal = data[:, 2]
        elif frmt == ".txt":
            shift = data[:, 0]
            signal = data[:, 1]
        else:
            raise Exception(f"Unknown Format: >>{frmt}")
        return shift, signal

    summary = []
    if fit_t == 1:
        for name, fle in files.items():
            print(f"Processing {fle}")
            try:
                shift, signal = parse_file(fle)
            except ValueError:
                raise Exception(f"File:>> {fle} provided cannot be parsed")
            summary_data, fit_report = fittingV_ver2.Multi_fit(
                shift, signal, fit_type, widths
            ).peaks_1()
            if generate_plots:
                plot_dir = data_path / "plots_one_peak"
                if not plot_dir.exists():
                    plot_dir.mkdir(parents=True, exist_ok=True)
                fig = plt.figure()
                plt.scatter(summary_data["raw_x"], summary_data["raw_y"], marker=".", s=1)
                plt.plot(
                    summary_data["raw_x"],
                    summary_data["best_fit"],
                    color="yellow",
                    label="best_fit",
                )
                plt.ylabel("Signal")
                plt.xlabel("Shift cm-3")
                plt.title(f"{name}_single_peak_fit")
                plt.legend(loc="best")
                print(f"Saving Plot {name} and fit_report in {plot_dir}")
                fig.savefig(plot_dir / f"{name}_plot.png")

                plot_data = pd.DataFrame(
                    dict(
                        raw_x=summary_data["raw_x"],
                        raw_y=summary_data["raw_y"],
                        best_fit_x=summary_data["raw_x"],
                        best_fit_y=summary_data["best_fit"],
                    )
                )
                plot_data.to_csv(plot_dir / f"{name}_plot_data.csv")

                with open(plot_dir / f"{name}_fit_result.txt", "w") as fh:
                    fh.write(fit_report)

            sum_data_to_write = pd.DataFrame(
                dict(
                    fwhm=summary_data["fwhm"],
                    peak_pos=summary_data["peak_pos"],
                    area=summary_data["area"],
                ),
                index=[name],
            )
            sum_data_to_write.index.name = "data_set_file_name"
            summary.append(sum_data_to_write)
        pd.concat(summary).to_csv(data_path / "fitting_summary.csv")
    elif fit_t == 2:
        for name, fle in files.items():
            print(f"Processing {fle}")
            try:
                shift, signal = parse_file(fle)
            except ValueError:
                raise Exception(f"File:>> {fle} provided cannot be parsed")
            summary_data, fit_report = fittingV_ver2.Multi_fit(
                shift, signal, fit_type, widths
            ).peaks_2()
            if generate_plots:
                plot_dir = data_path / "plots_two_peaks"
                if not plot_dir.exists():
                    plot_dir.mkdir(parents=True, exist_ok=True)
                fig = plt.figure()
                plt.scatter(summary_data["raw_x"], summary_data["raw_y"], marker=".", s=1)
                plt.plot(
                    summary_data["compound_x"],
                    summary_data["compound_y"],
                    color="yellow",
                    label="compund_fit",
                )
                plt.plot(
                    summary_data["peak1_x"],
                    summary_data["peak1_bestfit"],
                    color="red",
                    label="peak_1",
                )
                plt.plot(
                    summary_data["peak2_x"],
                    summary_data["peak2_bestfit"],
                    color="green",
                    label="peak_2",
                )
                plt.legend(loc="best")
                plt.ylabel("Signal")
                plt.xlabel("Shift cm-3")
                plt.title(f"{name}_double_peak_fit")
                print(f"Saving Plot {name} and fit_report in {plot_dir}")
                fig.savefig(plot_dir / f"{name}_plot.png")

                plot_data_raw = pd.DataFrame(
                    dict(
                        raw_x=summary_data["raw_x"],
                        raw_y=summary_data["raw_y"],
                        best_fit_x=summary_data["compound_x"],
                        best_fit_y=summary_data["compound_y"],
                    )
                )
                plot_data_1 = pd.DataFrame(
                    dict(
                        peak1_x=summary_data["peak1_x"],
                        peak1_bestfit=summary_data["peak1_bestfit"],
                    )
                )
                plot_data_2 = pd.DataFrame(
                    dict(
                        peak2_x=summary_data["peak2_x"],
                        peak2_bestfit=summary_data["peak2_bestfit"],
                    )
                )

                plot_data_raw.to_csv(plot_dir / f"{name}_plot_data.csv")
                plot_data_1.to_csv(plot_dir / f"{name}_peak_1_data.csv")
                plot_data_2.to_csv(plot_dir / f"{name}_peak_2_data.csv")
                with open(plot_dir / f"{name}_fit_result.txt", "w") as fh:
                    fh.write(fit_report)

            sum_data_to_write = pd.DataFrame(
                dict(
                    total_area=summary_data["total_area"],
                    area_1=summary_data["area_1"],
                    area_2=summary_data["area_2"],
                    ratio=summary_data["ratio"],
                    peak1_pos=summary_data["peak1_pos"],
                    peak2_pos=summary_data["peak2_pos"],
                    fwhm_1=summary_data["fwhm_1"],
                    fwhm_2=summary_data["fwhm_2"],
                ),
                index=[name],
            )
            sum_data_to_write.index.name = "data_set_file_name"
            summary.append(sum_data_to_write)
        pd.concat(summary).to_csv(data_path / "fitting_summary_double_peak_fit.csv")
    elif fit_t == 3:
        for name, fle in files.items():
            print(f"Processing {fle}")
            try:
                shift, signal = parse_file(fle)
            except ValueError:
                raise Exception(f"File:>> {fle} provided cannot be parsed")
            summary_data, fit_report = fittingV_ver2.Multi_fit(
                shift, signal, fit_type, widths
            ).peaks_3()
            if generate_plots:
                plot_dir = data_path / "plots_three_peaks"
                if not plot_dir.exists():
                    plot_dir.mkdir(parents=True, exist_ok=True)
                fig = plt.figure()
                plt.scatter(summary_data["raw_x"], summary_data["raw_y"], marker=".", s=1)
                plt.plot(
                    summary_data["compound_x"],
                    summary_data["compound_y"],
                    color="yellow",
                    label="compound_fit",
                )
                plt.plot(
                    summary_data["peak1_x"],
                    summary_data["peak1_bestfit"],
                    color="red",
                    label="peak_1",
                )
                plt.plot(
                    summary_data["peak2_x"],
                    summary_data["peak2_bestfit"],
                    color="green",
                    label="peak_2",
                )
                plt.plot(
                    summary_data["peak3_x"],
                    summary_data["peak3_bestfit"],
                    color="black",
                    label="peak_3",
                )
                plt.legend(loc="best")
                plt.ylabel("Signal")
                plt.xlabel("Shift cm-3")
                plt.title(f"{name}_triple_peak_fit")
                print(f"Saving Plot {name} and fit_report in {plot_dir}")
                fig.savefig(plot_dir / f"{name}_plot.png")

                plot_data_raw = pd.DataFrame(
                    dict(
                        raw_x=summary_data["raw_x"],
                        raw_y=summary_data["raw_y"],
                        best_fit_x=summary_data["compound_x"],
                        best_fit_y=summary_data["compound_y"],
                    )
                )
                plot_data_1 = pd.DataFrame(
                    dict(
                        peak1_x=summary_data["peak1_x"],
                        peak1_bestfit=summary_data["peak1_bestfit"],
                    )
                )
                plot_data_2 = pd.DataFrame(
                    dict(
                        peak2_x=summary_data["peak2_x"],
                        peak2_bestfit=summary_data["peak2_bestfit"],
                    )
                )
                plot_data_3 = pd.DataFrame(
                    dict(
                        peak3_x=summary_data["peak3_x"],
                        peak3_bestfit=summary_data["peak3_bestfit"],
                    )
                )
                plot_data_raw.to_csv(plot_dir / f"{name}_plot_data.csv")
                plot_data_1.to_csv(plot_dir / f"{name}_peak_1_data.csv")
                plot_data_2.to_csv(plot_dir / f"{name}_peak_2_data.csv")
                plot_data_3.to_csv(plot_dir / f"{name}_peak_3_data.csv")
                with open(plot_dir / f"{name}_fit_result.txt", "w") as fh:
                    fh.write(fit_report)

            sum_data_to_write = pd.DataFrame(
                dict(
                    total_area=summary_data["total_area"],
                    area_1=summary_data["area_1"],
                    area_2=summary_data["area_2"],
                    area_3=summary_data["area_3"],
                    peak1_pos=summary_data["peak1_pos"],
                    peak2_pos=summary_data["peak2_pos"],
                    peak3_pos=summary_data["peak3_pos"],
                    fwhm_1=summary_data["fwhm_1"],
                    fwhm_2=summary_data["fwhm_2"],
                    fwhm_3=summary_data["fwhm_3"],
                ),
                index=[name],
            )
            sum_data_to_write.index.name = "data_set_file_name"
            summary.append(sum_data_to_write)
        pd.concat(summary).to_csv(data_path / "fitting_summary_triple_peak_fit.csv")


if __name__ == "__main__":
    print(
        "Raman Spectra Fitting Program v0.0.1\n"
        "Based on the Lmfit fitting module\n"
        "https://lmfit.github.io/lmfit-py/\n\n"
        "written by: LTrindade\n\n"
    )
    while True:
        main()
        while True:
            try:
                rerun = str(input("All Files Processed Re-run? y/n\n")).lower()
                if rerun not in ["y", "n"]:
                    raise ValueError
                else:
                    rer = True if rerun == "y" else False
                    break
            except ValueError:
                print("please select a valid option")
        if not rer:
            break
