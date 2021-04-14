# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 23:28:51 2019

@author: leeva


"""
import pandas as pd
import scipy.signal
from lmfit.models import LinearModel
from lmfit.models import LorentzianModel
from lmfit.models import VoigtModel
from scipy.integrate import simps


class Multi_fit:
    def __init__(self, shift, signal, fit_t, peak_width=5):

        self.x = shift
        self.y = signal
        self.type = fit_t
        self.peak_width = peak_width
        self.data = pd.DataFrame(dict(Shift=shift, Signal=signal))

    def find_peaks(self):
        peak_indexes = scipy.signal.find_peaks(
            self.data["Signal"].to_numpy(), prominence=1, width=self.peak_width
        )[0]
        peak_df = self.data.iloc[peak_indexes]
        return peak_df.sort_values("Signal", ascending=False)

    def peaks_1(self):
        # Initialise Models
        if self.type == 1:
            peak_mod = LorentzianModel(prefix="v1_")
        elif self.type == 2:
            peak_mod = VoigtModel(prefix="v1_")
        else:
            raise Exception(f"{self.type} is not a selectable option\n")
        background = LinearModel(prefix="back_")  # Model to take into account background noise

        # Intial Fit with guessed parameters
        pars = peak_mod.guess(self.data["Signal"], x=self.data["Shift"])
        pars.update(background.make_params())  # add the background model
        model = peak_mod + background

        # Running linear regression fit on the data
        fit = model.fit(self.data["Signal"].values, pars, x=self.data["Shift"].values)
        fwhm_1 = fit.params["v1_fwhm"].value
        peak1_pos = fit.params["v1_center"].value
        area_1 = simps(fit.best_fit, dx=0.5)

        return (
            dict(
                raw_x=self.data["Shift"].values,
                raw_y=self.data["Signal"].values,
                best_fit=fit.best_fit,
                peak_pos=peak1_pos,
                area=area_1,
                fwhm=fwhm_1,
            ),
            fit.fit_report(min_correl=0.5),
        )

    def peaks_2(self):

        if self.type == 1:
            lor_fit1 = LorentzianModel(prefix="v1_")
            lor_fit2 = LorentzianModel(prefix="v2_")
        elif self.type == 2:
            lor_fit1 = VoigtModel(prefix="v1_")
            lor_fit2 = VoigtModel(prefix="v2_")
        else:
            raise Exception(f"{self.type} is not a selectable option\n")
        background = LinearModel(prefix="back_")

        # _find main peak_ #
        peaks = self.find_peaks()
        high_peak = peaks.iloc[0]["Shift"]
        second_peak = peaks.iloc[1]["Shift"]

        fit_range = 100  # Default to fit single peak in compound

        # _comopund fit_ #
        pars = lor_fit1.guess(self.data["Signal"].values, x=self.data["Shift"].values)
        pars["v1_center"].set(high_peak, min=high_peak - 2, max=high_peak + 2)
        pars.update(lor_fit2.make_params())
        pars["v2_center"].set(second_peak, min=second_peak - 2, max=second_peak + 2)
        pars.update(background.make_params())
        compound = lor_fit1 + lor_fit2 + background
        compound_fit = compound.fit(self.data["Signal"].values, pars, x=self.data["Shift"].values)
        total_area = simps(compound_fit.best_fit, dx=0.5)
        compound_x = self.data["Shift"]
        compound_y = compound_fit.best_fit

        # _single peak fits_ #
        fit_p1 = self.data[
            (
                (self.data["Shift"] > compound_fit.params["v1_center"].value - (fit_range / 2))
                & (self.data["Shift"] < compound_fit.params["v1_center"].value + (fit_range / 2))
            )
        ]
        fit_p2 = self.data[
            (
                (self.data["Shift"] > compound_fit.params["v2_center"].value - (fit_range / 2))
                & (self.data["Shift"] < compound_fit.params["v2_center"].value + (fit_range / 2))
            )
        ]

        fit_X_p1 = fit_p1["Shift"].values
        fit_Y_p1 = fit_p1["Signal"].values

        fit_X_p2 = fit_p2["Shift"].values
        fit_Y_p2 = fit_p2["Signal"].values

        peak1 = lor_fit1.guess(fit_Y_p1, x=fit_X_p1)
        peak1.update(background.make_params())
        peak2 = lor_fit2.guess(fit_Y_p2, x=fit_X_p2)
        peak2.update(background.make_params())

        peak1["v1_center"].set(
            compound_fit.params["v1_center"].value,
            min=compound_fit.params["v1_center"].value - 5,
            max=compound_fit.params["v1_center"].value + 5,
        )
        peak2["v2_center"].set(
            compound_fit.params["v2_center"].value,
            min=compound_fit.params["v2_center"].value - 5,
            max=compound_fit.params["v2_center"].value + 5,
        )
        peak2["v2_amplitude"].set(value=compound_fit.params["v2_amplitude"].value, vary=False)
        peak2["v2_sigma"].set(value=compound_fit.params["v2_sigma"].value, vary=False)
        mod_peak1 = lor_fit1 + background
        mod_peak2 = lor_fit2 + background
        peak1_fit = mod_peak1.fit(fit_Y_p1, peak1, x=fit_X_p1)
        peak2_fit = mod_peak2.fit(fit_Y_p2, peak2, x=fit_X_p2)

        area_1 = simps(peak1_fit.best_fit, dx=0.5)
        area_2 = simps(peak2_fit.best_fit, dx=0.5)

        ratio = area_1 / area_2

        peak1_pos = compound_fit.params["v1_center"].value
        peak2_pos = peak2_fit.params["v2_center"].value

        fwhm_1 = peak1_fit.params["v1_fwhm"].value
        fwhm_2 = peak2_fit.params["v2_fwhm"].value

        return (
            dict(
                raw_x=self.data["Shift"].values,
                raw_y=self.data["Signal"].values,
                compound_x=compound_x,
                compound_y=compound_y,
                peak1_x=fit_X_p1,
                peak1_bestfit=peak1_fit.best_fit,
                peak2_x=fit_X_p2,
                peak2_bestfit=peak2_fit.best_fit,
                total_area=total_area,
                area_1=area_1,
                area_2=area_2,
                ratio=ratio,
                peak1_pos=peak1_pos,
                peak2_pos=peak2_pos,
                fwhm_1=fwhm_1,
                fwhm_2=fwhm_2,
            ),
            compound_fit.fit_report(min_correl=0.5),
        )

    def peaks_3(self):

        """
        Add output functionality

        """

        if self.type == 1:
            peak_mod = LorentzianModel(prefix="v1_")
            peak_mod2 = LorentzianModel(prefix="v2_")
            peak_mod3 = LorentzianModel(prefix="v3_")
        elif self.type == 2:
            peak_mod = VoigtModel(prefix="v1_")
            peak_mod2 = VoigtModel(prefix="v2_")
            peak_mod3 = VoigtModel(prefix="v3_")
        else:
            raise Exception(f"{self.type} is not a selectable option\n")

        background = LinearModel(prefix="back_")  # Model to take into account background noise
        fit_range = 100
        peaks = self.find_peaks()

        fit_p1 = self.data[
            (
                (self.data["Shift"] > peaks.iloc[0]["Shift"] - (fit_range / 2))
                & (self.data["Shift"] < peaks.iloc[0]["Shift"] + (fit_range / 2))
            )
        ]
        fit_X_p1 = fit_p1["Shift"].values
        fit_Y_p1 = fit_p1["Signal"].values

        fit_p2 = self.data[
            (
                (self.data["Shift"] > peaks.iloc[1]["Shift"] - (fit_range / 2))
                & (self.data["Shift"] < peaks.iloc[1]["Shift"] + (fit_range / 2))
            )
        ]
        fit_X_p2 = fit_p2["Shift"].values
        fit_Y_p2 = fit_p2["Signal"].values

        fit_p3 = self.data[
            (
                (self.data["Shift"] > peaks.iloc[2]["Shift"] - (fit_range / 2))
                & (self.data["Shift"] < peaks.iloc[2]["Shift"] + (fit_range / 2))
            )
        ]
        fit_X_p3 = fit_p3["Shift"].values
        fit_Y_p3 = fit_p3["Signal"].values

        # Single peak fits
        pars_1 = pars = peak_mod.guess(fit_Y_p1, x=fit_X_p1)
        pars_1.update(background.make_params())
        mod_peak1 = peak_mod + background
        pars_1["v1_center"].set(
            peaks.iloc[0]["Shift"], min=peaks.iloc[0]["Shift"] - 2, max=peaks.iloc[0]["Shift"] + 2
        )
        out_1 = mod_peak1.fit(fit_Y_p1, pars_1, x=fit_X_p1)

        pars_2 = peak_mod2.guess(fit_Y_p2, x=fit_X_p2)
        pars_2.update(background.make_params())
        mod_peak2 = peak_mod2 + background
        pars_2["v2_center"].set(
            peaks.iloc[1]["Shift"], min=peaks.iloc[1]["Shift"] - 2, max=peaks.iloc[1]["Shift"] + 2
        )
        out_2 = mod_peak2.fit(fit_Y_p2, pars_2, x=fit_X_p2)

        pars_3 = peak_mod3.guess(fit_Y_p3, x=fit_X_p3)
        pars_3.update(background.make_params())
        mod_peak3 = peak_mod3 + background
        pars_3["v3_center"].set(
            peaks.iloc[2]["Shift"], min=peaks.iloc[2]["Shift"] - 2, max=peaks.iloc[2]["Shift"] + 2
        )
        out_3 = mod_peak3.fit(fit_Y_p3, pars_3, x=fit_X_p3)

        fit_area1 = simps(out_1.best_fit, dx=0.5)
        fit_area2 = simps(out_2.best_fit, dx=0.5)
        fit_area3 = simps(out_3.best_fit, dx=0.5)

        fit_xy_1 = (fit_X_p1, out_1.best_fit)
        fit_xy_2 = (fit_X_p2, out_2.best_fit)
        fit_xy_3 = (fit_X_p3, out_3.best_fit)
        fwhm_1 = out_1.params["v1_fwhm"].value
        fwhm_2 = out_2.params["v2_fwhm"].value
        fwhm_3 = out_3.params["v3_fwhm"].value

        peak_1_fit = out_1.params["v1_center"].value
        peak_2_fit = out_2.params["v2_center"].value
        peak_3_fit = out_3.params["v3_center"].value

        # compound fit
        pars["v1_center"].set(
            peaks.iloc[0]["Shift"], min=peaks.iloc[0]["Shift"] - 2, max=peaks.iloc[0]["Shift"] + 2
        )
        pars.update(peak_mod2.make_params())
        pars["v2_center"].set(
            peaks.iloc[1]["Shift"], min=peaks.iloc[1]["Shift"] - 2, max=peaks.iloc[1]["Shift"] + 2
        )
        pars.update(peak_mod3.make_params())
        pars["v3_center"].set(
            peaks.iloc[2]["Shift"], min=peaks.iloc[2]["Shift"] - 2, max=peaks.iloc[2]["Shift"] + 2
        )
        pars.update(background.make_params())
        compound = peak_mod + peak_mod2 + peak_mod3 + background
        compound_fit = compound.fit(self.data["Signal"].values, pars, x=self.data["Shift"].values)

        total_area = simps(compound_fit.best_fit, dx=0.5)
        compound_x = self.data["Shift"]
        compound_y = compound_fit.best_fit

        return (
            dict(
                raw_x=self.data["Shift"].values,
                raw_y=self.data["Signal"].values,
                compound_x=compound_x,
                compound_y=compound_y,
                peak1_x=fit_xy_1[0],
                peak1_bestfit=fit_xy_1[1],
                peak2_x=fit_xy_2[0],
                peak2_bestfit=fit_xy_2[1],
                peak3_x=fit_xy_3[0],
                peak3_bestfit=fit_xy_3[1],
                total_area=total_area,
                area_1=fit_area1,
                area_2=fit_area2,
                area_3=fit_area3,
                peak1_pos=peak_1_fit,
                peak2_pos=peak_2_fit,
                peak3_pos=peak_3_fit,
                fwhm_1=fwhm_1,
                fwhm_2=fwhm_2,
                fwhm_3=fwhm_3,
            ),
            compound_fit.fit_report(min_correl=0.5),
        )
