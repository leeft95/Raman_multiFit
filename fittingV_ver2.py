# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 23:28:51 2019

@author: leeva


ADD ERROR HANDLING

WHILE 
    TRY
    EXCEPT
    
change the peak find code to a function

change the fitting to a function call too...maybe



"""

import scipy
import scipy.signal
import scipy.interpolate as sci
import matplotlib.pyplot as plt
from scipy.integrate import simps
import numpy as np
from lmfit.models import (GaussianModel,LorentzianModel,PseudoVoigtModel,
VoigtModel)


class Multi_fit():
    
    def __init__(self,shift,signal,title,peaks,peak_width,fit_t):
        
        self.x = shift
        self.y = signal
        self.title = title
        self.num_peaks = peaks
        self.width = peak_width
        self.type = fit_t
    
    def peaks_2(self):
        
        if self.type == 1:
            print('l')
            lor_fit1 = LorentzianModel(prefix = 'v1_')
            lor_fit2 = LorentzianModel(prefix = 'v2_')
            
        elif self.type == 2:
            print('v')
            lor_fit1 = VoigtModel(prefix = 'v1_')
            lor_fit2 = VoigtModel(prefix = 'v2_')
        elif self.type > 3:
            print(str(self.type) + 'is not a selectable option\n')
            
            return 

        
        widths = np.arange(1,self.width,0.01)
        
        print(self.title)
        
        ####_find main peak_####
        
        smooth_indexes_scipy = scipy.signal.find_peaks_cwt(self.y, widths)
        peak_indexes_xs_ys = np.asarray([list(a) for a in list(zip(self.x[smooth_indexes_scipy], self.y[smooth_indexes_scipy]))])
        high_peak = sorted(peak_indexes_xs_ys, key=lambda pair: pair[1])[-self.num_peaks:]
        print('\n' + str(high_peak[-1][0]) +  '\n')
        
        ####_normalisation lists_####
        
        fit_X = []
        fit_Y = []
        
        fit_XX = []
        fit_YY = []
        
        fit_X_p1 = []
        fit_Y_p1 = []
        
        fit_X_p2 = []
        fit_Y_p2 = []
        
        ####_cut down peak search area_####
        
        fit_range = 100
        
        for i in range(len(self.x)):
            if self.x[i] > (high_peak[-1][0]-fit_range) and self.x[i] < (high_peak[-1][0]+fit_range):
                fit_X.append(self.x[i])
                fit_Y.append(self.y[i])
        
        
        fit_X = np.asarray(fit_X)
        fit_Y = np.asarray(fit_Y)
        
        ####_find second peak in cut down range and normalise data_####
        
        smooth_indexes_scipy = scipy.signal.find_peaks_cwt(fit_Y, widths)        
        peak_indexes_xs_ys = np.asarray([list(a) for a in list(zip(fit_X[smooth_indexes_scipy], fit_Y[smooth_indexes_scipy]))])
        high_peak2 = sorted(peak_indexes_xs_ys, key=lambda pair: pair[1])[-self.num_peaks:]
        print(high_peak2)
    
        for i in range(len(fit_X)):
            if fit_X[i] > high_peak[-1][0]-fit_range and fit_X[i] < high_peak[-1][0]+fit_range:
                x  = (fit_Y[i] - min(fit_Y))/(max(fit_Y) - min(fit_Y))
                fit_XX.append(fit_X[i])
                fit_YY.append(x)
                
        fit_XX = np.asarray(fit_XX)
        fit_YY = np.asarray(fit_YY)
        
        
        ####_comopund fit_####
        

        
        pars = lor_fit1.guess(fit_YY,x = fit_XX)
        
        pars['v1_center'].set(high_peak2[-1][0], min=high_peak[-1][0]-2, max=high_peak2[-1][0]+2)
    
        pars.update(lor_fit2.make_params())
        
        pars['v2_center'].set(high_peak2[0][0], min=high_peak2[0][0]-5, max=high_peak2[0][0]+5)
        
        compound = lor_fit1 + lor_fit2
        
        compound_fit = compound.fit(fit_YY,pars,x = fit_XX)
        
        total_area = simps(compound_fit.best_fit,dx =0.5)
        
        compound_x  = fit_XX
        compound_y  = compound_fit.best_fit
        
        plt.plot(compound_x,fit_YY,label = 'best fit')
        plt.plot(compound_x,compound_fit.best_fit,label = 'best fit')
        

        ####_single peak fits_####
        
        for i in range(len(fit_XX)):
            if fit_XX[i] >  compound_fit.params['v1_center'].value-(fit_range/2) and fit_XX[i] <  compound_fit.params['v1_center'].value+(fit_range/2):
                fit_X_p1.append(fit_XX[i])
                fit_Y_p1.append(fit_YY[i])
            if fit_XX[i] >  compound_fit.params['v2_center'].value-(fit_range/2) and fit_XX[i] <  compound_fit.params['v2_center'].value+(fit_range/2):
                fit_X_p2.append(fit_XX[i])
                fit_Y_p2.append(fit_YY[i])
                
        fit_X_p1 = np.asarray(fit_X_p1)
        fit_Y_p1 = np.asarray(fit_Y_p1)
        
        fit_X_p2 = np.asarray(fit_X_p2)
        fit_Y_p2 = np.asarray(fit_Y_p2)
        
        peak1 = lor_fit1.guess(fit_Y_p1,x = fit_X_p1)
        peak2 = lor_fit2.guess(fit_Y_p2,x = fit_X_p2)
        
        peak1['v1_center'].set(compound_fit.params['v1_center'].value, min=compound_fit.params['v1_center'].value-5, max=compound_fit.params['v1_center'].value+5)  
        
        peak2['v2_center'].set(compound_fit.params['v2_center'].value, min=compound_fit.params['v2_center'].value-5, max=compound_fit.params['v2_center'].value+5)
        peak2['v2_amplitude'].set(value =compound_fit.params['v2_amplitude'].value,vary = False)
        peak2['v2_sigma'].set(value = compound_fit.params['v2_sigma'].value, vary = False)
        
        peak1_fit = lor_fit1.fit(fit_Y_p1,peak1,x = fit_X_p1)
        peak2_fit = lor_fit2.fit(fit_Y_p2,peak2,x = fit_X_p2)
    
        area_1 = simps(peak1_fit.best_fit, dx=0.5)
        area_2 = simps(peak2_fit.best_fit, dx=0.5)
        
        plt.plot(fit_X_p1,peak1_fit.best_fit)
        plt.plot(fit_X_p2,peak2_fit.best_fit)
        plt.title(str(self.title))
        plt.savefig(str(self.title)+'.png', transparent=True)
        plt.close()
        ratio = area_1/area_2
        
        peak1_x = fit_X_p1
        peak2_x = fit_X_p2
        
        peak1_pos = compound_fit.params['v1_center'].value
        peak2_pos = peak2_fit.params['v2_center'].value

        fwhm_1 = peak1_fit.params['v1_fwhm'].value
        fwhm_2 = peak2_fit.params['v2_fwhm'].value

        print(len(peak2_fit.best_fit))
        print(len(peak1_x))


        return compound_x,compound_y,peak1_x,peak1_fit.best_fit,peak2_x,peak2_fit.best_fit,total_area,\
               area_1,area_2,ratio,peak1_pos,peak2_pos,fwhm_1,fwhm_2
            
           
    def peaks_3(self): 
        
        """
        Add output functionality
        
        """
        
        if self.type == 1:
            vigot = LorentzianModel(prefix='v1_')
            vigot2 = LorentzianModel(prefix='v2_')
            vigot3 = LorentzianModel(prefix='v3_')
        elif self.type == 2:
            vigot = VoigtModel(prefix='v1_')
            vigot2 = VoigtModel(prefix='v2_')
            vigot3 = VoigtModel(prefix='v3_')
        elif self.type > 3:
            print(str(self.type) + 'is not a selectable option\n')
            
            return
        
        widths=np.arange(1,self.width, 0.001)
        
        smooth_indexes_scipy = scipy.signal.find_peaks_cwt(self.y, widths)
        peak_indexes_xs_ys = np.asarray([list(a) for a in list(zip(self.x[smooth_indexes_scipy], self.y[smooth_indexes_scipy]))])
        high_peak = sorted(peak_indexes_xs_ys, key=lambda pair: pair[1])[-self.num_peaks:]
        print('\n' + str(high_peak[-1][0]) + '\n')
        print(high_peak)

        fit_X = []
        fit_Y = []
        fit_XX = []
        fit_YY = []
        
        fit_range1 = 500

        for i in range(len(self.x)):
            if self.x[i] > high_peak[-1][0]-fit_range1 and self.x[i] < high_peak[-1][0]+fit_range1:
                x  = (self.y[i] - min(self.y))/(max(self.y) - min(self.y))
                fit_XX.append(self.x[i])
                fit_YY.append(x)
                
        fit_XX = np.asarray(fit_XX)
        fit_YY = np.asarray(fit_YY)
        
        smooth_indexes_scipy = scipy.signal.find_peaks_cwt(fit_YY, widths)
        peak_indexes_xs_ys = np.asarray([list(a) for a in list(zip(fit_XX[smooth_indexes_scipy], fit_YY[smooth_indexes_scipy]))])
        high_peak1 = sorted(peak_indexes_xs_ys, key=lambda pair: pair[1])[-2:]
        print('\n' + str(high_peak[-1][0]) + '\n')
        print(peak_indexes_xs_ys)
        print(high_peak1)
                
        fit_range2 = 500
            
        for i in range(len(self.x)):
            if self.x[i] > high_peak[-1][0]-fit_range2 and self.x[i] < high_peak[-1][0]+fit_range2:
                x  = (self.y[i] - min(self.y))/(max(self.y) - min(self.y))
                fit_X.append(self.x[i])
                fit_Y.append(x)
        
        fit_X = np.asarray(fit_X)
        fit_Y = np.asarray(fit_Y)
        
        
        fit_X_p1 = []
        fit_Y_p1 = []
        
        fit_X_p2 = []
        fit_Y_p2 = []
        
        fit_X_p3 = []
        fit_Y_p3 = []                    
        
        for i in range(len(self.x)):
            if self.x[i] > high_peak1[-1][0]-15 and self.x[i] < high_peak1[-1][0]+15:
                x  = (self.y[i] - min(self.y))/(max(self.y) - min(self.y))
                fit_X_p1.append(self.x[i])
                fit_Y_p1.append(x)
            if self.x[i] > high_peak[0][0]-15 and self.x[i] < high_peak[0][0]+15:
                x  = (self.y[i] - min(self.y))/(max(self.y) - min(self.y))
                fit_X_p2.append(self.x[i])
                fit_Y_p2.append(x)
            if self.x[i] > high_peak1[0][0]-15 and self.x[i] < high_peak1[0][0]+15:
                x  = (self.y[i] - min(self.y))/(max(self.y) - min(self.y))
                fit_X_p3.append(self.x[i])
                fit_Y_p3.append(x)
            #return fit_X, fit _Y
        
        fit_X_p1 = np.asarray(fit_X_p1)
        fit_Y_p1 = np.asarray(fit_Y_p1)
        
        fit_X_p2 = np.asarray(fit_X_p2)
        fit_Y_p2 = np.asarray(fit_Y_p2)

        fit_X_p3 = np.asarray(fit_X_p3)
        fit_Y_p3 = np.asarray(fit_Y_p3)

        
        

        
        pars_1 = vigot.guess(fit_Y_p1,x = fit_X_p1)
        
        pars_1['v1_center'].set(high_peak[-1][0], min=high_peak[-1][0]-10, max=high_peak[-1][0]+10)
        
        out_1 = vigot.fit(fit_Y_p1,pars_1,x = fit_X_p1)
        
        
        pars_2 = vigot2.guess(fit_Y_p2,x = fit_X_p2)
        
        pars_2['v2_center'].set(high_peak[0][0], min=high_peak[0][0]-10, max=high_peak[0][0]+10)

        
        out_2 = vigot2.fit(fit_Y_p2,pars_2,x = fit_X_p2)
        
        
        pars_3 = vigot3.guess(fit_Y_p3,x = fit_X_p3)
        
        pars_3['v3_center'].set(high_peak1[0][0], min=high_peak1[0][0]-10, max=high_peak1[0][0]+10)

        
        out_3 = vigot3.fit(fit_Y_p3,pars_3,x = fit_X_p3)
        
        fit_area1 = simps(out_1.best_fit, dx=0.5)
        fit_area2 = simps(out_2.best_fit, dx=0.5)
        fit_area3 = simps(out_3.best_fit, dx=0.5)


        plt.close()
        plt.plot(fit_X,fit_Y,label = self.title  + '\n peaks = ' + str(high_peak[-1][0]) + ' , ' + str(high_peak[0][0]) )
        plt.plot(fit_X_p1,out_1.best_fit,label = 'peak 1 best fit')
        plt.plot(fit_X_p2,out_2.best_fit,label = 'peak 2 best fit')
        plt.plot(fit_X_p3,out_3.best_fit,label = 'peak 3 best fit')
        plt.savefig(str(self.title)+'.png', transparent=True)

        peak_1_fit = out_1.params['v1_center'].value
        peak_2_fit = out_2.params['v2_center'].value
        peak_3_fit = out_3.params['v3_center'].value

        compound_xy = (fit_X,fit_Y)
        fit_xy_1 = (fit_X_p1,out_1.best_fit)
        fit_xy_2 = (fit_X_p2,out_2.best_fit)
        fit_xy_3 = (fit_X_p3,out_3.best_fit)
        fwhm_1 = out_1.params['v1_fwhm'].value
        fwhm_2 = out_2.params['v2_fwhm'].value
        fwhm_3 = out_3.params['v3_fwhm'].value


        send_var = (fit_area1, fit_area2, fit_area3 , peak_1_fit, peak_2_fit, peak_3_fit,
                         compound_xy, fit_xy_1, fit_xy_2, fit_xy_3,fwhm_1,fwhm_2,fwhm_3)

        return send_var
    
    def peaks_1(self):
        
        return None
    
    def waterfall(self):
        
        return None
        
     
        
        
        
        
        