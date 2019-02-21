# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 00:15:54 2019

@author: leeva
"""

import os
from sys import path
path.append(os.getcwd() + "\\classes")
import numpy as np
import plot
import matplotlib.pyplot as plt

import fittingV_ver2
class Data:
    
    def __init__(self,shift,signal,outfile,name,filename):
        self.shift = shift
        self.signal = signal
        self.outfile = outfile
        self.name = name
        self.filename = filename
        
    
       

    
folders = []
files = []
name = []
spec = []
pre = []
rama = []
m = []

n = 10000
cond = 0
print('for all input 0 else input number of spectra: ')
num = int(input())
if num == 0:
   cond = 1
elif num == -1 or num == -11:
    pre.append(str('p1_'))
#    pre.append(str('p15_'))
#    pre.append(str('p16_'))
#    pre.append(str('p2_'))
#    pre.append(str('p3_'))
#    pre.append(str('p6_'))
#    pre.append(str('p5_'))
#    pre.append(str('p6_'))
#    pre.append(str('p29_'))
#    pre.append(str('p30_'))
#    pre.append(str('p32_'))
#    pre.append(str('p31_'))
#    pre.append(str('p2_'))
#    pre.append(str('p29_'))
    cond = 0
elif num !=0 or num != -1:
    print('input spectra prefixes ie \' p0\' + CR for the appropriate number of files')
    for i in range(num):
        pre.append(input() + '_')
        count = i
        print('number of files ' + str(count+1) + ' of ' + str(num) + '\n')
        if count+1 == num:
            print('done')

width = int(input('Input the max peak width:\n'))
n = int(input('Input the number of peaks to fit:\n'))
fit_type = int(input('Input the fit type:\n 1 = Lorentzian\n 2 = Voigt\n'))
#print(cond)
        
frmt = '1.txt'
'''
\\high pressure cell\\pressure 1
'''
if cond == 0:
    for i in range(len(pre)):
            for entry in os.scandir(path = os.getcwd() + "\\Mg(OH)2\\high pressure cell\\pressure 1"):
                        if entry.is_dir():
                            folders.append(entry.path)
                        if (entry.name.startswith(pre[i]) and (entry.name.endswith(frmt) == 1)):
                            files.append(entry.path)
                            name.append(entry.name)
if cond == 1:
    for entry in os.scandir(path = os.getcwd() + "\\Mg(OH)2\\high pressure cell\\pressure 1"):
                if entry.is_dir():
                    folders.append(entry.path)
                if entry.name.endswith(frmt) == 1:
                    files.append(entry.path)
                    name.append(entry.name)
    
#print(len(files))
for i in range(len(files)):
            infile = open(files[i], 'r')
            outfile = 'new_' + os.path.splitext(name[i])[0] + '.png'
            line = infile.readline()
            data = np.loadtxt(infile)
            shift = data[:,0]
            signal = data[:,1]
            
            spec.append(Data(shift,signal,outfile,name[i],os.path.splitext(name[i])[0]))
            
            
            rama.append(fittingV_ver2.Multi_fit(spec[i].shift,spec[i].signal,spec[i].name,n,width,fit_type))

            
            
if num == -1:
    sped = str(3600)
    print(sped)
    for i in range(len(spec)):
        v = sped in spec[i].name
        #d = 'ruby_after' in spec[i].name
        if v == True:
            m.append(i)
elif num == -11:
    sped = str(500)
    print(sped)
    for i in range(len(spec)):
        v = sped in spec[i].name
        #d = 'ruby_after' in spec[i].name
        if v == True:
            m.append(i)
    
else:
    print('\nInput keyword of spectra:\n')
    sped = str(input())
    for i in range(len(spec)):
        v = sped in spec[i].name
        #d = 'ruby_after' in spec[i].name
        if v == True:
            m.append(i)
        




#plot.Plotter.draw_multi(rama,m)
#for i in range(len(m)):


#if sped == 'ruby_after':
#    save_path = (os.getcwd() + '\\ruby_after\\' + sped + '_output data1.txt')
#    os.makedirs(os.path.dirname(save_path), exist_ok=True)
#    output = open(sped + '_output data.txt','w')
#    output.write('name\t shift\n')
#    area,ruby = fittingV.Voigtfit.fitV(spec,n,m,rama,sped)
#    for i in range(len(ruby)):
#        output.write(ruby[i])
if sped == '3600':
    save_path = (str(os.getcwd()) + '\\3600cm\\' + sped + '_output data.txt')
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    output = open(save_path,'w')
    output.write('name\t peak_1\t peak_2\t area_1\t area_2\t area_total\t ratio\t fwhm_1\t fwhm_2\n')
    for i in range(len(m)):
        if n == 2:
            compound_x,compound_y,peak1_x,peak1_fit,peak2_x,peak2_fit,total_area,area_1,area_2,ratio = fittingV_ver2.Multi_fit.peaks_2(rama[m[i]])
        elif n == 3:
            x,y = fittingV_ver2.Multi_fit.peaks_3(rama[m[i]])
#    plt.plot(rama[m[0]].x,rama[m[0]].y)
#    plt.show()
#    plt.plot(compound_x,compound_y)
#    plt.show()
#    plt.plot(peak2_x,peak2_fit)
#elif sped == '3300':
#    print('input number of peaks to fit: \n')
#    n = int(input())    
#    save_path = (str(os.getcwd()) + '\\3300cm\\' + sped + '_output data.txt')
#    os.makedirs(os.path.dirname(save_path), exist_ok=True)
#    output = open(save_path,'w')
#    output.write('name\t peak_1\t peak_2\t area_1\t area_2\t area_total\t ratio\t fwhm_1\t fwhm_2\n')
#    area,ruby = fittingV.Voigtfit.fitV(spec,n,m,rama,sped)
#    for i in range(len(area)):
#        output.write(area[i])
elif sped == '500':
#   n = 3#input('input number of peaks to fit: \n')
    save_path = (str(os.getcwd()) + '\\500cm\\' + sped + '_output data.txt')
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    output = open(save_path,'w')
    output.write('name\t area_1\t area_2\t area_total\t ratio\n')
    for i in range(len(m)):
        x,y,z = fittingV_ver2.Multi_fit.peaks_3(rama[m[i]])
#    area,ruby = fittingV.Voigtfit.fitV(spec,n,m,rama,sped)
#    for i in range(len(area)):
#        output.write(area[i])
output.close()
    #plt.plot(x,y)
    #plt.close()
#root = peak_find(rama[m[0]].x,rama[m[0]].y,n)

        
        

