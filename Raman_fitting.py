# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 00:15:54 2019

@author: leeva




"""

import os
from sys import path
import numpy as np
import fittingV_ver2

'''
Dependency: lmfit,numpy,scipy,matplotlib
'''

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

while True:
    try:
        data_path = input('Input data path: \n')
        if os.path.isdir(data_path) != True:
            raise ValueError
        else:
            break
    except ValueError:
        print('Please input valid path')

while True:
    try:
        frmt = input('Input file format (form ".txt" or ".dat"): \n')
        if frmt != '.dat' and frmt != '.txt':
            raise ValueError
        else:
            break
    except ValueError:
        print('please input a valid file path')


while True:    
    try:
        num = int(input('for all input 0 else input number of spectra: \n'))
        break
    except ValueError:
        print('Please input an integer')
        

if num !=0:
    print('input spectra prefixes ie \' p0\' + CR for the appropriate number of files')
    for i in range(num):
        pre.append(input() + '_')
        count = i
        print('number of files ' + str(count+1) + ' of ' + str(num) + '\n')
        if count+1 == num:
            print('done')

while True:
    try:
        width = int(input('Input the max peak width:\n'))
        break
    except ValueError:
        print('please input a valid number')

while True:
    try:
        fit_t = int(input('Input the number of peaks to fit \n 2 Peaks == 2 \n 3 Peaks == 3:\n'))
        if fit_t != 2 and fit_t != 3:
            raise ValueError
        else:
            break
    except ValueError:
        print('please select a valid option')

while True:
    try:
        fit_type = int(input('Input the fit type:\n 1 = Lorentzian\n 2 = Voigt\n'))
        if fit_type != 1 and fit_type != 2:
            raise ValueError
        else:
            break
    except ValueError:
        print('please select a valid option')

if num != 0:
    for i in range(len(pre)):
            for entry in os.scandir(path = data_path):
                        if entry.is_dir():
                            folders.append(entry.path)
                        if (entry.name.startswith(pre[i]) and (entry.name.endswith(frmt) == 1)):
                            files.append(entry.path)
                            name.append(entry.name)
                            print(entry.name)
else:
    for entry in os.scandir(path = data_path):
                if entry.is_dir():
                    folders.append(entry.path)
                if ((entry.name.endswith(frmt) == 1)):
                    files.append(entry.path)
                    name.append(entry.name)
                    print(entry.name)



for i in range(len(files)):
            infile = open(files[i], 'r')
            outfile = 'new_' + os.path.splitext(name[i])[0] + '.png'
            line = infile.readline()
            data = np.loadtxt(infile)
            if frmt == '.dat':
                shift = data[:, 0]
                signal = data[:, 2]
            elif frmt == '.txt':
                shift = data[:, 0]
                signal = data[:, 1]
            spec.append(Data(shift,signal,outfile,name[i],os.path.splitext(name[i])[0]))
            rama.append(fittingV_ver2.Multi_fit(spec[i].shift,spec[i].signal,spec[i].name,n,width,fit_type))

            


print('\nInput keyword of spectra ie: (Ruby_after) or (3600cm):\n')
sped = str(input())
for i in range(len(spec)):
    v = sped in spec[i].name
    if v == True:
        m.append(i)

'''
Change the sped == 'x' to where x is your own keyword for your file ie its central measurement wavenumber
'''

if fit_t == 2:
    print('here_2')
    save_path = (str(os.getcwd()) + '\\double_peak_fit\\' + sped + '_output data.txt')
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    output = open(save_path, 'w')
    output.write('name\t peak_1\t peak_2\t area_1\t area_2\t area_total\t ratio\tfwhm_1\t fwhm_2\n')

    for i in range(len(m)):
        save_path2 = (str(os.getcwd()) + '\\double_peak_fit\\' + rama[m[i]].title +'.txt')
        output2 = open(save_path2, 'w')
        output2.write('peak1_x\tpeak2_x\tpeak1_fit\tpeak2_fit\n')
        compound_x,compound_y,peak1_x,peak1_fit,peak2_x,peak2_fit,total_area,area_1,area_2,ratio,peak1_pos,\
        peak2_pos,fwhm_1,fwhm_2 = fittingV_ver2.Multi_fit.peaks_2(rama[m[i]])
        output.write(rama[m[i]].title+'\t'+str(peak1_pos)+'\t'+str(peak2_pos)+'\t'+str(area_1)+'\t'+str(area_2)+'\t'+str(total_area)
                     +'\t'+str(ratio)+'\t' + str(fwhm_1) + '\t' + str(fwhm_2) + '\n')

        output2.close()

elif fit_t == 3:
    print('here_3')
    save_path = (str(os.getcwd()) + '\\triple_peak_fit\\' + sped + '_output data.txt')
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    output = open(save_path, 'w')
    output.write('name\t peak_1\t peak_2\t peak_3\t area_1\t area_2\t area_3\t fwhm_1\t fwhm_2\t fwhm_3\n')

    for i in range(len(m)):
        save_path2 = (str(os.getcwd()) + '\\triple_peak_fit\\' + rama[m[i]].title +'.txt')
        output2 = open(save_path2, 'w')
        output2.write('peak1_x\tpeak2_x\tpeak3_x\tpeak1_fit\tpeak2_fit\tpeak3_fit\n')

        out = fittingV_ver2.Multi_fit.peaks_3(rama[m[i]])

        output.write(rama[m[i]].title+'\t'+str(out[3])+'\t'+str(out[4])+'\t'+str(out[5])+'\t'
                     +str(out[0])+'\t'+str(out[1])+'\t'+str(out[2])+'\t' +str(out[10])+ '\t' +
                     str(out[11]) + '\t' + str(out[12]) + '\n')


        output2.write(str(out[7][0]) + '\t' + str(out[8][0]) + '\t' + str(out[9][0]) + '\t'
                          + str(out[7][1]) + '\t' + str(out[8][1]) + '\t' + str(out[9][1]) + '\n')
        output2.close()



output.close()

'''
Output data is found in a folder with the keyword name
'''
        
        

