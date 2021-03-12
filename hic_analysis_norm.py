# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 11:23:03 2021

@author: dal1858
"""

#%%
'''
In this script I will attempt to normalize my allelic chromatin maps to each other. 
I will determine if the allelic bias changes with the distance between the contacts, 
and then downsample the two alleles in proportion to their bias at that contact distance. 
'''

#%%
#Load modules

import random 
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from statsmodels.nonparametric.smoothers_lowess import lowess

#%%
#Load data. Contacts are in form of medium files (see Juicer pre for details). Loading all data
#within my locus, without any downsampling of either allele. 

#mat = pd.read_csv('C1_S1_R1_2_001.hicup.G1_all.medium.sorted.filtered.subsampled.txt', sep = '\t', header = None)

mat_file = r'W:\Users\dloftus\tan_etal_2020_dipc_mouse_hippocampus\mat_pat_aggregated_contact_txt_files\GSE146397_aggregated_data.allele_maternal_maternal.contacts.pairs.filtered.txt'
pat_file = r'W:\Users\dloftus\tan_etal_2020_dipc_mouse_hippocampus\mat_pat_aggregated_contact_txt_files\GSE146397_aggregated_data.allele_paternal_paternal.contacts.pairs.filtered.txt'

pat = pd.read_csv(mat_file, sep = '\t', header = None)
mat = pd.read_csv(pat_file, sep = '\t', header = None)

mat = mat[mat[[3, 7]].notnull().all(1)]
pat = pat[pat[[3, 7]].notnull().all(1)]


#%%


def pair_dist(hic_df, bin_size = 1, juicer_format = 'Medium'):
    
    '''
    generates a numpy histogram of the distances between HiC contacts
    '''
    
    locus = [72460000, 73200000]
    locus_len = locus[1] - locus[0]
    
    if juicer_format == 'Medium':
        read1 = np.array(hic_df[3])
        read2 = np.array(hic_df[7])
    elif juicer_format == 'Short':
        read1 = np.array(hic_df[2])
        read2 = np.array(hic_df[4])        
    
    keep_contact = abs(read1 - read2) <= locus_len
    read1 = read1[keep_contact]
    read2 = read2[keep_contact]
        
    bins = round( (locus_len * 10**-3) / bin_size)
    
    dists = np.histogram(abs(read1 - read2), bins = bins)
    
    return dists
  
def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w      

def pair_dist_diff(hic_df1, hic_df2, bin_size = 1, avg_bin_size = 10, max_dist = None, average_only = False, label = None, legend_outside = False, legend = True, title = False, juicer_format = 'Medium'):
    

    dist1 = pair_dist(hic_df1, bin_size = bin_size, juicer_format = juicer_format)
    dist2 = pair_dist(hic_df2, bin_size = bin_size, juicer_format = juicer_format)

    diff = np.log2(dist1[0] / dist2[0])
    
    resolution = bin_size

    #plt.plot([0, len(diff[:max_dist])], [0, 0], color = 'black')
    
    if average_only == False:
        plt.plot(diff[:max_dist], label = f'{resolution} kb resolution')
    
    bins_to_average = int(avg_bin_size * bin_size**-1)
    if label:
        label = label
    elif label == None:
        label = f'{avg_bin_size} kb moving average'
        
    plt.plot(moving_average(diff[:max_dist], bins_to_average), label = label)
    if legend_outside:
        plt.legend(bbox_to_anchor=(0.999, 1.03), loc=2)
    elif legend:
        plt.legend()
   
    if title:
        plt.title(title)
    else:
        plt.title('Allelic bias with contact distance')
    plt.xlabel('Genomic distance (kb)')
    plt.ylabel('log2(Mat/Pat)')
    
    
    return diff
    

#my_bin_size = 1
#diff = pair_dist_diff(mat, pat, bin_size = my_bin_size, avg_bin_size = 10, legend_outside = True)
#%%

#curve fitting to get a formula for the allelic bias with contact distance (logarithmic)

limit = round(len(diff) * (600/734))
short = diff[:limit]

c, m, b= np.polyfit(np.log(np.linspace(1, len(short), len(short))), short, 2)

x = np.linspace(1, len(diff), len(diff))
y = c * np.log(x)**2 + m * np.log(x) + b 

#set a lower limit to how much low contacts can be downsampled maternally 


for i, value in enumerate(y):
    if value > max(short):
        y[i] = max(short)

plt.plot(diff, label = 'Data')
plt.plot(x, y, label = 'Fitted curve')
plt.title('Allelic bias vs contact distance, fitted')
plt.xlabel('Contact distance (kb)')
plt.ylabel('Mat/Pat')
plt.legend()

#%%

#plot residuals of logarithmic

residuals = diff - y
plt.plot(residuals)

#%%
#exponential 
def objective(x, a, b, c, d, e, f):
    return a * x**2 + b * x + c + d * x**3 + e * x**4 + f * x**5

#ignore nan
valid = ~(np.isnan(np.linspace(1, len(diff), len(diff))) | np.isnan(diff))
popt, pcov = curve_fit(objective, np.linspace(1, len(diff), len(diff))[valid], diff[valid])

a, b, c, d, e, f = popt
x = np.linspace(1, len(diff), len(diff))
y = a * x**2 + b * x + c + d * x**3 + e * x**4 + f * x**5

plt.plot(diff, label = 'Data')
plt.plot(x, y, label = 'Fitted curve')
plt.title('Allelic bias vs contact distance, fitted')
plt.xlabel('Contact distance (kb)')
plt.ylabel('Mat/Pat')
plt.legend()



#%%


def filter_pat(row, norm = 'logarithmic'):
    
    contact_dist = abs(int(row[3]) - int(row[7]))
    if norm == 'logarithmic':
        keep_prob = c * np.log(contact_dist * 10**-3)**2 + m * np.log(contact_dist * 10**-3) + b
    elif norm == 'exponential':
        keep_prob = a * (contact_dist * 10**-3)**2 + b * (contact_dist * 10**-3) + c + d * (contact_dist * 10**-3)**3 + e * (contact_dist * 10**-3)**4 + f * (contact_dist * 10**-3)**5
    
    if contact_dist < (my_bin_size * 10**3 * len(y[y > 1])):
        keep = True
    elif keep_prob > random.random():
        keep = True
    else:
        keep = False 
    
    return keep 

pat_filtered = pat[pat.apply(filter_pat, norm = 'exponential', axis = 1) == True]
        



#%%

c, m, b= np.polyfit(np.log(np.linspace(1, len(short), len(short))), short, 2)

x = np.linspace(1, len(diff), len(diff))
y = c * np.log(x)**2 + m * np.log(x) + b 
plt.plot(x, y)
#%%
def filter_mat(row, norm = 'logarithmic'):

    contact_dist = abs(int(row[3]) - int(row[7]))
    
    if norm == 'logarithmic':
       
        keep_prob = (c * np.log(contact_dist * 10**-3)**2 + m * np.log(contact_dist * 10**-3) + b)**-1
        if keep_prob < (max(short))**-1:
            keep_prob = (max(short))**-1
    elif norm == 'exponential':
        keep_prob = (a * (contact_dist * 10**-3)**2 + b * (contact_dist * 10**-3) + c + d * (contact_dist * 10**-3)**3 + e * (contact_dist * 10**-3)**4 + f * (contact_dist * 10**-3)**5)**-1
    
    if contact_dist > (my_bin_size * 10**3 * len(y[y > 1])):
        keep = True
    elif keep_prob > random.random():
        keep = True
    else:
        keep = False 

    return keep

mat_filtered = mat[mat.apply(filter_mat, norm = 'exponential', axis = 1) == True]

#%%


mat_dist_filt = pair_dist(mat_filtered)
pat_dist_filt = pair_dist(pat_filtered)

diff_filt = (mat_dist_filt[0] / pat_dist_filt[0])


plt.plot(diff_filt, label = '1 kb resolution')
plt.plot(moving_average(diff_filt, 10), label = '10 kb moving average')

plt.title('Allelic bias vs contact distance, corrected')
plt.xlabel('Contact distance (kb)')
plt.ylabel('Mat/Pat')
plt.legend()


#%%

#export medium files of filtered hic maps 

mat_filtered.to_csv('W:\\Users\\dloftus\\capture_hic\\results_run1\\hicup\\allelic_snpsplit\\output\\subsampling\\maternal_log_dist_corrected.medium.txt', sep = '\t', header = None, index = False)
pat_filtered.to_csv('W:\\Users\\dloftus\\capture_hic\\results_run1\\hicup\\allelic_snpsplit\\output\\subsampling\\paternal_log_dist_corrected.medium.txt', sep = '\t', header = None, index = False)
#%%

#repeat the correction process 

limit = round(len(diff_filt) * (600/734))
short = diff_filt[:limit]

c, m, b= np.polyfit(np.log(np.linspace(1, len(short), len(short))), short, 2)

x = np.linspace(1, len(diff_filt), len(diff_filt))
y = c * np.log(x)**2 + m * np.log(x) + b 

#set a lower limit to how much low contacts can be downsampled maternally 


for i, value in enumerate(y):
    if value > max(short):
        y[i] = max(short)

plt.plot(diff_filt, label = 'Data')
plt.plot(x, y, label = 'Fitted curve')
plt.legend()








#%%

#what about loess correction? 

from statsmodels.nonparametric.smoothers_lowess import lowess
lowess_normed = lowess(diff, range(len(diff)), return_sorted = False, frac = 0.2)

plt.plot([0, len(diff)], [0,0], color = 'black')
plt.plot(diff[:])
plt.plot(lowess_normed)
















#%%

#let's look at autosomes except chr15 

pat_auto = pd.read_csv('W:\\Users\\dloftus\\capture_hic\\results_run1\\hicup\\allelic_snpsplit\\output\\subsampling\\C1_S1_R1_2_001.hicup.G2_all.medium.sorted.autosomes_except_15.txt', sep = '\t', header = None)
mat_auto = pd.read_csv('W:\\Users\\dloftus\\capture_hic\\results_run1\\hicup\\allelic_snpsplit\\output\\subsampling\\C1_S1_R1_2_001.hicup.G1_all.medium.sorted.autosomes_except_15.txt', sep = '\t', header = None)
mat_auto = mat_auto[mat_auto[[3, 7]].notnull().all(1)]
pat_auto = pat_auto[pat_auto[[3, 7]].notnull().all(1)]



#%%

diff_auto = pair_dist_diff(mat_auto, pat_auto, bin_size = my_bin_size, avg_bin_size = 10)

#%%

def plot_diff(diff_array, max_dist = 300, log2 = True):
    
    if log2:
        plt.plot(np.log2(diff_array[:max_dist]))
        plt.plot(moving_average(np.log2(diff_array[:max_dist])))
    else:
        plt.plot(diff_array[:max_dist])
        plt.plot(moving_average(diff_array[:max_dist]))
    
#%%
   
#make a dictionary containing mat and pat contacts for chr 1, 2, and 3 separately (6 key:value pairs in total)
chr_123 = {}    
for chromosome, number in zip(['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr15'], [1, 2, 3, 4, 5, 6, 15]):
    
    mat_chr = pd.read_csv(f'W:\\Users\\dloftus\\capture_hic\\results_run1\\hicup\\allelic_snpsplit\\output\\subsampling\\C1_S1_R1_2_001.hicup.G1_all.medium.sorted.{chromosome}.txt', sep = '\t', header = None)
    pat_chr = pd.read_csv(f'W:\\Users\\dloftus\\capture_hic\\results_run1\\hicup\\allelic_snpsplit\\output\\subsampling\\C1_S1_R1_2_001.hicup.G2_all.medium.sorted.{chromosome}.txt', sep = '\t', header = None)
    mat_chr = mat_chr[mat_chr[[3, 7]].notnull().all(1)]
    pat_chr = pat_chr[pat_chr[[3, 7]].notnull().all(1)]
    
    chr_123[f'{chromosome}_mat'] = mat_chr
    chr_123[f'{chromosome}_pat'] = pat_chr
    
pat_locus = pd.read_csv('W:\\Users\\dloftus\\capture_hic\\results_run1\\hicup\\allelic_snpsplit\\output\\subsampling\\C1_S1_R1_2_001.hicup.G2_all.medium.sorted.filtered.txt', sep = '\t', header = None)
mat_locus = pd.read_csv('W:\\Users\\dloftus\\capture_hic\\results_run1\\hicup\\allelic_snpsplit\\output\\subsampling\\C1_S1_R1_2_001.hicup.G1_all.medium.sorted.filtered.txt', sep = '\t', header = None)
mat_locus = mat_locus[mat_locus[[3, 7]].notnull().all(1)]
pat_locus = pat_locus[pat_locus[[3, 7]].notnull().all(1)]

chr_123['locus_mat'] = mat_locus
chr_123['locus_pat'] = pat_locus

#%%
    
#loop through these and generate a log2 plot for all of them 

plt.plot([0, len(diff[:700])], [0, 0], color = 'black')

for chromosome in 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr15', 'locus': 
    
    diff_temp = pair_dist_diff(chr_123[f'{chromosome}_mat'], chr_123[f'{chromosome}_pat'], \
                               average_only = True, label = f'{chromosome}', max_dist = 700, \
                               legend_outside = True, avg_bin_size = 50)
    
    
    
#%%
    
#Subsample the mat and pat locus contacts to see what the log2(mat/pat) distribution is 
    
plt.plot([0, len(diff[:700])], [0, 0], color = 'black')

for i in range(10):
    mat_sub = mat_locus.sample(frac = 0.1)
    pat_sub = pat_locus.sample(frac = 0.1)
    print(len(pat_sub))
    
    pair_dist_diff(mat_sub, pat_sub,  average_only = True, max_dist = 700, legend = False, \
                   title = '10x 0.1 subsample of locus', avg_bin_size = 50)











#%%
    
#what about du et al, 2017? 
    
mat_du = pd.read_csv('W:\\Users\\dloftus\\due_etal_2017_allelic_hic\\8cell_rep123_maternal_allValidPairs_Peg13locus.txt', sep = '\t', header = None)
pat_du = pd.read_csv('W:\\Users\\dloftus\\due_etal_2017_allelic_hic\\8cell_rep123_paternal_allValidPairs_Peg13locus.txt', sep = '\t', header = None)

diff = pair_dist_diff(mat_du, pat_du, bin_size = my_bin_size, avg_bin_size = 10, juicer_format = 'Short', max_dist = 300, \
                      title = 'Du et al., 2017. Allelic HiC. Peg13 locus. 8-cell embryo.')
