# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 16:56:24 2021

@author: dal1858
"""

#%%

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

#%%

#mat = pd.read_csv('C1_S1_R1_2_001.hicup.G1_all.medium.sorted.locus.txt', sep = '\t', header = None)
mat = pd.read_csv('C1_S1_R1_2_001.hicup.G1_all.medium.sorted.filtered.subsampled.txt', sep = '\t', header = None)
pat = pd.read_csv('C1_S1_R1_2_001.hicup.G2_all.medium.sorted.filtered.txt', sep = '\t', header = None)

mat_all = pd.read_csv('C1_S1_R1_2_001.hicup.G1_all.medium.sorted.filtered.txt', sep = '\t', header = None)
mat = mat[mat[[3, 7]].notnull().all(1)]
pat = pat[pat[[3, 7]].notnull().all(1)]

#%%

biallelic = pd.read_csv('C1_S1_R1_2_001.hicup.no_header.medium.sorted.locus.txt', sep = '\t', header = None)

#%%

locus_len = 73200000 - 72460000

def plot_hic(medium_file, resolution = 10, plot = True, saturate = False, normalize = False):
    
    bin_count = locus_len/(resolution*10**3)
        
    heatmap, xedges, yedges = np.histogram2d(medium_file.iloc[:, 3], medium_file.iloc[:, 7], bins = bin_count)
    
    saturation_level = saturate 
    
    if normalize: 
        if saturate: 
            saturation_level = saturate/sum(sum(heatmap))
        heatmap = heatmap/sum(sum(heatmap))

    #set the maximum value for a given contact 
    if saturate: 
        heatmap = np.clip(heatmap, 0, saturation_level)
        
    if plot:
        plt.imshow(heatmap, cmap = 'Reds')
    
    return heatmap
    

def diff_hic(plot1, plot2, resolution = 20, saturate = False):
    
    observed = plot_hic(plot1, resolution = resolution, saturate = saturate, normalize = True)
    control = plot_hic(plot2, resolution = resolution, saturate = saturate, normalize = True)
    
    sub_plot = observed - control
    plt.imshow(sub_plot, cmap = 'seismic')
    
    


#%%
    
#I will try a logarithmic color plot 

import matplotlib.colors as clrs

fig, ax = plt.subplots()
c = ax.pcolor(bi_plot, norm=clrs.PowerNorm(gamma = 0.05))




#%%


#what about the ratios of Kcnk9-Peg13 vs Kcnk9-Tc9enh? 

def peg13_tc9enh_ratio(hic_df, resolution = 0.5):
    
    #define region coordinates 
    locus = [72460000, 73200000]
    kcnk9 = [72544000, 72553000]
    peg13 = [72805000, 72815000]
    tc9enh = [73017000, 73029000]
    
    locus_len = locus[1] - locus[0]

    hic_matrix = plot_hic(hic_df, resolution = resolution, plot = False)
    
    #change genomic coordinates to matrix coordinates 
    def adjust(region):
        
        adjusted_start = (region[0] - locus[0])/(locus[1] - locus[0]) * len(hic_matrix)
        adjusted_end = (region[1] - locus[0])/(locus[1] - locus[0]) * len(hic_matrix)
    
        adjusted_region = [adjusted_start, adjusted_end]
        
        return adjusted_region
    
    #make a submatrix of just the contact of interest
    def generate_submatrix(region1, region2):
        
        submatrix = hic_matrix[int(adjust(region1)[0]):int(adjust(region1)[1]), int(adjust(region2)[0]):int(adjust(region2)[1])]
        
        return submatrix
    
    peg13_contacts = generate_submatrix(kcnk9, peg13)
    tc9enh_contacts = generate_submatrix(kcnk9, tc9enh)
    
    ratio = sum(sum(peg13_contacts))/sum(sum(tc9enh_contacts))
    
    return ratio



#%%

#plot the mat and pat ratios 

res = 0.15

pat_ratio = peg13_tc9enh_ratio(pat, resolution = res)
mat_ratio = peg13_tc9enh_ratio(mat, resolution = res)


fig = plt.figure(figsize = (3, 5))

x_pos = [0.25, 0.75]
names = ['Maternal', 'Paternal']
plt.bar(x_pos, [mat_ratio, pat_ratio], width = 0.3, align = 'center', color = 'black')
plt.xticks(x_pos, names)
plt.xlim(0, 1)
plt.ylabel('Kcnk9 DMR/Tc9enh Contact Ratio')
plt.title('Kcnk9 Contact Ratio: DMR/Tc9enh')

#%%

ratios = []
resolutions = []
for res in np.linspace(0.05, 1, num = 30):
    
    resolutions.append(res)
        
    pat_ratio = peg13_tc9enh_ratio(pat, resolution = res)
    mat_ratio = peg13_tc9enh_ratio(mat, resolution = res)
    
    ratio = pat_ratio/mat_ratio
    
    ratios.append(ratio)
    
plt.plot(resolutions, ratios)

#%%




def insulation(hic_df, resolution = 0.5):
    
    locus = [72460000, 73200000]
    locus_len = locus[1] - locus[0]
    
    large_tad = [72564000, 73016000]
    subtad1 = [72564000, 72800000]
    subtad2 = [72821000, 73016000]
    
    intra_tad1 = [72567000, 72800000]
    intra_tad2 = [72820000, 73009000]
    
    hic_matrix = plot_hic(hic_df, resolution = resolution, plot = False)
    
    #change genomic coordinates to matrix coordinates 
    def adjust(region):
        
        adjusted_start = (region[0] - locus[0])/(locus[1] - locus[0]) * len(hic_matrix)
        adjusted_end = (region[1] - locus[0])/(locus[1] - locus[0]) * len(hic_matrix)
    
        adjusted_region = [adjusted_start, adjusted_end]
        
        return adjusted_region
        
    #extract just the tad of interest from the matrix 
    def generate_tad(region):
        
        adjusted_tad_region = adjust(region)
        tad_submatrix = hic_matrix[int(adjusted_tad_region[0]): int(adjusted_tad_region[1]), int(adjusted_tad_region[0]): int(adjusted_tad_region[1])]
        
        return tad_submatrix 
    
    #generate a submatrix between two regions 
    def generate_submatrix(region1, region2):
        
        submatrix = hic_matrix[int(adjust(region1)[0]):int(adjust(region1)[1]), int(adjust(region2)[0]):int(adjust(region2)[1])]
        
        return submatrix
    
    intra_tad_matrix = generate_submatrix(intra_tad1, intra_tad2)
    subtad1_matrix = generate_tad(subtad1)
    subtad2_matrix = generate_tad(subtad2)
    
    def matrix_sum(matrix):
        
        my_sum = sum(sum(matrix))
        return my_sum
    
    ratio = matrix_sum(intra_tad_matrix)/(matrix_sum(subtad1_matrix) + matrix_sum(subtad2_matrix)) * 2
    
    return ratio


pat_ins = insulation(pat)
mat_ins = insulation(mat)






#%%

def v4c(hic_df, resolution = 5, anchor_start = 72547000, anchor_end = 72552000, color = 'black', title = None, xlabel = None, ylabel = None, plot = True):

    locus = [72460000, 73200000]
    locus_len = locus[1] - locus[0]
    
    anchor1 = hic_df[(hic_df[3] > anchor_start) & (hic_df[3] < anchor_end)]
    anchor2 = hic_df[(hic_df[7] > anchor_start) & (hic_df[7] < anchor_end)]
    
    anchor1_contacts = anchor1[7]
    anchor2_contacts = anchor2[3]
    
    anchor_contacts = anchor1_contacts.append(anchor2_contacts)
    bins = round( (locus_len * 10**-3) / resolution)
    
    if plot:
        plt.hist(anchor_contacts, bins = bins, density = False, color = color)
        if title:
            plt.title(title)
        if xlabel:
            plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)
    
    return np.histogram(anchor_contacts, bins = bins)



    
#v4c(df, resolution = 3, anchor_start = 72547000, anchor_end = 72552000, color = color, title = f'{name} Virtual 4C: Anchor = Kcnk9', ylabel = 'Contacts with Anchor', xlabel = 'mm10 chr15 coords')

#%%
def v4c_diff(hic_df1, hic_df2, resolution = 5, anchor_start = 72547000, anchor_end = 72552000):
    
    import math 
    
    one = v4c(hic_df1, resolution = resolution, anchor_start = anchor_start, anchor_end = anchor_end, plot = False)
    two = v4c(hic_df2, resolution = resolution, anchor_start = anchor_start, anchor_end = anchor_end, plot = False)

    diff = one[0] / two[0]
    diff = np.log2(diff)
    
    plt.plot(diff)
    return one 

test = v4c_diff(mat, pat, resolution = 3)



#%%

def v4c_avg(hic_df, bin_res = 1, avg_rad = 5, anchor_start = 72547000, anchor_end = 72552000, color = 'black', title = None, xlabel = None, ylabel = None, plot = True, leg_label = None):

    '''
    Virtual 4C with a moving average 
    '''
    
    locus = [72460000, 73200000]
    locus_len = locus[1] - locus[0]
    
    anchor1 = hic_df[(hic_df[3] > anchor_start) & (hic_df[3] < anchor_end)]
    anchor2 = hic_df[(hic_df[7] > anchor_start) & (hic_df[7] < anchor_end)]
    
    anchor1_contacts = anchor1[7]
    anchor2_contacts = anchor2[3]
    
    anchor_contacts = anchor1_contacts.append(anchor2_contacts)
    
    hist_coords = np.linspace(locus[0] + avg_rad, locus[1] - avg_rad, round( (locus_len * 10**-3) / bin_res) - avg_rad * 2)
    
    bin_counts_moving_average = []
    for coord in hist_coords:
        
        contacts = [contact for contact in anchor_contacts if (contact > (coord - avg_rad * 10**3)) & (contact < (coord + avg_rad * 10**3))]
        
        bin_counts_moving_average.append(len(contacts))
            
    if plot:
        plt.plot(bin_counts_moving_average, color = color, label = leg_label)
        x = np.linspace(0, len(bin_counts_moving_average), len(bin_counts_moving_average)*10**-2 + 1)
        labels = np.linspace(locus[0] * 10**-6, locus[1] * 10**-6, 8).round(2)
        plt.xticks(x, labels)
        if title:
            plt.title(title)
        if xlabel:
            plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)

    
    return bin_counts_moving_average


#%%
    
temp_bin_res = 1
temp_avg_rad = 5

kcnk9 = [72544000, 72553000]
tc9enh = [73017000, 73029000]
peg13 = [72805000, 72815000]
locus = [72460000, 73200000]

def reg_len(reg):
    return reg[1] - reg[0]

anchor = kcnk9

mat_v4c = v4c_avg(mat_filtered, bin_res = temp_bin_res, avg_rad = temp_avg_rad, color = 'red', \
        anchor_start = anchor[0], anchor_end = anchor[1], leg_label = 'Maternal')
pat_v4c = v4c_avg(pat_filtered, bin_res = temp_bin_res, avg_rad = temp_avg_rad, color = 'blue', \
        anchor_start = anchor[0], anchor_end = anchor[1], leg_label = 'Paternal')

plt.title(f'Virtual 4C: Anchor = Kcnk9, Normalized')
plt.xlabel('mm10 chr15 coordinates (Mb)')
plt.xticks()
plt.legend()


#%%

diff = np.array(mat_v4c) - np.array(pat_v4c)

plt.plot(diff)
plt.plot(np.zeros(len(diff)), color = 'black')
plt.title('Ref - Alt')


#%%

div = np.array(mat_v4c) / np.array(pat_v4c)

div = np.log2(div)

plt.plot(div)
plt.plot(np.zeros(len(div)), color = 'black')
plt.title('log2(Ref/Alt)')




#%%

#Determine if the reference bias decays with read pair distance 

#my_bin_size = 1

def pair_dist(hic_df, bin_size = my_bin_size):
    
    locus = [72460000, 73200000]
    locus_len = locus[1] - locus[0]
    
    read1 = np.array(hic_df[3])
    read2 = np.array(hic_df[7])
        
    bins = round( (locus_len * 10**-3) / bin_size)
    
    dists = np.histogram(abs(read1 - read2), bins = bins)
    
    return dists
        

def pair_dist_diff(hic_df1, hic_df2, bin_size = 0.1):
    

    dist1 = pair_dist(hic_df1)
    dist2 = pair_dist(hic_df2)

    diff = (dist1[0] / dist2[0])
    plt.plot(diff, label = '1 kb resolution')

    plt.plot(moving_average(diff, 10), label = '10 kb moving average')
    plt.legend()
    plt.title('Allelic bias with contact distance')
    plt.xlabel('Genomic distance (kb)')
    plt.ylabel('log2(mat/mat)')
    
def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

pair_dist_diff(mat_all, pat)



#%%

#curve fitting to get a formula for the allelic bias with contact distance (logarithmic)

avg = moving_average(diff, 100)
limit = round(len(avg) * (600/734))
short = avg[:limit]

c, m, b= np.polyfit(np.log(np.linspace(1, len(short), len(short))), short, 2)

x = np.linspace(1, len(avg), len(avg))
y = c * np.log(x)**2 + m * np.log(x) + b 

#set a lower limit to how much low contacts can be downsampled maternally 

for i, value in enumerate(y):
    if value > max(short):
        y[i] = max(short)

plt.plot(x, y, label = 'Fitted curve')
plt.plot(avg, label = 'Data')
plt.legend()

#%%

#exponential 

from scipy.optimize import curve_fit
#%%
#exponential 
def objective(x, a, b, c, d, e, f):
    return a * x**2 + b * x + c + d * x**3 + e * x**4 + f * x**5

popt, pcov = curve_fit(objective, np.linspace(1, len(short), len(short)), short)

a, b, c, d, e, f = popt
x = np.linspace(1, len(short), len(short))
y = a * x**2 + b * x + c + d * x**3 + e * x**4 + f * x**5


plt.plot(x, y, label = 'Fitted curve')
plt.plot(short, label = 'Data')
plt.legend()



#%%


#subsample both the mat and pat alleles based on the contact distance 

import random 
    

my_bin_size = 0.1

def filter_pat(row):
    
    contact_dist = abs(int(row[3]) - int(row[7]))
    #keep_prob = c * np.log(contact_dist * 10**-2)**2 + m * np.log(contact_dist * 10**-2) + b
    keep_prob = a * (contact_dist * 10**-2)**2 + b * (contact_dist * 10**-2) + c + d * (contact_dist * 10**-2)**3 + e * (contact_dist * 10**-2)**4 + f * (contact_dist * 10**-2)**5
    
    if contact_dist < (my_bin_size * 10**3 * len(y[y > 1])):
        keep = True
    elif keep_prob > random.random():
        keep = True
    else:
        keep = False 
    
    return keep 

pat_filtered = pat[pat.apply(filter_pat, axis = 1) == True]
        



#%%

def filter_mat(row):


    contact_dist = abs(int(row[3]) - int(row[7]))
    #keep_prob = (c * np.log(contact_dist * 10**-2)**2 + m * np.log(contact_dist * 10**-2) + b)**-1
    keep_prob = (a * (contact_dist * 10**-2)**2 + b * (contact_dist * 10**-2) + c + d * (contact_dist * 10**-2)**3 + e * (contact_dist * 10**-2)**4 + f * (contact_dist * 10**-2)**5)**-1
    
    if contact_dist > (my_bin_size * 10**3 * len(y[y > 1])):
        keep = True
    elif keep_prob > random.random():
        keep = True
    else:
        keep = False 

    return keep

mat_filtered = mat_all[mat_all.apply(filter_mat, axis = 1) == True]

#%%


mat_dist_filt = pair_dist(mat_filtered)
pat_dist_filt = pair_dist(pat_filtered)

diff_filt = (mat_dist_filt[0] / pat_dist_filt[0])
plt.plot(diff_filt, label = '1 kb resolution')