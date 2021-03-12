#!/usr/bin/env python
# -*- coding: utf-8 -*-


###
# This script takes in a medium file and filters it based on 
# coordinates in a bed file, and subsamples the resulting 
# medium file based on user input
###


import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-m', 
                    '--medium', 
                    type = argparse.FileType('r'), 
                    help = 'Input file in Medium format', 
                    required = True)
parser.add_argument('-b', 
                    '--bed', 
                    type = argparse.FileType('r'), 
                    help = 'Bed file containing regions to filter for. Currently only works when chromosome name is an integer, without chr',
                    required = False)
parser.add_argument('-s', 
                    '--subsample', 
                    type = float, 
                    default = 1, 
                    help = 'Subsampling fraction to keep', 
                    required = False)
parser.add_argument('-o', 
                    '--output',
                    type = str, 
                    help = 'Output file name', 
                    required = True)


args = parser.parse_args()


def main(med = args.medium, 
         bed = args.bed, 
         s = args.subsample, 
         out = args.output):
    
    import pandas as pd 
    
    print("Loading medium file")
    med_df = pd.read_csv(med, sep = '\t', header = None)
    print("Medium file loaded")
    
    
    #remove chr from medium file
    def remove_chr(med_file):
        try: 
            med_file = int(med_file.replace('chr', ''))
        except ValueError:
            med_file = med_file.replace('chr', '')
        
        return med_file 
        
    try:
        if 'chr' in med_df.iloc[0, 2]:
            med_df[2] = med_df[2].apply(remove_chr)
            med_df[6] = med_df[6].apply(remove_chr)
    except TypeError:
        pass
        
    #filter medium file for coordinates in bed file 
    if bed:
        print("Filtering from bed")
        bed_df = pd.read_csv(bed, sep = '\t', header = None)
        filtered_df_list = []
        for index, coordinate in bed_df.iterrows():
            coordinate_med_df = med_df[(med_df[2] == coordinate[0]) \
                            & (med_df[3] > coordinate[1]) \
                            & (med_df[3] < coordinate[2]) \
                            & (med_df[6] == coordinate[0]) \
                            & (med_df[7] > coordinate[1]) \
                            & (med_df[7] < coordinate[2])]
            filtered_df_list.append(coordinate_med_df)
            
        filtered_df = pd.concat(filtered_df_list)
            
        subsampled = filtered_df.sample(frac = s)
    else:
        subsampled = med_df.sample(frac = s)
        
    #add 'chr' back into chromosome name for juicer pre
    def add_chr(chromosome):
        chromosome = 'chr' + str(chromosome)
        return chromosome
    
    subsampled[2] = subsampled[2].apply(add_chr)
    subsampled[6] = subsampled[6].apply(add_chr)
    
    subsampled.to_csv(out, sep = '\t', header = None, index = False)
                    
        

    
    
    
    
if __name__ == '__main__': 
    main()
