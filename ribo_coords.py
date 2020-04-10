import csv
from Bio.Seq import Seq
import os 
import textwrap
import pandas as pd
from datetime import datetime
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
import cPickle as pickle
import numpy as np
import ribo_util
from IPython.display import display
from collections import Counter



def pausescore_g_coords(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    multi     = inputs['multiprocess']
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started pause score analysis at " + str(datetime.now())

    for fname in files:
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/'
        if not os.path.exists(path_pausescore):
            os.makedirs(path_pausescore)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
       
        if not multi == 'yes':
            run_pausescore_g_coords(fname, settings, plus, minus, gff_dict, path_pausescore)
        else:     
            argument = [fname, settings, plus, minus, gff_dict, path_pausescore]
            arguments.append(argument)
    
    if multi == 'yes':
        ribo_util.multiprocess(run_pausescore_g_coords, arguments, threads)
    
    print "Finished pause score analysis at " + str(datetime.now())
    
    return


def run_pausescore_g_coords(fname, settings, plus, minus, gff, path_pausescore):
    
    '''define variables'''
    
    minlength   = settings['minlength']
    maxlength   = settings['maxlength']
    alignment   = settings['alignment']
    threshold   = settings['threshold']


    lengthindex = range(minlength, maxlength + 1)
    
    A_site = settings['A_site shift']
    P_site = A_site - 3
    E_site = A_site - 6
    two_site = A_site + 6
    one_site = A_site + 3
    mone_site = A_site - 9
    mtwo_site = A_site - 12
    
    
    frameshift      = settings['frameshift']

    plot_upstream   = settings['plot_upstream'] / 3 * 3        #change window to interval of 3
    plot_downstream = settings['plot_downstream'] / 3 * 3
    
    next_codon = settings['next_codon']
    
    
    # define plot length
    plotlength = plot_upstream + plot_downstream + 1
    positionindex = range(0, plotlength)
    
    # load density files
    density_plus  = plus
    density_minus = minus 
    
    # load annotation
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    seq_list    = gff_dict['Sequence'] 
 
    gff_extra  = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    start_trim = settings['start_trim'] / 3 * 3
    stop_trim  = settings['stop_trim'] / 3 * 3   
    
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    plot_upstream_1   = str(plot_upstream)  +'_'
    plot_downstream_1 = str(plot_downstream)+'_'
    start_trim_1      = str(start_trim)     +'_'
    stop_trim_1       = str(stop_trim)      +'_'
    frameshift_1      = str(frameshift)     +'_'
    a_site_1          = str(A_site)         +'_'
    
    name_settings = minlength_1+maxlength_1+plot_upstream_1+plot_downstream_1
    name_settings += start_trim_1+stop_trim_1+frameshift_1
        
    # import genetic code
    aa_code, codon_code = ribo_util.get_genetic_code()
    
    
    '''output data structure'''
    
    # aa/codon avgplot      = { aa: {length: [values]}} 
    # aa/codon count        = { aa: N}
    # aa/codon score        = { aa: [aa_list], A site score: [A_list]...
                                
    aa_avgplot  = {}
    aa_count    = {}
    aa_score    = {}
    aa_score_df = {}
    aa_list     = []
    aa_A_list   = []
    aa_P_list   = []
    aa_E_list   = []
    aa_1_list   = []
    aa_2_list   = []
    aa_m1_list   = []
    aa_m2_list   = []
    
    codon_avgplot  = {}
    codon_count    = {}
    codon_score    = {}
    codon_score_df = {}
    codon_list     = []
    codon_A_list   = []
    codon_P_list   = []
    codon_E_list   = []
    codon_1_list   = []
    codon_2_list   = []
    codon_m1_list   = []
    codon_m2_list   = []
    
    # create empty datastructures to store density info:
    for i in range(1, 214):
        aa_avgplot[i] = {length : [0]*(plotlength) for length in lengthindex}
        aa_count[i]   = 0
        aa_score[i]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
        
    for codon in aa_code['G']:
        for i in range(1, 214):
            codon_avgplot[codon + '_' + str(i)] = {length : [0]*(plotlength) for length in lengthindex}
            codon_count[codon + '_' + str(i)]   = 0
            codon_score[codon + '_' + str(i)]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
    
    '''genes in data''' 
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']         = [0, []]
    excluded_genes['low_density']   = [0, []]
    excluded_genes['not_divisible'] = [0, []]
    # count included genes
    included_genes = [0, []]
    
    
    
    #list_location = '/Volumes/HDD/Ribo_seq/libraries/analysis/reference_information/Workbook1.csv'
    
    #infile = pd.read_csv(list_location)
    #TE_list = infile['yeaR'].tolist()
        
    '''iterate through every annotated gene to get codon density info:'''
    
    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,stop_list, strand_list, seq_list):
        
        #if not alias in TE_list:
            #continue
        
        
        ''' define start and stop positions for codons to analyze:
        
        # = codon 
        #################################### = GFF sequence  (50 extra nt)
           ##############################    = AUG to UGA
             ##########################      = density to analyze : remove start and stop peaks
                 ##################          = codons to analyze : remove plot window
        
        '''
        
        # First, define density without start and stop peaks:
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim + frameshift
            density_stop  = stop  - stop_trim + frameshift
            
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim - frameshift
            density_stop  = stop  + stop_trim - 3 + frameshift
            
            period = -1
        
        # GFF seq has 50 extra nucleotides, so remove:
        # Also remove several codons from start and stop positions:
        
        codon_seq_start = gff_extra + start_trim + plot_upstream + frameshift
        codon_seq_stop  = -gff_extra - stop_trim - plot_downstream + frameshift
            
        seq       = sequence[codon_seq_start : codon_seq_stop]
        seqlength = len(seq)
        
        # exclude genes that are not divisable by 3
        if seqlength % 3 != 0:
            excluded_genes['not_divisible'][0] += 1
            excluded_genes['not_divisible'][1].append(alias)
            continue
        
        # exclude genes shorter than plot
        if seqlength < plotlength + 1:
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        elif codon_seq_start - codon_seq_stop > len(sequence):
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
    
        # make a list of codons in the sequence:
        
        G_codons_count = {x : 0 for x in aa_code['G']}
        codons_seq = []
        codon_indexes = []

        for i in range(0, seqlength, 3):
            if seq[i:i+3] in G_codons_count:
                G_codons_count[seq[i:i+3]] += 1
                codons_seq.append(seq[i:i+3] + '_' + str(G_codons_count[seq[i:i+3]]))
                codon_indexes.append(i)
        
        aa_seq     = [i for i in range(1, len(codons_seq) + 1)]
        
        
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        # fill density dict with density info
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            length_density_float = [float(i) for i in length_density]
            gene_density[length] = length_density
            total_density        = [x + y for x, y in itertools.izip(total_density, length_density)]
            
        gene_avgreads = float(sum(total_density)) / float(genelength)

            
        # normalize density by total gene density
        if gene_avgreads * 3 < threshold:
            excluded_genes['low_density'][0] += 1
            excluded_genes['low_density'][1].append(alias)
            continue
            
        else: 
            relative_density = {}
            
            for length in lengthindex:
                relative_density[length] = [reads / gene_avgreads for reads in gene_density[length]]
        
        # add data to dataframe
        #codon_index = 0
        for codon, aa, codon_position in itertools.izip(codons_seq, aa_seq, codon_indexes):
            
            #codon_position = (codon_index * 3)
            #aa = aa_seq[codon_index]
            
            #codon_index        += 1
            '''
            if codon_index > 50:
                continue
            '''
            codon_count[codon] += 1
            aa_count[aa]       += 1
            
            for length in lengthindex:
                for position in positionindex:
                    
                    density_position = codon_position + position
                    density          = relative_density[length][density_position]
                    
                    codon_avgplot[codon][length][position] += density
                    aa_avgplot[aa][length][position]       += density
                    
        included_genes[0] += 1
        included_genes[1].append(alias)
    
    #divide data by total instances of the codon or aa
    for codon in codon_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if codon_count[codon] == 0:
                    continue 
                else: 
                    codon_avgplot[codon][length][position] = codon_avgplot[codon][length][position] / codon_count[codon]
    
    for aa in aa_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if aa_count[aa] == 0:
                    continue
                else: 
                    aa_avgplot[aa][length][position] = aa_avgplot[aa][length][position] / aa_count[aa]
                
   #Convert data for plotting and csv
    
    codon_data     = {}
    codon_data_sum = {}
    aa_data     = {}
    aa_data_sum = {}
    
    A_site_shift = plot_upstream - A_site 
    P_site_shift = plot_upstream - P_site 
    E_site_shift = plot_upstream - E_site
    one_site_shift = plot_upstream - one_site 
    two_site_shift = plot_upstream - two_site 
    mone_site_shift = plot_upstream - mone_site 
    mtwo_site_shift = plot_upstream - mtwo_site 
    
    for codon in codon_avgplot.keys():
        df = pd.DataFrame(codon_avgplot[codon])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        codon_data[codon]     = df
        codon_data_sum[codon] = df.sum(0)
        
        codon_score[codon][0] = sum(codon_data_sum[codon][A_site_shift: A_site_shift + 3:1])/3
        codon_score[codon][1] = sum(codon_data_sum[codon][P_site_shift: P_site_shift + 3:1])/3                            
        codon_score[codon][2] = sum(codon_data_sum[codon][E_site_shift: E_site_shift + 3:1])/3
        codon_score[codon][3] = sum(codon_data_sum[codon][one_site_shift: one_site_shift + 3:1])/3
        codon_score[codon][4] = sum(codon_data_sum[codon][two_site_shift: two_site_shift + 3:1])/3                            
        codon_score[codon][5] = sum(codon_data_sum[codon][mone_site_shift: mone_site_shift + 3:1])/3
        codon_score[codon][6] = sum(codon_data_sum[codon][mtwo_site_shift: mtwo_site_shift + 3:1])/3
        
        codon_list.append(codon)
        codon_A_list.append(codon_score[codon][0])
        codon_P_list.append(codon_score[codon][1])
        codon_E_list.append(codon_score[codon][2])
        codon_1_list.append(codon_score[codon][3])
        codon_2_list.append(codon_score[codon][4])
        codon_m1_list.append(codon_score[codon][5])
        codon_m2_list.append(codon_score[codon][6])

    for aa in aa_avgplot.keys():
        df = pd.DataFrame(aa_avgplot[aa])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        aa_data[aa]     = df
        aa_data_sum[aa] = df.sum(0) 
                                    
        aa_score[aa][0] = sum(aa_data_sum[aa][A_site_shift: A_site_shift + 3:1])/3    
        aa_score[aa][1] = sum(aa_data_sum[aa][P_site_shift: P_site_shift + 3:1])/3
        aa_score[aa][2] = sum(aa_data_sum[aa][E_site_shift: E_site_shift + 3:1])/3
        aa_score[aa][3] = sum(aa_data_sum[aa][one_site_shift: one_site_shift + 3:1])/3    
        aa_score[aa][4] = sum(aa_data_sum[aa][two_site_shift: two_site_shift + 3:1])/3
        aa_score[aa][5] = sum(aa_data_sum[aa][mone_site_shift: mone_site_shift + 3:1])/3 
        aa_score[aa][6] = sum(aa_data_sum[aa][mtwo_site_shift: mtwo_site_shift + 3:1])/3 

        aa_list.append(aa)
        aa_A_list.append(aa_score[aa][0])
        aa_P_list.append(aa_score[aa][1])
        aa_E_list.append(aa_score[aa][2])
        aa_1_list.append(aa_score[aa][3])
        aa_2_list.append(aa_score[aa][4])
        aa_m1_list.append(aa_score[aa][5])
        aa_m2_list.append(aa_score[aa][6])
        
    codon_score_df['Codon']  = codon_list
    codon_score_df['A_site'] = codon_A_list
    codon_score_df['P_site'] = codon_P_list
    codon_score_df['E_site'] = codon_E_list
    codon_score_df['1_site'] = codon_1_list
    codon_score_df['2_site'] = codon_2_list
    codon_score_df['-1_site'] = codon_m1_list
    codon_score_df['-2_site'] = codon_m2_list

    aa_score_df['Amino Acid'] = aa_list
    aa_score_df['A_site']     = aa_A_list
    aa_score_df['P_site']     = aa_P_list
    aa_score_df['E_site']     = aa_E_list
    aa_score_df['1_site']     = aa_1_list
    aa_score_df['2_site']     = aa_2_list
    aa_score_df['-1_site']     = aa_m1_list
    aa_score_df['-2_site']     = aa_m2_list

    codon_df = pd.DataFrame(codon_score_df)
    codon_df.to_csv(path_pausescore + 'codon_scores.csv')
    
    aa_df = pd.DataFrame(aa_score_df)
    aa_df.to_csv(path_pausescore + 'aa_scores.csv')
    
    ribo_util.makePickle(codon_data,     path_pausescore +'codon_HM_data'  + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(codon_data_sum, path_pausescore +'codon_plot_data'+ name_settings , protocol=pickle.HIGHEST_PROTOCOL)  
    ribo_util.makePickle(codon_score_df, path_pausescore +'codon_scores'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_data,        path_pausescore +'aa_HM_data'     + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(aa_data_sum,    path_pausescore +'aa_plot_data'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_score_df,    path_pausescore +'aa_scores'      + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    
    #print excluded_genes['short'][0]    
    #print excluded_genes['low_density'][0]
    #print excluded_genes['not_divisible'][0]
    #print included_genes[0]
    return


def run_avg_dist(fname, settings, plus, minus, gff, path_pausescore):
    
    '''define variables'''
    
    minlength   = settings['minlength']
    maxlength   = settings['maxlength']
    alignment   = settings['alignment']
    threshold   = settings['threshold']


    lengthindex = range(minlength, maxlength + 1)
    

    frameshift      = settings['frameshift']

    plot_upstream   = settings['plot_upstream'] / 3 * 3        #change window to interval of 3
    plot_downstream = settings['plot_downstream'] / 3 * 3
    
    next_codon = settings['next_codon']
    
    
    # define plot length
    plotlength = plot_upstream + plot_downstream + 1
    
    # load density files
    density_plus  = plus
    density_minus = minus 
    
    # load annotation
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    seq_list    = gff_dict['Sequence'] 
 
    gff_extra  = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    start_trim = settings['start_trim'] / 3 * 3
    stop_trim  = settings['stop_trim'] / 3 * 3   


    '''genes in data''' 
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']         = [0, []]
    excluded_genes['low_density']   = [0, []]
    excluded_genes['not_divisible'] = [0, []]
    # count included genes
    included_genes = [0, []]
    
     
    '''iterate through every annotated gene to get codon density info:'''
    
    avgreads_list = []

    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,stop_list, strand_list, seq_list):

        # First, define density without start and stop peaks:
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim + frameshift
            density_stop  = stop  - stop_trim + frameshift
            
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim - frameshift
            density_stop  = stop  + stop_trim - 3 + frameshift
            
            period = -1
        
        # GFF seq has 50 extra nucleotides, so remove:
        # Also remove several codons from start and stop positions:
        
        codon_seq_start = gff_extra + start_trim + plot_upstream + frameshift
        codon_seq_stop  = -gff_extra - stop_trim - plot_downstream + frameshift
        
        seq       = sequence[codon_seq_start : codon_seq_stop]
        seqlength = len(seq)
        
        # exclude genes that are not divisable by 3
        if seqlength % 3 != 0:
            excluded_genes['not_divisible'][0] += 1
            excluded_genes['not_divisible'][1].append(alias)
            continue
        
        # exclude genes shorter than plot
        if seqlength < plotlength + 1:
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        elif codon_seq_start - codon_seq_stop > len(sequence):
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
             
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        # fill density dict with density info
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            length_density_float = [float(i) for i in length_density]
            gene_density[length] = length_density
            total_density        = [x + y for x, y in itertools.izip(total_density, length_density)]
            
        gene_avgreads = float(sum(total_density)) / float(genelength)

        avgreads_list.append(gene_avgreads) 

    return avgreads_list

def avg_dist(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    multi     = inputs['multiprocess']
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started pause score analysis at " + str(datetime.now())

    avg_dict = {}

    for fname in files:
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/'
        if not os.path.exists(path_pausescore):
            os.makedirs(path_pausescore)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
       
        avg_dict[fname.split('.')[0]] = run_avg_dist(fname, settings, plus, minus, gff_dict, path_pausescore)
        
    print "Finished pause score analysis at " + str(datetime.now())
    
    return avg_dict

def run_pausescore_upper_limit(fname, settings, plus, minus, gff, path_pausescore):
    
    '''define variables'''
    
    minlength   = settings['minlength']
    maxlength   = settings['maxlength']
    alignment   = settings['alignment']
    threshold   = settings['threshold']
    upper_limit = settings['upper_limit']

    lengthindex = range(minlength, maxlength + 1)
    
    A_site = settings['A_site shift']
    P_site = A_site - 3
    E_site = A_site - 6
    two_site = A_site + 6
    one_site = A_site + 3
    mone_site = A_site - 9
    mtwo_site = A_site - 12
    
    
    frameshift      = settings['frameshift']

    plot_upstream   = settings['plot_upstream'] / 3 * 3        #change window to interval of 3
    plot_downstream = settings['plot_downstream'] / 3 * 3
    
    next_codon = settings['next_codon']
    
    
    # define plot length
    plotlength = plot_upstream + plot_downstream + 1
    positionindex = range(0, plotlength)
    
    # load density files
    density_plus  = plus
    density_minus = minus 
    
    # load annotation
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    seq_list    = gff_dict['Sequence'] 
 
    gff_extra  = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    start_trim = settings['start_trim'] / 3 * 3
    stop_trim  = settings['stop_trim'] / 3 * 3   
    
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    plot_upstream_1   = str(plot_upstream)  +'_'
    plot_downstream_1 = str(plot_downstream)+'_'
    start_trim_1      = str(start_trim)     +'_'
    stop_trim_1       = str(stop_trim)      +'_'
    frameshift_1      = str(frameshift)     +'_'
    a_site_1          = str(A_site)         +'_'
    
    name_settings = minlength_1+maxlength_1+plot_upstream_1+plot_downstream_1
    name_settings += start_trim_1+stop_trim_1+frameshift_1
        
    # import genetic code
    aa_code, codon_code = ribo_util.get_genetic_code()
    
    
    '''output data structure'''
    
    # aa/codon avgplot      = { aa: {length: [values]}} 
    # aa/codon count        = { aa: N}
    # aa/codon score        = { aa: [aa_list], A site score: [A_list]...
                                
    aa_avgplot  = {}
    aa_count    = {}
    aa_score    = {}
    aa_score_df = {}
    aa_list     = []
    aa_A_list   = []
    aa_P_list   = []
    aa_E_list   = []
    aa_1_list   = []
    aa_2_list   = []
    aa_m1_list   = []
    aa_m2_list   = []
    
    codon_avgplot  = {}
    codon_count    = {}
    codon_score    = {}
    codon_score_df = {}
    codon_list     = []
    codon_A_list   = []
    codon_P_list   = []
    codon_E_list   = []
    codon_1_list   = []
    codon_2_list   = []
    codon_m1_list   = []
    codon_m2_list   = []
    
    # create empty datastructures to store density info:
    for aa in aa_code.keys():
        aa_avgplot[aa] = {length : [0]*(plotlength) for length in lengthindex}
        aa_count[aa]   = 0
        aa_score[aa]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
        
    for codon in codon_code.keys():
        codon_avgplot[codon] = {length : [0]*(plotlength) for length in lengthindex}
        codon_count[codon]   = 0
        codon_score[codon]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
    
    '''genes in data''' 
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']         = [0, []]
    excluded_genes['low_density']   = [0, []]
    excluded_genes['not_divisible'] = [0, []]
    # count included genes
    included_genes = [0, []]
    
    
    
    #list_location = '/Volumes/HDD/Ribo_seq/libraries/analysis/reference_information/Workbook1.csv'
    
    #infile = pd.read_csv(list_location)
    #TE_list = infile['yeaR'].tolist()
        
    '''iterate through every annotated gene to get codon density info:'''
    
    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,stop_list, strand_list, seq_list):
        
        #if not alias in TE_list:
            #continue
        
        
        ''' define start and stop positions for codons to analyze:
        
        # = codon 
        #################################### = GFF sequence  (50 extra nt)
           ##############################    = AUG to UGA
             ##########################      = density to analyze : remove start and stop peaks
                 ##################          = codons to analyze : remove plot window
        
        '''
        
        # First, define density without start and stop peaks:
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim + frameshift
            density_stop  = stop  - stop_trim + frameshift
            
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim - frameshift
            density_stop  = stop  + stop_trim - 3 + frameshift
            
            period = -1
        
        # GFF seq has 50 extra nucleotides, so remove:
        # Also remove several codons from start and stop positions:
        
        codon_seq_start = gff_extra + start_trim + plot_upstream + frameshift
        codon_seq_stop  = -gff_extra - stop_trim - plot_downstream + frameshift
            
        seq       = sequence[codon_seq_start : codon_seq_stop]
        seqlength = len(seq)
        
        # exclude genes that are not divisable by 3
        if seqlength % 3 != 0:
            excluded_genes['not_divisible'][0] += 1
            excluded_genes['not_divisible'][1].append(alias)
            continue
        
        # exclude genes shorter than plot
        if seqlength < plotlength + 1:
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        elif codon_seq_start - codon_seq_stop > len(sequence):
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
    
        # make a list of codons in the sequence:
        codons_seq = [seq[i:i+3] for i in range(0, seqlength, 3)]
        aa_seq     = [codon_code[codon] for codon in codons_seq]
        
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        # fill density dict with density info
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            length_density_float = [float(i) for i in length_density]
            gene_density[length] = length_density
            total_density        = [x + y for x, y in itertools.izip(total_density, length_density)]
            
        gene_avgreads = float(sum(total_density)) / float(genelength)

            
        # normalize density by total gene density
        if (gene_avgreads * 3 < threshold) or (gene_avgreads * 3 > upper_limit):
            excluded_genes['low_density'][0] += 1
            excluded_genes['low_density'][1].append(alias)
            continue
            
        else: 
            relative_density = {}
            
            for length in lengthindex:
                relative_density[length] = [reads / gene_avgreads for reads in gene_density[length]]
        
        # add data to dataframe
        codon_index = 0
        for codon in codons_seq:
            #if codon == 'TAA':
            #    print "kekekeek"
            codon_position = (codon_index * 3)
            aa = aa_seq[codon_index]
            if aa == '_':
                print 'Stop codon found in'
                print alias
            codon_index        += 1
            
            '''if codon_index > 50:
                continue
            '''
            codon_count[codon] += 1
            aa_count[aa]       += 1
            
            for length in lengthindex:
                for position in positionindex:
                    
                    density_position = codon_position + position
                    density          = relative_density[length][density_position]
                    
                    codon_avgplot[codon][length][position] += density
                    aa_avgplot[aa][length][position]       += density
                    
        included_genes[0] += 1
        included_genes[1].append(alias)
    
    #divide data by total instances of the codon or aa
    for codon in codon_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if codon_count[codon] == 0:
                    continue 
                else: 
                    codon_avgplot[codon][length][position] = codon_avgplot[codon][length][position] / codon_count[codon]
    
    for aa in aa_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if aa_count[aa] == 0:
                    continue
                else: 
                    aa_avgplot[aa][length][position] = aa_avgplot[aa][length][position] / aa_count[aa]
                
   #Convert data for plotting and csv
    
    codon_data     = {}
    codon_data_sum = {}
    aa_data     = {}
    aa_data_sum = {}
    
    A_site_shift = plot_upstream - A_site 
    P_site_shift = plot_upstream - P_site 
    E_site_shift = plot_upstream - E_site
    one_site_shift = plot_upstream - one_site 
    two_site_shift = plot_upstream - two_site 
    mone_site_shift = plot_upstream - mone_site 
    mtwo_site_shift = plot_upstream - mtwo_site 
    
    for codon in codon_avgplot.keys():
        df = pd.DataFrame(codon_avgplot[codon])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        codon_data[codon]     = df
        codon_data_sum[codon] = df.sum(0)
        
        codon_score[codon][0] = sum(codon_data_sum[codon][A_site_shift: A_site_shift + 3:1])/3
        codon_score[codon][1] = sum(codon_data_sum[codon][P_site_shift: P_site_shift + 3:1])/3                            
        codon_score[codon][2] = sum(codon_data_sum[codon][E_site_shift: E_site_shift + 3:1])/3
        codon_score[codon][3] = sum(codon_data_sum[codon][one_site_shift: one_site_shift + 3:1])/3
        codon_score[codon][4] = sum(codon_data_sum[codon][two_site_shift: two_site_shift + 3:1])/3                            
        codon_score[codon][5] = sum(codon_data_sum[codon][mone_site_shift: mone_site_shift + 3:1])/3
        codon_score[codon][6] = sum(codon_data_sum[codon][mtwo_site_shift: mtwo_site_shift + 3:1])/3
        
        codon_list.append(codon)
        codon_A_list.append(codon_score[codon][0])
        codon_P_list.append(codon_score[codon][1])
        codon_E_list.append(codon_score[codon][2])
        codon_1_list.append(codon_score[codon][3])
        codon_2_list.append(codon_score[codon][4])
        codon_m1_list.append(codon_score[codon][5])
        codon_m2_list.append(codon_score[codon][6])

    for aa in aa_avgplot.keys():
        df = pd.DataFrame(aa_avgplot[aa])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        aa_data[aa]     = df
        aa_data_sum[aa] = df.sum(0) 
                                    
        aa_score[aa][0] = sum(aa_data_sum[aa][A_site_shift: A_site_shift + 3:1])/3    
        aa_score[aa][1] = sum(aa_data_sum[aa][P_site_shift: P_site_shift + 3:1])/3
        aa_score[aa][2] = sum(aa_data_sum[aa][E_site_shift: E_site_shift + 3:1])/3
        aa_score[aa][3] = sum(aa_data_sum[aa][one_site_shift: one_site_shift + 3:1])/3    
        aa_score[aa][4] = sum(aa_data_sum[aa][two_site_shift: two_site_shift + 3:1])/3
        aa_score[aa][5] = sum(aa_data_sum[aa][mone_site_shift: mone_site_shift + 3:1])/3 
        aa_score[aa][6] = sum(aa_data_sum[aa][mtwo_site_shift: mtwo_site_shift + 3:1])/3 

        aa_list.append(aa)
        aa_A_list.append(aa_score[aa][0])
        aa_P_list.append(aa_score[aa][1])
        aa_E_list.append(aa_score[aa][2])
        aa_1_list.append(aa_score[aa][3])
        aa_2_list.append(aa_score[aa][4])
        aa_m1_list.append(aa_score[aa][5])
        aa_m2_list.append(aa_score[aa][6])
        
    codon_score_df['Codon']  = codon_list
    codon_score_df['A_site'] = codon_A_list
    codon_score_df['P_site'] = codon_P_list
    codon_score_df['E_site'] = codon_E_list
    codon_score_df['1_site'] = codon_1_list
    codon_score_df['2_site'] = codon_2_list
    codon_score_df['-1_site'] = codon_m1_list
    codon_score_df['-2_site'] = codon_m2_list

    aa_score_df['Amino Acid'] = aa_list
    aa_score_df['A_site']     = aa_A_list
    aa_score_df['P_site']     = aa_P_list
    aa_score_df['E_site']     = aa_E_list
    aa_score_df['1_site']     = aa_1_list
    aa_score_df['2_site']     = aa_2_list
    aa_score_df['-1_site']     = aa_m1_list
    aa_score_df['-2_site']     = aa_m2_list

    codon_df = pd.DataFrame(codon_score_df)
    codon_df.to_csv(path_pausescore + 'codon_scores.csv')
    
    aa_df = pd.DataFrame(aa_score_df)
    aa_df.to_csv(path_pausescore + 'aa_scores.csv')
    
    ribo_util.makePickle(codon_data,     path_pausescore +'codon_HM_data'  + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(codon_data_sum, path_pausescore +'codon_plot_data'+ name_settings , protocol=pickle.HIGHEST_PROTOCOL)  
    ribo_util.makePickle(codon_score_df, path_pausescore +'codon_scores'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_data,        path_pausescore +'aa_HM_data'     + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(aa_data_sum,    path_pausescore +'aa_plot_data'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_score_df,    path_pausescore +'aa_scores'      + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    
    #print excluded_genes['short'][0]    
    #print excluded_genes['low_density'][0]
    #print excluded_genes['not_divisible'][0]
    #print included_genes[0]
    return 


def pausescore_upper_limit(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    multi     = inputs['multiprocess']
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started pause score analysis at " + str(datetime.now())

    for fname in files:
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/'
        if not os.path.exists(path_pausescore):
            os.makedirs(path_pausescore)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
       
        if not multi == 'yes':
            run_pausescore_upper_limit(fname, settings, plus, minus, gff_dict, path_pausescore)
        else:     
            argument = [fname, settings, plus, minus, gff_dict, path_pausescore]
            arguments.append(argument)
    
    if multi == 'yes':
        ribo_util.multiprocess(run_pausescore_upper_limit, arguments, threads)
    
    print "Finished pause score analysis at " + str(datetime.now())
    
    return

def run_pausescore_inter_GG(fname, settings, plus, minus, gff, path_pausescore):
    
    '''define variables'''
    
    minlength   = settings['minlength']
    maxlength   = settings['maxlength']
    alignment   = settings['alignment']
    threshold   = settings['threshold']
    #upper_limit = settings['upper_limit']
    interval_threshold = settings['interval_threshold']


    lengthindex = range(minlength, maxlength + 1)
    
    A_site = settings['A_site shift']
    P_site = A_site - 3
    E_site = A_site - 6
    two_site = A_site + 6
    one_site = A_site + 3
    mone_site = A_site - 9
    mtwo_site = A_site - 12
    
    
    frameshift      = settings['frameshift']

    plot_upstream   = settings['plot_upstream'] / 3 * 3        #change window to interval of 3
    plot_downstream = settings['plot_downstream'] / 3 * 3
    
    next_codon = settings['next_codon']
    
    
    # define plot length
    plotlength = plot_upstream + plot_downstream + 1
    positionindex = range(0, plotlength)
    
    # load density files
    density_plus  = plus
    density_minus = minus 
    
    # load annotation
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    seq_list    = gff_dict['Sequence'] 
 
    gff_extra  = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    start_trim = settings['start_trim'] / 3 * 3
    stop_trim  = settings['stop_trim'] / 3 * 3   
    
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    plot_upstream_1   = str(plot_upstream)  +'_'
    plot_downstream_1 = str(plot_downstream)+'_'
    start_trim_1      = str(start_trim)     +'_'
    stop_trim_1       = str(stop_trim)      +'_'
    frameshift_1      = str(frameshift)     +'_'
    a_site_1          = str(A_site)         +'_'
    
    name_settings = minlength_1+maxlength_1+plot_upstream_1+plot_downstream_1
    name_settings += start_trim_1+stop_trim_1+frameshift_1
        
    # import genetic code
    aa_code, codon_code = ribo_util.get_genetic_code()
    
    
    '''output data structure'''
    
    # aa/codon avgplot      = { aa: {length: [values]}} 
    # aa/codon count        = { aa: N}
    # aa/codon score        = { aa: [aa_list], A site score: [A_list]...
                                
    aa_avgplot  = {}
    aa_count    = {}
    aa_score    = {}
    aa_score_df = {}
    aa_list     = []
    aa_A_list   = []
    aa_P_list   = []
    aa_E_list   = []
    aa_1_list   = []
    aa_2_list   = []
    aa_m1_list   = []
    aa_m2_list   = []
    
    codon_avgplot  = {}
    codon_count    = {}
    codon_score    = {}
    codon_score_df = {}
    codon_list     = []
    codon_A_list   = []
    codon_P_list   = []
    codon_E_list   = []
    codon_1_list   = []
    codon_2_list   = []
    codon_m1_list   = []
    codon_m2_list   = []
    
    # create empty datastructures to store density info:

    for i in range(1, 214):
        aa_avgplot[i] = {length : [0]*(plotlength) for length in lengthindex}
        aa_count[i]   = 0
        aa_score[i]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
        
    for codon in aa_code['G']:
        for i in range(1, 214):
            codon_avgplot[codon + '_' + str(i)] = {length : 0 for length in lengthindex}
            codon_count[codon + '_' + str(i)]   = 0
            codon_score[codon + '_' + str(i)]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 

    
    '''genes in data''' 
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']         = [0, []]
    excluded_genes['low_density']   = [0, []]
    excluded_genes['not_divisible'] = [0, []]
    excluded_genes['short_intervals'] = [0, []]
    # count included genes
    included_genes = [0, []]
    
    
    
    #list_location = '/Volumes/HDD/Ribo_seq/libraries/analysis/reference_information/Workbook1.csv'
    
    #infile = pd.read_csv(list_location)
    #TE_list = infile['yeaR'].tolist()
        
    '''iterate through every annotated gene to get codon density info:'''
    
    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,stop_list, strand_list, seq_list):
        
        #if not alias in TE_list:
            #continue
        
        
        ''' define start and stop positions for codons to analyze:
        
        # = codon 
        #################################### = GFF sequence  (50 extra nt)
           ##############################    = AUG to UGA
             ##########################      = density to analyze : remove start and stop peaks
                 ##################          = codons to analyze : remove plot window
        
        '''
        
        # First, define density without start and stop peaks:
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim + frameshift
            density_stop  = stop  - stop_trim + frameshift
            
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim - frameshift
            density_stop  = stop  + stop_trim - 3 + frameshift
            
            period = -1
        
        # GFF seq has 50 extra nucleotides, so remove:
        # Also remove several codons from start and stop positions:
        
        codon_seq_start = gff_extra + start_trim + plot_upstream + frameshift
        codon_seq_stop  = -gff_extra - stop_trim - plot_downstream + frameshift
            
        seq       = sequence[codon_seq_start : codon_seq_stop]
        seqlength = len(seq)
        
        bp_seq = Seq(seq)
        translation = str(bp_seq.translate())
        intervals = translation.split('G')
        
        for interval in intervals:
            if len(interval) < interval_threshold:
                excluded_genes['short_intervals'][0] += 1
                excluded_genes['short_intervals'][1].append(alias)
                continue

        # exclude genes that are not divisable by 3
        if seqlength % 3 != 0:
            excluded_genes['not_divisible'][0] += 1
            excluded_genes['not_divisible'][1].append(alias)
            continue
        
        # exclude genes shorter than plot
        if seqlength < plotlength + 1:
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        elif codon_seq_start - codon_seq_stop > len(sequence):
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
    
        # make a list of codons in the sequence:

        G_codons_count = {x : 0 for x in aa_code['G']}
        codons_seq = []
        codon_indexes = []

        for i in range(0, seqlength, 3):
            if seq[i:i+3] in G_codons_count:
                G_codons_count[seq[i:i+3]] += 1
                codons_seq.append(seq[i:i+3] + '_' + str(G_codons_count[seq[i:i+3]]))
                codon_indexes.append(i)
        
        aa_seq     = [i for i in range(1, len(codons_seq) + 1)]
        

        
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        # fill density dict with density info
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            length_density_float = [float(i) for i in length_density]
            gene_density[length] = length_density
            total_density        = [x + y for x, y in itertools.izip(total_density, length_density)]
            
        gene_avgreads = float(sum(total_density)) / float(genelength)

            
        # normalize density by total gene density
        if gene_avgreads * 3 < threshold: #or (gene_avgreads * 3 > upper_limit):
            excluded_genes['low_density'][0] += 1
            excluded_genes['low_density'][1].append(alias)
            continue
            
        else: 
            relative_density = {}
            
            for length in lengthindex:
                relative_density[length] = [reads / gene_avgreads for reads in gene_density[length]]
        
        # add data to dataframe
        
        prev_position = 0
        for codon, aa, codon_position in itertools.izip(codons_seq, aa_seq, codon_indexes):
        
            # doesn't count last intervals!
            codon_count[codon] += 1
            aa_count[aa]       += 1
            
            for length in lengthindex:
                
                    density = sum(relative_density[length][prev_position + 35:codon_position - 35])
                    
                    codon_avgplot[codon][length] += density/float(codon_position - prev_position - 70)
                    
                    

                    prev_position = codon_position

        included_genes[0] += 1
        included_genes[1].append(alias)
    
    intervals_avg = {}
    for codon in codon_avgplot:
        intervals_avg[codon] = 0
        if codon_count[codon] == 0:
            continue 
        else:
            for length in codon_avgplot[codon]:
                intervals_avg[codon] += codon_avgplot[codon][length]  / codon_count[codon]

  

    return intervals_avg


def pausescore_inter_GG(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    multi     = inputs['multiprocess']
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started pause score analysis at " + str(datetime.now())

    avg_dict = {}

    for fname in files:
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/'
        if not os.path.exists(path_pausescore):
            os.makedirs(path_pausescore)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
       
        avg_dict[fname.split('.')[0]] = run_pausescore_inter_GG(fname, settings, plus, minus, gff_dict, path_pausescore)
        
    print "Finished pause score analysis at " + str(datetime.now())
    
    return avg_dict



def run_pausescore_no_pseudo(fname, settings, plus, minus, gff, path_pausescore):
    
    '''define variables'''
    
    minlength   = settings['minlength']
    maxlength   = settings['maxlength']
    alignment   = settings['alignment']
    threshold   = settings['threshold']


    lengthindex = range(minlength, maxlength + 1)
    
    A_site = settings['A_site shift']
    P_site = A_site - 3
    E_site = A_site - 6
    two_site = A_site + 6
    one_site = A_site + 3
    mone_site = A_site - 9
    mtwo_site = A_site - 12
    
    
    frameshift      = settings['frameshift']

    plot_upstream   = settings['plot_upstream'] / 3 * 3        #change window to interval of 3
    plot_downstream = settings['plot_downstream'] / 3 * 3
    
    next_codon = settings['next_codon']
    
    
    # define plot length
    plotlength = plot_upstream + plot_downstream + 1
    positionindex = range(0, plotlength)
    
    # load density files
    density_plus  = plus
    density_minus = minus 
    
    # load annotation
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    seq_list    = gff_dict['Sequence'] 
 
    gff_extra  = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    start_trim = settings['start_trim'] / 3 * 3
    stop_trim  = settings['stop_trim'] / 3 * 3   
    
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    plot_upstream_1   = str(plot_upstream)  +'_'
    plot_downstream_1 = str(plot_downstream)+'_'
    start_trim_1      = str(start_trim)     +'_'
    stop_trim_1       = str(stop_trim)      +'_'
    frameshift_1      = str(frameshift)     +'_'
    a_site_1          = str(A_site)         +'_'
    
    name_settings = minlength_1+maxlength_1+plot_upstream_1+plot_downstream_1
    name_settings += start_trim_1+stop_trim_1+frameshift_1
        
    # import genetic code
    aa_code, codon_code = ribo_util.get_genetic_code()
    
    
    '''output data structure'''
    
    # aa/codon avgplot      = { aa: {length: [values]}} 
    # aa/codon count        = { aa: N}
    # aa/codon score        = { aa: [aa_list], A site score: [A_list]...
                                
    aa_avgplot  = {}
    aa_count    = {}
    aa_score    = {}
    aa_score_df = {}
    aa_list     = []
    aa_A_list   = []
    aa_P_list   = []
    aa_E_list   = []
    aa_1_list   = []
    aa_2_list   = []
    aa_m1_list   = []
    aa_m2_list   = []
    
    codon_avgplot  = {}
    codon_count    = {}
    codon_score    = {}
    codon_score_df = {}
    codon_list     = []
    codon_A_list   = []
    codon_P_list   = []
    codon_E_list   = []
    codon_1_list   = []
    codon_2_list   = []
    codon_m1_list   = []
    codon_m2_list   = []
    
    # create empty datastructures to store density info:
    for aa in aa_code.keys():
        aa_avgplot[aa] = {length : [0]*(plotlength) for length in lengthindex}
        aa_count[aa]   = 0
        aa_score[aa]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
        
    for codon in codon_code.keys():
        codon_avgplot[codon] = {length : [0]*(plotlength) for length in lengthindex}
        codon_count[codon]   = 0
        codon_score[codon]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
    
    '''genes in data''' 
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']         = [0, []]
    excluded_genes['low_density']   = [0, []]
    excluded_genes['not_divisible'] = [0, []]
    # count included genes
    included_genes = [0, []]
    
    
    
    list_location = 'pseudogenes.csv'
    
    infile = pd.read_csv(list_location)
    pseudo_list = infile['Alias'].tolist()
        
    '''iterate through every annotated gene to get codon density info:'''
    
    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,stop_list, strand_list, seq_list):
        
        if alias in pseudo_list:
            continue
        
        
        ''' define start and stop positions for codons to analyze:
        
        # = codon 
        #################################### = GFF sequence  (50 extra nt)
           ##############################    = AUG to UGA
             ##########################      = density to analyze : remove start and stop peaks
                 ##################          = codons to analyze : remove plot window
        
        '''
        
        # First, define density without start and stop peaks:
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim + frameshift
            density_stop  = stop  - stop_trim + frameshift
            
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim - frameshift
            density_stop  = stop  + stop_trim - 3 + frameshift
            
            period = -1
        
        # GFF seq has 50 extra nucleotides, so remove:
        # Also remove several codons from start and stop positions:
        
        codon_seq_start = gff_extra + start_trim + plot_upstream + frameshift
        codon_seq_stop  = -gff_extra - stop_trim - plot_downstream + frameshift
            
        seq       = sequence[codon_seq_start : codon_seq_stop]
        seqlength = len(seq)
        
        # exclude genes that are not divisable by 3
        if seqlength % 3 != 0:
            excluded_genes['not_divisible'][0] += 1
            excluded_genes['not_divisible'][1].append(alias)
            continue
        
        # exclude genes shorter than plot
        if seqlength < plotlength + 1:
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        elif codon_seq_start - codon_seq_stop > len(sequence):
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
    
        # make a list of codons in the sequence:
        codons_seq = [seq[i:i+3] for i in range(0, seqlength, 3)]
        aa_seq     = [codon_code[codon] for codon in codons_seq]
        
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        # fill density dict with density info
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            length_density_float = [float(i) for i in length_density]
            gene_density[length] = length_density
            total_density        = [x + y for x, y in itertools.izip(total_density, length_density)]
            
        gene_avgreads = float(sum(total_density)) / float(genelength)

            
        # normalize density by total gene density
        if gene_avgreads * 3 < threshold:
            excluded_genes['low_density'][0] += 1
            excluded_genes['low_density'][1].append(alias)
            continue
            
        else: 
            relative_density = {}
            
            for length in lengthindex:
                relative_density[length] = [reads / gene_avgreads for reads in gene_density[length]]
        
        # add data to dataframe
        codon_index = 0
        for codon in codons_seq:
            #if codon == 'TAA':
            #    print "kekekeek"
            codon_position = (codon_index * 3)
            aa = aa_seq[codon_index]
            if aa == '_':
                print 'Stop codon found in'
                print alias
            codon_index        += 1
            
            '''if codon_index > 50:
                continue
            '''
            codon_count[codon] += 1
            aa_count[aa]       += 1
            
            for length in lengthindex:
                for position in positionindex:
                    
                    density_position = codon_position + position
                    density          = relative_density[length][density_position]
                    
                    codon_avgplot[codon][length][position] += density
                    aa_avgplot[aa][length][position]       += density
                    
        included_genes[0] += 1
        included_genes[1].append(alias)
    
    #divide data by total instances of the codon or aa
    for codon in codon_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if codon_count[codon] == 0:
                    continue 
                else: 
                    codon_avgplot[codon][length][position] = codon_avgplot[codon][length][position] / codon_count[codon]
    
    for aa in aa_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if aa_count[aa] == 0:
                    continue
                else: 
                    aa_avgplot[aa][length][position] = aa_avgplot[aa][length][position] / aa_count[aa]
                
   #Convert data for plotting and csv
    
    codon_data     = {}
    codon_data_sum = {}
    aa_data     = {}
    aa_data_sum = {}
    
    A_site_shift = plot_upstream - A_site 
    P_site_shift = plot_upstream - P_site 
    E_site_shift = plot_upstream - E_site
    one_site_shift = plot_upstream - one_site 
    two_site_shift = plot_upstream - two_site 
    mone_site_shift = plot_upstream - mone_site 
    mtwo_site_shift = plot_upstream - mtwo_site 
    
    for codon in codon_avgplot.keys():
        df = pd.DataFrame(codon_avgplot[codon])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        codon_data[codon]     = df
        codon_data_sum[codon] = df.sum(0)
        
        codon_score[codon][0] = sum(codon_data_sum[codon][A_site_shift: A_site_shift + 3:1])/3
        codon_score[codon][1] = sum(codon_data_sum[codon][P_site_shift: P_site_shift + 3:1])/3                            
        codon_score[codon][2] = sum(codon_data_sum[codon][E_site_shift: E_site_shift + 3:1])/3
        codon_score[codon][3] = sum(codon_data_sum[codon][one_site_shift: one_site_shift + 3:1])/3
        codon_score[codon][4] = sum(codon_data_sum[codon][two_site_shift: two_site_shift + 3:1])/3                            
        codon_score[codon][5] = sum(codon_data_sum[codon][mone_site_shift: mone_site_shift + 3:1])/3
        codon_score[codon][6] = sum(codon_data_sum[codon][mtwo_site_shift: mtwo_site_shift + 3:1])/3
        
        codon_list.append(codon)
        codon_A_list.append(codon_score[codon][0])
        codon_P_list.append(codon_score[codon][1])
        codon_E_list.append(codon_score[codon][2])
        codon_1_list.append(codon_score[codon][3])
        codon_2_list.append(codon_score[codon][4])
        codon_m1_list.append(codon_score[codon][5])
        codon_m2_list.append(codon_score[codon][6])

    for aa in aa_avgplot.keys():
        df = pd.DataFrame(aa_avgplot[aa])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        aa_data[aa]     = df
        aa_data_sum[aa] = df.sum(0) 
                                    
        aa_score[aa][0] = sum(aa_data_sum[aa][A_site_shift: A_site_shift + 3:1])/3    
        aa_score[aa][1] = sum(aa_data_sum[aa][P_site_shift: P_site_shift + 3:1])/3
        aa_score[aa][2] = sum(aa_data_sum[aa][E_site_shift: E_site_shift + 3:1])/3
        aa_score[aa][3] = sum(aa_data_sum[aa][one_site_shift: one_site_shift + 3:1])/3    
        aa_score[aa][4] = sum(aa_data_sum[aa][two_site_shift: two_site_shift + 3:1])/3
        aa_score[aa][5] = sum(aa_data_sum[aa][mone_site_shift: mone_site_shift + 3:1])/3 
        aa_score[aa][6] = sum(aa_data_sum[aa][mtwo_site_shift: mtwo_site_shift + 3:1])/3 

        aa_list.append(aa)
        aa_A_list.append(aa_score[aa][0])
        aa_P_list.append(aa_score[aa][1])
        aa_E_list.append(aa_score[aa][2])
        aa_1_list.append(aa_score[aa][3])
        aa_2_list.append(aa_score[aa][4])
        aa_m1_list.append(aa_score[aa][5])
        aa_m2_list.append(aa_score[aa][6])
        
    codon_score_df['Codon']  = codon_list
    codon_score_df['A_site'] = codon_A_list
    codon_score_df['P_site'] = codon_P_list
    codon_score_df['E_site'] = codon_E_list
    codon_score_df['1_site'] = codon_1_list
    codon_score_df['2_site'] = codon_2_list
    codon_score_df['-1_site'] = codon_m1_list
    codon_score_df['-2_site'] = codon_m2_list

    aa_score_df['Amino Acid'] = aa_list
    aa_score_df['A_site']     = aa_A_list
    aa_score_df['P_site']     = aa_P_list
    aa_score_df['E_site']     = aa_E_list
    aa_score_df['1_site']     = aa_1_list
    aa_score_df['2_site']     = aa_2_list
    aa_score_df['-1_site']     = aa_m1_list
    aa_score_df['-2_site']     = aa_m2_list

    codon_df = pd.DataFrame(codon_score_df)
    codon_df.to_csv(path_pausescore + 'codon_scores.csv')
    
    aa_df = pd.DataFrame(aa_score_df)
    aa_df.to_csv(path_pausescore + 'aa_scores.csv')
    
    ribo_util.makePickle(codon_data,     path_pausescore +'codon_HM_data'  + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(codon_data_sum, path_pausescore +'codon_plot_data'+ name_settings , protocol=pickle.HIGHEST_PROTOCOL)  
    ribo_util.makePickle(codon_score_df, path_pausescore +'codon_scores'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_data,        path_pausescore +'aa_HM_data'     + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(aa_data_sum,    path_pausescore +'aa_plot_data'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_score_df,    path_pausescore +'aa_scores'      + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    
    #print excluded_genes['short'][0]    
    #print excluded_genes['low_density'][0]
    #print excluded_genes['not_divisible'][0]
    #print included_genes[0]
    return 


def pausescore_no_pseudo(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    multi     = inputs['multiprocess']
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started pause score analysis at " + str(datetime.now())

    for fname in files:
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/'
        if not os.path.exists(path_pausescore):
            os.makedirs(path_pausescore)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
       
        if not multi == 'yes':
            run_pausescore_no_pseudo(fname, settings, plus, minus, gff_dict, path_pausescore)
        else:     
            argument = [fname, settings, plus, minus, gff_dict, path_pausescore]
            arguments.append(argument)
    
    if multi == 'yes':
        ribo_util.multiprocess(run_pausescore_no_pseudo, arguments, threads)
    
    print "Finished pause score analysis at " + str(datetime.now())
    
    return

def run_pausescore_CDS_first_half(fname, settings, plus, minus, gff, path_pausescore):
    
       
    '''define variables'''
    
    minlength   = settings['minlength']
    maxlength   = settings['maxlength']
    alignment   = settings['alignment']
    threshold   = settings['threshold']


    lengthindex = range(minlength, maxlength + 1)
    
    A_site = settings['A_site shift']
    P_site = A_site - 3
    E_site = A_site - 6
    two_site = A_site + 6
    one_site = A_site + 3
    mone_site = A_site - 9
    mtwo_site = A_site - 12
    
    
    frameshift      = settings['frameshift']

    plot_upstream   = settings['plot_upstream'] / 3 * 3        #change window to interval of 3
    plot_downstream = settings['plot_downstream'] / 3 * 3
    
    next_codon = settings['next_codon']
    
    
    # define plot length
    plotlength = plot_upstream + plot_downstream + 1
    positionindex = range(0, plotlength)
    
    # load density files
    density_plus  = plus
    density_minus = minus 
    
    # load annotation
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    seq_list    = gff_dict['Sequence'] 
 
    gff_extra  = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    start_trim = settings['start_trim'] / 3 * 3
    stop_trim  = settings['stop_trim'] / 3 * 3   
    
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    plot_upstream_1   = str(plot_upstream)  +'_'
    plot_downstream_1 = str(plot_downstream)+'_'
    start_trim_1      = str(start_trim)     +'_'
    stop_trim_1       = str(stop_trim)      +'_'
    frameshift_1      = str(frameshift)     +'_'
    a_site_1          = str(A_site)         +'_'
    
    name_settings = minlength_1+maxlength_1+plot_upstream_1+plot_downstream_1
    name_settings += start_trim_1+stop_trim_1+frameshift_1
        
    # import genetic code
    aa_code, codon_code = ribo_util.get_genetic_code()
    
    
    '''output data structure'''
    
    # aa/codon avgplot      = { aa: {length: [values]}} 
    # aa/codon count        = { aa: N}
    # aa/codon score        = { aa: [aa_list], A site score: [A_list]...
                                
    aa_avgplot  = {}
    aa_count    = {}
    aa_score    = {}
    aa_score_df = {}
    aa_list     = []
    aa_A_list   = []
    aa_P_list   = []
    aa_E_list   = []
    aa_1_list   = []
    aa_2_list   = []
    aa_m1_list   = []
    aa_m2_list   = []
    
    codon_avgplot  = {}
    codon_count    = {}
    codon_score    = {}
    codon_score_df = {}
    codon_list     = []
    codon_A_list   = []
    codon_P_list   = []
    codon_E_list   = []
    codon_1_list   = []
    codon_2_list   = []
    codon_m1_list   = []
    codon_m2_list   = []
    
    # create empty datastructures to store density info:
    for aa in aa_code.keys():
        aa_avgplot[aa] = {length : [0]*(plotlength) for length in lengthindex}
        aa_count[aa]   = 0
        aa_score[aa]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
        
    for codon in codon_code.keys():
        codon_avgplot[codon] = {length : [0]*(plotlength) for length in lengthindex}
        codon_count[codon]   = 0
        codon_score[codon]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
    
    '''genes in data''' 
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']         = [0, []]
    excluded_genes['low_density']   = [0, []]
    excluded_genes['not_divisible'] = [0, []]
    # count included genes
    included_genes = [0, []]
    
    
    
    #list_location = '/Volumes/HDD/Ribo_seq/libraries/analysis/reference_information/Workbook1.csv'
    
    #infile = pd.read_csv(list_location)
    #TE_list = infile['yeaR'].tolist()
        
    '''iterate through every annotated gene to get codon density info:'''
    
    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,stop_list, strand_list, seq_list):
        
        #if not alias in TE_list:
            #continue
        
        
        ''' define start and stop positions for codons to analyze:
        
        # = codon 
        #################################### = GFF sequence  (50 extra nt)
           ##############################    = AUG to UGA
             ##########################      = density to analyze : remove start and stop peaks
                 ##################          = codons to analyze : remove plot window
        
        '''
        
        # First, define density without start and stop peaks:
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim + frameshift
            density_stop  = stop  - stop_trim + frameshift
            
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim - frameshift
            density_stop  = stop  + stop_trim - 3 + frameshift
            
            period = -1
        
        # GFF seq has 50 extra nucleotides, so remove:
        # Also remove several codons from start and stop positions:
        
        codon_seq_start = gff_extra + start_trim + plot_upstream + frameshift
        codon_seq_stop  = -gff_extra - stop_trim - plot_downstream + frameshift
            
        seq       = sequence[codon_seq_start : codon_seq_stop]
        seqlength = len(seq)
        
        # exclude genes that are not divisable by 3
        if seqlength % 3 != 0:
            excluded_genes['not_divisible'][0] += 1
            excluded_genes['not_divisible'][1].append(alias)
            continue
        
        # exclude genes shorter than plot
        if seqlength < plotlength + 1:
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        elif codon_seq_start - codon_seq_stop > len(sequence):
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
    
        # make a list of codons in the sequence:
        codons_seq = [seq[i:i+3] for i in range(0, seqlength/6*3, 3)]
        aa_seq     = [codon_code[codon] for codon in codons_seq]
        
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        # fill density dict with density info
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            length_density_float = [float(i) for i in length_density]
            gene_density[length] = length_density
            total_density        = [x + y for x, y in itertools.izip(total_density, length_density)]
            
        gene_avgreads = float(sum(total_density)) / float(genelength)

            
        # normalize density by total gene density
        if gene_avgreads * 3 < threshold:
            excluded_genes['low_density'][0] += 1
            excluded_genes['low_density'][1].append(alias)
            continue
            
        else: 
            relative_density = {}
            
            for length in lengthindex:
                relative_density[length] = [reads / gene_avgreads for reads in gene_density[length]]
        
        # add data to dataframe
        codon_index = 0
        for codon in codons_seq:
            #if codon == 'TAA':
            #    print "kekekeek"
            codon_position = (codon_index * 3)
            aa = aa_seq[codon_index]
            if aa == '_':
                print 'Stop codon found in'
                print alias
            codon_index        += 1
            
            '''if codon_index > 50:
                continue
            '''
            codon_count[codon] += 1
            aa_count[aa]       += 1
            
            for length in lengthindex:
                for position in positionindex:
                    
                    density_position = codon_position + position
                    density          = relative_density[length][density_position]
                    
                    codon_avgplot[codon][length][position] += density
                    aa_avgplot[aa][length][position]       += density
                    
        included_genes[0] += 1
        included_genes[1].append(alias)
    
    #divide data by total instances of the codon or aa
    for codon in codon_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if codon_count[codon] == 0:
                    continue 
                else: 
                    codon_avgplot[codon][length][position] = codon_avgplot[codon][length][position] / codon_count[codon]
    
    for aa in aa_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if aa_count[aa] == 0:
                    continue
                else: 
                    aa_avgplot[aa][length][position] = aa_avgplot[aa][length][position] / aa_count[aa]
                
   #Convert data for plotting and csv
    
    codon_data     = {}
    codon_data_sum = {}
    aa_data     = {}
    aa_data_sum = {}
    
    A_site_shift = plot_upstream - A_site 
    P_site_shift = plot_upstream - P_site 
    E_site_shift = plot_upstream - E_site
    one_site_shift = plot_upstream - one_site 
    two_site_shift = plot_upstream - two_site 
    mone_site_shift = plot_upstream - mone_site 
    mtwo_site_shift = plot_upstream - mtwo_site 
    
    for codon in codon_avgplot.keys():
        df = pd.DataFrame(codon_avgplot[codon])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        codon_data[codon]     = df
        codon_data_sum[codon] = df.sum(0)
        
        codon_score[codon][0] = sum(codon_data_sum[codon][A_site_shift: A_site_shift + 3:1])/3
        codon_score[codon][1] = sum(codon_data_sum[codon][P_site_shift: P_site_shift + 3:1])/3                            
        codon_score[codon][2] = sum(codon_data_sum[codon][E_site_shift: E_site_shift + 3:1])/3
        codon_score[codon][3] = sum(codon_data_sum[codon][one_site_shift: one_site_shift + 3:1])/3
        codon_score[codon][4] = sum(codon_data_sum[codon][two_site_shift: two_site_shift + 3:1])/3                            
        codon_score[codon][5] = sum(codon_data_sum[codon][mone_site_shift: mone_site_shift + 3:1])/3
        codon_score[codon][6] = sum(codon_data_sum[codon][mtwo_site_shift: mtwo_site_shift + 3:1])/3
        
        codon_list.append(codon)
        codon_A_list.append(codon_score[codon][0])
        codon_P_list.append(codon_score[codon][1])
        codon_E_list.append(codon_score[codon][2])
        codon_1_list.append(codon_score[codon][3])
        codon_2_list.append(codon_score[codon][4])
        codon_m1_list.append(codon_score[codon][5])
        codon_m2_list.append(codon_score[codon][6])

    for aa in aa_avgplot.keys():
        df = pd.DataFrame(aa_avgplot[aa])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        aa_data[aa]     = df
        aa_data_sum[aa] = df.sum(0) 
                                    
        aa_score[aa][0] = sum(aa_data_sum[aa][A_site_shift: A_site_shift + 3:1])/3    
        aa_score[aa][1] = sum(aa_data_sum[aa][P_site_shift: P_site_shift + 3:1])/3
        aa_score[aa][2] = sum(aa_data_sum[aa][E_site_shift: E_site_shift + 3:1])/3
        aa_score[aa][3] = sum(aa_data_sum[aa][one_site_shift: one_site_shift + 3:1])/3    
        aa_score[aa][4] = sum(aa_data_sum[aa][two_site_shift: two_site_shift + 3:1])/3
        aa_score[aa][5] = sum(aa_data_sum[aa][mone_site_shift: mone_site_shift + 3:1])/3 
        aa_score[aa][6] = sum(aa_data_sum[aa][mtwo_site_shift: mtwo_site_shift + 3:1])/3 

        aa_list.append(aa)
        aa_A_list.append(aa_score[aa][0])
        aa_P_list.append(aa_score[aa][1])
        aa_E_list.append(aa_score[aa][2])
        aa_1_list.append(aa_score[aa][3])
        aa_2_list.append(aa_score[aa][4])
        aa_m1_list.append(aa_score[aa][5])
        aa_m2_list.append(aa_score[aa][6])
        
    codon_score_df['Codon']  = codon_list
    codon_score_df['A_site'] = codon_A_list
    codon_score_df['P_site'] = codon_P_list
    codon_score_df['E_site'] = codon_E_list
    codon_score_df['1_site'] = codon_1_list
    codon_score_df['2_site'] = codon_2_list
    codon_score_df['-1_site'] = codon_m1_list
    codon_score_df['-2_site'] = codon_m2_list

    aa_score_df['Amino Acid'] = aa_list
    aa_score_df['A_site']     = aa_A_list
    aa_score_df['P_site']     = aa_P_list
    aa_score_df['E_site']     = aa_E_list
    aa_score_df['1_site']     = aa_1_list
    aa_score_df['2_site']     = aa_2_list
    aa_score_df['-1_site']     = aa_m1_list
    aa_score_df['-2_site']     = aa_m2_list

    codon_df = pd.DataFrame(codon_score_df)
    codon_df.to_csv(path_pausescore + 'codon_scores.csv')
    
    aa_df = pd.DataFrame(aa_score_df)
    aa_df.to_csv(path_pausescore + 'aa_scores.csv')
    
    ribo_util.makePickle(codon_data,     path_pausescore +'codon_HM_data'  + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(codon_data_sum, path_pausescore +'codon_plot_data'+ name_settings , protocol=pickle.HIGHEST_PROTOCOL)  
    ribo_util.makePickle(codon_score_df, path_pausescore +'codon_scores'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_data,        path_pausescore +'aa_HM_data'     + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(aa_data_sum,    path_pausescore +'aa_plot_data'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_score_df,    path_pausescore +'aa_scores'      + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    
    #print excluded_genes['short'][0]    
    #print excluded_genes['low_density'][0]
    #print excluded_genes['not_divisible'][0]
    #print included_genes[0]
    return 



def pausescore_first_half(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    multi     = inputs['multiprocess']
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started pause score analysis at " + str(datetime.now())

    for fname in files:
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/'
        if not os.path.exists(path_pausescore):
            os.makedirs(path_pausescore)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
       
        if not multi == 'yes':
            run_pausescore_CDS_first_half(fname, settings, plus, minus, gff_dict, path_pausescore)
        else:     
            argument = [fname, settings, plus, minus, gff_dict, path_pausescore]
            arguments.append(argument)
    
    if multi == 'yes':
        ribo_util.multiprocess(run_pausescore_CDS_first_half, arguments, threads)
    
    print "Finished pause score analysis at " + str(datetime.now())
    
    return


def run_pausescore_first_pos(fname, settings, plus, minus, gff, path_pausescore):
    
    '''define variables'''
    
    minlength   = settings['minlength']
    maxlength   = settings['maxlength']
    alignment   = settings['alignment']
    threshold   = settings['threshold']


    lengthindex = range(minlength, maxlength + 1)
    
    A_site = settings['A_site shift']
    P_site = A_site - 3
    E_site = A_site - 6
    two_site = A_site + 6
    one_site = A_site + 3
    mone_site = A_site - 9
    mtwo_site = A_site - 12
    
    
    frameshift      = settings['frameshift']

    plot_upstream   = settings['plot_upstream'] / 3 * 3        #change window to interval of 3
    plot_downstream = settings['plot_downstream'] / 3 * 3
    
    next_codon = settings['next_codon']
    
    
    # define plot length
    plotlength = plot_upstream + plot_downstream + 1
    positionindex = range(0, plotlength)
    
    # load density files
    density_plus  = plus
    density_minus = minus 
    
    # load annotation
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    seq_list    = gff_dict['Sequence'] 
 
    gff_extra  = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    start_trim = settings['start_trim'] / 3 * 3
    stop_trim  = settings['stop_trim'] / 3 * 3   
    
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    plot_upstream_1   = str(plot_upstream)  +'_'
    plot_downstream_1 = str(plot_downstream)+'_'
    start_trim_1      = str(start_trim)     +'_'
    stop_trim_1       = str(stop_trim)      +'_'
    frameshift_1      = str(frameshift)     +'_'
    a_site_1          = str(A_site)         +'_'
    
    name_settings = minlength_1+maxlength_1+plot_upstream_1+plot_downstream_1
    name_settings += start_trim_1+stop_trim_1+frameshift_1
        
    # import genetic code
    aa_code, codon_code = ribo_util.get_genetic_code()
    
    
    '''output data structure'''
    
    # aa/codon avgplot      = { aa: {length: [values]}} 
    # aa/codon count        = { aa: N}
    # aa/codon score        = { aa: [aa_list], A site score: [A_list]...
                                
    aa_avgplot  = {}
    aa_count    = {}
    aa_score    = {}
    aa_score_df = {}
    aa_list     = []
    aa_A_list   = []
    aa_P_list   = []
    aa_E_list   = []
    aa_1_list   = []
    aa_2_list   = []
    aa_m1_list   = []
    aa_m2_list   = []
    
    codon_avgplot  = {}
    codon_count    = {}
    codon_score    = {}
    codon_score_df = {}
    codon_list     = []
    codon_A_list   = []
    codon_P_list   = []
    codon_E_list   = []
    codon_1_list   = []
    codon_2_list   = []
    codon_m1_list   = []
    codon_m2_list   = []
    
    # create empty datastructures to store density info:
    for aa in aa_code.keys():
        aa_avgplot[aa] = {length : [0]*(plotlength) for length in lengthindex}
        aa_count[aa]   = 0
        aa_score[aa]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
        
    for codon in codon_code.keys():
        codon_avgplot[codon] = {length : [0]*(plotlength) for length in lengthindex}
        codon_count[codon]   = 0
        codon_score[codon]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
    
    '''genes in data''' 
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']         = [0, []]
    excluded_genes['low_density']   = [0, []]
    excluded_genes['not_divisible'] = [0, []]
    # count included genes
    included_genes = [0, []]
    
    
    
    #list_location = '/Volumes/HDD/Ribo_seq/libraries/analysis/reference_information/Workbook1.csv'
    
    #infile = pd.read_csv(list_location)
    #TE_list = infile['yeaR'].tolist()
        
    '''iterate through every annotated gene to get codon density info:'''
    
    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,stop_list, strand_list, seq_list):
        
        #if not alias in TE_list:
            #continue
        
        
        ''' define start and stop positions for codons to analyze:
        
        # = codon 
        #################################### = GFF sequence  (50 extra nt)
           ##############################    = AUG to UGA
             ##########################      = density to analyze : remove start and stop peaks
                 ##################          = codons to analyze : remove plot window
        
        '''
        
        # First, define density without start and stop peaks:
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim + frameshift
            density_stop  = stop  - stop_trim + frameshift
            
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim - frameshift
            density_stop  = stop  + stop_trim - 3 + frameshift
            
            period = -1
        
        # GFF seq has 50 extra nucleotides, so remove:
        # Also remove several codons from start and stop positions:
        
        codon_seq_start = gff_extra + start_trim + plot_upstream + frameshift
        codon_seq_stop  = -gff_extra - stop_trim - plot_downstream + frameshift
            
        seq       = sequence[codon_seq_start : codon_seq_stop]
        seqlength = len(seq)
        
        # exclude genes that are not divisable by 3
        if seqlength % 3 != 0:
            excluded_genes['not_divisible'][0] += 1
            excluded_genes['not_divisible'][1].append(alias)
            continue
        
        # exclude genes shorter than plot
        if seqlength < plotlength + 1:
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        elif codon_seq_start - codon_seq_stop > len(sequence):
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
    
        # make a list of codons in the sequence:
        codons_seq = [seq[i:i+3] for i in range(0, seqlength, 3)]
        aa_seq     = [codon_code[codon] for codon in codons_seq]
        
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        # fill density dict with density info
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            length_density_float = [float(i) for i in length_density]
            gene_density[length] = length_density
            total_density        = [x + y for x, y in itertools.izip(total_density, length_density)]
            
        gene_avgreads = float(sum(total_density)) / float(genelength)

            
        # normalize density by total gene density
        if gene_avgreads * 3 < threshold:
            excluded_genes['low_density'][0] += 1
            excluded_genes['low_density'][1].append(alias)
            continue
            
        else: 
            relative_density = {}
            
            for length in lengthindex:
                relative_density[length] = [reads / gene_avgreads for reads in gene_density[length]]
        
        # add data to dataframe
        codon_index = 0
        for codon in codons_seq:
            #if codon == 'TAA':
            #    print "kekekeek"
            codon_position = (codon_index * 3)
            aa = aa_seq[codon_index]
            if aa == '_':
                print 'Stop codon found in'
                print alias
            codon_index        += 1
            
            '''if codon_index > 50:
                continue
            '''
            codon_count[codon] += 1
            aa_count[aa]       += 1
            
            for length in lengthindex:
                for position in positionindex:
                    
                    density_position = codon_position + position
                    density          = relative_density[length][density_position]
                    
                    codon_avgplot[codon][length][position] += density
                    aa_avgplot[aa][length][position]       += density
                    
        included_genes[0] += 1
        included_genes[1].append(alias)
    
    #divide data by total instances of the codon or aa
    for codon in codon_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if codon_count[codon] == 0:
                    continue 
                else: 
                    codon_avgplot[codon][length][position] = codon_avgplot[codon][length][position] / codon_count[codon]
    
    for aa in aa_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if aa_count[aa] == 0:
                    continue
                else: 
                    aa_avgplot[aa][length][position] = aa_avgplot[aa][length][position] / aa_count[aa]
                
   #Convert data for plotting and csv
    
    codon_data     = {}
    codon_data_sum = {}
    aa_data     = {}
    aa_data_sum = {}
    
    A_site_shift = plot_upstream - A_site 
    P_site_shift = plot_upstream - P_site 
    E_site_shift = plot_upstream - E_site
    one_site_shift = plot_upstream - one_site 
    two_site_shift = plot_upstream - two_site 
    mone_site_shift = plot_upstream - mone_site 
    mtwo_site_shift = plot_upstream - mtwo_site 
    
    for codon in codon_avgplot.keys():
        df = pd.DataFrame(codon_avgplot[codon])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        codon_data[codon]     = df
        codon_data_sum[codon] = df.sum(0)
        
        codon_score[codon][0] = codon_data_sum[codon][A_site_shift]
        codon_score[codon][1] = codon_data_sum[codon][P_site_shift]                          
        codon_score[codon][2] = codon_data_sum[codon][E_site_shift]
        codon_score[codon][3] = codon_data_sum[codon][one_site_shift]
        codon_score[codon][4] = codon_data_sum[codon][two_site_shift]                          
        codon_score[codon][5] = codon_data_sum[codon][mone_site_shift]
        codon_score[codon][6] = codon_data_sum[codon][mtwo_site_shift]
        
        codon_list.append(codon)
        codon_A_list.append(codon_score[codon][0])
        codon_P_list.append(codon_score[codon][1])
        codon_E_list.append(codon_score[codon][2])
        codon_1_list.append(codon_score[codon][3])
        codon_2_list.append(codon_score[codon][4])
        codon_m1_list.append(codon_score[codon][5])
        codon_m2_list.append(codon_score[codon][6])

    for aa in aa_avgplot.keys():
        df = pd.DataFrame(aa_avgplot[aa])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        aa_data[aa]     = df
        aa_data_sum[aa] = df.sum(0) 
                                    
        aa_score[aa][0] = aa_data_sum[aa][A_site_shift]   
        aa_score[aa][1] = aa_data_sum[aa][P_site_shift]
        aa_score[aa][2] = aa_data_sum[aa][E_site_shift]
        aa_score[aa][3] = aa_data_sum[aa][one_site_shift]   
        aa_score[aa][4] = aa_data_sum[aa][two_site_shift]
        aa_score[aa][5] = aa_data_sum[aa][mone_site_shift] 
        aa_score[aa][6] = aa_data_sum[aa][mtwo_site_shift]

        aa_list.append(aa)
        aa_A_list.append(aa_score[aa][0])
        aa_P_list.append(aa_score[aa][1])
        aa_E_list.append(aa_score[aa][2])
        aa_1_list.append(aa_score[aa][3])
        aa_2_list.append(aa_score[aa][4])
        aa_m1_list.append(aa_score[aa][5])
        aa_m2_list.append(aa_score[aa][6])
        
    codon_score_df['Codon']  = codon_list
    codon_score_df['A_site'] = codon_A_list
    codon_score_df['P_site'] = codon_P_list
    codon_score_df['E_site'] = codon_E_list
    codon_score_df['1_site'] = codon_1_list
    codon_score_df['2_site'] = codon_2_list
    codon_score_df['-1_site'] = codon_m1_list
    codon_score_df['-2_site'] = codon_m2_list

    aa_score_df['Amino Acid'] = aa_list
    aa_score_df['A_site']     = aa_A_list
    aa_score_df['P_site']     = aa_P_list
    aa_score_df['E_site']     = aa_E_list
    aa_score_df['1_site']     = aa_1_list
    aa_score_df['2_site']     = aa_2_list
    aa_score_df['-1_site']     = aa_m1_list
    aa_score_df['-2_site']     = aa_m2_list

    codon_df = pd.DataFrame(codon_score_df)
    codon_df.to_csv(path_pausescore + 'codon_scores.csv')
    
    aa_df = pd.DataFrame(aa_score_df)
    aa_df.to_csv(path_pausescore + 'aa_scores.csv')
    
    ribo_util.makePickle(codon_data,     path_pausescore +'codon_HM_data'  + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(codon_data_sum, path_pausescore +'codon_plot_data'+ name_settings , protocol=pickle.HIGHEST_PROTOCOL)  
    ribo_util.makePickle(codon_score_df, path_pausescore +'codon_scores'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_data,        path_pausescore +'aa_HM_data'     + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(aa_data_sum,    path_pausescore +'aa_plot_data'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_score_df,    path_pausescore +'aa_scores'      + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    
    #print excluded_genes['short'][0]    
    #print excluded_genes['low_density'][0]
    #print excluded_genes['not_divisible'][0]
    #print included_genes[0]
    return 


def pausescore_first_pos(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    multi     = inputs['multiprocess']
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started pause score analysis at " + str(datetime.now())

    for fname in files:
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/'
        if not os.path.exists(path_pausescore):
            os.makedirs(path_pausescore)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
       
        if not multi == 'yes':
            run_pausescore_first_pos(fname, settings, plus, minus, gff_dict, path_pausescore)
        else:     
            argument = [fname, settings, plus, minus, gff_dict, path_pausescore]
            arguments.append(argument)
    
    if multi == 'yes':
        ribo_util.multiprocess(run_pausescore_first_pos, arguments, threads)
    
    print "Finished pause score analysis at " + str(datetime.now())
    
    return

def run_pausescore_high(fname, settings, plus, minus, gff, path_pausescore):
    
    '''define variables'''
    
    minlength   = settings['minlength']
    maxlength   = settings['maxlength']
    alignment   = settings['alignment']
    threshold   = settings['threshold']

    

    lengthindex = range(minlength, maxlength + 1)
    
    A_site = settings['A_site shift']
    P_site = A_site - 3
    E_site = A_site - 6
    two_site = A_site + 6
    one_site = A_site + 3
    mone_site = A_site - 9
    mtwo_site = A_site - 12
    
    
    frameshift      = settings['frameshift']

    plot_upstream   = settings['plot_upstream'] / 3 * 3        #change window to interval of 3
    plot_downstream = settings['plot_downstream'] / 3 * 3
    
    next_codon = settings['next_codon']
    
    
    # define plot length
    plotlength = plot_upstream + plot_downstream + 1
    positionindex = range(0, plotlength)
    
    # load density files
    density_plus  = plus
    density_minus = minus 
    
    # load annotation
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    seq_list    = gff_dict['Sequence'] 
 
    gff_extra  = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    start_trim = settings['start_trim'] / 3 * 3
    stop_trim  = settings['stop_trim'] / 3 * 3   
    
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    plot_upstream_1   = str(plot_upstream)  +'_'
    plot_downstream_1 = str(plot_downstream)+'_'
    start_trim_1      = str(start_trim)     +'_'
    stop_trim_1       = str(stop_trim)      +'_'
    frameshift_1      = str(frameshift)     +'_'
    a_site_1          = str(A_site)         +'_'
    
    name_settings = minlength_1+maxlength_1+plot_upstream_1+plot_downstream_1
    name_settings += start_trim_1+stop_trim_1+frameshift_1
        
    # import genetic code
    aa_code, codon_code = ribo_util.get_genetic_code()
    
    
    '''output data structure'''
    
    # aa/codon avgplot      = { aa: {length: [values]}} 
    # aa/codon count        = { aa: N}
    # aa/codon score        = { aa: [aa_list], A site score: [A_list]...
                                
    aa_avgplot  = {}
    aa_count    = {}
    aa_score    = {}
    aa_score_df = {}
    aa_list     = []
    aa_A_list   = []
    aa_P_list   = []
    aa_E_list   = []
    aa_1_list   = []
    aa_2_list   = []
    aa_m1_list   = []
    aa_m2_list   = []
    
    codon_avgplot  = {}
    codon_count    = {}
    codon_score    = {}
    codon_score_df = {}
    codon_list     = []
    codon_A_list   = []
    codon_P_list   = []
    codon_E_list   = []
    codon_1_list   = []
    codon_2_list   = []
    codon_m1_list   = []
    codon_m2_list   = []
    
    # create empty datastructures to store density info:
    for aa in aa_code.keys():
        aa_avgplot[aa] = {length : [0]*(plotlength) for length in lengthindex}
        aa_count[aa]   = 0
        aa_score[aa]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
        
    for codon in codon_code.keys():
        codon_avgplot[codon] = {length : [0]*(plotlength) for length in lengthindex}
        codon_count[codon]   = 0
        codon_score[codon]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
    
    '''genes in data''' 
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']         = [0, []]
    excluded_genes['low_density']   = [0, []]
    excluded_genes['not_divisible'] = [0, []]
    # count included genes
    included_genes = [0, []]
    
    
    
    #list_location = '/Volumes/HDD/Ribo_seq/libraries/analysis/reference_information/Workbook1.csv'
    
    #infile = pd.read_csv(list_location)
    #TE_list = infile['yeaR'].tolist()
        
    '''iterate through every annotated gene to get codon density info:'''
    
    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,stop_list, strand_list, seq_list):
        
        #if not alias in TE_list:
            #continue
        
        
        ''' define start and stop positions for codons to analyze:
        
        # = codon 
        #################################### = GFF sequence  (50 extra nt)
           ##############################    = AUG to UGA
             ##########################      = density to analyze : remove start and stop peaks
                 ##################          = codons to analyze : remove plot window
        
        '''
        
        # First, define density without start and stop peaks:
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim + frameshift
            density_stop  = stop  - stop_trim + frameshift
            
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim - frameshift
            density_stop  = stop  + stop_trim - 3 + frameshift
            
            period = -1
        
        # GFF seq has 50 extra nucleotides, so remove:
        # Also remove several codons from start and stop positions:
        
        codon_seq_start = gff_extra + start_trim + plot_upstream + frameshift
        codon_seq_stop  = -gff_extra - stop_trim - plot_downstream + frameshift
            
        seq       = sequence[codon_seq_start : codon_seq_stop]
        seqlength = len(seq)
        
        # exclude genes that are not divisable by 3
        if seqlength % 3 != 0:
            excluded_genes['not_divisible'][0] += 1
            excluded_genes['not_divisible'][1].append(alias)
            continue
        
        # exclude genes shorter than plot
        if seqlength < plotlength + 1:
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        elif codon_seq_start - codon_seq_stop > len(sequence):
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
    
        # make a list of codons in the sequence:
        codons_seq = [seq[i:i+3] for i in range(0, seqlength, 3)]
        aa_seq     = [codon_code[codon] for codon in codons_seq]
        
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        # fill density dict with density info
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            length_density_float = [float(i) for i in length_density]
            gene_density[length] = length_density
            total_density        = [x + y for x, y in itertools.izip(total_density, length_density)]
            
        gene_avgreads = float(sum(total_density)) / float(genelength)

            
        # normalize density by total gene density
        if gene_avgreads * 3 < threshold:
            excluded_genes['low_density'][0] += 1
            excluded_genes['low_density'][1].append(alias)
            continue
            
        else: 
            relative_density = {}
            
            for length in lengthindex:
                relative_density[length] = [reads / gene_avgreads for reads in gene_density[length]]
        
        # add data to dataframe
        codon_index = 0
        for codon in codons_seq:
            #if codon == 'TAA':
            #    print "kekekeek"
            codon_position = (codon_index * 3)
            aa = aa_seq[codon_index]
            if aa == '_':
                print 'Stop codon found in'
                print alias
            codon_index        += 1
            
            '''if codon_index > 50:
                continue
            '''
            codon_count[codon] += 1
            aa_count[aa]       += 1
            
            for length in lengthindex:
                for position in positionindex:
                    
                    density_position = codon_position + position
                    density          = relative_density[length][density_position]
                    
                    codon_avgplot[codon][length][position] += density
                    aa_avgplot[aa][length][position]       += density
                    
        included_genes[0] += 1
        included_genes[1].append(alias)
    
    #divide data by total instances of the codon or aa
    for codon in codon_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if codon_count[codon] == 0:
                    continue 
                else: 
                    codon_avgplot[codon][length][position] = codon_avgplot[codon][length][position] / codon_count[codon]
    
    for aa in aa_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if aa_count[aa] == 0:
                    continue
                else: 
                    aa_avgplot[aa][length][position] = aa_avgplot[aa][length][position] / aa_count[aa]
                
   #Convert data for plotting and csv
    
    codon_data     = {}
    codon_data_sum = {}
    aa_data     = {}
    aa_data_sum = {}
    
    A_site_shift = plot_upstream - A_site 
    P_site_shift = plot_upstream - P_site 
    E_site_shift = plot_upstream - E_site
    one_site_shift = plot_upstream - one_site 
    two_site_shift = plot_upstream - two_site 
    mone_site_shift = plot_upstream - mone_site 
    mtwo_site_shift = plot_upstream - mtwo_site 
    
    for codon in codon_avgplot.keys():
        df = pd.DataFrame(codon_avgplot[codon])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        codon_data[codon]     = df
        codon_data_sum[codon] = df.sum(0)
        
        codon_score[codon][0] = sum(codon_data_sum[codon][A_site_shift: A_site_shift + 3:1])/3
        codon_score[codon][1] = sum(codon_data_sum[codon][P_site_shift: P_site_shift + 3:1])/3                            
        codon_score[codon][2] = sum(codon_data_sum[codon][E_site_shift: E_site_shift + 3:1])/3
        codon_score[codon][3] = sum(codon_data_sum[codon][one_site_shift: one_site_shift + 3:1])/3
        codon_score[codon][4] = sum(codon_data_sum[codon][two_site_shift: two_site_shift + 3:1])/3                            
        codon_score[codon][5] = sum(codon_data_sum[codon][mone_site_shift: mone_site_shift + 3:1])/3
        codon_score[codon][6] = sum(codon_data_sum[codon][mtwo_site_shift: mtwo_site_shift + 3:1])/3
        
        codon_list.append(codon)
        codon_A_list.append(codon_score[codon][0])
        codon_P_list.append(codon_score[codon][1])
        codon_E_list.append(codon_score[codon][2])
        codon_1_list.append(codon_score[codon][3])
        codon_2_list.append(codon_score[codon][4])
        codon_m1_list.append(codon_score[codon][5])
        codon_m2_list.append(codon_score[codon][6])

    for aa in aa_avgplot.keys():
        df = pd.DataFrame(aa_avgplot[aa])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        aa_data[aa]     = df
        aa_data_sum[aa] = df.sum(0) 
                                    
        aa_score[aa][0] = sum(aa_data_sum[aa][A_site_shift: A_site_shift + 3:1])/3    
        aa_score[aa][1] = sum(aa_data_sum[aa][P_site_shift: P_site_shift + 3:1])/3
        aa_score[aa][2] = sum(aa_data_sum[aa][E_site_shift: E_site_shift + 3:1])/3
        aa_score[aa][3] = sum(aa_data_sum[aa][one_site_shift: one_site_shift + 3:1])/3    
        aa_score[aa][4] = sum(aa_data_sum[aa][two_site_shift: two_site_shift + 3:1])/3
        aa_score[aa][5] = sum(aa_data_sum[aa][mone_site_shift: mone_site_shift + 3:1])/3 
        aa_score[aa][6] = sum(aa_data_sum[aa][mtwo_site_shift: mtwo_site_shift + 3:1])/3 

        aa_list.append(aa)
        aa_A_list.append(aa_score[aa][0])
        aa_P_list.append(aa_score[aa][1])
        aa_E_list.append(aa_score[aa][2])
        aa_1_list.append(aa_score[aa][3])
        aa_2_list.append(aa_score[aa][4])
        aa_m1_list.append(aa_score[aa][5])
        aa_m2_list.append(aa_score[aa][6])
        
    codon_score_df['Codon']  = codon_list
    codon_score_df['A_site'] = codon_A_list
    codon_score_df['P_site'] = codon_P_list
    codon_score_df['E_site'] = codon_E_list
    codon_score_df['1_site'] = codon_1_list
    codon_score_df['2_site'] = codon_2_list
    codon_score_df['-1_site'] = codon_m1_list
    codon_score_df['-2_site'] = codon_m2_list

    aa_score_df['Amino Acid'] = aa_list
    aa_score_df['A_site']     = aa_A_list
    aa_score_df['P_site']     = aa_P_list
    aa_score_df['E_site']     = aa_E_list
    aa_score_df['1_site']     = aa_1_list
    aa_score_df['2_site']     = aa_2_list
    aa_score_df['-1_site']     = aa_m1_list
    aa_score_df['-2_site']     = aa_m2_list

    codon_df = pd.DataFrame(codon_score_df)
    codon_df.to_csv(path_pausescore + 'codon_scores.csv')
    
    aa_df = pd.DataFrame(aa_score_df)
    aa_df.to_csv(path_pausescore + 'aa_scores.csv')
    
    ribo_util.makePickle(codon_data,     path_pausescore +'codon_HM_data'  + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(codon_data_sum, path_pausescore +'codon_plot_data'+ name_settings , protocol=pickle.HIGHEST_PROTOCOL)  
    ribo_util.makePickle(codon_score_df, path_pausescore +'codon_scores'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_data,        path_pausescore +'aa_HM_data'     + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(aa_data_sum,    path_pausescore +'aa_plot_data'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_score_df,    path_pausescore +'aa_scores'      + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    
    #print excluded_genes['short'][0]    
    #print excluded_genes['low_density'][0]
    #print excluded_genes['not_divisible'][0]
    #print included_genes[0]
    return 


def pausescore_high(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    multi     = inputs['multiprocess']
    arguments = []
    penalty     = settings['penalty']

    if not files:
        print("There are no files")
        return
    
    print "Started pause score analysis at " + str(datetime.now())
    plus_dict, minus_dict = high_dict(plus_dict, minus_dict, penalty)
    for fname in files:
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/'
        if not os.path.exists(path_pausescore):
            os.makedirs(path_pausescore)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
       
        if not multi == 'yes':
            run_pausescore_high(fname, settings, plus, minus, gff_dict, path_pausescore)
        else:     
            argument = [fname, settings, plus, minus, gff_dict, path_pausescore]
            arguments.append(argument)
    
    if multi == 'yes':
        ribo_util.multiprocess(run_pausescore_high, arguments, threads)
    
    print "Finished pause score analysis at " + str(datetime.now())
    
    return

def high_dict(plus_dict, minus_dict, penalty):
    high_plus = {}
    high_minus = {}
    
    for sample in plus_dict:
        
        sample_plus = pd.DataFrame.from_dict(plus_dict[sample]) - penalty
        sample_minus = pd.DataFrame.from_dict(minus_dict[sample]) - penalty

        sample_plus[sample_plus < 0] = 0
        sample_minus[sample_minus < 0] = 0

        high_plus[sample] = sample_plus.to_dict('list')
        high_minus[sample] = sample_minus.to_dict('list')

    return high_plus, high_minus


def pausescore_waves_1st_g(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    multi     = inputs['multiprocess']
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started pause score analysis at " + str(datetime.now())

    for fname in files:
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/waves/'
        if not os.path.exists(path_pausescore):
            os.makedirs(path_pausescore)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
       
        if not multi == 'yes':
            run_pausescore_waves_1st_g(fname, settings, plus, minus, gff_dict, path_pausescore)
        else:     
            argument = [fname, settings, plus, minus, gff_dict, path_pausescore]
            arguments.append(argument)
    
    if multi == 'yes':
        ribo_util.multiprocess(run_pausescore_waves_1st_g, arguments, threads)
    
    print "Finished pause score analysis at " + str(datetime.now())
    
    return



def run_pausescore_waves_1st_g(fname, settings, plus, minus, gff, path_pausescore):
    
    '''define variables'''
    
    minlength   = settings['minlength']
    maxlength   = settings['maxlength']
    lengthindex = range(minlength, maxlength + 1)
    
    A_site = settings['A_site shift']
    P_site = A_site - 3
    E_site = A_site - 6
    
    frameshift      = settings['frameshift']

    plot_upstream   = settings['plot_upstream_wave'] / 3 * 3        #change window to interval of 3
    plot_downstream = settings['plot_downstream_wave'] / 3 * 3
    
    next_codon = settings['next_codon']
    
    # define plot length
    plotlength = plot_upstream + plot_downstream + 1
    positionindex = range(0, plotlength)
    
    # load density files
    density_plus  = plus
    density_minus = minus 
    
    # load annotation
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    seq_list    = gff_dict['Sequence'] 
 
    gff_extra  = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    start_trim = settings['start_trim'] / 3 * 3
    stop_trim  = settings['stop_trim'] / 3 * 3   
    
    #proximity_threshold = settings['proximity_threshold']

    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    plot_upstream_1   = str(plot_upstream)  +'_'
    plot_downstream_1 = str(plot_downstream)+'_'
    start_trim_1      = str(start_trim)     +'_'
    stop_trim_1       = str(stop_trim)      +'_'
    frameshift_1      = str(frameshift)     +'_'
    next_codon_1      =  str(next_codon)    +'_'
    
    name_settings = plot_upstream_1+plot_downstream_1+start_trim_1+stop_trim_1
    name_settings += minlength_1+maxlength_1+frameshift_1+next_codon_1
    
    
    
    # import genetic code
    aa_code, codon_code = ribo_util.get_genetic_code()
    
    
    '''output data structure'''
    
    # aa/codon avgplot      = { aa: {length: [values]}} 
    # aa/codon count        = { aa: N}
    # aa/codon score        = { aa: [aa_list], A site score: [A_list]...
                                
    aa_avgplot  = {}
    aa_count    = {}
    aa_score    = {}
    aa_score_df = {}
    aa_list     = []
    aa_A_list   = []
    aa_P_list   = []
    aa_E_list   = []
    
    codon_avgplot  = {}
    codon_count    = {}
    codon_score    = {}
    codon_score_df = {}
    codon_list     = []
    codon_A_list   = []
    codon_P_list   = []
    codon_E_list   = []
    
    # create empty datastructures to store density info:
    for i in range(1, 214):
        aa_avgplot[i] = {length : [0]*(plotlength) for length in lengthindex}
        aa_count[i]   = 0
        aa_score[i]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
        
    for codon in aa_code['G']:
        for i in range(1, 214):
            codon_avgplot[codon + '_' + str(i)] = {length : [0]*(plotlength) for length in lengthindex}
            codon_count[codon + '_' + str(i)]   = 0
            codon_score[codon + '_' + str(i)]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
    '''
    for aa in aa_code.keys():
        aa_avgplot[aa] = {length : [0]*(plotlength) for length in lengthindex}
        aa_count[aa]   = 0
        aa_score[aa]   = [0, 0, 0] # [Asite, P site, E site] 
        
    for codon in codon_code.keys():
        codon_avgplot[codon] = {length : [0]*(plotlength) for length in lengthindex}
        codon_count[codon]   = 0
        codon_score[codon]   = [0, 0, 0] # [Asite, P site, E site] '''
    
    '''genes in data''' 
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']         = [0, []]
    excluded_genes['low_density']   = [0, []]
    excluded_genes['not_divisible'] = [0, []]
    # count included genes
    included_genes = [0, []]
        
        
    '''iterate through every annotated gene to get codon density info:'''
    
    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,stop_list, strand_list, seq_list):  
                
        ''' define start and stop positions for codons to analyze:
        
        # = codon 
        #################################### = GFF sequence  (50 extra nt)
           ##############################    = AUG to UGA
             ##########################      = density to analyze : remove start and stop peaks
                 ##################          = codons to analyze : remove plot window
        
        '''
        
        # First, define density without start and stop peaks:
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim + frameshift
            density_stop  = stop  - stop_trim + frameshift
            
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim - frameshift
            density_stop  = stop  + stop_trim - 3 + frameshift
            
            period = -1
        
        # GFF seq has 50 extra nucleotides, so remove:
        # Also remove several codons from start and stop positions:
        
        codon_seq_start = gff_extra + start_trim + plot_upstream + frameshift
        codon_seq_stop  = -gff_extra - stop_trim - plot_downstream + frameshift
            
        seq       = sequence[codon_seq_start : codon_seq_stop]
        seqlength = len(seq)
    
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        # fill density dict with density info
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            length_density_float = [float(i) for i in length_density]
            gene_density[length] = length_density
            total_density        = [x + y for x, y in itertools.izip(total_density, length_density)]
            
        gene_avgreads = float(sum(total_density)) / float(genelength)
        # exclude genes that are not divisable by 3
        if seqlength % 3 != 0:
            excluded_genes['not_divisible'][0] += 1
            excluded_genes['not_divisible'][1].append(alias)
            continue
        
        # exclude genes shorter than plot
        if seqlength < plotlength + 1:
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        elif codon_seq_start - codon_seq_stop > len(sequence):
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        # normalize density by total gene density
        if gene_avgreads < 0.1:
            excluded_genes['low_density'][0] += 1
            excluded_genes['low_density'][1].append(alias)
            continue
            
        else: 
            relative_density = {}
            
            for length in lengthindex:
                relative_density[length] = [reads / gene_avgreads for reads in gene_density[length]]
        
        # add data to dataframe
        codon_index = 0
        
        # make a list of codons in the sequence:
        '''
        codons_seq = [seq[i:i+3] for i in range(0, seqlength, 3)]
        aa_seq     = [codon_code[codon] for codon in codons_seq]'''


        G_codons_count = {x : 0 for x in aa_code['G']}
        codons_seq = []
        codon_indexes = []

        for i in range(0, seqlength, 3):
            if seq[i:i+3] in G_codons_count:
                G_codons_count[seq[i:i+3]] += 1
                codons_seq.append(seq[i:i+3] + '_' + str(G_codons_count[seq[i:i+3]]))
                codon_indexes.append(i)
        
        aa_seq     = [i for i in range(1, len(codons_seq) + 1)]
        
        skip_mask = []

        for i in range(len(codon_indexes)):
            if i == 0:
                interval = codon_indexes[i]
            else:
                interval = codon_indexes[i] - codon_indexes[i-1]
            
            if interval < plot_upstream:
                skip_mask.append(i)
        
        codon_indexes = [codon_indexes[i] for i in range(len(codon_indexes)) if i not in skip_mask]
        aa_seq = [aa_seq[i] for i in range(len(codon_indexes)) if i not in skip_mask]
        codons_seq = [codons_seq[i] for i in range(len(codon_indexes)) if i not in skip_mask]
        
        
        '''
        if len(codon_indexes) > 0:
            intervals = [codon_indexes[i+1] - codon_indexes[i] for i in range(len(codon_indexes) - 1)]
        #    print codon_indexes
            intervals.append(codon_indexes[0])
            intervals = [x for x in intervals if x < plot_upstream]
            if len(intervals) > 0:
                continue
        '''
        #for codon in codons_seq:
            #codon_position = (codon_index * 3)
            #aa = aa_seq[codon_index]

        for codon, aa, codon_position in itertools.izip(codons_seq, aa_seq, codon_indexes):   
            if next_codon == 'yes':
                if aa in aa_seq[codon_index+1:plot_downstream + 1]:
                    codon_index += 1
                    continue
            
            codon_index        += 1
            codon_count[codon] += 1
            aa_count[aa]       += 1
            
            for length in lengthindex:
                for position in positionindex:
                    
                    density_position = codon_position + position
                    density          = relative_density[length][density_position]
                    
                    codon_avgplot[codon][length][position] += density
                    aa_avgplot[aa][length][position]       += density
                    
        included_genes[0] += 1
        included_genes[1].append(alias)
    
    #divide data by total instances of the codon or aa
    for codon in codon_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if codon_count[codon] == 0:
                    continue 
                else: 
                    codon_avgplot[codon][length][position] = codon_avgplot[codon][length][position] / codon_count[codon]
    
    for aa in aa_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if aa_count[aa] == 0:
                    continue
                else: 
                    aa_avgplot[aa][length][position] = aa_avgplot[aa][length][position] / aa_count[aa]
                
    '''Convert data for plotting and csv'''
    
    codon_data     = {}
    codon_data_sum = {}
    aa_data     = {}
    aa_data_sum = {}
    
    A_site_shift = plot_upstream - A_site 
    P_site_shift = plot_upstream - P_site 
    E_site_shift = plot_upstream - E_site 
    
    for codon in codon_avgplot.keys():
        df = pd.DataFrame(codon_avgplot[codon])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        codon_data[codon]     = df
        codon_data_sum[codon] = df.sum(0)
        
        codon_score[codon][0] = sum(codon_data_sum[codon][A_site_shift: A_site_shift + 3:1])/3
        codon_score[codon][1] = sum(codon_data_sum[codon][P_site_shift: P_site_shift + 3:1])/3                            
        codon_score[codon][2] = sum(codon_data_sum[codon][E_site_shift: E_site_shift + 3:1])/3
        
        codon_list.append(codon)
        codon_A_list.append(codon_score[codon][0])
        codon_P_list.append(codon_score[codon][1])
        codon_E_list.append(codon_score[codon][2])
        
    for aa in aa_avgplot.keys():
        df = pd.DataFrame(aa_avgplot[aa])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        aa_data[aa]     = df
        aa_data_sum[aa] = df.sum(0) 
                                    
        aa_score[aa][0] = sum(aa_data_sum[aa][A_site_shift: A_site_shift + 3:1])/3    
        aa_score[aa][1] = sum(aa_data_sum[aa][P_site_shift: P_site_shift + 3:1])/3
        aa_score[aa][2] = sum(aa_data_sum[aa][E_site_shift: E_site_shift + 3:1])/3  
        
        aa_list.append(aa)
        aa_A_list.append(aa_score[aa][0])
        aa_P_list.append(aa_score[aa][1])
        aa_E_list.append(aa_score[aa][2])
    
    codon_score_df['Codon']  = codon_list
    codon_score_df['A_site'] = codon_A_list
    codon_score_df['P_site'] = codon_P_list
    codon_score_df['E_site'] = codon_E_list
    
    aa_score_df['Amino Acid'] = aa_list
    aa_score_df['A_site']     = aa_A_list
    aa_score_df['P_site']     = aa_P_list
    aa_score_df['E_site']     = aa_E_list
    
    
    codon_df = pd.DataFrame(codon_score_df)
    codon_df.to_csv(path_pausescore + 'codon_scores.csv')
    
    aa_df = pd.DataFrame(aa_score_df)
    aa_df.to_csv(path_pausescore + 'aa_scores.csv')
    
    ribo_util.makePickle(codon_data,     path_pausescore + name_settings +'codon_HM_data'  , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(codon_data_sum, path_pausescore + name_settings +'codon_plot_data', protocol=pickle.HIGHEST_PROTOCOL)  
    ribo_util.makePickle(codon_score_df, path_pausescore + name_settings +'codon_scores'   , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_data,        path_pausescore + name_settings +'aa_HM_data'     , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(aa_data_sum,    path_pausescore + name_settings +'aa_plot_data'   , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_score_df,    path_pausescore + name_settings +'aa_scores'      , protocol=pickle.HIGHEST_PROTOCOL)
    
    print excluded_genes['short'][0]    
    print excluded_genes['low_density'][0]
    print excluded_genes['not_divisible'][0]
    print included_genes[0]
    return 




def pausescore_g_swarm(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    multi     = inputs['multiprocess']
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started pause score analysis at " + str(datetime.now())

    for fname in files:
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/'
        if not os.path.exists(path_pausescore):
            os.makedirs(path_pausescore)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
       
        if not multi == 'yes':
            run_pausescore_g_swarm(fname, settings, plus, minus, gff_dict, path_pausescore)
        else:     
            argument = [fname, settings, plus, minus, gff_dict, path_pausescore]
            arguments.append(argument)
    
    if multi == 'yes':
        ribo_util.multiprocess(run_pausescore_g_swarm, arguments, threads)
    
    print "Finished pause score analysis at " + str(datetime.now())
    
    return


def run_pausescore_g_swarm(fname, settings, plus, minus, gff, path_pausescore):
    
    '''define variables'''
    
    minlength   = settings['minlength']
    maxlength   = settings['maxlength']
    alignment   = settings['alignment']
    threshold   = settings['threshold']


    lengthindex = range(minlength, maxlength + 1)
    
    A_site = settings['A_site shift']
    P_site = A_site - 3
    E_site = A_site - 6
    two_site = A_site + 6
    one_site = A_site + 3
    mone_site = A_site - 9
    mtwo_site = A_site - 12
    
    
    frameshift      = settings['frameshift']

    plot_upstream   = settings['plot_upstream'] / 3 * 3        #change window to interval of 3
    plot_downstream = settings['plot_downstream'] / 3 * 3
    
    next_codon = settings['next_codon']
    
    
    # define plot length
    plotlength = plot_upstream + plot_downstream + 1
    positionindex = range(0, plotlength)
    
    # load density files
    density_plus  = plus
    density_minus = minus 
    
    # load annotation
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    seq_list    = gff_dict['Sequence'] 
 
    gff_extra  = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    start_trim = settings['start_trim'] / 3 * 3
    stop_trim  = settings['stop_trim'] / 3 * 3   
    
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    plot_upstream_1   = str(plot_upstream)  +'_'
    plot_downstream_1 = str(plot_downstream)+'_'
    start_trim_1      = str(start_trim)     +'_'
    stop_trim_1       = str(stop_trim)      +'_'
    frameshift_1      = str(frameshift)     +'_'
    a_site_1          = str(A_site)         +'_'
    
    name_settings = minlength_1+maxlength_1+plot_upstream_1+plot_downstream_1
    name_settings += start_trim_1+stop_trim_1+frameshift_1
        
    # import genetic code
    aa_code, codon_code = ribo_util.get_genetic_code()
    
    
    '''output data structure'''
    
    # aa/codon avgplot      = { aa: {length: [values]}} 
    # aa/codon count        = { aa: N}
    # aa/codon score        = { aa: [aa_list], A site score: [A_list]...
                                
    aa_avgplot  = {}
    aa_count    = {}
    aa_score    = {}
    aa_score_df = {}
    aa_list     = []
    aa_A_list   = []
    aa_P_list   = []
    aa_E_list   = []
    aa_1_list   = []
    aa_2_list   = []
    aa_m1_list   = []
    aa_m2_list   = []
    
    codon_avgplot  = {}
    codon_count    = {}
    codon_score    = {}
    codon_score_df = {}
    codon_list     = []
    codon_A_list   = []
    codon_P_list   = []
    codon_E_list   = []
    codon_1_list   = []
    codon_2_list   = []
    codon_m1_list   = []
    codon_m2_list   = []
    
    # create empty datastructures to store density info:
    for i in range(1, 71):
        aa_avgplot[i] = {}
        #aa_count[i]   = 0
        #aa_score[i]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
        
    for codon in aa_code['G']:
        for i in range(1, 71):
            codon_avgplot[codon + '_' + str(i)] = {}
            #codon_count[codon + '_' + str(i)]   = 0
            #codon_score[codon + '_' + str(i)]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
    
    '''genes in data''' 
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']         = [0, []]
    excluded_genes['low_density']   = [0, []]
    excluded_genes['not_divisible'] = [0, []]
    # count included genes
    included_genes = [0, []]
    
    
    
    #list_location = '/Volumes/HDD/Ribo_seq/libraries/analysis/reference_information/Workbook1.csv'
    
    #infile = pd.read_csv(list_location)
    #TE_list = infile['yeaR'].tolist()
        
    '''iterate through every annotated gene to get codon density info:'''
    
    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,stop_list, strand_list, seq_list):
        
        #if not alias in TE_list:
            #continue
        
        
        ''' define start and stop positions for codons to analyze:
        
        # = codon 
        #################################### = GFF sequence  (50 extra nt)
           ##############################    = AUG to UGA
             ##########################      = density to analyze : remove start and stop peaks
                 ##################          = codons to analyze : remove plot window
        
        '''
        
        # First, define density without start and stop peaks:
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim + frameshift
            density_stop  = stop  - stop_trim + frameshift
            
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim - frameshift
            density_stop  = stop  + stop_trim - 3 + frameshift
            
            period = -1
        
        # GFF seq has 50 extra nucleotides, so remove:
        # Also remove several codons from start and stop positions:
        
        codon_seq_start = gff_extra + start_trim + plot_upstream + frameshift
        codon_seq_stop  = -gff_extra - stop_trim - plot_downstream + frameshift
            
        seq       = sequence[codon_seq_start : codon_seq_stop]
        seqlength = len(seq)
        
        # exclude genes that are not divisable by 3
        if seqlength % 3 != 0:
            excluded_genes['not_divisible'][0] += 1
            excluded_genes['not_divisible'][1].append(alias)
            continue
        
        # exclude genes shorter than plot
        if seqlength < plotlength + 1:
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        elif codon_seq_start - codon_seq_stop > len(sequence):
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
    
        # make a list of codons in the sequence:
        
        G_codons_count = {x : 0 for x in aa_code['G']}
        codons_seq = []
        codon_indexes = []

        for i in range(0, seqlength, 3):
            if seq[i:i+3] in G_codons_count:
                G_codons_count[seq[i:i+3]] += 1
                codons_seq.append(seq[i:i+3] + '_' + str(G_codons_count[seq[i:i+3]]))
                codon_indexes.append(i)
        
        aa_seq     = [i for i in range(1, len(codons_seq) + 1)]
        
        
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        # fill density dict with density info
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            length_density_float = [float(i) for i in length_density]
            gene_density[length] = length_density
            total_density        = [x + y for x, y in itertools.izip(total_density, length_density)]
            
        gene_avgreads = float(sum(total_density)) / float(genelength)

            
        # normalize density by total gene density
        if gene_avgreads * 3 < threshold:
            excluded_genes['low_density'][0] += 1
            excluded_genes['low_density'][1].append(alias)
            continue
            
        else: 
            relative_density = {}
            
            for length in lengthindex:
                relative_density[length] = [reads / gene_avgreads for reads in gene_density[length]]
        
        # add data to dataframe
        #codon_index = 0
        for codon, aa, codon_position in itertools.izip(codons_seq, aa_seq, codon_indexes):
            
            #codon_position = (codon_index * 3)
            #aa = aa_seq[codon_index]
            
            #codon_index        += 1
            '''
            if codon_index > 50:
                continue
            '''
            #codon_count[codon] += 1
            #aa_count[aa]       += 1

            if aa > 70:
                break
            
            codon_avgplot[codon][alias + '-' + str(codon_position)] = {length : [0]*(plotlength) for length in lengthindex}
            aa_avgplot[aa][alias + '-' + str(codon_position)]       = {length : [0]*(plotlength) for length in lengthindex}

            for length in lengthindex:
                for position in positionindex:
                    
                    density_position = codon_position + position
                    density          = relative_density[length][density_position]
                    
                    codon_avgplot[codon][alias + '-' + str(codon_position)][length][position] += density
                    aa_avgplot[aa][alias + '-' + str(codon_position)][length][position]       += density
                    
        included_genes[0] += 1
        included_genes[1].append(alias)
   
    #divide data by total instances of the codon or aa
    '''
    for codon in codon_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if codon_count[codon] == 0:
                    continue 
                else: 
                    codon_avgplot[codon][length][position] = codon_avgplot[codon][length][position] / codon_count[codon]
    
    for aa in aa_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if aa_count[aa] == 0:
                    continue
                else: 
                    aa_avgplot[aa][length][][position] = aa_avgplot[aa][length][position] / aa_count[aa]
    '''            
   #Convert data for plotting and csv
    
    codon_data     = {}
    codon_data_sum = {}
    aa_data     = {}
    aa_data_sum = {}
    
    A_site_shift = plot_upstream - A_site 
    P_site_shift = plot_upstream - P_site 
    E_site_shift = plot_upstream - E_site
    one_site_shift = plot_upstream - one_site 
    two_site_shift = plot_upstream - two_site 
    mone_site_shift = plot_upstream - mone_site 
    mtwo_site_shift = plot_upstream - mtwo_site 
    
    #print aa_avgplot[7]

    path_swarm = path_pausescore + '/swarm/'

    try:
        os.mkdir(path_swarm)
    except:
        pass
    #print codon_avgplot.keys()
    for gly_number in codon_avgplot.keys():
        
        codon_list  = []
        codon_A_list  = []
        codon_P_list  = []
        codon_E_list  = []
        codon_1_list  = []
        codon_2_list  = []
        codon_m1_list  = []
        codon_m2_list  = []
        codon_score = {}
        codon_score_df = {}
        #print codon_avgplot[gly_number].keys()
        for codon in codon_avgplot[gly_number].keys():
            # codon : SecM-56
            df = pd.DataFrame(codon_avgplot[gly_number][codon])
            df = df.T
            df = df.reindex(index=df.index[::-1])
            #print df
            #break
        
            plot_range  = len(df.columns)
            shift_start = - plot_upstream
            shift_stop  = plot_range - plot_upstream

            df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

            codon_data[codon]     = df
            codon_data_sum[codon] = df.sum(0)
            codon_score[codon]   = [0, 0, 0, 0, 0, 0, 0]
            codon_score[codon][0] = sum(codon_data_sum[codon][A_site_shift: A_site_shift + 3:1])/3
            codon_score[codon][1] = sum(codon_data_sum[codon][P_site_shift: P_site_shift + 3:1])/3                            
            codon_score[codon][2] = sum(codon_data_sum[codon][E_site_shift: E_site_shift + 3:1])/3
            codon_score[codon][3] = sum(codon_data_sum[codon][one_site_shift: one_site_shift + 3:1])/3
            codon_score[codon][4] = sum(codon_data_sum[codon][two_site_shift: two_site_shift + 3:1])/3                            
            codon_score[codon][5] = sum(codon_data_sum[codon][mone_site_shift: mone_site_shift + 3:1])/3
            codon_score[codon][6] = sum(codon_data_sum[codon][mtwo_site_shift: mtwo_site_shift + 3:1])/3
            
            codon_list.append(codon)
            codon_A_list.append(codon_score[codon][0])
            codon_P_list.append(codon_score[codon][1])
            codon_E_list.append(codon_score[codon][2])
            codon_1_list.append(codon_score[codon][3])
            codon_2_list.append(codon_score[codon][4])
            codon_m1_list.append(codon_score[codon][5])
            codon_m2_list.append(codon_score[codon][6])
        
        codon_score_df['Codon']  = codon_list
        codon_score_df['A_site'] = codon_A_list
        codon_score_df['P_site'] = codon_P_list
        codon_score_df['E_site'] = codon_E_list
        codon_score_df['1_site'] = codon_1_list
        codon_score_df['2_site'] = codon_2_list
        codon_score_df['-1_site'] = codon_m1_list
        codon_score_df['-2_site'] = codon_m2_list
        #print codon_score_df
        codon_df = pd.DataFrame(codon_score_df)
        codon_df.to_csv(path_swarm + gly_number + '_codon_scores.csv')
        
        #ribo_util.makePickle(codon_data,     path_swarm + gly_number + '_codon_HM_data'  + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
        #ribo_util.makePickle(codon_data_sum, path_swarm + gly_number + '_codon_plot_data'+ name_settings , protocol=pickle.HIGHEST_PROTOCOL)  
        #ribo_util.makePickle(codon_score_df, path_swarm + gly_number + '_codon_scores'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)

    for gly_number in aa_avgplot.keys():

        aa_list  = []
        aa_score = {}
        aa_score_df = {}
        aa_A_list   = []
        aa_P_list   = []
        aa_E_list   = []
        aa_1_list   = []
        aa_2_list   = []
        aa_m1_list   = []
        aa_m2_list   = []
        for aa in aa_avgplot[gly_number].keys():
            df = pd.DataFrame(aa_avgplot[gly_number][aa])
            df = df.T
            df = df.reindex(index=df.index[::-1])

            plot_range  = len(df.columns)
            shift_start = - plot_upstream
            shift_stop  = plot_range - plot_upstream

            df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

            aa_data[aa]     = df
            aa_data_sum[aa] = df.sum(0) 
            aa_score[aa]   = [0, 0, 0, 0, 0, 0, 0]                            
            aa_score[aa][0] = sum(aa_data_sum[aa][A_site_shift: A_site_shift + 3:1])/3    
            aa_score[aa][1] = sum(aa_data_sum[aa][P_site_shift: P_site_shift + 3:1])/3
            aa_score[aa][2] = sum(aa_data_sum[aa][E_site_shift: E_site_shift + 3:1])/3
            aa_score[aa][3] = sum(aa_data_sum[aa][one_site_shift: one_site_shift + 3:1])/3    
            aa_score[aa][4] = sum(aa_data_sum[aa][two_site_shift: two_site_shift + 3:1])/3
            aa_score[aa][5] = sum(aa_data_sum[aa][mone_site_shift: mone_site_shift + 3:1])/3 
            aa_score[aa][6] = sum(aa_data_sum[aa][mtwo_site_shift: mtwo_site_shift + 3:1])/3 

            aa_list.append(aa)
            aa_A_list.append(aa_score[aa][0])
            aa_P_list.append(aa_score[aa][1])
            aa_E_list.append(aa_score[aa][2])
            aa_1_list.append(aa_score[aa][3])
            aa_2_list.append(aa_score[aa][4])
            aa_m1_list.append(aa_score[aa][5])
            aa_m2_list.append(aa_score[aa][6])
        
    

        aa_score_df['Amino Acid'] = aa_list
        aa_score_df['A_site']     = aa_A_list
        aa_score_df['P_site']     = aa_P_list
        aa_score_df['E_site']     = aa_E_list
        aa_score_df['1_site']     = aa_1_list
        aa_score_df['2_site']     = aa_2_list
        aa_score_df['-1_site']     = aa_m1_list
        aa_score_df['-2_site']     = aa_m2_list
        
        aa_df = pd.DataFrame(aa_score_df)
        aa_df.to_csv(path_swarm + str(gly_number) + '_aa_scores.csv')
        
        #ribo_util.makePickle(aa_data,        path_swarm + gly_number + '_aa_HM_data'     + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
        #ribo_util.makePickle(aa_data_sum,    path_swarm + gly_number + '_aa_plot_data'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
        #ribo_util.makePickle(aa_score_df,    path_swarm + gly_number + '_aa_scores'      + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    
    #print excluded_genes['short'][0]    
    #print excluded_genes['low_density'][0]
    #print excluded_genes['not_divisible'][0]
    #print included_genes[0]
    print path_swarm
    return


def run_logscore(fname, settings, plus, minus, gff, path_pausescore):
    
    '''define variables'''
    
    minlength   = settings['minlength']
    maxlength   = settings['maxlength']
    alignment   = settings['alignment']
    threshold   = settings['threshold']


    #lengthindex = range(minlength, maxlength + 1)
    lengthindex = ['gi|49175990|ref|NC_000913.2|']
    
    A_site = settings['A_site shift']
    P_site = A_site - 3
    E_site = A_site - 6
    two_site = A_site + 6
    one_site = A_site + 3
    mone_site = A_site - 9
    mtwo_site = A_site - 12
    
    
    frameshift      = settings['frameshift']

    plot_upstream   = settings['plot_upstream'] / 3 * 3        #change window to interval of 3
    plot_downstream = settings['plot_downstream'] / 3 * 3
    
    next_codon = settings['next_codon']
    
    
    # define plot length
    plotlength = plot_upstream + plot_downstream + 1
    positionindex = range(0, plotlength)
    
    # load density files
    density_plus  = plus
    density_minus = minus 
    
    # load annotation
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    seq_list    = gff_dict['Sequence'] 
 
    gff_extra  = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    start_trim = settings['start_trim'] / 3 * 3
    stop_trim  = settings['stop_trim'] / 3 * 3   
    
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    plot_upstream_1   = str(plot_upstream)  +'_'
    plot_downstream_1 = str(plot_downstream)+'_'
    start_trim_1      = str(start_trim)     +'_'
    stop_trim_1       = str(stop_trim)      +'_'
    frameshift_1      = str(frameshift)     +'_'
    a_site_1          = str(A_site)         +'_'
    
    name_settings = minlength_1+maxlength_1+plot_upstream_1+plot_downstream_1
    name_settings += start_trim_1+stop_trim_1+frameshift_1
        
    # import genetic code
    aa_code, codon_code = ribo_util.get_genetic_code()
    
    
    '''output data structure'''
    
    # aa/codon avgplot      = { aa: {length: [values]}} 
    # aa/codon count        = { aa: N}
    # aa/codon score        = { aa: [aa_list], A site score: [A_list]...
                                
    aa_avgplot  = {}
    aa_count    = {}
    aa_score    = {}
    aa_score_df = {}
    aa_list     = []
    aa_A_list   = []
    aa_P_list   = []
    aa_E_list   = []
    aa_1_list   = []
    aa_2_list   = []
    aa_m1_list   = []
    aa_m2_list   = []
    
    codon_avgplot  = {}
    codon_count    = {}
    codon_score    = {}
    codon_score_df = {}
    codon_list     = []
    codon_A_list   = []
    codon_P_list   = []
    codon_E_list   = []
    codon_1_list   = []
    codon_2_list   = []
    codon_m1_list   = []
    codon_m2_list   = []
    
    # create empty datastructures to store density info:
    for aa in aa_code.keys():
        aa_avgplot[aa] = {length : [0]*(plotlength) for length in lengthindex}
        aa_count[aa]   = 0
        aa_score[aa]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
        
    for codon in codon_code.keys():
        codon_avgplot[codon] = {length : [0]*(plotlength) for length in lengthindex}
        codon_count[codon]   = 0
        codon_score[codon]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
    
    '''genes in data''' 
    
    
    # count included genes
    with open(path_pausescore + 'included_genes', 'r') as inf:
        included_genes = pickle.load(inf)[1]
    
    
    
    #list_location = '/Volumes/HDD/Ribo_seq/libraries/analysis/reference_information/Workbook1.csv'
    
    #infile = pd.read_csv(list_location)
    #TE_list = infile['yeaR'].tolist()
        
    '''iterate through every annotated gene to get codon density info:'''
    
    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,stop_list, strand_list, seq_list):
        
        if not alias in included_genes:
            continue
        
        
        ''' define start and stop positions for codons to analyze:
        
        # = codon 
        #################################### = GFF sequence  (50 extra nt)
           ##############################    = AUG to UGA
             ##########################      = density to analyze : remove start and stop peaks
                 ##################          = codons to analyze : remove plot window
        
        '''
        
        # First, define density without start and stop peaks:
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim + frameshift
            density_stop  = stop  - stop_trim + frameshift
            
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim - frameshift
            density_stop  = stop  + stop_trim - 3 + frameshift
            
            period = -1
        
        # GFF seq has 50 extra nucleotides, so remove:
        # Also remove several codons from start and stop positions:
        
        codon_seq_start = gff_extra + start_trim + plot_upstream + frameshift
        codon_seq_stop  = -gff_extra - stop_trim - plot_downstream + frameshift
            
        seq       = sequence[codon_seq_start : codon_seq_stop]
        seqlength = len(seq)
        
    
        # make a list of codons in the sequence:
        codons_seq = [seq[i:i+3] for i in range(0, seqlength, 3)]
        aa_seq     = [codon_code[codon] for codon in codons_seq]
        
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        # fill density dict with density info
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            length_density_float = [float(i) for i in length_density]
            gene_density[length] = length_density
            
        
        relative_density = gene_density
            
        
        # add data to dataframe
        codon_index = 0
        for codon in codons_seq:
            #if codon == 'TAA':
            #    print "kekekeek"
            codon_position = (codon_index * 3)
            aa = aa_seq[codon_index]
            if aa == '_':
                print 'Stop codon found in'
                print alias
            codon_index        += 1
            
            '''if codon_index > 50:
                continue
            '''
            codon_count[codon] += 1
            aa_count[aa]       += 1
            
            for length in lengthindex:
                for position in positionindex:
                    
                    density_position = codon_position + position
                    density          = relative_density[length][density_position]
                    
                    codon_avgplot[codon][length][position] += density
                    aa_avgplot[aa][length][position]       += density
                    
        
    
    #divide data by total instances of the codon or aa
    for codon in codon_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if codon_count[codon] == 0:
                    continue 
                else: 
                    codon_avgplot[codon][length][position] = codon_avgplot[codon][length][position] / codon_count[codon]
    
    for aa in aa_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if aa_count[aa] == 0:
                    continue
                else: 
                    aa_avgplot[aa][length][position] = aa_avgplot[aa][length][position] / aa_count[aa]
                
   #Convert data for plotting and csv
    
    codon_data     = {}
    codon_data_sum = {}
    aa_data     = {}
    aa_data_sum = {}
    
    A_site_shift = plot_upstream - A_site 
    P_site_shift = plot_upstream - P_site 
    E_site_shift = plot_upstream - E_site
    one_site_shift = plot_upstream - one_site 
    two_site_shift = plot_upstream - two_site 
    mone_site_shift = plot_upstream - mone_site 
    mtwo_site_shift = plot_upstream - mtwo_site 
    
    for codon in codon_avgplot.keys():
        df = pd.DataFrame(codon_avgplot[codon])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        codon_data[codon]     = df
        codon_data_sum[codon] = df.sum(0)
        
        codon_score[codon][0] = sum(codon_data_sum[codon][A_site_shift: A_site_shift + 3:1])/3
        codon_score[codon][1] = sum(codon_data_sum[codon][P_site_shift: P_site_shift + 3:1])/3                            
        codon_score[codon][2] = sum(codon_data_sum[codon][E_site_shift: E_site_shift + 3:1])/3
        codon_score[codon][3] = sum(codon_data_sum[codon][one_site_shift: one_site_shift + 3:1])/3
        codon_score[codon][4] = sum(codon_data_sum[codon][two_site_shift: two_site_shift + 3:1])/3                            
        codon_score[codon][5] = sum(codon_data_sum[codon][mone_site_shift: mone_site_shift + 3:1])/3
        codon_score[codon][6] = sum(codon_data_sum[codon][mtwo_site_shift: mtwo_site_shift + 3:1])/3
        
        codon_list.append(codon)
        codon_A_list.append(codon_score[codon][0])
        codon_P_list.append(codon_score[codon][1])
        codon_E_list.append(codon_score[codon][2])
        codon_1_list.append(codon_score[codon][3])
        codon_2_list.append(codon_score[codon][4])
        codon_m1_list.append(codon_score[codon][5])
        codon_m2_list.append(codon_score[codon][6])

    for aa in aa_avgplot.keys():
        df = pd.DataFrame(aa_avgplot[aa])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        aa_data[aa]     = df
        aa_data_sum[aa] = df.sum(0) 
                                    
        aa_score[aa][0] = sum(aa_data_sum[aa][A_site_shift: A_site_shift + 3:1])/3    
        aa_score[aa][1] = sum(aa_data_sum[aa][P_site_shift: P_site_shift + 3:1])/3
        aa_score[aa][2] = sum(aa_data_sum[aa][E_site_shift: E_site_shift + 3:1])/3
        aa_score[aa][3] = sum(aa_data_sum[aa][one_site_shift: one_site_shift + 3:1])/3    
        aa_score[aa][4] = sum(aa_data_sum[aa][two_site_shift: two_site_shift + 3:1])/3
        aa_score[aa][5] = sum(aa_data_sum[aa][mone_site_shift: mone_site_shift + 3:1])/3 
        aa_score[aa][6] = sum(aa_data_sum[aa][mtwo_site_shift: mtwo_site_shift + 3:1])/3 

        aa_list.append(aa)
        aa_A_list.append(aa_score[aa][0])
        aa_P_list.append(aa_score[aa][1])
        aa_E_list.append(aa_score[aa][2])
        aa_1_list.append(aa_score[aa][3])
        aa_2_list.append(aa_score[aa][4])
        aa_m1_list.append(aa_score[aa][5])
        aa_m2_list.append(aa_score[aa][6])
        
    codon_score_df['Codon']  = codon_list
    codon_score_df['A_site'] = codon_A_list
    codon_score_df['P_site'] = codon_P_list
    codon_score_df['E_site'] = codon_E_list
    codon_score_df['1_site'] = codon_1_list
    codon_score_df['2_site'] = codon_2_list
    codon_score_df['-1_site'] = codon_m1_list
    codon_score_df['-2_site'] = codon_m2_list

    aa_score_df['Amino Acid'] = aa_list
    aa_score_df['A_site']     = aa_A_list
    aa_score_df['P_site']     = aa_P_list
    aa_score_df['E_site']     = aa_E_list
    aa_score_df['1_site']     = aa_1_list
    aa_score_df['2_site']     = aa_2_list
    aa_score_df['-1_site']     = aa_m1_list
    aa_score_df['-2_site']     = aa_m2_list

    codon_df = pd.DataFrame(codon_score_df)
    codon_df.to_csv(path_pausescore + 'codon_scores.csv')
    
    aa_df = pd.DataFrame(aa_score_df)
    aa_df.to_csv(path_pausescore + 'aa_scores.csv')
    
    ribo_util.makePickle(codon_data,     path_pausescore +'codon_HM_data'  + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(codon_data_sum, path_pausescore +'codon_plot_data'+ name_settings , protocol=pickle.HIGHEST_PROTOCOL)  
    ribo_util.makePickle(codon_score_df, path_pausescore +'codon_scores'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_data,        path_pausescore +'aa_HM_data'     + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(aa_data_sum,    path_pausescore +'aa_plot_data'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_score_df,    path_pausescore +'aa_scores'      + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
   
    return 


def logscore(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    multi     = inputs['multiprocess']
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started pause score analysis at " + str(datetime.now())

    for fname in files:
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/'
        if not os.path.exists(path_pausescore):
            os.makedirs(path_pausescore)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
       
        if not multi == 'yes':
            run_logscore(fname, settings, plus, minus, gff_dict, path_pausescore)
        else:     
            argument = [fname, settings, plus, minus, gff_dict, path_pausescore]
            arguments.append(argument)
    
    if multi == 'yes':
        ribo_util.multiprocess(run_logscore, arguments, threads)
    
    print "Finished pause score analysis at " + str(datetime.now())
    
    return


def loadlargePickles(inputs, settings, paths_in, paths_out):
    '''load GFF and density dictionaries into memory to speed up analysis
    WARNING :::: do not run to many samples at once!!!!!'''
    
    files     = inputs['files']
    alignment = settings['alignment']
    minlength = settings['minlength']
    maxlength = settings['maxlength']
    
    
    #load gff
    path_gff_dict = paths_in['path_gff_dict']
    gff_dict = unPickle(path_gff_dict) 
    
    #load density as a larger dict with fname as keys
    path_density = paths_out['path_density']
    plus_dict  = {}
    minus_dict = {}
    
    for fname in files: 
        path_den_plus  = path_density + fname + "/plus_lc"
        path_den_minus = path_density + fname + "/minus_lc"
        density_plus   = unPickle(path_den_plus)
        density_minus  = unPickle(path_den_minus)
        
        
        plus_dict[fname]  = density_plus
        minus_dict[fname] = density_minus

    return gff_dict, plus_dict, minus_dict

def unPickle(path_pickle):
    #returns the pickled object stored in a pickle file
    f = open(path_pickle, 'r')
    data = pickle.load(f)
    f.close()
    return data



def run_included(fname, settings, plus, minus, gff, path_pausescore):
    
    '''define variables'''
    
    minlength   = settings['minlength']
    maxlength   = settings['maxlength']
    alignment   = settings['alignment']
    threshold   = settings['threshold']


    lengthindex = range(minlength, maxlength + 1)
    
    A_site = settings['A_site shift']
    P_site = A_site - 3
    E_site = A_site - 6
    two_site = A_site + 6
    one_site = A_site + 3
    mone_site = A_site - 9
    mtwo_site = A_site - 12
    
    
    frameshift      = settings['frameshift']

    plot_upstream   = settings['plot_upstream'] / 3 * 3        #change window to interval of 3
    plot_downstream = settings['plot_downstream'] / 3 * 3
    
    next_codon = settings['next_codon']
    
    
    # define plot length
    plotlength = plot_upstream + plot_downstream + 1
    positionindex = range(0, plotlength)
    
    # load density files
    density_plus  = plus
    density_minus = minus 
    
    # load annotation
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    seq_list    = gff_dict['Sequence'] 
 
    gff_extra  = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    start_trim = settings['start_trim'] / 3 * 3
    stop_trim  = settings['stop_trim'] / 3 * 3   
    
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    plot_upstream_1   = str(plot_upstream)  +'_'
    plot_downstream_1 = str(plot_downstream)+'_'
    start_trim_1      = str(start_trim)     +'_'
    stop_trim_1       = str(stop_trim)      +'_'
    frameshift_1      = str(frameshift)     +'_'
    a_site_1          = str(A_site)         +'_'
    
    name_settings = minlength_1+maxlength_1+plot_upstream_1+plot_downstream_1
    name_settings += start_trim_1+stop_trim_1+frameshift_1
        
    # import genetic code
    aa_code, codon_code = ribo_util.get_genetic_code()
    
    
    '''output data structure'''
    
    # aa/codon avgplot      = { aa: {length: [values]}} 
    # aa/codon count        = { aa: N}
    # aa/codon score        = { aa: [aa_list], A site score: [A_list]...
                                
    aa_avgplot  = {}
    aa_count    = {}
    aa_score    = {}
    aa_score_df = {}
    aa_list     = []
    aa_A_list   = []
    aa_P_list   = []
    aa_E_list   = []
    aa_1_list   = []
    aa_2_list   = []
    aa_m1_list   = []
    aa_m2_list   = []
    
    codon_avgplot  = {}
    codon_count    = {}
    codon_score    = {}
    codon_score_df = {}
    codon_list     = []
    codon_A_list   = []
    codon_P_list   = []
    codon_E_list   = []
    codon_1_list   = []
    codon_2_list   = []
    codon_m1_list   = []
    codon_m2_list   = []
    
    # create empty datastructures to store density info:
    for aa in aa_code.keys():
        aa_avgplot[aa] = {length : [0]*(plotlength) for length in lengthindex}
        aa_count[aa]   = 0
        aa_score[aa]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
        
    for codon in codon_code.keys():
        codon_avgplot[codon] = {length : [0]*(plotlength) for length in lengthindex}
        codon_count[codon]   = 0
        codon_score[codon]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
    
    '''genes in data''' 
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']         = [0, []]
    excluded_genes['low_density']   = [0, []]
    excluded_genes['not_divisible'] = [0, []]
    # count included genes
    included_genes = [0, []]
    
    
    
    #list_location = '/Volumes/HDD/Ribo_seq/libraries/analysis/reference_information/Workbook1.csv'
    
    #infile = pd.read_csv(list_location)
    #TE_list = infile['yeaR'].tolist()
        
    '''iterate through every annotated gene to get codon density info:'''
    
    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,stop_list, strand_list, seq_list):
        
        #if not alias in TE_list:
            #continue
        
        
        ''' define start and stop positions for codons to analyze:
        
        # = codon 
        #################################### = GFF sequence  (50 extra nt)
           ##############################    = AUG to UGA
             ##########################      = density to analyze : remove start and stop peaks
                 ##################          = codons to analyze : remove plot window
        
        '''
        
        # First, define density without start and stop peaks:
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim + frameshift
            density_stop  = stop  - stop_trim + frameshift
            
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim - frameshift
            density_stop  = stop  + stop_trim - 3 + frameshift
            
            period = -1
        
        # GFF seq has 50 extra nucleotides, so remove:
        # Also remove several codons from start and stop positions:
        
        codon_seq_start = gff_extra + start_trim + plot_upstream + frameshift
        codon_seq_stop  = -gff_extra - stop_trim - plot_downstream + frameshift
            
        seq       = sequence[codon_seq_start : codon_seq_stop]
        seqlength = len(seq)
        
        # exclude genes that are not divisable by 3
        if seqlength % 3 != 0:
            excluded_genes['not_divisible'][0] += 1
            excluded_genes['not_divisible'][1].append(alias)
            continue
        
        # exclude genes shorter than plot
        if seqlength < plotlength + 1:
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        elif codon_seq_start - codon_seq_stop > len(sequence):
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
    
        # make a list of codons in the sequence:
        codons_seq = [seq[i:i+3] for i in range(0, seqlength, 3)]
        aa_seq     = [codon_code[codon] for codon in codons_seq]
        
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        # fill density dict with density info
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            length_density_float = [float(i) for i in length_density]
            gene_density[length] = length_density
            total_density        = [x + y for x, y in itertools.izip(total_density, length_density)]
            
        gene_avgreads = float(sum(total_density)) / float(genelength)

            
        # normalize density by total gene density
        if gene_avgreads * 3 < threshold:
            excluded_genes['low_density'][0] += 1
            excluded_genes['low_density'][1].append(alias)
            continue
            
        else: 
            included_genes[0] += 1
            included_genes[1].append(alias)
        
    #print excluded_genes['short'][0]    
    #print excluded_genes['low_density'][0]
    #print excluded_genes['not_divisible'][0]
    #print included_genes[0]
    print path_pausescore + 'included_genes'
    ribo_util.makePickle(included_genes, path_pausescore + 'included_genes')
    return 


def included(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    multi     = inputs['multiprocess']
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started pause score analysis at " + str(datetime.now())

    for fname in files:
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/'
        if not os.path.exists(path_pausescore):
            os.makedirs(path_pausescore)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
       
        if not multi == 'yes':
            run_included(fname, settings, plus, minus, gff_dict, path_pausescore)
        else:     
            argument = [fname, settings, plus, minus, gff_dict, path_pausescore]
            arguments.append(argument)
    
    if multi == 'yes':
        ribo_util.multiprocess(run_included, arguments, threads)
    
    print "Finished pause score analysis at " + str(datetime.now())
    
    return

def asymmetry(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started asymmetry analysis at " + str(datetime.now())

    for fname in files:
        path_analysis = paths_out['path_analysis'] + fname + '/'
        if not os.path.exists(path_analysis):
            os.makedirs(path_analysis)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
        
        run_asymmetry(fname, settings, plus, minus, gff_dict, path_analysis)
        argument = [fname, settings, plus, minus, gff_dict, path_analysis]
        arguments.append(argument)
    
    #ribo_util.multiprocess(run_asymmetry, arguments, threads)
    print "Finished asymmetry analysis at " + str(datetime.now())

def run_asymmetry(fname, settings, plus, minus, gff, path_analysis): 
    
    minlength    = settings['minlength']
    maxlength    = settings['maxlength']
    threshold    = settings['threshold']
    subgroup_csv = settings['subgroup']
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    
    name_settings = minlength_1+maxlength_1
    
    density_plus  = plus
    density_minus = minus 
    gff_dict      = gff
    
    start_trim = settings['start_trim']# / 3 * 3
    stop_trim  = settings['stop_trim']# / 3 * 3

    if not subgroup_csv == 'none':
        subgroup       = pd.read_csv(subgroup_csv)
        subgroup       = subgroup.to_dict(orient='list')
        subgroup_alias = subgroup['Alias']
    else:
        subgroup_alias = ['none']

    lengthindex = range(minlength, maxlength+1)
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    
    asymmetry_alias      = []
    asymmetry_subgroup   = []
    asymmetry_score      = []
    asymmetry_reads      = []
    asymmetry_genelength = []
    asymmetry_dict  = {}
    genes_excluded  = []
    
    #bins = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
    #total_reads = 0

    result = {}

    import pickle
    with open(path_analysis + 'pause_score/included_genes.csv', 'r') as inf:
        included = pickle.load(inf)[1]
    print included
    
    for alias, start, stop, strand in itertools.izip(alias_list, start_list,stop_list, strand_list):  
        
        
        if alias not in included:
            continue  
        
        
        if strand == '+':
            density_dict = density_plus
            startp = start + start_trim
            stopp  = stop  - stop_trim            
            period = 1
            genelength = abs(startp - stopp) + 1 
            tenth_gene = genelength / 10
            p10 = startp + tenth_gene
            p20 = p10 + tenth_gene
            p30 = p20 + tenth_gene
            p40 = p30 + tenth_gene
            p50 = p40 + tenth_gene
            p60 = p50 + tenth_gene
            p70 = p60 + tenth_gene
            p80 = p70 + tenth_gene
            p90 = p80 + tenth_gene
            
        elif strand == '-':
            density_dict = density_minus
            startp = start - start_trim
            stopp  = stop  + stop_trim
            period = -1
            genelength = abs(startp - stopp) + 1
            tenth_gene = genelength / 10
            p10 = startp - tenth_gene
            p20 = p10 - tenth_gene
            p30 = p20 - tenth_gene
            p40 = p30 - tenth_gene
            p50 = p40 - tenth_gene
            p60 = p50 - tenth_gene
            p70 = p60 - tenth_gene
            p80 = p70 - tenth_gene
            p90 = p80 - tenth_gene
        

        #if genelength < 100: 
        #   genes_excluded.append(alias)
        #    continue
        
        b10 = 0.
        b20 = 0.
        b30 = 0.
        b40 = 0.
        b50 = 0.
        b60 = 0.
        b70 = 0.
        b80 = 0.
        b90 = 0.
        b100 = 0.   
        
        reads = 0.

        for length in lengthindex:
            
            #first  = float(sum(density_dict[length][startp: midpoint + 1: period]))
            #second = density_dict[length][midpoint: stopp  + 1: period]
            genelength = len(density_dict[length][start: stop + 1: period])
            
            #first  = float(sum(first))
            #second = float(sum(second))
            
            b10  += float(sum(density_dict[length][startp: p10 + 1: period]))
            b20  += float(sum(density_dict[length][p10: p20 + 1: period]))
            b30  += float(sum(density_dict[length][p20: p30 + 1: period]))
            b40  += float(sum(density_dict[length][p30: p40 + 1: period]))
            b50  += float(sum(density_dict[length][p40: p50 + 1: period]))
            b60  += float(sum(density_dict[length][p50: p60 + 1: period]))
            b70  += float(sum(density_dict[length][p60: p70 + 1: period]))
            b80  += float(sum(density_dict[length][p70: p80 + 1: period]))
            b90  += float(sum(density_dict[length][p80: p90 + 1: period]))
            b100  += float(sum(density_dict[length][p90: stopp + 1: period]))
                                    
        reads = b10 + b20 + b30 + b40 + b50 + b60 + b70 + b80 + b90 + b100
        #print reads
                
        #if reads < 150:
        #    genes_excluded.append(alias)
        #    continue
                                
       
        #bins += np.array([b10, b20, b30, b40, b50, b60, b70, b80, b90, b100])
        #total_reads += reads

        result[alias] = [np.array([b10, b20, b30, b40, b50, b60, b70, b80, b90, b100]), reads, genelength]
    #print total_reads
    #print bins
    #bins = bins/float(total_reads)*100
    #ribo_util.makePickle(bins, path_analysis + 'asymmetry_per_cent', protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(result, path_analysis + 'asymmetry_per_cent_solo', protocol=pickle.HIGHEST_PROTOCOL)
