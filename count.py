#!/HPCTMP_NOBKUP/yl1217/local/anaconda2/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 16:11:45 2018

@author: Yu Liu
"""

#%%Config
import os

#Read in config file
parent_dir = os.getcwd()
conf = os.path.join(parent_dir, 'config')
configfile = open(conf, 'r')
config0 = configfile.readlines()
configfile.close()
del configfile
config = config0[:-1]
config1 = config0[-1]
del config0

if config1[-1] != '\n':
    config1 = config1 + '\n'
else:
    config1 = config1
config.append(config1)
del config1

#Parse config parameters
envdic = {}
for line in config:
    line = line[:-1]
    name = line.strip().split('=')[0].strip()
    env = line.strip().split('=')[1].strip()
    envdic[name] = env
    
if not envdic['countreads'].upper() == 'YES':
    exit
else:
    project = envdic['project']
    samtools = envdic['samtools']
    pos_f = envdic['pos_f']
    result_tag = envdic['result_tag']
    width1 = float(envdic['width1'])
    width2 = float(envdic['width2'])
    
del envdic

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import HTSeq
import subprocess
import re

#Prepare folders
initial_fastq_dir = os.path.join(parent_dir, project)
work_dir = os.path.join(initial_fastq_dir, 'data')

result_dir = os.path.join(initial_fastq_dir, 'result')
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
    
#Define sample list
sample_names = []
for each in os.listdir(work_dir):
    if re.search('\.bam$', each):
        sample_names.append(each)

#%% annotate_coordinates
os.chdir(work_dir)

symbols = set()
symbols_dic = {}
pos_file = open(pos_f, 'r')
pos_content = pos_file.readlines()[1:]

for each in pos_content:
    eachList = each.strip().split('\t')
    name = eachList[0]
    chrom = eachList[1]
    strand = eachList[2]
    start = int(eachList[3])
    end = int(eachList[4])
    iv = HTSeq.GenomicInterval(chrom, start, end, strand)
    symbols.add(iv)
    symbols_dic[name] = iv
    
tssnum = len(symbols)
ttsnum = len(symbols)


#%% read depth

totalreads = []
for each in sample_names:
    n = subprocess.check_output(samtools + ' view -c ' + \
                                each, shell=True)[:-1]
    totalreads.append(int(n))

data_dict_up_fwd = {}
data_dict_up_rev = {}
data_dict_dn_fwd = {}
data_dict_dn_rev = {}

#%% change artifical reads to true reads
for index in range(0, len(sample_names)):
    
    profile_up_fwd = np.zeros(2*width1, dtype = 'd')
    profile_dn_fwd = np.zeros(2*width2, dtype = 'd')
    
    profile_up_rev = np.zeros(2*width1, dtype = 'd')
    profile_dn_rev = np.zeros(2*width2, dtype = 'd')
    
    sampleid = sample_names[index].strip().split('.')[0]
    fwdfilename = sampleid + '.FWD'
    revfilename = sampleid + '.REV'
    
    filename = sample_names[index]
    N = totalreads[index]
    factor = (10.0**6)/N
    bamfile = HTSeq.BAM_Reader(filename)
    coverage = HTSeq.GenomicArray('auto', stranded = True, typecode = 'd')
    for almnt in bamfile:
        if almnt.pe_which == 'first':
            newiv = almnt.iv
            try:
                coverage[newiv] += factor
            except IndexError:
                trimmed_newiv = HTSeq.GenomicInterval(newiv.chrom, 0, \
                                                      newiv.end, \
                                                      newiv.strand)
                coverage[trimmed_newiv] += factor
        else:
            chrom = almnt.iv.chrom
            start = almnt.iv.start
            end = almnt.iv.end
            strandori = almnt.iv.strand
            if strandori == '+':
                strand = '-'
            else:
                strand = '+'
            newiv = HTSeq.GenomicInterval(chrom, start, end, strand)
            try:
                coverage[newiv] += factor
            except IndexError:
                trimmed_newiv = HTSeq.GenomicInterval(newiv.chrom, 0, \
                                                      newiv.end, \
                                                      newiv.strand)
                coverage[trimmed_newiv] += factor

    
    tsspos = set()
    ttspos = set()

###define count regions    
    for iv in symbols:
        chrom = iv.chrom
        startpos = iv.start_d_as_pos.pos
        endpos = iv.end_d_as_pos.pos
        strand = iv.strand

        p = HTSeq.GenomicPosition(chrom, startpos, strand)
        q = HTSeq.GenomicPosition(chrom, endpos, strand)

        if strand == '+':
            window = HTSeq.GenomicInterval(p.chrom, p.pos - width1, \
                                           p.pos + width1, '+')
            window_rev = HTSeq.GenomicInterval(p.chrom, p.pos - width1, \
                                               p.pos + width1, '-')
            except_window = HTSeq.GenomicInterval(p.chrom, 0, \
                                                  p.pos + width1, '+')
            except_window_rev = HTSeq.GenomicInterval(p.chrom, 0, \
                                                      p.pos + width1, '-')
            
            windowq = HTSeq.GenomicInterval(q.chrom, q.pos - 1 - width2, \
                                            q.pos - 1 + width2, '+')
            windowq_rev = HTSeq.GenomicInterval(q.chrom, q.pos - 1 - width2, \
                                                q.pos - 1 + width2, '-')
            except_windowq = HTSeq.GenomicInterval(q.chrom, 0, \
                                                   q.pos - 1 + width2, '+')
            except_windowq_rev = HTSeq.GenomicInterval(q.chrom, 0, \
                                                       q.pos - 1 + width2, '-')
            
            
            
        else:
            window = HTSeq.GenomicInterval(p.chrom, p.pos - width1, \
                                           p.pos + width1, '-')
            window_rev = HTSeq.GenomicInterval(p.chrom, p.pos - width1, \
                                               p.pos + width1, '+')
            except_window = HTSeq.GenomicInterval(p.chrom, 0, \
                                                  p.pos + width1, '-')
            except_window_rev = HTSeq.GenomicInterval(p.chrom, 0, \
                                                      p.pos + width1, '+')
            
            windowq = HTSeq.GenomicInterval(q.chrom, q.pos + 1 - width2, \
                                            q.pos + 1 + width2, '-')
            windowq_rev = HTSeq.GenomicInterval(q.chrom, q.pos + 1 - width2, \
                                            q.pos + 1 + width2, '+')
            except_windowq = HTSeq.GenomicInterval(q.chrom, 0, \
                                                   q.pos + 1 + width2, '-')
            except_windowq_rev = HTSeq.GenomicInterval(q.chrom, 0, \
                                                       q.pos + 1 + width2, '+')
            
#######################################TSS region 

###fwd file
        try:
            wincvg_fwd = np.fromiter(coverage[window], dtype = 'd', \
                                     count = 2*width1)
        except IndexError:
            try:
                wincvg_fwd = np.fromiter(coverage[except_window], dtype = 'd', \
                                         count = 2*width1)
            except ValueError:
                tssnum = tssnum - 1
                pass
###rev file            
        try:
            wincvg_rev = np.fromiter(coverage[window_rev], dtype = 'd', \
                                     count = 2*width1)
        except IndexError:
            try:
                wincvg_rev = np.fromiter(coverage[except_window_rev], dtype = 'd', \
                                         count = 2*width1)
            except ValueError:
                pass
            
###################################TTS region
###fwd file
        try:
            wincvgq_fwd = np.fromiter(coverage[windowq], dtype = 'd', \
                                      count = 2*width2)
        except IndexError:
            try:
                wincvgq_fwd = np.fromiter(coverage[except_windowq], \
                                          dtype = 'd', \
                                          count = 2*width2)
            except ValueError:
                ttsnum = ttsnum - 1
                pass
###rev file
        try:
            wincvgq_rev = np.fromiter(coverage[windowq_rev], dtype = 'd', \
                                      count = 2*width2)
        except IndexError:
            try:
                wincvgq_rev = np.fromiter(coverage[except_windowq_rev], \
                                          dtype = 'd', \
                                          count = 2*width2)
            except ValueError:
                pass

##################integrate
        if strand == '+':
            tssfwd = wincvg_fwd
            ttsfwd = wincvgq_fwd
            tssrev = wincvg_rev
            ttsrev = wincvgq_rev
        else:
            tssfwd = wincvg_fwd[::-1]
            ttsfwd = wincvgq_fwd[::-1]
            tssrev = wincvg_rev[::-1]
            ttsrev = wincvgq_rev[::-1]
            
        profile_up_fwd += tssfwd
        profile_up_rev += tssrev
        profile_dn_fwd += ttsfwd
        profile_dn_rev += ttsrev
###        
    profile_up_fwd_mean = profile_up_fwd/tssnum
    profile_up_rev_mean = profile_up_rev/tssnum
    profile_dn_fwd_mean = profile_dn_fwd/ttsnum
    profile_dn_rev_mean = profile_dn_rev/ttsnum
    
    up_fwd_s = Series(profile_up_fwd_mean)
    up_rev_s = Series(profile_up_rev_mean)
    dn_fwd_s = Series(profile_dn_fwd_mean)
    dn_rev_s = Series(profile_dn_rev_mean)


    data_dict_up_fwd[fwdfilename] = up_fwd_s
    data_dict_up_rev[revfilename] = up_rev_s
    data_dict_dn_fwd[fwdfilename] = dn_fwd_s
    data_dict_dn_rev[revfilename] = dn_rev_s
    
    temp_up_fwd = DataFrame(data_dict_up_fwd[fwdfilename])
    temp_up_fwd['sample'] = fwdfilename
    if index == 0:
        final_up_fwd = temp_up_fwd
    else:
        final_up_fwd = pd.concat([final_up_fwd, temp_up_fwd], axis = 1)
    
    temp_up_rev = DataFrame(data_dict_up_rev[revfilename])
    temp_up_rev['sample'] = revfilename
    if index == 0:
        final_up_rev = temp_up_rev
    else:
        final_up_rev = pd.concat([final_up_rev, temp_up_rev], axis = 1)
    
    temp_dn_fwd = DataFrame(data_dict_dn_fwd[fwdfilename])
    temp_dn_fwd['sample'] = fwdfilename
    if index == 0:
        final_dn_fwd = temp_dn_fwd
    else:
        final_dn_fwd = pd.concat([final_dn_fwd, temp_dn_fwd], axis = 1)
    
    temp_dn_rev = DataFrame(data_dict_dn_rev[revfilename])
    temp_dn_rev['sample'] = revfilename
    if index == 0:
        final_dn_rev = temp_dn_rev
    else:
        final_dn_rev = pd.concat([final_dn_rev, temp_dn_rev], axis = 1)

os.chdir(result_dir)

final_up_fwd.to_csv('sense' + result_tag + '_tss.htseq.txt', 
                    sep = '\t', index = False)
final_up_rev.to_csv('antisense' + result_tag + '_tss.htseq.txt', 
                    sep = '\t', index = False)
final_dn_fwd.to_csv('sense' + result_tag + '_tts.htseq.txt', 
                    sep = '\t', index = False)
final_dn_rev.to_csv('antisense' + result_tag + '_tts.htseq.txt', 
                    sep = '\t', index = False)

