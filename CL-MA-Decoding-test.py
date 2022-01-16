import random
import copy
import time
import os
import re
import operator
import deBruijnGraph

from utils import *
from test_utils import *
from DNAdroplet import DNADroplet
from DNAfountain import DNAFountain
from glass import Glass


work_dir = r'/work/1.DeBruijnGraphDecoding/6MePCR/1.rawdata/YC10_5_4_BDDP210000406-1A/'  #
sim_file = r'/work/1.DeBruijnGraphDecoding/6M/1.rawdata/' + r'6.5MB.DNAs.newids.tab.nohead.fa.noPP.tab' # The original strand sequences for calculating of Strand recovery rate.
input_file = work_dir + r'out.extendedFrags.fastq.cl.d4'  # Clustered sequences by Starcode.

muscle = r'muscle '  #Location of muscle program
mafft = r'mafft '


max_seq_copy = 50 # Maximal number of sequence copies for alignments.
min_seq_copy = 3  # Minimal number of sequence copies for alignments.
switch_align_seqnum = 28 # Swithch aligner

os.system('mkdir ' + work_dir + 'muscle')

dps_seqs = read_sim(sim_file)
res = {}


print('Reading starcode clusters....')
clu_seqs = read_starcode_clusters(input_file, min_seq_copy)

for clu in clu_seqs:
    clu_seqs_file = work_dir + 'muscle/' + str(clu)
    fo = open(clu_seqs_file, "tw")
    i=1
    for seq in clu_seqs[clu]:
        if i > max_seq_copy:
            break
        fo.write(">" + str(i) + "\n")
        fo.write(seq)
        fo.write("\n")
        i = i + 1
    fo.close()


res['align'] = []
res['cons'] = []


print('Running Alignment 1 ....')
a = time.perf_counter()
for clu in clu_seqs:
    clu_seqs_file = work_dir + 'muscle/' + str(clu)
    if len(clu_seqs[clu]) >= min_seq_copy:
        if len(clu_seqs[clu]) < switch_align_seqnum:
            #os.system(muscle + r' -align ' + clu_seqs_file + r' -output ' + clu_seqs_file + r'.aln')
            os.system(muscle + r' -in ' + clu_seqs_file + r' -out ' + clu_seqs_file +  r'.aln')
        else:
            os.system(mafft + clu_seqs_file + r' > ' + clu_seqs_file + r'.aln')
res['align'].append(time.perf_counter() - a)


a = time.perf_counter()
clu_cons_seqs = []
for clu in clu_seqs:
    if len(clu_seqs[clu]) >= min_seq_copy:
        clu_seqs_file = work_dir + 'muscle/' + str(clu)
        clu_aln_file = clu_seqs_file +  r'.aln'
        reads = read_fasta(clu_aln_file)
        cons = majority_merge(reads)
        clu_cons_seqs.append(cons)
        
res['cons'].append(time.perf_counter() - a)


dec_num = check_cons(clu_cons_seqs, dps_seqs)
res['Sr'] = dec_num

