import random
import copy
import time
import os
import re
import operator

from utils import *
from test_utils import *
from DNAdroplet import DNADroplet
from DNAfountain import DNAFountain
from glass import Glass
from deBruijnGraph import DeBruijnGraph




check_index = True
number_of_droplets = 1
# exp_seq_copy_num = 10000
err_rate = 0.01
kmer_length = 18
data_block_length = 35
fountain_seed = 3
max_repeat_num_kmer = 5
reps = 1000
seq_copies = 20

primerF = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
primerE = 'CTGACACTGATGCATCCG'  # complement seq of P2

work_dir = r'/work/0.Simulations/3.DecodingSpeedCompare/lig_test/run' + str(fountain_seed) +  "/"

input_file = r'input_files/6.8MB.zip'

muscle = r'/home/lifu/Downloads/software/muscle3.8.31/muscle3.8.31_i86linux64 '

run = 'noCUt'

DBGPS  = r'/home/lifu/CodingProjects/switch-codes-switch-codes-kmer-cnt-DBG-master/DBGPS '

run_name = 'Speed Comparison'
file1 = open(input_file, 'rb')


fountain_init_index = 1
filebytes1 = file1.read()
file1.close()


# fdna1 = DNAFountain(filebytes1,data_block_length)
fdna1 = DNAFountain(filebytes1, data_block_length, fountain_init_index, fountain_seed)
fdna1.gen_degrees()

# dna_file = input_file + run_name + r'.eDNAs'
# file2 = open(dna_file,'tw')
#
# file3 = open(dna_file + r'.sim','tw')
#
#
# # def get_droplets_check_repeat_kmer(num, fdna1, kmer_len=21):
#
# # adrop = None
# j = 0
dps_seqs = {}
seq_arr = []

print('Generating Droplets')



i = 0
while i < number_of_droplets:
    dps = fdna1.DNAdroplet()
    dps_seqs[dps.head_index] = primerF + dps.to_DNA_CRC() + primerE
    seq_arr.append(primerF + dps.to_DNA_CRC() + primerE)
    i = i + 1

dna_hd = DNAHandler()

res= {}
res_mus = {}
for ligation_num in range(1, 11):
# for exp_seq_copy_num in [10]:
    #, 30, 40, 50, 60, 70, 80, 90, 100]:
#for exp_seq_copy_num in [10, 20, 30]:
    # exp_seq_copy_num = 15
    # res[exp_seq_copy_num] = {}
    cov_cut_off = cut_num(seq_copies)
    res[ligation_num] = {}

    suc = 0
    suc_mus = 0
    for a in range(0, reps):
        print('Generating eDNAs')
        eDNAs_all = []
        for id in dps_seqs:
            eDNAs = dna_hd.copy_randnum([primerF + dps_seqs[id] + primerE ], seq_copies, seq_copies)
            eDNAs = dna_hd.add_rand_sub_new(eDNAs, err_rate / 2)
            eDNAs = dna_hd.add_rand_ins_new(eDNAs, err_rate / 4)
            eDNAs = dna_hd.add_rand_del_new(eDNAs, err_rate / 4)
            eDNAs_all.extend(eDNAs)

            eDNAs = dna_hd.break_ligation(eDNAs, err_rate, ligation_num)

            fo = open(work_dir + 'muscle/' + str(id), "tw")
            i = 1
            for ee in eDNAs:
                fo.write(">" + str(i) + "\n")
                fo.write(ee + "\n")
                i = i + 1
            fo.close()



        print('Recovering with de Bruijn Graph')
        a = time.perf_counter()
        deG = DeBruijnGraph()
        deG.kmer_len = kmer_length
        deG.max_path_num = 1000
        print('Adding seqs')
        deG.add_seqs(eDNAs_all)
        deG.remove_low_cov_kmers(cov_cut_off)
        #deG.remove_low_cov_kmers(1)
        print('Recovery seqs')
        if test_deG1(dps_seqs, deG, data_block_length) > 0:
            suc = suc + 1



        for id in dps_seqs:
            id_seq_file = work_dir + 'muscle/' + str(id)
            os.system(muscle + r' -in ' + id_seq_file + r' -out ' + id_seq_file +  r'.aln')

        # res['mus']['muti-align'][exp_seq_copy_num].append(time.perf_counter() - a)

        i = 0
        for id in dps_seqs:
            id_seq_file = work_dir + 'muscle/' + str(id)
            id_aln_file = id_seq_file +  r'.aln'
            reads = read_fasta(id_aln_file)
            cons = majority_merge(reads)
            if dps_seqs[id] in cons:
                suc_mus = suc_mus + 1

    res[ligation_num][cov_cut_off] = suc
    res_mus[ligation_num] = suc_mus

