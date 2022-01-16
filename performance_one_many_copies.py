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
import deBruijnGraph




check_index = True
number_of_droplets = 1
exp_seq_copy_num = 1
err_rate = 0.03
kmer_length = 12
data_block_length = 35
fountain_seed = 1

reps = 1000

primerF = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
primerE = 'CTGACACTGATGCATCCG'  # complement seq of P2

# file1 = open(r'../Apollo program.kux', 'rb')
work_dir = r'/work/0.Simulations/3.DecodingSpeedCompare/run_one_many' + str(fountain_seed) + r'/'


# work_dir = r'E:\work\dec_speed' + '\\'
input_file = r'input_files/6.8MB.zip'



muscle = r'/home/lifu/Downloads/software/muscle3.8.31/muscle3.8.31_i86linux64 '

DBGPS  = r'/usr/local/bin/jellyfish '

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

res1 = {}
res2 = {}
res3 = {}

res = res1
res['deG'] = {}
res['mus'] = {}


res['deG']['add_seqs'] = {}
res['deG']['Sr'] = {}
res['deG']['Strand_decoding'] = {}
res['deG']['counting'] = {}
res['deG']['add_dump'] = {}
res['deG']['cov_cut_off'] = {}


res['mus']['muti-align'] = {}
res['mus']['consensus'] = {}
res['mus']['Sr'] = {}


# for exp_seq_copy_num in range(1:101):
for exp_seq_copy_num in range(5, 8, 1):
# for exp_seq_copy_num in [10, 15, 20, 30]:
    # exp_seq_copy_num = 15
    # res[exp_seq_copy_num] = {}

    res['deG']['add_seqs'][exp_seq_copy_num] = []
    res['deG']['Sr'][exp_seq_copy_num] = []
    res['deG']['Strand_decoding'][exp_seq_copy_num] = []
    res['deG']['counting'][exp_seq_copy_num] = []
    res['deG']['add_dump'][exp_seq_copy_num] = []
    res['deG']['cov_cut_off'][exp_seq_copy_num] = []


    res['mus']['muti-align'][exp_seq_copy_num] = []
    res['mus']['consensus'][exp_seq_copy_num] = []
    res['mus']['Sr'][exp_seq_copy_num] = []

    cov_cut_off = cut_num(exp_seq_copy_num)

    for rep in range(0, reps):
        print('Generating eDNAs')
        eDNAs_all = []
        for id in dps_seqs:
            eDNAs = dna_hd.copy_randnum([primerF + dps_seqs[id] + primerE ], exp_seq_copy_num, exp_seq_copy_num)
            eDNAs = dna_hd.add_rand_sub_new(eDNAs, err_rate / 2)
            eDNAs = dna_hd.add_rand_ins_new(eDNAs, err_rate / 4)
            eDNAs = dna_hd.add_rand_del_new(eDNAs, err_rate / 4)
            eDNAs_all.extend(eDNAs)
            fo = open(work_dir + 'muscle/' + str(id), "tw")
            i = 1
            for ee in eDNAs:
                fo.write(">" + str(i) + "\n")
                fo.write(ee + "\n")
                i = i + 1
            fo.close()


        fo = open(work_dir + 'deBru/all_eDNAs.fa',  "tw")
        i = 1
        for ee in eDNAs_all:
            fo.write(">" + str(i) + "\n")
            fo.write(ee + "\n")
            i = i + 1
        fo.close()


        # print('Recovering with de Bruijn Graph')
        #
        #
        # a = time.perf_counter()
        # deG = deBruijnGraph.DeBruijnGraph()
        # deG.kmer_len = kmer_length
        # deG.max_path_num = 2000
        # deG.add_seqs(eDNAs_all)
        # deG.remove_low_cov_kmers(cov_cut_off)
        #
        # res['deG']['cov_cut_off'][exp_seq_copy_num].append(cov_cut_off)
        # res['deG']['add_seqs'][exp_seq_copy_num].append(time.perf_counter() - a)
        #
        # a = time.perf_counter()
        # res['deG']['Sr'][exp_seq_copy_num].append(test_deG1(dps_seqs, deG, data_block_length))
        # res['deG']['Strand_decoding'][exp_seq_copy_num].append(time.perf_counter() - a)
        #
        #
        # a = time.perf_counter()
        #
        # eDNAs_file_loc = work_dir + 'deBru/all_eDNAs.fa'
        #
        # # os.system(DBGPS + r' -k ' + str(kmer_length) + r' -a ' + str() + eDNAs_file_loc + r' -o ' + eDNAs_file_loc + r'.count')
        # # res['deG']['counting'][exp_seq_copy_num].append(time.perf_counter() - a)


        print('Run multiple alignments')

        a = time.perf_counter()
        for id in dps_seqs:
            id_seq_file = work_dir + 'muscle/' + str(id)
            os.system(muscle + r' -in ' + id_seq_file + r' -out ' + id_seq_file +  r'.aln')

        res['mus']['muti-align'][exp_seq_copy_num].append(time.perf_counter() - a)


        a = time.perf_counter()
        i = 0
        for id in dps_seqs:
            id_seq_file = work_dir + 'muscle/' + str(id)
            id_aln_file = id_seq_file +  r'.aln'
            reads = read_fasta(id_aln_file)
            cons = majority_merge(reads)
            if dps_seqs[id] in cons:
                i = i + 1

        res['mus']['consensus'][exp_seq_copy_num].append(time.perf_counter() - a)
        res['mus']['Sr'][exp_seq_copy_num].append(i)


#
# res = res2
# res['deG'] = {}
# res['mus'] = {}
#
#
# res['deG']['add_seqs'] = {}
# res['deG']['Sr'] = {}
# res['deG']['Strand_decoding'] = {}
# res['deG']['counting'] = {}
# res['deG']['add_dump'] = {}
# res['deG']['cov_cut_off'] = {}
#
#
# res['mus']['muti-align'] = {}
# res['mus']['consensus'] = {}
# res['mus']['Sr'] = {}
#
#
# # for exp_seq_copy_num in range(1:101):
# for exp_seq_copy_num in range(3, 105, 2):
# # for exp_seq_copy_num in [10, 15, 20, 30]:
#     # exp_seq_copy_num = 15
#     # res[exp_seq_copy_num] = {}
#
#     res['deG']['add_seqs'][exp_seq_copy_num] = []
#     res['deG']['Sr'][exp_seq_copy_num] = []
#     res['deG']['Strand_decoding'][exp_seq_copy_num] = []
#     res['deG']['counting'][exp_seq_copy_num] = []
#     res['deG']['add_dump'][exp_seq_copy_num] = []
#     res['deG']['cov_cut_off'][exp_seq_copy_num] = []
#
#
#     res['mus']['muti-align'][exp_seq_copy_num] = []
#     res['mus']['consensus'][exp_seq_copy_num] = []
#     res['mus']['Sr'][exp_seq_copy_num] = []
#
#     cov_cut_off = cut_num(exp_seq_copy_num)
#
#     for rep in range(0, reps):
#         print('Generating eDNAs')
#         eDNAs_all = []
#         for id in dps_seqs:
#             eDNAs = dna_hd.copy_randnum([primerF + dps_seqs[id] + primerE ], exp_seq_copy_num, exp_seq_copy_num)
#             eDNAs = dna_hd.add_rand_sub_new(eDNAs, err_rate / 2)
#             eDNAs = dna_hd.add_rand_ins_new(eDNAs, err_rate / 4)
#             eDNAs = dna_hd.add_rand_del_new(eDNAs, err_rate / 4)
#             eDNAs_all.extend(eDNAs)
#             fo = open(work_dir + 'muscle/' + str(id), "tw")
#             i = 1
#             for ee in eDNAs:
#                 fo.write(">" + str(i) + "\n")
#                 fo.write(ee + "\n")
#                 i = i + 1
#             fo.close()
#
#
#         fo = open(work_dir + 'deBru/all_eDNAs.fa',  "tw")
#         i = 1
#         for ee in eDNAs_all:
#             fo.write(">" + str(i) + "\n")
#             fo.write(ee + "\n")
#             i = i + 1
#         fo.close()
#
#
#         print('Recovering with de Bruijn Graph')
#
#
#         a = time.perf_counter()
#         deG = deBruijnGraph.DeBruijnGraph()
#         deG.kmer_len = 18
#         deG.max_path_num = 2000
#         deG.add_seqs(eDNAs_all)
#         deG.remove_low_cov_kmers(cov_cut_off)
#
#         res['deG']['cov_cut_off'][exp_seq_copy_num].append(cov_cut_off)
#         res['deG']['add_seqs'][exp_seq_copy_num].append(time.perf_counter() - a)
#
#         a = time.perf_counter()
#         res['deG']['Sr'][exp_seq_copy_num].append(test_deG1(dps_seqs, deG, data_block_length))
#         res['deG']['Strand_decoding'][exp_seq_copy_num].append(time.perf_counter() - a)
#
#
#         a = time.perf_counter()
#
#         eDNAs_file_loc = work_dir + 'deBru/all_eDNAs.fa'
#
#         os.system(DBGPS + r' -k ' + str(kmer_length) + r' -a ' + str() + eDNAs_file_loc + r' -o ' + eDNAs_file_loc + r'.count')
#         res['deG']['counting'][exp_seq_copy_num].append(time.perf_counter() - a)
#
#
#         print('Run multiple alignments')
#
#         a = time.perf_counter()
#         for id in dps_seqs:
#             id_seq_file = work_dir + 'muscle/' + str(id)
#             os.system(muscle + r' -in ' + id_seq_file + r' -out ' + id_seq_file +  r'.aln')
#
#         res['mus']['muti-align'][exp_seq_copy_num].append(time.perf_counter() - a)
#
#
#         a = time.perf_counter()
#         i = 0
#         for id in dps_seqs:
#             id_seq_file = work_dir + 'muscle/' + str(id)
#             id_aln_file = id_seq_file +  r'.aln'
#             reads = read_fasta(id_aln_file)
#             cons = majority_merge(reads)
#             if dps_seqs[id] in cons:
#                 i = i + 1
#
#         res['mus']['consensus'][exp_seq_copy_num].append(time.perf_counter() - a)
#         res['mus']['Sr'][exp_seq_copy_num].append(i)
#
#
#
# res = res3
# res['deG'] = {}
# res['mus'] = {}
#
#
# res['deG']['add_seqs'] = {}
# res['deG']['Sr'] = {}
# res['deG']['Strand_decoding'] = {}
# res['deG']['counting'] = {}
# res['deG']['add_dump'] = {}
# res['deG']['cov_cut_off'] = {}
#
#
# res['mus']['muti-align'] = {}
# res['mus']['consensus'] = {}
# res['mus']['Sr'] = {}
#
#
# # for exp_seq_copy_num in range(1:101):
# for exp_seq_copy_num in range(3, 105, 2):
# # for exp_seq_copy_num in [10, 15, 20, 30]:
#     # exp_seq_copy_num = 15
#     # res[exp_seq_copy_num] = {}
#
#     res['deG']['add_seqs'][exp_seq_copy_num] = []
#     res['deG']['Sr'][exp_seq_copy_num] = []
#     res['deG']['Strand_decoding'][exp_seq_copy_num] = []
#     res['deG']['counting'][exp_seq_copy_num] = []
#     res['deG']['add_dump'][exp_seq_copy_num] = []
#     res['deG']['cov_cut_off'][exp_seq_copy_num] = []
#
#
#     res['mus']['muti-align'][exp_seq_copy_num] = []
#     res['mus']['consensus'][exp_seq_copy_num] = []
#     res['mus']['Sr'][exp_seq_copy_num] = []
#
#     cov_cut_off = cut_num(exp_seq_copy_num)
#
#     for rep in range(0, reps):
#         print('Generating eDNAs')
#         eDNAs_all = []
#         for id in dps_seqs:
#             eDNAs = dna_hd.copy_randnum([primerF + dps_seqs[id] + primerE ], exp_seq_copy_num, exp_seq_copy_num)
#             eDNAs = dna_hd.add_rand_sub_new(eDNAs, err_rate / 2)
#             eDNAs = dna_hd.add_rand_ins_new(eDNAs, err_rate / 4)
#             eDNAs = dna_hd.add_rand_del_new(eDNAs, err_rate / 4)
#             eDNAs_all.extend(eDNAs)
#             fo = open(work_dir + 'muscle/' + str(id), "tw")
#             i = 1
#             for ee in eDNAs:
#                 fo.write(">" + str(i) + "\n")
#                 fo.write(ee + "\n")
#                 i = i + 1
#             fo.close()
#
#
#         fo = open(work_dir + 'deBru/all_eDNAs.fa',  "tw")
#         i = 1
#         for ee in eDNAs_all:
#             fo.write(">" + str(i) + "\n")
#             fo.write(ee + "\n")
#             i = i + 1
#         fo.close()
#
#
#         print('Recovering with de Bruijn Graph')
#
#
#         # a = time.perf_counter()
#         # deG = deBruijnGraph.DeBruijnGraph()
#         # deG.kmer_len = 18
#         # deG.max_path_num = 100000
#         # deG.add_seqs(eDNAs_all)
#         # deG.remove_low_cov_kmers(cov_cut_off)
#         #
#         # res['deG']['cov_cut_off'][exp_seq_copy_num].append(cov_cut_off)
#         # res['deG']['add_seqs'][exp_seq_copy_num].append(time.perf_counter() - a)
#         #
#         # a = time.perf_counter()
#         # res['deG']['Sr'][exp_seq_copy_num].append(test_deG1(dps_seqs, deG, data_block_length))
#         # res['deG']['Strand_decoding'][exp_seq_copy_num].append(time.perf_counter() - a)
#
#
#         a = time.perf_counter()
#
#         eDNAs_file_loc = work_dir + 'deBru/all_eDNAs.fa'
#
#         os.system(DBGPS + r' -k ' + str(kmer_length) + r' -a ' + str() + eDNAs_file_loc + r' -o ' + eDNAs_file_loc + r'.count')
#         res['deG']['counting'][exp_seq_copy_num].append(time.perf_counter() - a)
#
#
#         print('Run multiple alignments')
#
#         a = time.perf_counter()
#         for id in dps_seqs:
#             id_seq_file = work_dir + 'muscle/' + str(id)
#             os.system(muscle + r' -in ' + id_seq_file + r' -out ' + id_seq_file +  r'.aln')
#
#         res['mus']['muti-align'][exp_seq_copy_num].append(time.perf_counter() - a)
#
#
#         a = time.perf_counter()
#         i = 0
#         for id in dps_seqs:
#             id_seq_file = work_dir + 'muscle/' + str(id)
#             id_aln_file = id_seq_file +  r'.aln'
#             reads = read_fasta(id_aln_file)
#             cons = majority_merge(reads)
#             if dps_seqs[id] in cons:
#                 i = i + 1
#
#         res['mus']['consensus'][exp_seq_copy_num].append(time.perf_counter() - a)
#         res['mus']['Sr'][exp_seq_copy_num].append(i)


# io_time = []
#
# a = time.perf_counter()
# for id in range(0, 10000):
#     id_seq_file = work_dir + '1seq'
#     os.system(muscle + r' -in ' + id_seq_file + r' -out ' + id_seq_file +  r'.aln')
#
# io_time.append(time.perf_counter() - a)
#
#
# a = time.perf_counter()
# for id in range(0, 10000):
#     id_seq_file = work_dir + '1seq'
#     os.system(muscle + r' -in ' + id_seq_file + r' -out ' + id_seq_file +  r'.aln')
#
# io_time.append(time.perf_counter() - a)
#
#
# a = time.perf_counter()
# for id in range(0, 10000):
#     id_seq_file = work_dir + '1seq'
#     os.system(muscle + r' -in ' + id_seq_file + r' -out ' + id_seq_file +  r'.aln')
#
# io_time.append(time.perf_counter() - a)
#
# #Scripts for output the results
# for cc in range(3, 101):
#     sr = 0
#     dec_tms = []
#     for it in range(0,1000):
#         if res1['deG']['Sr'][cc][it] > 0:
#             sr = sr + 1
#             dec_tms.append(  res1['deG']['Strand_decoding'][cc][it] )
#     print(cc, end="\t")
#     print(sr, end = "\t")
#     print(np.average(dec_tms), end="\n")

for cc in range(5, 8):
    sr = 0
    dec_tms = []
    dec_cons = []
    for it in range(0,1000):
        if res1['mus']['Sr'][cc][it] > 0:
            sr = sr + 1
            dec_tms.append(  res1['mus']['muti-align'][cc][it] )
            dec_cons.append(  res1['mus']['consensus'][cc][it])
    print(cc, end="\t")
    print(sr, end = "\t")
    print(np.average(dec_tms), end="\t")
    print(np.average(dec_cons), end="\n")