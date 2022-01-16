import random
import copy
import time

from utils import *
from test_utils import *
from DNAdroplet import DNADroplet
from DNAfountain import DNAFountain
from glass import Glass
from deBruijnGraph import DeBruijnGraph
from DNAHandler import DNAHandler


check_index = True
number_of_droplets = 1
droplet_rep_num = 1000
rep_num = 3
kmer_length = 10
data_block_length = 35
max_seq_copies = 25

err_rates = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]

primerF = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
primerE = 'CTGACACTGATGCATCCG'  # complement seq of P2

file1 = open('input_files/Dunhuang.6.8MB.zip', 'rb')


filebytes1 = file1.read()
file1.close()
fountain_init_index = 1

# file1 = open('2.pdf', 'rb')
fountain_seed = 1
fdna1 = DNAFountain(filebytes1, data_block_length, fountain_init_index, fountain_seed)
fdna1.gen_degrees()

droplet_all = get_droplets(number_of_droplets, fdna1)

droplet_IDs = []
droplet_ID_DNA = {}
droplet_DNAs = []
maxID = 1

for dps in droplet_all:
    droplet_ID_DNA[dps.head_index] = dps.to_DNA_CRC()
    droplet_IDs.append(dps.head_index)
    droplet_DNAs.append(primerF + dps.to_DNA_CRC() + primerE)
    if dps.head_index > maxID:
        maxID = dps.head_index
    # deG.addSeq(dps.to_DNA_CRC())

dna_hd = DNAHandler()

# -----------------------------------

res = {}
res['sub'] = {}
res['ins'] = {}
res['del'] = {}
res['brk'] = {}
res['reg'] = {}
res['mix'] = {}

succ_num_sub = res['sub']
succ_num_ins = res['ins']
succ_num_del = res['del']
succ_num_brk = res['brk']
succ_num_reg = res['reg']
succ_num_mix = res['mix']

for seqnum in range(2, max_seq_copies + 1):
    succ_num_sub[seqnum] = {}
    succ_num_ins[seqnum] = {}
    succ_num_del[seqnum] = {}
    succ_num_brk[seqnum] = {}
    succ_num_reg[seqnum] = {}
    succ_num_mix[seqnum] = {}
    print(seqnum)

    #Testing mix errors
    succ_res = succ_num_mix
    for err_rate in err_rates:
        # for errrate in [0.01, 0.02]:
        succ_res[seqnum][err_rate] = 0
        suc_nums = []
        for rep in range(0, rep_num):
            suc_num = 0
            for test in range(0, droplet_rep_num):
                for id in droplet_IDs:
                    eDNAs = dna_hd.copy_randnum([droplet_ID_DNA[id]], seqnum, seqnum)

                    eDNAs = dna_hd.add_rand_sub_new(eDNAs, err_rate/2, True)
                    eDNAs = dna_hd.add_rand_ins_new(eDNAs, err_rate/8, True)
                    eDNAs = dna_hd.add_rand_del_new(eDNAs, err_rate/8, True)
                    eDNAs = dna_hd.add_rand_brk_new(eDNAs, err_rate/4, True)
                    eDNAs = dna_hd.random_ligation(eDNAs)
                    if strand_theory_recover_rate(droplet_all, eDNAs, kmer_length, 0) > 0:
                        suc_num = suc_num + 1
            suc_nums.append(suc_num)

        succ_res[seqnum][err_rate] = {}
        succ_res[seqnum][err_rate]['Sr'] = np.average(suc_nums)
        succ_res[seqnum][err_rate]['std'] = np.std(suc_nums, ddof=1)

    # Testing substitution errors
    succ_res = succ_num_sub
    for err_rate in err_rates:
        # for errrate in [0.01, 0.02]:
        succ_res[seqnum][err_rate] = 0
        suc_nums = []
        for rep in range(0, rep_num):
            suc_num = 0
            for test in range(0, droplet_rep_num):
                for id in droplet_IDs:
                    eDNAs = dna_hd.copy_randnum([droplet_ID_DNA[id]], seqnum, seqnum)

                    eDNAs = dna_hd.add_rand_sub_new(eDNAs, err_rate, True)
                    if strand_theory_recover_rate(droplet_all, eDNAs, kmer_length, 0) > 0:
                        suc_num = suc_num + 1
            suc_nums.append(suc_num)

        succ_res[seqnum][err_rate] = {}
        succ_res[seqnum][err_rate]['Sr'] = np.average(suc_nums)
        succ_res[seqnum][err_rate]['std'] = np.std(suc_nums, ddof=1)

    # Testing insertion errors
    succ_res = succ_num_ins
    for err_rate in err_rates:
        # for errrate in [0.01, 0.02]:
        succ_res[seqnum][err_rate] = 0
        suc_nums = []
        for rep in range(0, rep_num):
            suc_num = 0
            for test in range(0, droplet_rep_num):
                for id in droplet_IDs:
                    eDNAs = dna_hd.copy_randnum([droplet_ID_DNA[id]], seqnum, seqnum)

                    eDNAs = dna_hd.add_rand_ins_new(eDNAs, err_rate, True)
                    if strand_theory_recover_rate(droplet_all, eDNAs, kmer_length, 0) > 0:
                        suc_num = suc_num + 1
            suc_nums.append(suc_num)

        succ_res[seqnum][err_rate] = {}
        succ_res[seqnum][err_rate]['Sr'] = np.average(suc_nums)
        succ_res[seqnum][err_rate]['std'] = np.std(suc_nums, ddof=1)

    # Testing deletion errors
    succ_res = succ_num_del
    for err_rate in err_rates:
        succ_res[seqnum][err_rate] = 0
        suc_nums = []
        for rep in range(0, rep_num):
            suc_num = 0
            for test in range(0, droplet_rep_num):
                for id in droplet_IDs:
                    eDNAs = dna_hd.copy_randnum([droplet_ID_DNA[id]], seqnum, seqnum)

                    eDNAs = dna_hd.add_rand_del_new(eDNAs, err_rate, True)
                    if strand_theory_recover_rate(droplet_all, eDNAs, kmer_length, 0) > 0:
                        suc_num = suc_num + 1
            suc_nums.append(suc_num)

        succ_res[seqnum][err_rate] = {}
        succ_res[seqnum][err_rate]['Sr'] = np.average(suc_nums)
        succ_res[seqnum][err_rate]['std'] = np.std(suc_nums, ddof=1)

    # Testing DNA break and rearrangement errors
    succ_res = succ_num_brk
    for err_rate in err_rates:
        # for errrate in [0.01, 0.02]:
        succ_res[seqnum][err_rate] = 0
        suc_nums = []
        suc_reg_nums = []
        for rep in range(0, rep_num):
            suc_num = 0
            reg_suc_num = 0
            for test in range(0, droplet_rep_num):
                for id in droplet_IDs:
                    eDNAs = dna_hd.copy_randnum([droplet_ID_DNA[id]], seqnum, seqnum)

                    eDNAs = dna_hd.add_rand_brk_new(eDNAs, err_rate, True)
                    if strand_theory_recover_rate(droplet_all, eDNAs, kmer_length, 0) > 0:
                        suc_num = suc_num + 1
                    eDNAs = dna_hd.random_ligation(eDNAs)
                    if strand_theory_recover_rate(droplet_all, eDNAs, kmer_length, 0) > 0:
                        reg_suc_num = reg_suc_num + 1

            suc_nums.append(suc_num)
            suc_reg_nums.append(reg_suc_num)


        succ_res[seqnum][err_rate] = {}
        succ_res[seqnum][err_rate]['Sr'] = np.average(suc_nums)
        succ_res[seqnum][err_rate]['std'] = np.std(suc_nums, ddof=1)

        succ_num_reg[seqnum][err_rate] = {}
        succ_num_reg[seqnum][err_rate]['Sr'] = np.average(suc_reg_nums)
        succ_num_reg[seqnum][err_rate]['std'] = np.std(suc_reg_nums, ddof=1)





