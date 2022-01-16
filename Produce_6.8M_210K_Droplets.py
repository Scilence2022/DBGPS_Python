import random

from utils import *
from DNAdroplet import DNADroplet
#from fountain import Fountain
from DNAfountain import DNAFountain
from glass import Glass
from deBruijnGraph import DeBruijnGraph
from crc16pure import *
from test_utils import *
import copy


kmer_length = 21
data_block_length = 35

fountain_seed = 2

work_dir = r'input_files'  + '/'
file1 = open(work_dir + r'Dunhuang.6.8MB.zip', 'rb')

fountain_init_index = 101010101

filebytes1 = file1.read()
file1.close()
failed_gts = []
suc_num = 0

fdna1 = DNAFountain(filebytes1, data_block_length, fountain_init_index, fountain_seed, 16, 16, 8)
fdna1.degree_table_folds = 1000
fdna1.ROBUST_FAILURE_PROBABILITY = 0.01
fdna1.c_value = 0.01
fdna1.gen_degrees()

total_size = 210000
droplet_all = get_droplets_check_repeat_kmer(total_size, fdna1, kmer_length)


p1 = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
p2 = 'CTGACACTGATGCATCCG'  # complement seq of P2

file2 = open(work_dir + r'6.8MB.DNAs.1','tw')
# file2.write('Head Index\tData\tDNA\tDNA-Primers\tDegree\tChunk Nums\n')

i = 0
for dps in droplet_all:
    file2.write(str(dps.head_index))
    # file2.write('\t')
    # file2.write(str(dps.data))
    file2.write('\t')
    file2.write(dps.to_DNA_CRC_sIndex())
    # file2.write('\t')
    # file2.write(p1)
    # file2.write(dps.to_DNA_CRC_sIndex())
    # file2.write(p2)
    # file2.write('\t')
    # file2.write(str(dps.degree))
    # file2.write('\t')
    # file2.write(str(dps.get_chunk_nums()))
    # file2.write('\t')
    # file2.write(str(dps.tail_index))
    file2.write('\n')
    i = i + 1
file2.close()

print("Encode finished!")











