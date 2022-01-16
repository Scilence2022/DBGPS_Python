import sys; print('Python %s on %s' % (sys.version, sys.platform))
from utils import *
from test_utils import *
from deBruijnGraph import DeBruijnGraph
from glass import Glass
from DNAdroplet import DNADroplet
import time
import os

kmer_len = 18
kmer_cut_off = 3

work_dir = r'/work/Grass_work_test/'
seqs = read_fasta(work_dir + "File1_ODNA.fasta")
jellyfish = r'jellyfish '

deG = DeBruijnGraph()
deG.kmer_len = kmer_len
input_file = work_dir + 'I16_S2_R1_001.fastq'
output_file = work_dir + 'I16_S2_R1_001.fastq.km' + str(kmer_len)
dump_file = output_file + r'.dump.L' + str(kmer_cut_off)

print(jellyfish + ' count -s 1000000 -m ' + str(kmer_len) + r' ' + input_file + r' -o ' + output_file)
os.system(jellyfish + ' count -s 1000000 -m ' + str(kmer_len) + r' ' + input_file + r' -o ' + output_file)
print(jellyfish + ' dump -c -t ' + output_file + r' -o ' + dump_file + r' -L ' + str(kmer_cut_off) )
os.system(jellyfish + ' dump -c -t ' + output_file + r' -o ' + dump_file + r' -L ' + str(kmer_cut_off) )

deG.open_dump(dump_file)


res = {}
for block_len in range(3, 16):
    find_path = 0
    for seq in seqs:
        if deG.test_find_path_grass(seq, block_len):
            find_path = find_path + 1
    res[block_len] = find_path


