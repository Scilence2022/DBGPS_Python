import sys; print('Python %s on %s' % (sys.version, sys.platform))

import random
from utils import *
from DNAfountain import DNAFountain
from test_utils import *

import copy
import time
import os.path
import sys
import getopt


chunk_size = 35
fountain_seed = 1
droplet_num = 210000
double_index = False
start_index = 1
work_dir = r''
input_file = ''
output_file = ''
index_bytes = 4
anchor_bytes = 4
ec_bytes = 2
redundancy_rate = 0.05


p1 = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
p2 = 'CTGACACTGATGCATCCG'  # complement seq of P2


opts,args = getopt.getopt(sys.argv[1:],'-h-i:-o:-c:-n:-s:-r:-d-l:',
                          ['help','input=','output=',  'chunk_size=', 'redundancy_rate=', 'droplet_num=', 'seed=', 'initial_index=', 'index_bytes=', 'ec_bytes='])

usage = 'Usage:\n' + r'      python encode.py -i input_file -n number_of_droplets -o output.fasta [Options]'
options = 'Options:\n'
options = options + r'      -h, --help                              Show help information' + '\n'
options = options + r'      -i, --input   <input file>              Input file' + '\n'
options = options + r'      -o, --output  <output file>             Output file' + '\n'
options = options + r'      -r, --redundancy_rate  <number>         Redundancy rate, default 0.05 ' + '\n'
options = options + r'      -n, --droplet_num  <number>             Number of droplets, default 210,000 ' + '\n'
options = options + r'      -c, --chunk_size  <size>                Chunk size, default 35 (bytes)' + '\n'
options = options + r'      -s, --seed    <seed>                    Fountain random seed, default 1' + '\n'
options = options + r'      -l, --initial_index  <initial index>    Initial index, default 1' + '\n'
options = options + r'      --index_bytes  <number>                 Length of index and anchor codes, default 4 (bytes)' + '\n'
options = options + r'      --ec_bytes  <number>                    Length of ec codes, default 2 (bytes)' + '\n'


for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print(usage)
        print(options)
        sys.exit()
    if opt_name in ('-i','--input'):
        input_file = opt_value
    if opt_name in ('-o','--output'):
        output_file = opt_value
    if opt_name in ('-c','--chunk_size'):
        chunk_size = int(opt_value)
    if opt_name in ('-n', '--droplet_num'):
        droplet_num = int(opt_value)

    if opt_name in ('-r', '--redundancy_rate'):
        redundancy_rate = float(opt_value)

    if opt_name in ('-s', '--seed'):
        fountain_seed = int(opt_value)

    if opt_name in ('-l', '--low_index'):
        start_index = int(opt_value)
    if opt_name in ('--index_bytes', '--notmatch'):
        index_bytes = int(opt_value)
        anchor_bytes = index_bytes
    if opt_name in ('--ec_bytes', '--notmatch'):
        ec_bytes = int(opt_value)

start = time.perf_counter()

file1 = open(work_dir + input_file, 'rb')
filebytes1 = file1.read()
file1.close()


fdna1 = DNAFountain(filebytes1, chunk_size, start_index, fountain_seed, index_bytes*4, anchor_bytes*4, ec_bytes*4)
fdna1.degree_table_folds = 5
fdna1.gen_degrees()

print('Chunk number: ')
print(fdna1.num_of_chunks)

if redundancy_rate > 0:
    droplet_num = int(fdna1.num_of_chunks * (1 + redundancy_rate)) + 1

print('Generating droplets' + str(droplet_num) + '............')

droplet_all = get_droplets(droplet_num, fdna1)

print('Max index: ')
print(droplet_all[-1].head_index)

file2 = open(output_file + '.log','tw')
# file2.write('Head Index\tData\tDNA\tDegree\tChunk Nums\n')

for dps in droplet_all:
    file2.write(str(dps.head_index))
    file2.write('\t')
    file2.write(str(dps.data))
    file2.write('\t')
    if double_index:
        file2.write(dps.to_DNA_CRC())
    else:
        file2.write(dps.to_DNA_CRC_sIndex())
    file2.write('\t')
    file2.write(str(dps.degree))
    file2.write('\t')
    file2.write(str(dps.get_chunk_nums()))
    file2.write('\t')
    file2.write(str(dps.tail_index))
    file2.write('\n')
file2.close()

file2 = open(output_file,'tw')
for dps in droplet_all:
    file2.write(">")
    file2.write(str(dps.head_index))
    file2.write("\n")
    file2.write(p1)
    if double_index:
        file2.write(dps.to_DNA_CRC())
    else:
        file2.write(dps.to_DNA_CRC_sIndex())
    file2.write(p2)
    file2.write('\n')
file2.close()


print('The DNA sequences in fasta format: ', end='')
print(output_file)

print('Details about the encoding: ', end='')
print(output_file + '.log')

end = time.perf_counter() - start
print('Encoding time: ', end ='')
print(end)
