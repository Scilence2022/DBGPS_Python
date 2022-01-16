import sys; print('Python %s on %s' % (sys.version, sys.platform))

from utils import *
from test_utils import *
from deBruijnGraph import DeBruijnGraph
from glass import Glass
from DNAdroplet import DNADroplet
import time
import os.path
import getopt



input_file = ''
file_type = 'dump_kmers'
output_file = 'output.rar'
kmer_size = 31
kmer_cut_off = 0
chunk_size = 35
chunk_num = 194818
f_seed = 2
index_bytes = 4
crc_bytes = 2
double_index = False
both_way = False
max_path_num = 1050
index_l = 1
index_u = (chunk_num * 2) + index_l

opts,args = getopt.getopt(sys.argv[1:],'-h-i:-o:t:-k:-c:-n:-s:-a-b',
                          ['help','input=','output=', 'file_type=', 'kmer_size=', 'chunk_size=', 'chunk_num=', 'seed=', 'anchor', 'both_way', 'min_index=',  'max_index=', 'index_bytes=', 'ec_bytes='])

usage = 'Usage:\n' + r'      python decode.py -i input_file -t type_of_seqs -o outfile [Options]'
options = 'Options:\n'
options = options + r'      -h, --help                              Show help information' + '\n'
options = options + r'      -i, --input   <input file>              Input file' + '\n'
options = options + r'      -t, --file_type   <file type>           Input file type: FastQ, Fasta or Jellyfish dumped k-mers (default)' + '\n'
options = options + r'      -o, --output  <output file>             Output file' + '\n'
options = options + r'      -k, --kmer_size  <number>               k-mer size, default = 21 ' + '\n'
options = options + r'      -c, --chunk_size  <size>                Chunk size, default = 35 (bytes)' + '\n'
options = options + r'      -n, --chunk_num  <number>               Chunk number, default = 194,818 (for testing only)' + '\n'
options = options + r'      --cut  <number>                         Cut_off for elimination of low coverage k-mers, default = 0 ' + '\n'
options = options + r'      -s, --seed    <seed>                    Fountain random seed, default 1' + '\n'
options = options + r'      -a, --anchor                            Anchor codes, default On' + '\n'
options = options + r'      -b, --both_way                          Both-way search mode, default On' + '\n'
options = options + r'      --min_index  <initial index>            Initial index, default = 1' + '\n'
options = options + r'      --max_index  <max index>                Max index, default = 20000' + '\n'
options = options + r'      --index_bytes  <number>                 Length of index and anchor codes, default = 4 (bytes)' + '\n'
options = options + r'      --ec_bytes  <number>                    Length of ec codes, default = 2 (bytes)' + '\n'


for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print(usage)
        print(options)
        sys.exit()
    if opt_name in ('-i','--input'):
        input_file = opt_value
    if opt_name in ('-o','--output'):
        output_file = opt_value
    if opt_name in ('-t','--file_type'):
        file_type = opt_value
    if opt_name in ('-k','--kmer_size'):
        kmer_size = int(opt_value)
    if opt_name in ('-c','--chunk_size'):
        chunk_size = int(opt_value)
    if opt_name in ('-n', '--chunk_num'):
        chunk_num = int(opt_value)
    if opt_name in ('-s', '--seed'):
        f_seed = int(opt_value)
    if opt_name in ('-a', '--anchor'):
        double_index = True
    if opt_name in ('--min_index'):
        index_l = int(opt_value)

    if opt_name in ('--max_index'):
        index_u = int(opt_value)

    if opt_name in ('--index_bytes'):
        index_bytes = int(opt_value)

    if opt_name in ('--crc_bytes'):
        crc_bytes = int(opt_value)


start = time.perf_counter()


deG = DeBruijnGraph()
deG.kmer_len = kmer_size
deG.veri_kmers = False
deG.max_path_num = max_path_num

if file_type == 'FastQ' or file_type == 'fastq':
    deG.open_fastq(input_file)
else:
    if file_type == 'fasta' or file_type == 'Fasta':
        deG.open_fasta(input_file)
    else:
        deG.open_dump(input_file)

strands_recovered = {}
strand_search_details = {}
strand_search_details['crc_fail'] = 0

print('\nRemoving low coverage k-mers ......')
deG.remove_low_cov_kmers(kmer_cut_off)

print('\nReconstructing DNA strands ......')

a = time.perf_counter()
for index in range(index_l, index_u+1):
    deG.find_droplet_DNA(index, chunk_size)
    if len(deG.crc_paths) > 1:
        strand_search_details['crc_fail']  = strand_search_details['crc_fail']  + 1
    else:
        if len(deG.crc_paths) > 0:
            strands_recovered[index] = deG.crc_paths[0]

a = time.perf_counter() - a
# print('CRC fails:')
# print(strand_search_details['crc_fail'])
print('Recovered: ', end='')
print(len(strands_recovered), end= '')
print(' strands used ', end='')
print(a, end=' seconds\n')


print('Rebuilding DNA droplets........')
degree_table = get_degrees(chunk_num, chunk_num * 1000, f_seed)
cup = Glass(chunk_num)
for id in strands_recovered:
    adp = DNADroplet(bytes(chunk_size))
    adp.num_of_chunks = chunk_num
    adp.set_head_index(id)
    adp.set_droplet_from_DNA_CRC(strands_recovered[id])
    adp.degree = degree_table[id-1]
    adp.update()
    cup.addDroplet(adp)




print('Decoding by fountain codes .........')

cup.decode()
cup.writeToFile(output_file)

end = time.perf_counter() - start
print('Decoding time: ', end ='')
print(end)


