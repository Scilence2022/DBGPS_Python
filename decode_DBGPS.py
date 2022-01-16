import sys; print('Python %s on %s' % (sys.version, sys.platform))

from utils import *
from test_utils import *
from deBruijnGraph import DeBruijnGraph
from glass import Glass
from DNAdroplet import DNADroplet
import time
import os.path
import getopt


input_file = r'input_files/6.8MB.DNAs.1'
output_file = input_file + '.dec.new.zip'

chunk_size = 35
chunk_num = 194818
f_seed = 2
index_bytes = 4
crc_bytes = 2
total_bytes = 6818623


opts,args = getopt.getopt(sys.argv[1:],'-h-i:-o:-d:-n:',
                          ['help','input=','output=',  'chunk_size=', 'chunk_num=', 'seed=',  'index_bytes=', 'ec_bytes='])

usage = 'Usage:\n' + r'      python decode_DBGPS.py -i input_file -o outfile [Options]'
options = 'Options:\n'
options = options + r'      -h, --help                              Show help information' + '\n'
options = options + r'      -i, --input   <input file>              The decoded strands by DBGPS' + '\n'
options = options + r'      -o, --output  <output file>             Output file' + '\n'

options = options + r'      -d, --chunk_size  <size>                Chunk size, default = 35 (bytes)' + '\n'
options = options + r'      -n, --chunk_num  <number>               Chunk number, default = 194,818 (for testing only)' + '\n'
options = options + r'      --seed    <seed>                        Fountain random seed, default 1' + '\n'

options = options + r'      --index_bytes  <number>                 Bytes of index codes, default = 4' + '\n'
options = options + r'      --ec_bytes  <number>                    Bytes of ec codes, default = 2' + '\n'


for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print(usage)
        print(options)
        sys.exit()
    if opt_name in ('-i','--input'):
        input_file = opt_value
    if opt_name in ('-o','--output'):
        output_file = opt_value
    if opt_name in ('-d','--chunk_size'):
        chunk_size = int(opt_value)
    if opt_name in ('-n', '--chunk_num'):
        chunk_num = int(opt_value)
    if opt_name in ('-s', '--seed'):
        f_seed = int(opt_value)
    # if opt_name in ('--index_bytes'):
    #     index_bytes = int(opt_value)
    if opt_name in ('--crc_bytes'):
        crc_bytes = int(opt_value)

start = time.perf_counter()

strands_recovered = read_sim(input_file)

print('Rebuilding degree table ........')
degree_table = get_degrees(chunk_num, chunk_num*1000, f_seed, 0.01, 0.01)

print('Rebuilding DNA droplets ........')
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

zero_bytes = chunk_num * chunk_size - total_bytes
if cup.isDone():
    #Fix the file terminal issue
    cup.chunks[-1] = cup.chunks[-1][:-zero_bytes]
    cup.writeToFile(output_file)
    print('Successfully decoded by the Fountain codes. ')
else:
    print('Decoding failed')

end = time.perf_counter() - start
print('Decoding time: ', end ='')
print(end)


