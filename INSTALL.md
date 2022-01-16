

### This is a pure Python implementation which only requires Numpy package. 

### To install Numpy, run the following command:

	pip install numpy

After the installation of Numpy, the scripts are ready to run.

Optional. Cythonize the code by running the following command in this folder:
python setup.py build_ext --inplace


### To encode a digital file into DNA sequences:

	python encode.py -i a_digital_file_name -n number_of_droplets -o encode.fasta

### For example: 
	python encode.py -i input_files/314kb.rar -n 12000 -o encode.fasta

    Generating droplets............
    The DNA sequences in fasta format: encode.fasta
    Details about the encoding: encode.fasta.log

### This will generate a Fasta file, i.e. encode.fasta, containing the DNA sequences with data encoded.
Additionally, a log file (encode.fasta.log) containing details of the encoding process was also generated.

### For additional parameters and options to run encode.py,
    python encode.py -h

Usage:

      python encode.py -i input_file -n number_of_droplets -o output.fasta [Options]
Options:

      -h, --help                              Show help information
      -i, --input   <input file>              Input file
      -o, --output  <output file>             Output file
      -n, --droplet_num  <number>             Number of droplets, default 12,000 
      -c, --chunk_size  <size>                Chunk size, default 30 (bytes)
      -s, --seed    <seed>                    Fountain random seed, default 1
      -d, --double_index                      Double index mode, default On
      -l, --initial_index  <initial index>    Initial index, default 1
      --index_bytes  <number>                 Length of index and anchor codes, default 4 (bytes)
      --ec_bytes  <number>                    Length of ec codes, default 2 (bytes)


### To decode the DNA sequences into digital file:

    python decode.py -i <a DNA sequences file in Fasta, FastQ or Jellyfish dumped format> -t fasta -o <decoded file name>

### For example:

	python decode.py -i encode.fasta -t fasta -o decoded.rar

	Opening fasta file
	.......................
	Removing low coverage k-mers ......
	Reconstructing DNA strands ......
	Recovered: 12000 strands used 45.3998139 seconds
	Rebuilding DNA droplets........
	Decoding by fountain codes .........

### For additional parameters and options to run decode.py,
	python decode.py -h

Usage:

      python decode.py -i input_file -t type_of_seqs -o outfile [Options]
	  
Options:

      -h, --help                              Show help information
      -i, --input   <input file>              Input file
      -t, --file_type   <file type>           Input file type: FastQ, Fasta or Jellyfish dumped k-mers (default)
      -o, --output  <output file>             Output file
      -k, --kmer_size  <number>               k-mer size, default = 21 
      -c, --chunk_size  <size>                Chunk size, default = 30 (bytes)
      -n, --chunk_num  <number>               Chunk number, default = 10,740 (for testing only)
      --cut  <number>                         Cut_off for elimination of low coverage k-mers, default = 0 
      -s, --seed    <seed>                    Fountain random seed, default 1
      -a, --anchor                            Anchor codes, default On
      -b, --both_way                          Both-way search mode, default On
      --min_index  <initial index>            Initial index, default = 1
      --max_index  <max index>                Max index, default = 20000
      --index_bytes  <number>                 Length of index and anchor codes, default = 4 (bytes)
      --ec_bytes  <number>                    Length of ec codes, default = 2 (bytes)


