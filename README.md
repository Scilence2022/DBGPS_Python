

## De Bruijn Graph-based de novo strand assembly for reliable and fast data retrieval in DNA data storage



The Fountain codes is modified from [https://github.com/dbieber/fountaincode](https://github.com/dbieber/fountaincode)




#### Source Files
[deBruijnGraph.py](deBruijnGraph.py) defines the DeBruijnGraph object, use de Bruijn path searching for decoding of DNA strands.

[fountain.py](fountain.py) defines the Fountain, produces Droplets according to the Fountain Code implementation.

[DNAfountain.py](DNAfountain.py) defines the DNA Fountain which is inherited from [fountain.py]

[droplet.py](droplet.py) defines the Droplet, containing a seed and the XOR'd data.

[DNAdroplet.py](DNAdroplet.py) defines the DNA Droplet which is inherited from [droplet.py]

[crc16pure.py](crc16pure.py) contains the function implements the CRC16 codes.

[glass.py](glass.py) defines a Glass, used to collect DNADroplets and/or Droplets and reconstruct the original data. The decoding functions have been modified to enhance decoding performance especially with large files. 

[DNAHandler.py](DNAHandler.py) defines the DNAHandler, used for DNA sequences manupulations, e.g. replication and introduction of errors in specified rates.

[utils.py](utils.py) contains helper functions used to make the (DNA)Fountain and (DNA)Droplets compatible.

[test_utils.py](test_utils.py) contains helper functions for performance testing.

[seqFountain.py](seqFountain.py) defines the seqFountain object which can open a sequencing file and generate specific number of random reads.

