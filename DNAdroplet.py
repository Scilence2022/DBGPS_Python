import crc16pure
from utils import *
from droplet import Droplet
# from crc16pure import crc16xmodem


import json
import random

class DNADroplet(Droplet):
    def __init__(self, data=b'', index=-1, num_of_chunks=100, head_index_len=16, tail_index_len=16, seed=1, crc_len=8):
        super(DNADroplet, self).__init__(data, index, num_of_chunks)
        self.seed = seed
        self.head_index_len = head_index_len
        self.tail_index_len = tail_index_len
        
        self.crc_len = crc_len

        self.head_index = index
        self.tail_index = self.get_tail_index()
        self.anchor = False

        self.checksum = "CRC"
        self.crc = 0
        self.allbytes = b''
        
        if data:
            self.data_len = len(data) * 4
        else:
            self.data_len = len(self.data) * 4


    #2020.02.19 Deal with the index-->chunks problem
    def set_head_index(self, index_num):
        self.head_index = index_num
        self.symbol_index = index_num
        self.update()

    # 2020.02.19 Deal with the index-->chunks problem
    def update(self):
        #self.get_tail_index()
        self.get_chunk_nums()

    def get_head_index_dna(self):
        aDNA = self.to_DNA_CRC()
        return aDNA[0:self.head_index_len]

    def get_tail_index_dna(self):
        aDNA = self.to_DNA_CRC()
        return aDNA[len(aDNA) - self.tail_index_len:len(aDNA)]

    def get_head_index(self):
        return self.head_index

    def get_tail_index(self):
        random.seed(self.head_index)
        self.tail_index = random.randint(0, 4 ** (int(self.tail_index_len/4)*4) - 1)
        return self.tail_index

    def _crc(self):

        databytes = self.head_index.to_bytes(int(self.head_index_len/4), 'big',
                                            signed=False) + self.data
        self.crc = crc16pure.crc16xmodem(databytes)

    def get_crc(self):
        if self.crc <= 0:
            self._crc()
        return self.crc

    def to_DNA_CRC(self):
        if self.anchor:
            return self.to_DNA_CRC_dIndex()
        else:
            return self.to_DNA_CRC_sIndex()

    def to_DNA_CRC_dIndex(self):
        self.get_crc()
        allbytes = self.head_index.to_bytes(int(self.head_index_len/4), 'big', signed=False)\
                   + self.data \
                   + self.crc.to_bytes(int(self.crc_len/4), 'big', signed=False) \
                   + self.tail_index.to_bytes(int(self.tail_index_len/4), 'big', signed=False)
        self.allbytes = allbytes
        return bytesToDNA(allbytes)

    def to_DNA_CRC_sIndex(self):
        self.get_crc()
        allbytes = self.head_index.to_bytes(int(self.head_index_len/4), 'big', signed=False)\
                   + self.data \
                   + self.crc.to_bytes(int(self.crc_len/4), 'big', signed=False)
        self.allbytes = allbytes
        return bytesToDNA(allbytes)

    def set_droplet_from_DNA_CRC(self, dnastr):

        data_bytes = DNAToBytes(dnastr[-(self.data_len + self.crc_len):-self.crc_len])
        self.data = data_bytes
        self.update()
        return True

    def set_droplet_from_DNA_CRC_old(self, dnastr):
  
        data_bytes = DNAToBytes(dnastr[self.head_index_len:self.head_index_len + self.data_len])
        self.data = data_bytes
        self.update()
        return True


