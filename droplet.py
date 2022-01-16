from utils import randChunkNums, generate_chunk_nums, get_degrees

import json
import random

class Droplet:
    def __init__(self, data=b'', symbol_index=-1, num_of_chunks=100, seed=1):
        self.data = data
        self.symbol_index = symbol_index
        self.num_of_chunks = num_of_chunks
        self.degree = 0
        self.chunk_num_list = self.get_chunk_nums()

    def set_symbol_index(self, symbol_index):
        self.symbol_index = symbol_index
        #self.num_of_chunks = num_of_chunks
        #self.degree = self.get_degree()
        self.chunk_num_list = self.get_chunk_nums()

    def get_chunk_nums(self):
        #random.seed(self.seed)
        if self.symbol_index > 0:
            self.chunk_num_list = generate_chunk_nums(self.num_of_chunks, [self.degree], self.symbol_index)
        else:
            self.chunk_num_list = None
        return self.chunk_num_list

    def get_degree(self):
#        if self.symbol_index > 0:
#            self.degree = get_degrees(self.num_of_chunks, 1, self.symbol_index)[0]
#        else:
#            self.degree = None
        return self.degree

    def to_string(self):
        return json.dumps(
            {
                'seed': self.symbol_index,
                'num_chunks': self.num_of_chunks,
                'data': self.data
            })

