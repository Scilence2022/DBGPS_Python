from droplet import Droplet
from math import ceil
from utils import *

import random

class Fountain:
    def __init__(self, data=b'', chunk_size=35, index=1, seed=1):
        self.seed = seed
        self.index = index
        self.data = data
        self.chunk_size = chunk_size
        self.num_of_chunks = int(ceil(len(data) / float(chunk_size)))
        #self.ROBUST_FAILURE_PROBABILITY = 0.1   #Large-Scale simulations
        #self.c_value = 0.2  #Large-Scale simulations
        self.ROBUST_FAILURE_PROBABILITY = 0.1 # 0.01
        self.c_value = 0.2 #0.01

        self.degree_table_folds = 3
        #self.degree_table_folds = 1000# For 6M data
        #self.degrees = get_degrees(self.num_of_chunks, int(self.num_of_chunks * self.degree_table_folds), self.seed, self.ROBUST_FAILURE_PROBABILITY, self.c_value, 'robust')

        self.VERBOSE = False

       # self.random_degrees = get_degrees_from("robust", blocks_n, k=drops_quantity)
        #random.seed(seed)
    #replaced by DNAdroplet in DNAFountain
    # def droplet(self):
    #     #self.updateSeed()
    #     aDroplet = Droplet(b'', self.seed, self.num_of_chunks)
    #     #return aDroplet
    #     chunk_num_list = aDroplet.chunk_num_list
    #     data = None
    #     for num in chunk_num_list:
    #         if data is None:
    #             data = self.chunk(num)
    #         else:
    #             data = xor(data, self.chunk(num))
    #     aDroplet.data = data
    #     return aDroplet

    def chunk(self, num):
        start = self.chunk_size * num
        end = min(self.chunk_size * (num+1), len(self.data))
        return self.data[start:end]

    def updateIndex(self):
        self.index = self.index + 1

    def updateSeed(self):
        self.seed = self.seed + 1
       # self.seed = random.randint(0,2**31-1)
        #random.seed(self.seed)

    def gen_degrees(self):
        self.degrees = get_degrees(self.num_of_chunks, int(self.num_of_chunks * self.degree_table_folds), self.seed, self.ROBUST_FAILURE_PROBABILITY, self.c_value, 'robust')
 #   def openfile(self, str):
