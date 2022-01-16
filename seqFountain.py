import re
import numpy as np
import random
from test_utils import sim_file_seqs
from test_utils import read_fasta

class SeqFountain:

    def __init__(self):
        self.seqs = []
        self.max_seq_num = 1000000000
        self.total_bases = 0

    def av_len(self):
        return int(self.total_length()/len(self.seqs))

    def total_length(self):
        tl = 0
        for a in self.seqs:
            tl = tl + len(a)
        return tl

    def total_bases_calc(self):
        for s in self.seqs:
            self.total_bases += len(s)
        return self.total_bases

    def read_FQ(self, file):

        f = open(file, "r")
        # matchLineA = re.compile('^(@)')
        matchLineB = re.compile('^(\+)')
        # matchLineC = re.compile('(F)')

        last_line = f.readline()
        line = f.readline()
        reads = 0
        while (line.strip()):
            if (matchLineB.match(line)):
                self.add_seq(last_line.strip())
                reads = reads + 1
                if reads > 100000:
                    print('.', end='')
                    reads = 0
            else:
                last_line = line
            line = f.readline()

    def read_sim(self, file):
        self.seqs.extend(sim_file_seqs(file))

    def read_fasta(self, file):
        self.seqs.extend(read_fasta(file))

    def add_seq(self, seq):
        self.seqs.append(seq)

    def rd_seq(self):
        return self.seqs[random.randint(0, len(self.seqs)-1)]

    def rd_seq_num(self, num):
        rd_nums = random.sample(range(0, len(self.seqs)-1), num)
        seqs = []
        for num in rd_nums:
            seqs.append(self.seqs[num])
        return seqs

    def rd_seq_base(self, bases):
        num = int(bases/self.av_len()) + 1
        # rd_nums = random.sample(range(0, len(self.seqs)-1), num)
        rd_nums = np.random.choice(range(0, len(self.seqs) - 1), num, False)
        seqs = []
        for num in rd_nums:
            seqs.append(self.seqs[num])
        return seqs