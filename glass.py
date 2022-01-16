import sys
from utils import xor

sys.setrecursionlimit(100000)


class Glass:
    def __init__(self, num_chunks):
        self.entries = []
        self.degree_one_entries = []
        self.droplets = []
        self.num_chunks = num_chunks
        self.chunks = [None] * num_chunks
        self.chunk_entries = [None] * num_chunks
        self.decoded_chunk_num = 0

    def addDroplet(self, d):
        self.droplets.append(d)
        ck_nums = d.get_chunk_nums()
        ck_nums_hs = {}
        for nm in ck_nums:
            ck_nums_hs[nm] = nm

        entry = [ck_nums_hs, d.data]
        if len(ck_nums) > 1:
            self.entries.append(entry)
        else:
            self.degree_one_entries.append(entry)

    def chunk_entry_rel(self):
        print(r'caculating chunk entry rels')
        i = 0
        while i < len(self.chunk_entries):
            self.chunk_entries[i] = {}
            i = i + 1
        for i, entry in enumerate(self.entries):
            # print(r'.')
            for chunk_num in entry[0].keys():
                self.chunk_entries[chunk_num][i] = 1

    def decode(self):
        self.chunk_entry_rel()
        entry = self.degree_one_entries.pop()
        while entry:
            chunk_num = list(entry[0].keys())[0]
            data = entry[1]
            self.detach_chunk(chunk_num, data)
            self.decoded_chunk_num = self.decoded_chunk_num + 1
            if len(self.degree_one_entries) > 0:
                entry = self.degree_one_entries.pop()
            else:
                entry = None

    def detach_chunk(self, chunk_num, data):
        self.chunks[chunk_num] = data
        for entry_id in self.chunk_entries[chunk_num]:

            if len(self.entries[entry_id][0]) < 2: # degree 1
                self.degree_one_entries.append(self.entries[entry_id])
            else:
                self.entries[entry_id][1] = xor(self.entries[entry_id][1], data)
                self.entries[entry_id][0].pop(chunk_num)
            if len(self.entries[entry_id][0]) < 2: # degree 1
                self.degree_one_entries.append(self.entries[entry_id])
        self.chunk_entries[chunk_num] = []

    def decode1(self):
        if len(self.droplets) > 10:
            solved_blocks_count = 0
            iteration_solved_count = 0
            while iteration_solved_count > 0 or solved_blocks_count == 0:
                iteration_solved_count = 0
                for i, entry in enumerate(self.entries):
                    if len(entry[0]) == 1:
                        iteration_solved_count += 1
                        entry_index = next(iter(entry[0]))
                        # if self.chunks[entry_index] is not None:
                        #     assert self.chunks[entry_index] == entry[1], "Wrong droplet detected!"
                        self.chunks[entry_index] = entry[1]

                        self.entries.pop(i)

                        solved_blocks_count += 1
                        self.detach_entry(entry)
        else:
            print('Too few droplets!')

    def detach_entry1(self, entry):
        for other_entry in self.entries:
            if len(other_entry[0]) > 1 and next(iter(entry[0])) in other_entry[0]:
                other_entry[1] = xor(other_entry[1], entry[1])
                other_entry[0].remove(next(iter(entry[0])))

    def updateEntry1(self, entry):
        for chunk_num in entry[0]:
            if self.chunks[chunk_num] is not None:
                entry[1] = xor(entry[1], self.chunks[chunk_num])
                entry[0].remove(chunk_num)
        if len(entry[0]) == 1:
            self.chunks[entry[0][0]] = entry[1]
            self.entries.remove(entry)
            for d in self.entries:
                if entry[0][0] in d[0]:
                    self.updateEntry(d)
                    
    def getString(self):
        return ''.join(x or ' _ ' for x in self.chunks)
        
    def isDone(self):
        return None not in self.chunks

    def chunksDone(self):
        count = 0
        for c in self.chunks:
            if c is not None:
                count+=1
        return count

    def toDNA(self):
        text = ''
        for drops in self.droplets:
            text = text + drops.toDNA() + "\n"
        return text

    def writeToFile(self, output_file):
        OUT = open(output_file, 'wb')
        i = 0
        while i < self.num_chunks:
            OUT.write(self.chunks[i])
            i = i + 1
        OUT.close()
        return True
