import re
import crc16pure
from utils import *
from DNAdroplet import DNADroplet
# from crc16pure import crc16xmodem


class DeBruijnGraph:
    def __init__(self, data=bytes(), kmer_len=21):
        self.data = data
        self.kmers = {}
        self.kmer_len = kmer_len
        # self.veri_kmers = False
        self.pKmer = {}
        self.cRule = {}
        self.nextBlocks = {}
        self.prevBlocks = {}

        self.nextBases = {}
        self.prevBases = {}

        self.both_search = False
        self.pathsA = []
        self.pathsB = []
        self.pathAB = []

        self.crc_paths = []

        # self.overlap = 20
        self.pathABLen = 0
        # self.dataEncodingDnaLength = 150
        self.ratio_tolerance = 200


        self.set_max_kmer_num = False
        self.max_path_num = 1000
        self.max_num_of_repeat_kmer = 4
        # self.set_compress_converter()

        # Just for test only, will be removed later
        self.primerF = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
        self.primerE = 'CTGACACTGATGCATCCG'  # complement seq of P2



    # 20200807 modified to fix a bug with quanlity scores
    def open_fastq(self, file):
        f = open(file, "r")
        matchLineB = re.compile('^(\+)')
        last_line = f.readline()
        line = f.readline()
        reads = 0
        while line.strip():
            if matchLineB.match(line):
                self.add_seq(last_line.strip())
                reads = reads + 1
                if reads > 1000:
                    #print('.', end='')
                    reads = 0
            else:
                last_line = line
            line = f.readline()

    #open specific number of FastQ sequences
    def open_fastq_num(self, file, num =12000 * 100):
        f = open(file, "r")
        matchLineB = re.compile('^(\+)')
        last_line = f.readline()
        line = f.readline()
        reads = 0
        read_num = 1
        while line.strip() and read_num <= num:
            if matchLineB.match(line):
                self.add_seq(last_line.strip())
                reads = reads + 1
                if reads > 1000:
#                    print('.', end='')
                    reads = 0
            else:
                last_line = line
            line = f.readline()
            read_num = read_num + 1

    def open_dump(self, file):
        f = open(file, "r")
        line = f.readline()
        kmer_num = 0

        #Check K-mer length
        while (line.strip()):
            arr = line.split('\t')
            # print(arr)
            self.add_kmer(arr[0], int(arr[1]))
            rev_seq = DNA_rev_complement(arr[0])
            if rev_seq != arr[0]:
                self.add_kmer(rev_seq, int(arr[1]))
            kmer_num = kmer_num + 1
            if kmer_num > 10000000:
#                print('.', end='')
                kmer_num = 0
            line = f.readline()

    # 20200502
    def open_fasta(self, file):
        print('Opening fasta file')
        f = open(file, "r")
        matchLineA = re.compile('^(>)')
        line = f.readline()
        seq_num = 0
        while (line.strip()):
            if not matchLineA.match(line):
                self.add_seq(line.strip())
            seq_num = seq_num + 1
            if seq_num > 1000:
#                print('.', end='')
                seq_num = 0
            line = f.readline()

    # 20201101
    def open_row_seq_file(self, file):
        f = open(file, "r")
        line = f.readline()
        while (line.strip()):
            self.add_seq(line.strip())
            line = f.readline()

    def add_seqs(self, seqs):
        for seq in seqs:
            self.add_seq(seq)

    # Modified to handle NNNs 2019/11/21 Lifu Song
    # Modified to handle complement Seq
    # Modified to handle degenerate bases 2020/01/16
    def add_seq(self, str):
        if len(str) >= self.kmer_len:
            i = 0
            kmstr = ''
            while i <= len(str) - self.kmer_len:
                kmstr = str[i:i + self.kmer_len]
                if not ("N" in kmstr):
                    self.add_kmer(kmstr)
                    self.add_kmer(DNA_rev_complement(kmstr))
                else:
                    if kmstr.count("N") == 1:
                        self.add_kmer(kmstr.replace("N", "A"))
                        self.add_kmer(kmstr.replace("N", "T"))
                        self.add_kmer(kmstr.replace("N", "G"))
                        self.add_kmer(kmstr.replace("N", "C"))
                        self.add_kmer(DNA_rev_complement(kmstr.replace("N", "A")))
                        self.add_kmer(DNA_rev_complement(kmstr.replace("N", "T")))
                        self.add_kmer(DNA_rev_complement(kmstr.replace("N", "G")))
                        self.add_kmer(DNA_rev_complement(kmstr.replace("N", "C")))
                i = i + 1


    # 2020/05/02
    # 2021/06/15
    def add_kmer(self, kmer, num=1):
        # comp_kmer = self.compressDNA(kmer)
        if kmer in self.kmers:
            self.kmers[kmer] = self.kmers[kmer] + num
        else:
            self.kmers[kmer] = num


    def kmer_freq_in_path(self, kmer, path):
        kmers = {}
        kmer_len = len(kmer)
        if len(path) >= kmer_len:
            i = 0
            kmstr = ''
            while i <= len(path) - kmer_len:
                kmstr = path[i:i + kmer_len]
                if kmstr in kmers.keys():
                    kmers[kmstr] = kmers[kmstr] + 1
                else:
                    kmers[kmstr] = 1
                i = i + 1
        if kmer in kmers.keys():
            return kmers[kmer]
        else:
            return 0

    # This function is for test only
    def highst_km_freq(self, str):
        maxFreq = 0
        if len(str) >= self.kmer_len:
            i = 0
            kmstr = ''
            while i <= len(str) - self.kmer_len:
                kmstr = str[i:i + self.kmer_len]
                if self.kmer_freq(kmstr) > maxFreq:
                    maxFreq = self.kmer_freq(kmstr)
                i = i + 1
        return maxFreq


    def set_max_kmer_value(self):
        if not self.set_max_kmer_num:
            avg_value = self.get_avg_kmer_value()
            for aKmer in self.kmers.keys():
                if self.kmers[aKmer] > avg_value:
                    self.kmers[aKmer] = avg_value
            self.set_max_kmer_num = True

    def set_min_kmer_value(self):
        # avg_value = self.get_avg_kmer_value()
        for aKmer in self.kmers.keys():
            if self.kmers[aKmer] < 2:
                self.kmers[aKmer] = 0

    def get_avg_kmer_value(self):
        allvalue = 0
        num_of_kmers = len(self.kmers)
        for aKmer in self.kmers.keys():
            if self.kmers[aKmer] < 2:
                num_of_kmers = num_of_kmers - 1
            else:
                allvalue = allvalue + self.kmers[aKmer]

        return int(allvalue / num_of_kmers) + 1


    def kmer_freq(self, kmstr):
        if kmstr in self.kmers.keys():
            return self.kmers[kmstr]
        else:
            return 0

    #20210114
    #20210614
    def has_kmer(self, kmer):
        if kmer in self.kmers:
            return True
        else:
            return False


    def test_find_path_grass(self, seq, block_len=10):
        path_len = len(seq) - self.kmer_len
        self.pathsA = [seq[:self.kmer_len]]

        i = 1
        b_pos = 1
        while i < path_len+1 and len(self.pathsA) <= self.max_path_num:
            tmpPaths = []
            for path in self.pathsA:
                self.pKmer['a'] = path[len(path) - self.kmer_len:len(path)]
                maxScore = self.score_next_bases()
                for block in self.nextBases.keys():
                    if self.nextBases[block] >= maxScore / self.ratio_tolerance and self.nextBases[block] > 0:
                        # acceptable block
                        tmpPaths.append(path + block)
            self.pathsA = tmpPaths

            #Removing the error detected paths. 3/4 of the paths with two and more error bases were removed randomly.
            if b_pos == block_len:
                c_kmer = seq[i:i+self.kmer_len]
                #print(c_kmer)
                tmpPaths = []
                path_kmer_dict = {}
                for path in self.pathsA:
                    path_tail_kmer = path[-self.kmer_len:]
                    if c_kmer in path_tail_kmer:
                        tmpPaths.append(path)
                    else:
                        if self.kmer_distance(path_tail_kmer, c_kmer) > 1:
                            if randomATGC() == "A":
                                tmpPaths.append(path)

            if len(tmpPaths) > 0:
                self.pathsA = tmpPaths
            else:
                return False

            i = i + 1
            b_pos = b_pos + 1
            if b_pos > block_len:
                b_pos = 1

        if len(self.pathsA) > 0 and len(self.pathsA[0]) == len(seq):
            for path in self.pathsA:
                if seq in path:
                    return True
        else:
            return False

    def kmer_distance(self, k1, k2):
        d = 0
        for i in range(0, len(k1)):
            if k1[i:i + 1] != k2[i:i + 1]:
                d = d + 1
        return d

    def find_droplet_DNA(self, index, data_bytes_len):
        adrop = DNADroplet(bytes(data_bytes_len))
        adrop.head_index = index
        # adrop.get_tail_index()
        indexa = adrop.get_head_index_dna()
        # indexb = adrop.get_tail_index_dna()
        self.get_crc_path(indexa, data_bytes_len * 4 + adrop.crc_len, adrop.crc_len)

    def get_crc_path(self, indexa, path_len, crc_len=8):
        self.crc_paths = []
        # if self.both_search:
        #     self.findSimPaths(indexa, indexb, path_len)
        # else:
        self.find_paths_se(indexa, path_len)
        for path in self.pathAB:
            pF_len = len(self.primerF)
            pathTrimPrimer = path[pF_len:]

            if self.crc_path_check(pathTrimPrimer, len(indexa), path_len, crc_len):
                self.crc_paths.append(pathTrimPrimer)
            # else:
            #     print(len(self.pathAB))
            #     print("Wrong CRC Path detected")

    def crc_path_check(self, dnastr, indexa_dna_length, data_dna_length, crc_dna_length=8):
        # print('hello')
        #print(dnastr)
        # print(dnastr[0:indexa_dna_length + data_dna_length - crc_dna_length])
        # print(dnastr[0:indexa_dna_length + data_dna_length - crc_dna_length])
        # exit()
        data_bytes = DNAToBytes(dnastr[0:indexa_dna_length + data_dna_length - crc_dna_length])
        crc_bytes = DNAToBytes(
            dnastr[indexa_dna_length + data_dna_length - crc_dna_length:indexa_dna_length + data_dna_length])

        if crc16pure.crc16xmodem(data_bytes) == int.from_bytes(crc_bytes, byteorder='big', signed=False):
            return True
        else:
            return False


    def find_paths_se(self, indexa, path_len):

        self.pathAB = []
        self.pathsA = [self.primerF + indexa]
        self.pathABLen = len(indexa) + path_len
        # print(self.pathABLen)
        i = 0
        while i < path_len and len(self.pathsA) < self.max_path_num:
            tmpPaths = []
            for path in self.pathsA:
                self.pKmer['a'] = path[len(path) - self.kmer_len:len(path)]
                maxScore = self.score_next_bases()
                for block in self.nextBases.keys():
                    if self.nextBases[block] >= maxScore / self.ratio_tolerance and self.nextBases[block] > 0:
                        # acceptable block
                        tmpPaths.append(path + block)

            if len(tmpPaths) > 0:
                self.pathsA = tmpPaths
            else:
                return 'No path found'
            i = i + 1
        if len(self.pathsA) > 0 and len(self.pathsA[0]) == self.pathABLen + len(self.primerF):
        #print(len(self.pathsA[0]))
            self.pathAB = self.pathsA




    def score_next_bases(self):
        self.nextBases = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        maxScore = 0
        for m in self.nextBases.keys():
            self.nextBases[m] = self.score_next_base(m)
            if self.nextBases[m] > maxScore:
                maxScore = self.nextBases[m]
        return maxScore

    def score_prev_bases(self):
        self.prevBases = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        maxScore = 0
        for m in self.prevBases.keys():
            self.prevBases[m] = self.score_prev_base(m)
            if self.prevBases[m] > maxScore:
                maxScore = self.prevBases[m]
        return maxScore

    def score_next_base(self, aBase):
        pKmer = self.pKmer['a']
        nKmer1 = pKmer[1:len(pKmer)] + aBase[0:1]
        if not self.kmers.get(nKmer1):
            return 0
        return self.kmers.get(nKmer1)

    def score_prev_base(self, aBase):
        pKmer = self.pKmer['b']
        nKmer1 = aBase[0:1] + pKmer[0:len(pKmer) - 1]
        if not self.kmers.get(nKmer1):
            return 0
        return self.kmers.get(nKmer1)


    # 2019-05-14
    # 2021-06-15
    def max_kmer_freq(self, str):
        max_freq = 0
        if len(str) >= self.kmer_len:
            i = 0
            while i <= len(str) - self.kmer_len:
                kmstr = str[i:i + self.kmer_len]
                if self.kmer_freq(kmstr) > max_freq:
                    max_freq = self.kmer_freq(kmstr)
                i = i + 1
        return max_freq


    # 2020-06-24
    def remove_low_cov_kmers(self, min_cov=10):
        if min_cov > 0:
            nKmlist = {}
            for a in self.kmers:
                if (self.kmers[a] > min_cov):
                    nKmlist[a] = self.kmers[a]
            self.kmers = nKmlist



    # Modified to make sure the pathAB length is correct 2020/02/05
    def join_paths(self):
        self.pathAB = []
        pF_len = len(self.primerF)
        pE_len = len(self.primerE)
        for leftpath in self.pathsA:
            for rightpath in self.pathsB:
                overlap = len(leftpath) - pF_len + len(rightpath) - pE_len - self.pathABLen
                if overlap > 0:
                    if leftpath[len(leftpath) - overlap:] == rightpath[:overlap]:
                        self.pathAB.append(leftpath + rightpath[overlap:])

    def nodes_for_cyto(self):
        aDict = {}
        nodes_text = ""
        for kk in self.kmers.keys():
            # if kk not in aDict.keys():
            aa = kk[1:len(kk)]
            if kk not in self.pathAB[0]:
                if (aa + "A") in self.kmers.keys():
                    nodes_text = nodes_text + kk + "\t-\t" + aa + "A\n"

                if (aa + "T") in self.kmers.keys():
                    nodes_text = nodes_text + kk + "\t-\t" + aa + "T\n"

                if (aa + "G") in self.kmers.keys():
                    nodes_text = nodes_text + kk + "\t-\t" + aa + "G\n"

                if (aa + "C") in self.kmers.keys():
                    nodes_text = nodes_text + kk + "\t-\t" + aa + "C\n"
                    # aDict[aa + "C"] = self.kmerList[aa + "C"]
            else:
                if (aa + "A") in self.kmers.keys():
                    if (aa + "A") in self.pathAB[0]:
                        nodes_text = nodes_text + kk + "\t+\t" + aa + "A\n"
                    else:
                        nodes_text = nodes_text + kk + "\t-\t" + aa + "A\n"

                if (aa + "T") in self.kmers.keys():
                    if (aa + "T") in self.pathAB[0]:
                        nodes_text = nodes_text + kk + "\t+\t" + aa + "T\n"
                    else:
                        nodes_text = nodes_text + kk + "\t-\t" + aa + "T\n"

                if (aa + "G") in self.kmers.keys():
                    if (aa + "G") in self.pathAB[0]:
                        nodes_text = nodes_text + kk + "\t+\t" + aa + "G\n"
                    else:
                        nodes_text = nodes_text + kk + "\t-\t" + aa + "G\n"

                if (aa + "C") in self.kmers.keys():
                    if (aa + "C") in self.pathAB[0]:
                        nodes_text = nodes_text + kk + "\t+\t" + aa + "C\n"
                    else:
                        nodes_text = nodes_text + kk + "\t-\t" + aa + "C\n"
                    # aDict[aa + "C"] = self.kmerList[aa + "C"]

        return nodes_text

    def path_kmers(self):
        nodes_text = ""
        tt = self.pathAB[0]
        i = 0
        while i < len(tt) - 12:
            nodes_text = nodes_text + tt[i:i + 12] + "\n"
            i = i + 1
        return nodes_text


    #
    # def scoreNextBlocks(self):
    #     # pKmer = dnastr[len(dnastr)-self.kmer_len:len(dnastr)]
    #     # pKmer = self.pKmer
    #     cRule = self.cRule['a']
    #     # pHandle = pKmer[3:len(pKmer)]
    #     # allValidBlocks = {}
    #     self.nextBlocks = triBlocksOfRule(cRule)
    #     maxScore = 0
    #     for m in self.nextBlocks.keys():
    #         self.nextBlocks[m] = self.scoreNextBlock(m)
    #         if self.nextBlocks[m] > maxScore:
    #             maxScore = self.nextBlocks[m]
    #     return maxScore
    #
    # def scorePrevBlocks(self):
    #     # pKmer = dnastr[len(dnastr)-self.kmer_len:len(dnastr)]
    #     # pKmer = self.pKmer
    #     cRule = self.cRule['b']
    #     # pHandle = pKmer[3:len(pKmer)]
    #     # allValidBlocks = {}
    #     self.prevBlocks = triBlocksOfRule(cRule)
    #     maxScore = 0
    #     for m in self.prevBlocks.keys():
    #         self.prevBlocks[m] = self.scorePrevBlock(m)
    #         if self.prevBlocks[m] > maxScore:
    #             maxScore = self.prevBlocks[m]
    #     return maxScore
    #
    # def scoreNextBlock(self, someBases):
    #     pKmer = self.pKmer['a']
    #     nKmer1 = pKmer[1:len(pKmer)] + someBases[0:1]
    #     nKmer2 = pKmer[2:len(pKmer)] + someBases[0:2]
    #     nKmer3 = pKmer[3:len(pKmer)] + someBases[0:3]
    #     cnKmer1 = self.compressDNA(nKmer1)
    #     cnKmer2 = self.compressDNA(nKmer2)
    #     cnKmer3 = self.compressDNA(nKmer3)
    #     if not self.kmerList.get(cnKmer1):
    #         return 0
    #     if not self.kmerList.get(cnKmer2):
    #         return 0
    #     if not self.kmerList.get(cnKmer3):
    #         return 0
    #     return self.kmerList.get(cnKmer1) * self.kmerList.get(cnKmer2) * self.kmerList.get(cnKmer3)
    #
    # def scorePrevBlock(self, someBases):
    #     pKmer = self.pKmer['b']
    #     nKmer1 = someBases[2:3] + pKmer[0:len(pKmer) - 1]
    #     nKmer2 = someBases[1:3] + pKmer[0:len(pKmer) - 2]
    #     nKmer3 = someBases[0:3] + pKmer[0:len(pKmer) - 3]
    #     cnKmer1 = self.compressDNA(nKmer1)
    #     cnKmer2 = self.compressDNA(nKmer2)
    #     cnKmer3 = self.compressDNA(nKmer3)
    #     if not self.kmerList.get(cnKmer1):
    #         return 0
    #     if not self.kmerList.get(cnKmer2):
    #         return 0
    #     if not self.kmerList.get(cnKmer3):
    #         return 0
    #     return self.kmerList.get(cnKmer1) * self.kmerList.get(cnKmer2) * self.kmerList.get(cnKmer3)

    # def set_compress_converter(self):
    #
    #     converter = ('A', 'T', 'G', 'C')
    #     # self.bytesConverter = {}
    #     # self.bytesConverter['A'] = 0
    #     # self.bytesConverter['T'] = 1
    #     # self.bytesConverter['G'] = 2
    #     # self.bytesConverter['C'] = 3
    #
    #     i = 0
    #     while i < 4:
    #         j = 0
    #         while j < 4:
    #             k = 0
    #             while k < 4:
    #                 m = 0
    #                 while m < 4:
    #                     fourBases = converter[i] + converter[j] + converter[k] + converter[m]
    #                     self.bytesConverter[fourBases] = bytes([i * 4 * 4 * 4 + j * 4 * 4 + k * 4 + m])
    #                     self.ATGConverter[bytes([i * 4 * 4 * 4 + j * 4 * 4 + k * 4 + m])] = fourBases
    #                     self.ATGConverter[i * 4 * 4 * 4 + j * 4 * 4 + k * 4 + m] = fourBases
    #                     m = m + 1
    #                 k = k + 1
    #             j = j + 1
    #         i = i + 1
    #
    # def compressDNA(self, dnaStr):
    #     # return dnaStr
    #     if self.compress:
    #         i = 0
    #         nBytes = b''
    #         while i < len(dnaStr):
    #             nBytes = nBytes + self.bytesConverter[dnaStr[i:i + 4]]
    #             i = i + 4
    #         return nBytes
    #     else:
    #         return dnaStr
    #
    # def compress_kmer(self, kmer):
    #     # return dnaStr
    #     a = len(kmer) % 4
    #     akmer = 'A' * a
    #     akmer = akmer + kmer
    #
    #     if self.compress:
    #         i = 0
    #         nBytes = b''
    #         while i < len(akmer):
    #             nBytes = nBytes + self.bytesConverter[akmer[i:i + 4]]
    #             i = i + 4
    #         return nBytes
    #     else:
    #         return kmer
    #
    # def depress_kmer(self, manyBytes):
    #     # define the converter of four bits to DNA chars
    #     # converter = ('A', 'T', 'G', 'C')
    #     dnastr = ''
    #     for aByte in manyBytes:
    #              dnastr = dnastr + self.ATGConverter[aByte]
    #     return dnastr

    # def depressDNA(self, manyBytes):
    #     # define the converter of four bits to DNA chars
    #     converter = ('A', 'T', 'G', 'C')
    #     dnastr = ''
    #     for aByte in manyBytes:
    #         fourBase = ''
    #         fourBase = converter[aByte % 4] + fourBase
    #         aByte = aByte // 4
    #         aByte = aByte // 4
    #         fourBase = converter[aByte % 4] + fourBase
    #         aByte = aByte // 4
    #         fourBase = converter[aByte % 4] + fourBase
    #         aByte = aByte // 4
    #         fourBase = converter[aByte % 4] + fourBase
    #         dnastr = dnastr + fourBase
    #     return dnastr