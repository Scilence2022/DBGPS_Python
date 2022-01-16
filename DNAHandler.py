from utils import DNA_rev_complement
#version=20200922
import numpy as np
import random



#20200829 Rename Functions
class DNAHandler:
    def __init__(self):
        self.initSubProb = 0.001
        self.subProbSeed = 0.001
        #self.delProbSeed = self.subProbSeed * 10
        #self.insProbSeed = self.delProbSeed/10
        # self.termProb = 0.0005
        self.termProb = 0
        self.zero_copy_seqs = 0
        self.one_copy_seqs = 0

    def synthesis(self, dnastr):
        dnastr3 = dnastr[::-1]
        newSeq = ''
        i = 0
        while i < len(dnastr):
            aBase = dnastr3[i:i + 1]
            #if not deleted
            if not self.randomProbTrue(self.delProb(i)):
                if self.randomProbTrue(self.subProb(i)):
                    newSeq = newSeq + self.randomATGC()
                else:
                    newSeq = newSeq + aBase
            if self.randomProbTrue(self.insProb(i)):
                newSeq = newSeq + self.randomATGC()
            #if synthesis is terminated
            if self.randomProbTrue(self.termProb):
                return newSeq[::-1]
            i = i + 1

        return newSeq[::-1]




    def copy_randnum_with_sub(self, DNAs,  min_copies=0, max_copies=100, err_rate=0.01):
        nDNAs = self.copy_randnum(DNAs, min_copies, max_copies)
        self.add_rand_sub(nDNAs, err_rate=err_rate)
        return nDNAs

#20200911
    def random_int(self, exp_num, size=1, type="poisson"):
        if type == "poisson":
            return np.random.poisson(exp_num, size)
        if type == "random":
            return np.random.randint(0, exp_num*2+1, size)

#20201020
    def random_err_nums(self, exp_num, dna_len, size=1, type="poisson"):

        nums = self.random_int(exp_num, size, type)
        if type == 'random':
            add_error_base_seq_num = int((exp_num*2 - int(exp_num*2)) * size/2)
            e_seqs = np.random.randint(1, size, add_error_base_seq_num)
            for a in e_seqs:
                nums[a] += 1
        if type == 'poisson':
            add_error_base_seq_num = int((exp_num*2 - int(exp_num*2)) * size/2)
            e_seqs = np.random.randint(1, size, add_error_base_seq_num)
            for a in e_seqs:
                nums[a] += 1
        return nums

    def add_rand_brk_new(self, DNAs, err_rate=0.01, fix_num=False):
   # NOT tested yet 20210225
        total_dnas_len = 0
        for a in DNAs:
            total_dnas_len = total_dnas_len + len(a)

        if fix_num:
            site_num = int(total_dnas_len * err_rate)
        else:
            site_num = np.random.poisson(int(total_dnas_len * err_rate), 1)
        # print(site_num)
        sites = np.random.choice(range(1, total_dnas_len + 1), site_num, False)
        sites.sort()
        # print(sites)

        sites_list = sites.tolist()
        sites_list.reverse()
        # print(sites_list)

        added_length = 0
        add_seq_num = 0
        n_seqs = []
        anum = 0

        for seq in DNAs:
            # print('ori seq here-------------')
            # print(seq)
            seq_len = len(seq)
            if len(sites_list) < 1:
                anum = -1
            else:
                anum = sites_list.pop()

            o_seq = seq
            n_seq = ''
            while anum < len(o_seq) + added_length + 1 and anum > 0:
                #
                seq_pos = anum - added_length
                n_seqs.append(o_seq[0:seq_pos])
                o_seq = o_seq[seq_pos:]
                added_length = anum
                if len(sites_list) < 1:
                    anum = -1
                else:
                    anum = sites_list.pop()

            if len(o_seq) > 0:
                n_seqs.append(o_seq)
                added_length = added_length + len(o_seq)
        return n_seqs

    # 20201208
    def add_rand_sub_new(self, DNAs, err_rate=0.01, fix_num = False):

        total_dnas_len = 0
        for a in DNAs:
            total_dnas_len = total_dnas_len + len(a)
        if fix_num:
            site_num = int(total_dnas_len * err_rate)
        else:
            site_num = np.random.poisson(int(total_dnas_len * err_rate), 1)
        # print(site_num)
        sites = np.random.choice(range(1, total_dnas_len + 1), site_num, False)
        sites.sort()
        # print(sites)

        sites_list = sites.tolist()
        sites_list.reverse()
        # print(sites_list)

        added_length = 0
        add_seq_num = 0
        n_seqs = []
        anum = 0
        if len(sites_list) < 1:
            anum = -1
        else:
            anum = sites_list.pop()

        for seq in DNAs:
            # print('ori seq here-------------')
            # print(seq)
            seq_len = len(seq)

            o_seq = seq
            n_seq = ''
            while anum < len(o_seq) + added_length + 1 and anum > 0:
                #
                seq_pos = anum - added_length
                n_seq = n_seq + o_seq[0:seq_pos-1] + self.randomATGC(o_seq[seq_pos-1:seq_pos], 1)
                o_seq = o_seq[seq_pos:]
                added_length = anum

                if len(sites_list) < 1:
                    anum = -1
                else:
                    anum = sites_list.pop()

                # print(n_seq)

            if len(o_seq) > 0:
                n_seq = n_seq + o_seq
                added_length = added_length + len(o_seq)


            n_seqs.append(n_seq)
        return n_seqs

        # 20201208

    def add_rand_ins_new(self, DNAs, err_rate=0.01, fix_num=False):

        total_dnas_len = 0
        for a in DNAs:
            total_dnas_len = total_dnas_len + len(a)

        if fix_num:
            site_num = int(total_dnas_len * err_rate)
        else:
            site_num = np.random.poisson(int(total_dnas_len * err_rate), 1)
        # print(site_num)
        sites = np.random.choice(range(1, total_dnas_len + 1), site_num, False)
        sites.sort()
        # print(sites)

        sites_list = sites.tolist()
        sites_list.reverse()
        # print(sites_list)

        added_length = 0
        add_seq_num = 0
        n_seqs = []
        anum = 0
        if len(sites_list) < 1:
            anum = -1
        else:
            anum = sites_list.pop()

        for seq in DNAs:
            # print('ori seq here-------------')
            # print(seq)
            seq_len = len(seq)

            o_seq = seq
            n_seq = ''
            while anum < len(o_seq) + added_length + 1 and anum > 0:
                #
                seq_pos = anum - added_length

                # print(seq_pos)

                n_seq = n_seq + o_seq[0:seq_pos] + self.randomATGC()
                o_seq = o_seq[seq_pos:]
                added_length = anum
                if len(sites_list) < 1:
                    anum = -1
                else:
                    anum = sites_list.pop()
                # print(n_seq)

            if len(o_seq) > 0:
                n_seq = n_seq + o_seq
                added_length = added_length + len(o_seq)
            # # seq = n_seq
            # print('one seq')
            # print(seq)
            # print(n_seq)
            # print('\n')

            n_seqs.append(n_seq)
        return n_seqs

        # 20201208

    def add_rand_del_new(self, DNAs, err_rate=0.01, fix_num=False):

        total_dnas_len = 0
        for a in DNAs:
            total_dnas_len = total_dnas_len + len(a)
        site_num = np.random.poisson(int(total_dnas_len * err_rate), 1)
        # print(site_num)
        if fix_num:
            site_num = int(total_dnas_len * err_rate)
        else:
            site_num = np.random.poisson(int(total_dnas_len * err_rate), 1)

        sites = np.random.choice(range(1, total_dnas_len + 1), site_num, False)
        sites.sort()
        # print(sites)

        sites_list = sites.tolist()
        sites_list.reverse()
        # print(sites_list)

        added_length = 0
        add_seq_num = 0
        n_seqs = []
        anum = 0
        if len(sites_list) < 1:
            anum = -1
        else:
            anum = sites_list.pop()

        for seq in DNAs:
            # print('ori seq here-------------')
            # print(seq)
            seq_len = len(seq)
            o_seq = seq
            n_seq = ''
            while anum < len(o_seq) + added_length + 1 and anum > 0:
                #
                seq_pos = anum - added_length
                # print(seq_pos)
                n_seq = n_seq + o_seq[0:seq_pos - 1]
                o_seq = o_seq[seq_pos:]
                added_length = anum
                if len(sites_list) < 1:
                    anum = -1
                else:
                    anum = sites_list.pop()

                # print(n_seq)

            if len(o_seq) > 0:
                n_seq = n_seq + o_seq
                added_length = added_length + len(o_seq)
            # seq = n_seq
            n_seqs.append(n_seq)
        return n_seqs

    def add_rand_indel_new(self, DNAs, ins_err_rate=0.01, del_err_rate=0.01, fix_num=False):
        total_bases = 0
        for a in DNAs:
            total_bases = total_bases + len(a)

        if fix_num:
            site_num = int(total_bases * (ins_err_rate + del_err_rate))
        else:
            site_num = np.random.poisson(int(total_bases * (ins_err_rate + del_err_rate)), 1)
        # print(site_num)
        sites = np.random.choice(range(1, total_bases + 1), site_num, False)
        sites.sort()
        # print(sites)

        sites_list = sites.tolist()
        sites_list.reverse()
        # print(sites_list)

        added_length = 0
        add_seq_num = 0
        n_seqs = []
        anum = 0

        if len(sites_list) < 1:
            anum = -1
        else:
            anum = sites_list.pop()

        for seq in DNAs:
            # print('ori seq here-------------')
            # print(seq)
            seq_len = len(seq)
            o_seq = seq
            n_seq = ''
            while anum < len(o_seq) + added_length + 1 and anum > 0:
                #
                seq_pos = anum - added_length
                # print(seq_pos)
                if self.ins_or_del(ins_err_rate, del_err_rate):
                    n_seq = n_seq + o_seq[0:seq_pos] + self.randomATGC()

                else:
                    n_seq = n_seq + o_seq[0:seq_pos - 1]

                o_seq = o_seq[seq_pos:]
                added_length = anum
                if len(sites_list) < 1:
                    anum = -1
                else:
                    anum = sites_list.pop()

            if len(o_seq) > 0:
                n_seq = n_seq + o_seq
                added_length = added_length + len(o_seq)

            n_seqs.append(n_seq)
        return n_seqs

    def ins_or_del(self, ins_err_rate, del_err_rate):
        ins_num = int(ins_err_rate * 10000)
        del_num = int(del_err_rate * 10000)
        a = random.randint(1, ins_num + del_num)
        if a < ins_num:
            return True
        else:
            return False



    def copy_randnum(self, DNAs, min_copies=0, max_copies=100):
        nDNAs = []
        if min_copies == max_copies:
            #
            i = 0
            while i < len(DNAs):
                nDNAs.extend(self.copy_one_to_randnum(DNAs[i], min_copies))
                i = i + 1
        else:
            rd_nums = np.random.randint(min_copies, max_copies, size=len(DNAs))
            i = 0
            while i < len(DNAs):
                nDNAs.extend(self.copy_one_to_randnum(DNAs[i], rd_nums[i]))
                i = i + 1
        return nDNAs

    def copy_seqs(self, DNAs, copies=100):
        nDNAs = []
        for seq in DNAs:
            for n in range(0,copies):
                nDNAs.append(seq)
        return nDNAs

    def copy_seq_poisson(self, DNAs, exp_num=50):
        nDNAs = []
        rd_nums = np.random.poisson(exp_num, len(DNAs))
        i = 0
        while i < len(DNAs):
            # rd_num = poisson_seq_num(10, 0.5, rep_times)
            nDNAs.extend(self.copy_one_to_randnum(DNAs[i], rd_nums[i]))
            i = i + 1
        return nDNAs

    def copy_one_to_randnum(self, DNA, copies):
        aa = []
        for a in range(0,copies):
            aa.append(DNA)
        return aa

    def random_indel_num(self, err_rate):
        num = 1
        while(self.randomProbTrue(err_rate)):
            num = num + 1
        return num

    def randomDNA(self, len):
        dna = ["A", "G", "C", "T"]
        random_sequence = ''
        for i in range(0, len):
            random_sequence += random.choice(dna)
        return random_sequence

    def randomATGC(self, noATGC='', len=1):
        all_bases = ["A", "G", "C", "T"]
        allow_ATGC = []

        for aBase in all_bases:
            if not aBase in noATGC:
                allow_ATGC.append(aBase)

        random_sequence = ''
        for i in range(0, len):
            random_sequence += random.choice(allow_ATGC)
        return random_sequence

    def randomDNAs(self, len, num):
        arr1 = []
        for i in range(1, num):
            arr1.append(self.randomDNA(len))
        return arr1

    def degenerate_DNA(self, dnastr, rate):

        num_break_points = int(np.random.exponential(rate, 1) *len(dnastr))
        if(num_break_points >= len(dnastr) - 1):
            return []
        if(num_break_points > 0):
            return self.break_DNA(dnastr, num_break_points)
        else:
            return [dnastr]

    def degenerate_DNAs_random_copies(self, dnas, rate, min_copies, max_copies):
        self.arr_deg_dnas = []
        self.zero_copy_seqs = 0
        self.one_copy_seqs = 0

        for dnastr in dnas:
            if(min_copies == max_copies):
                copies = min_copies
            else:
                copies = random.randint(min_copies, max_copies)
                if(copies == 1):
                    self.one_copy_seqs = self.one_copy_seqs + 1
                if (copies == 0):
                    self.zero_copy_seqs = self.zero_copy_seqs + 1
            j = 0
            arr_break_points_nums = []
            arr_break_points_nums = np.random.poisson(int(len(dnas[0]) * rate), size=copies)
            # print(arr_break_points_nums)
            self.arr_break_points_nums = arr_break_points_nums
            while(j < copies):
                break_point_num = arr_break_points_nums[j]
                if(break_point_num < len(dnastr)):
                    self.arr_deg_dnas.extend(self.break_DNA(dnastr,break_point_num))
                j = j + 1
        return self.arr_deg_dnas


    def break_DNA(self, dnastr, num_break_points):
        # print("DNAstr len: " + str(len(dnastr)) + "num_break_points: " + str(num_break_points))
        rand_ints = random.sample(range(0, len(dnastr)), num_break_points)
        rand_ints.sort()
        pieces_of_dna = []
        b_pos = 0
        for nums in rand_ints:
            if(nums - b_pos > 0):
                pieces_of_dna.append(dnastr[b_pos:nums])
            b_pos = nums + 1
        pieces_of_dna.append(dnastr[b_pos:len(dnastr)])
        return pieces_of_dna

    def break_ligation(self, dnas, rate, itr):
        i = 0
        n_dnas = []
        n_dnas.extend(dnas)

        while i < itr:
            nn_dnas = []
            for a in n_dnas:
                nn_dnas.extend(self.break_DNA(a, round(int(len(a) * rate))))
            n_dnas = self.random_ligation(nn_dnas)
            i = i + 1
        return n_dnas


    def random_ligation(self, dnas):
        random.shuffle(dnas)
        n_dnas = []
        i = int(len(dnas)/2)
        j=0

        while j < i:
            if random.choice([True, False]):
                n_dnas.append(dnas[j*2] + dnas[j*2+1])
            else:
                #random reverse ligation
                n_dnas.append(dnas[j * 2] + DNA_rev_complement(dnas[j * 2 + 1]))
            j = j + 1
        return n_dnas

    # def synthesisDNAs(self, dnastr, copy_num):
    #     seqs = []
    #     i = 0
    #     while i< copy_num:
    #         seqs.append(self.synthesisDNA(dnastr))
    #         i = i + 1
    #     return seqs

    def randomProbTrue(self, prob):
        rd = random.random()
        if rd <= prob:
            return True
        else:
            return False


    def delProb(self, pos):
        return self.subProb(pos)/2

    def subProb(self, pos):
        return pos * self.subProbSeed + self.initSubProb

    def insProb(self, pos):

        return self.subProb(pos)/10

    # def randomATGC(self):
    #     a = "ATGC"
    #     ar = random.randint(0,3)
    #     return a[ar:ar+1]

    def randNums(self, minNum, maxNum, num):
        arr1 = []
        i = 0
        while(len(arr1) < num):
            arr1.append(random.randint(minNum,maxNum))
        return arr1



""" Old functions, to be deleted. 
    def copy_one_with_sub(self, dnastr, copies, err_rate):
        dna_len = len(dnastr)
        seq_copies = []
        site_num_list = np.random.poisson(int(dna_len * err_rate), size=copies)
        for site_num in site_num_list:
            sites = random.sample(range(1, dna_len + 1), site_num)
            sites.sort()
            aSeqCopy = ''
            lastSite = 0
            for site in sites:
                if site - lastSite > 1:
                    aSeqCopy = aSeqCopy + dnastr[lastSite:site-1]
                    aSeqCopy = aSeqCopy + self.randomATGC(dnastr[site-1:site],1)

                else:
                    aSeqCopy = aSeqCopy + self.randomATGC(dnastr[site - 1:site], 1)
                lastSite = site

            if lastSite < dna_len:
                aSeqCopy = aSeqCopy + dnastr[lastSite:]
            seq_copies.append(aSeqCopy)
        return seq_copies


    def copy_randnum_with_ins(self, DNAs,  min_copies=0, max_copies=100, err_rate=0.01):
        arr_copy_seqs = []
        for dnastr in DNAs:
            copies = 0
            # print(dnastr)
            if(min_copies == max_copies):
                copies = min_copies
            else:
                copies = random.randint(min_copies, max_copies)
            if(copies == 1):
                self.one_copy_seqs = self.one_copy_seqs + 1
            else:
                if(copies == 0):
                    self.zero_copy_seqs = self.zero_copy_seqs + 1
                else:
                    arr_copy_seqs.extend(self.copy_one_with_ins(dnastr,copies,err_rate))
        return arr_copy_seqs

    def copy_one_with_ins(self, dnastr, copies, err_rate):
        dna_len = len(dnastr)
        seq_copies = []
        site_num_list = np.random.poisson(int(dna_len * err_rate)+1, size=copies)
        for site_num in site_num_list:
            # print(dna_len)
            # print(site_num)
            # print("")
            if site_num > dna_len:
                site_num = dna_len
            sites = random.sample(range(0, dna_len), site_num)
            sites.sort()
            aSeqCopy = ''
            lastSite = 0
            for site in sites:
                # if site - lastSite > 1:
                aSeqCopy = aSeqCopy + dnastr[lastSite:site]
                aSeqCopy = aSeqCopy + self.randomDNA(self.random_indel_num(err_rate))

                # else:
                #
                #     aSeqCopy = aSeqCopy + self.randomDNA(self.random_indel_num(err_rate))
                lastSite = site

            if lastSite < dna_len:
                aSeqCopy = aSeqCopy + dnastr[lastSite:]
            seq_copies.append(aSeqCopy)
        return seq_copies

    def copy_randnum_with_del(self, DNAs,  min_copies=0, max_copies=100, err_rate=0.01):
        arr_copy_seqs = []
        for dnastr in DNAs:
            copies = 0
            # print(dnastr)
            if(min_copies == max_copies):
                copies = min_copies
            else:
                copies = random.randint(min_copies, max_copies)
            if(copies == 1):
                self.one_copy_seqs = self.one_copy_seqs + 1
            else:
                if(copies == 0):
                    self.zero_copy_seqs = self.zero_copy_seqs + 1
                else:
                    arr_copy_seqs.extend(self.copy_one_with_del(dnastr,copies,err_rate))
        return arr_copy_seqs


    def copy_one_with_del(self, dnastr, copies, err_rate):
        dna_len = len(dnastr)
        seq_copies = []
        site_num_list = np.random.poisson(int(dna_len * err_rate)+1, size=copies)
        for site_num in site_num_list:
            # print(dna_len)
            # print(site_num)
            # print("")
            if site_num > dna_len:
                site_num = dna_len
            sites = random.sample(range(0, dna_len), site_num)
            sites.sort()
            aSeqCopy = ''
            lastSite = 0
            for site in sites:
                # if site - lastSite > 1:
                if site >= lastSite:
                    aSeqCopy = aSeqCopy + dnastr[lastSite:site]
                    lastSite = site + self.random_indel_num(err_rate)
            if lastSite < dna_len:
                aSeqCopy = aSeqCopy + dnastr[lastSite:]
            seq_copies.append(aSeqCopy)
        return seq_copies



    def copyDNAsWithIndelErr(self, DNAs,  min_copies=0, max_copies=100, ins_err_rate=0.01, del_err_rate=0.01):
        arr_copy_seqs = []
        for dnastr in DNAs:
            copies = 0
            # print(dnastr)
            if(min_copies == max_copies):
                copies = min_copies
            else:
                copies = random.randint(min_copies, max_copies)
            if(copies == 1):
                self.one_copy_seqs = self.one_copy_seqs + 1
            else:
                if(copies == 0):
                    self.zero_copy_seqs = self.zero_copy_seqs + 1
                else:
                    arr_copy_seqs.extend(self.copyDNAWithIndelErr(dnastr,copies,err_rate))
        return arr_copy_seqs


    def copyDNAWithIndelErr(self, dnastr, copies, ins_err_rate=0.01, del_err_rate=0.01):
        dna_len = len(dnastr)
        seq_copies = []
        site_num_list = np.random.poisson(int(dna_len * err_rate)+1, size=copies)
        for site_num in site_num_list:
            # print(dna_len)
            # print(site_num)
            # print("")
            if site_num > dna_len:
                site_num = dna_len
            sites = random.sample(range(0, dna_len), site_num)
            sites.sort()
            aSeqCopy = ''
            lastSite = 0
            for site in sites:
                # if site - lastSite > 1:
                if site >= lastSite:
                    aSeqCopy = aSeqCopy + dnastr[lastSite:site]
                    lastSite = site + self.random_indel_num(err_rate)
            if lastSite < dna_len:
                aSeqCopy = aSeqCopy + dnastr[lastSite:]
            seq_copies.append(aSeqCopy)
        return seq_copies
"""




'''
    #20201208
    def add_rand_indel(self, DNAs, err_rate=0.01, type="poisson"):

        dna_len = len(DNAs[0])
        exp_num = dna_len * err_rate

        site_num_list = self.random_err_nums(exp_num, dna_len, len(DNAs), type)


        i = 0
        while i < len(site_num_list):
            site_num = site_num_list[i]
            dnastr = DNAs[i]
            dna_len = len(dnastr)
            if site_num > dna_len - 1:
                site_num = dna_len - 1
            #Modified to adapt to various length of DNAs
            sites = random.sample(range(1, dna_len + 1), site_num)
            sites.sort()
            aseq_copy = ''
            last_site = 0
            for site in sites:
                if site >= last_site:
                    aseq_copy = aseq_copy + dnastr[last_site:site]
                    last_site = site + self.random_indel_num(err_rate)
            if last_site < dna_len:
                aseq_copy = aseq_copy + dnastr[last_site:]
            DNAs[i] = aseq_copy
            i = i + 1

    def add_rand_sub(self, DNAs, err_rate=0.01, type = "poisson"):


        dna_len = len(DNAs[0])
        site_num_list = self.random_err_nums(int(dna_len * err_rate), dna_len, len(DNAs), type)

        i = 0
        while i < len(site_num_list):
            site_num = site_num_list[i]
            dnastr = DNAs[i]
            sites = random.sample(range(1, dna_len + 1), site_num)
            sites.sort()
            aseq_copy = ''
            last_site = 0
            for site in sites:
                if site - last_site > 1:
                    aseq_copy = aseq_copy + dnastr[last_site:site-1]
                    aseq_copy = aseq_copy + self.randomATGC(dnastr[site-1:site],1)

                else:
                    aseq_copy = aseq_copy + self.randomATGC(dnastr[site - 1:site], 1)
                last_site = site

            if last_site < dna_len:
                aseq_copy = aseq_copy + dnastr[last_site:]
            DNAs[i] = aseq_copy
            i = i + 1


    def add_rand_del(self, DNAs, err_rate=0.01, type="poisson"):

        dna_len = len(DNAs[0])
        exp_num = dna_len * err_rate

        site_num_list = self.random_err_nums(exp_num, dna_len, len(DNAs), type)


        i = 0
        while i < len(site_num_list):
            site_num = site_num_list[i]
            dnastr = DNAs[i]
            dna_len = len(dnastr)
            if site_num > dna_len - 1:
                site_num = dna_len - 1
            #Modified to adapt to various length of DNAs
            sites = random.sample(range(1, dna_len + 1), site_num)
            sites.sort()
            aseq_copy = ''
            last_site = 0
            for site in sites:
                if site >= last_site:
                    aseq_copy = aseq_copy + dnastr[last_site:site]
                    last_site = site + self.random_indel_num(err_rate)
            if last_site < dna_len:
                aseq_copy = aseq_copy + dnastr[last_site:]
            DNAs[i] = aseq_copy
            i = i + 1



    def add_rand_ins(self, DNAs, err_rate=0.01, type="poisson"):

        dna_len = len(DNAs[0])

        exp_num = dna_len * err_rate

        site_num_list = self.random_err_nums(exp_num, dna_len, len(DNAs), type)

        i = 0
        while i < len(site_num_list):
            site_num = site_num_list[i]
            dnastr = DNAs[i]

            dna_len = len(dnastr)
            if site_num > dna_len - 1:
                site_num = dna_len - 1

            sites = random.sample(range(1, dna_len + 1), site_num)
            sites.sort()
            aseq_copy = ''
            last_site = 0
            for site in sites:
                # if site - last_site > 1:
                aseq_copy = aseq_copy + dnastr[last_site:site]
                aseq_copy = aseq_copy + self.randomDNA(self.random_indel_num(err_rate))

                # else:
                #     aseq_copy = aseq_copy + self.randomATGC(dnastr[site - 1:site], 1)
                last_site = site

            if last_site < dna_len:
                aseq_copy = aseq_copy + dnastr[last_site:]
            DNAs[i] = aseq_copy
            i = i + 1

'''