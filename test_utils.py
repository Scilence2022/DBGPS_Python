# import random
# import fountain
# import droplet
import copy
from DNAdroplet import DNADroplet
from utils import *
from glass import Glass
from deBruijnGraph import DeBruijnGraph
from DNAHandler import DNAHandler
from DNAfountain import DNAFountain
import re
import operator
import math
import glass

#20190701
def get_droplets(num, fdna1):
    i = 0
    droplets = []
    gc_drop_num = 0
    adrop = None
    print(i)
    print(num)
    while i < num:
        adrop = fdna1.DNAdroplet()
        if check_dna(adrop.to_DNA_CRC_sIndex()):
            droplets.append(adrop)
            # print(str(i)+"     Good Droplet detected!")
            i = i + 1
        else:
            gc_drop_num = gc_drop_num + 1
        #     print("Bad Droplet detected!")
        #     print(i)
    print('GC drop num:',end='\t')
    print(gc_drop_num)
    return droplets

def get_droplets_check_repeat_kmer(num, fdna1, kmer_len=21):
    i = 0
    droplets = []
    #kmer_lenth = kmer_len - 1
    kms = {}
    gc_drop_num = 0
    km_drop_num = 0
    # adrop = None
    j = 0
    while i < num:
        j = j + 1
        adrop = fdna1.DNAdroplet()
        if check_dna(adrop.to_DNA_CRC()):
            # if deGraph.highst_km_freq(adrop.to_DNA_CRC()) < 5:
            dps_kms = kmers_of_str(adrop.to_DNA_CRC_sIndex(), kmer_len-1)
            # kmers_in_dict()
            if not any_kmers_in_dict(dps_kms, kms):
                droplets.append(adrop)
                i = i + 1
                # Adding droplet kmers into kms
                for km in dps_kms:
                    kms[km] = 1
            else:
                km_drop_num = km_drop_num + 1
        else:
            gc_drop_num = gc_drop_num + 1

    #     if j > 10000:
    #         j = 0
    #         print(gc_drop_num, end='\t')
    #         print(km_drop_num, end='\t')
    #         print(num)
    # print('GC drop num:', end='\t')
    # print(gc_drop_num, end='\t')
    # print('Km drop num:', end='\t')
    # print(km_drop_num, end='\t')
    # print(num)
    return droplets



def get_droplets_check_repeat_kmer_multi_ft(num, fdna1, deGraph):

    kmer_len = deGraph.kmer_len
    i = 0
    droplets = []
    gc_drop_num = 0
    km_drop_num = 0
    # adrop = None
    j = 0
    while i < num:
        j = j + 1
        adrop = fdna1.DNAdroplet()
        if check_dna(adrop.to_DNA_CRC()):

            dps_kms = kmers_of_str(adrop.to_DNA_CRC_sIndex(), kmer_len - 1)

            if not any_kmers_in_dict(dps_kms, deGraph.kmers):
                droplets.append(adrop)
                i = i + 1
                for km in dps_kms:
                    deGraph.kmers[km] = 1
            else:
                km_drop_num = km_drop_num + 1
        else:
            gc_drop_num = gc_drop_num + 1

        if j > 10000:
            j = 0
            print(len(droplets), end='\t')
            print(gc_drop_num, end='\t')
            print(km_drop_num, end='\t')
            print(num)
    print(len(droplets), end='\t')
    print(gc_drop_num, end='\t')
    print(km_drop_num, end='\t')
    print(num)
    return droplets


def get_droplets_check_repeat_kmer_4deG(num, fdna1, deGraph, max_kmer_repeat_num=5):

    kmer_len = deGraph.kmer_len
    i = 0
    droplets = []
    gc_drop_num = 0
    km_drop_num = 0
    # adrop = None
    j = 0
    while i < num:
        j = j + 1
        adrop = fdna1.DNAdroplet()
        if check_dna(adrop.to_DNA_CRC()):
            if deGraph.highst_km_freq(adrop.to_DNA_CRC()) <= max_kmer_repeat_num:
                droplets.append(adrop)
                i = i + 1
                deGraph.add_seq(adrop.to_DNA_CRC())
            else:
                km_drop_num = km_drop_num + 1
        else:
            gc_drop_num = gc_drop_num + 1

        if j > 10000:
            j = 0
            print(gc_drop_num, end='\t')
            print(km_drop_num, end='\t')
            print(num)
    print(gc_drop_num, end='\t')
    print(km_drop_num, end='\t')
    print(num)
    return droplets

def add_droplets_to_deG_dist(num, fdna1, deGraph, min_kmer_dist=1):
    # kmer_len = deGraph.kmer_len
    i = 0
    droplets = []
    gc_drop_num = 0
    km_drop_num = 0
    # adrop = None
    j = 0
    while i < num:
        j = j + 1
        adrop = fdna1.DNAdroplet()
        if check_dna(adrop.to_DNA_CRC()):
            if kmer_editing_dist_seq_deG(adrop.to_DNA_CRC(), deGraph) >= min_kmer_dist:
                droplets.append(adrop)
                i = i + 1
                deGraph.add_seq(adrop.to_DNA_CRC())
            else:
                km_drop_num = km_drop_num + 1
        else:
            gc_drop_num = gc_drop_num + 1

        if j > 10000:
            j = 0
            print(gc_drop_num, end='\t')
            print(km_drop_num, end='\t')
            print(num)
    print(gc_drop_num, end='\t')
    print(km_drop_num, end='\t')
    print(num)
    return droplets


def dropout_rate(num, fdna1, min_gc=0.45, max_gc=0.55, max_homo=5, max_repeat_kmer=5):
    i = 0
    # droplets = []
    # deGraph = DeBruijnGraph()
    drop_result = {}
    drop_result['gc'] = 0
    drop_result['homo'] = 0

    while i < num:
        adrop = fdna1.DNAdroplet()
        dnastr = adrop.to_DNA_CRC()
        gc_rate = calc_gc(dnastr)
        if gc_rate > max_gc:
            drop_result['gc'] = drop_result['gc'] + 1
        if gc_rate < min_gc:
            drop_result['gc'] = drop_result['gc'] + 1
        homo_poly_len = max_homo_len(dnastr)
        if homo_poly_len > max_homo:
            drop_result['homo'] = drop_result['homo'] + 1
#         if check_dna(adrop.to_DNA_CRC()):
#             # print(str(i)+"     Good Droplet detected!")
# #           print(deGraph.highstKmFreq(adrop.to_DNA_CRC()))
#             droplets.append(adrop)
#             # deGraph.addSeq(adrop.to_DNA_CRC())
#             #print("Good Droplet detected!")
#
#         else:
#             drop_num = drop_num + 1
        i = i + 1
    return drop_result

def get_degree_droplet(degree, fdna1):
    adrop = fdna1.DNAdroplet()
    while adrop.degree != degree:
        adrop = fdna1.DNAdroplet()

    return adrop

def test_glass(gt, fdna1):
    i = 0
    while i < gt.num_chunks - 1:
        if gt.chunks[i] != fdna1.chunk(i):
            # print(i)
            return False
        i = i + 1
    gap = len(gt.chunks[i]) - len(fdna1.chunk(i))
    if gt.chunks[i] != fdna1.chunk(i) + b'\x00' * gap:
        # print(i)
        return False
    return True


def test_droplets(droplets, fdna1):
    gt = glass.Glass(fdna1.num_of_chunks)
    i = 0
    for droplet in droplets:
        # print(str(i) + "  adding droplets")
        # print(droplet.degree)
        gt.addDroplet(droplet)
        i = i + 1
        #print(i)
    gt.decode()
    i = 0
    while i < gt.num_chunks-1:
        if gt.chunks[i] != fdna1.chunk(i):
            #print(i)
            return False
        i = i + 1
    gap = len(gt.chunks[i]) - len(fdna1.chunk(i))
    if gt.chunks[i] != fdna1.chunk(i) + b'\x00' * gap:
        #print(i)
        return False
    return True

def test_DNA_collection(DNAs, fdna1, kmer_len=15 , min_index=1, max_index=100000):

    degree_table = get_degrees(fdna1.num_of_chunks, int(fdna1.num_of_chunks * 10), fdna1.seed)
    deG = DeBruijnGraph()
    deG.kmer_len = kmer_len
    deG.add_seqs(DNAs)
    test_results = {}

    print("Recovering data...")
    i = min_index
    find_path_num = 0
    find_multi_path_num = 0
    full_recover = False
    aGlass = Glass(fdna1.num_of_chunks)
    recov_droplets = []
    while (i <= max_index):
        # print("Finding index " + str(i) + " path")
        deG.find_droplet_DNA(i, 19)
        if (len(deG.crc_paths) > 0):
            find_path_num = find_path_num + 1
            # print("findpathnum  " + str(find_path_num))
            if (len(deG.crc_paths) > 1):
                find_multi_path_num = find_multi_path_num + 1
            else:
                adrop = DNADroplet()
                adrop.num_of_chunks = fdna1.num_of_chunks
                adrop.degree = degree_table[i - 1]
                adrop.set_droplet_from_DNA_CRC(deG.crc_paths[0])
                recov_droplets.append(copy.deepcopy(adrop))
                aGlass.addDroplet(adrop)
        print("Searched index: " + str(i) + " Found Droplets: " + str(find_path_num) + " Found multiple Paths: " + str(find_multi_path_num))
        i = i + 1
    if (test_glass(aGlass, fdna1)):
        full_recover = True

    test_results['Droplet_found'] = find_path_num
    test_results['Multi_path_found'] = find_multi_path_num
    test_results['Full_recover'] = full_recover

    return test_results


def recover_file(droplets, num_of_chunks,output_file,fdna1):
    OUT = open(output_file, 'wb')
    error_chunk_num = 0
    gt = Glass(num_of_chunks)
    i = 0
    for droplet in droplets:
        gt.addDroplet(droplet)
        i = i + 1
        #print(i)
    i = 0
    while i < gt.num_chunks-1:
        OUT.write(gt.chunks[i])
        if gt.chunks[i] != fdna1.chunk(i):
            #print(i)
            error_chunk_num = error_chunk_num + 1
        i = i + 1

    gap = len(gt.chunks[i]) - len(fdna1.chunk(i))
    OUT.write(gt.chunks[i])
    if gt.chunks[i] != fdna1.chunk(i) + b'\x00' * gap:
        error_chunk_num = error_chunk_num + 1

    return error_chunk_num


def get_droplets_noCheck(num, fdna1):
    i = 0
    droplets = []
    adrop = None
    while i < num:
        adrop = fdna1.DNAdroplet()
        #if check_dna(adrop.to_DNA_CRC()):
        droplets.append(adrop)
            #print("Good Droplet detected!")
        i = i + 1
        #else:
            #print("Bad Droplet detected!")
         #print(i)
    return droplets

def count_length_DNAs(dnas, max_len = 162):
    lengths = {}
    for j in range(0, max_len+1):
        lengths[j] = 0
    for dna in dnas:
        lengths[len(dna)] = lengths[len(dna)] + 1
    return lengths

def count_avg_length_DNAs(dnas):
    num_of_DNAs = len(dnas)
    length_all_DNAs = 0
    for dna in dnas:
        length_all_DNAs = length_all_DNAs + len(dna)
    return int(length_all_DNAs/num_of_DNAs)


# 2020/06/30
def hashToFile(hs, file):
    out = open(file, 'tw')
    for a in hs:
        out.write(str(a))
        out.write("\t")
        out.write(str(hs[a]))
        out.write("\n")
    out.close()


# Max value of hash
def max_value_hash(hs):
    max = 0
    for aKey in hs:
        if hs[aKey] > max:
            max = hs[aKey]
    return max

def poisson_seq_num(ranmd, rep_f_exp, rep_times=20, size=1):
     init_seq_nums = np.random.poisson(ranmd)
     rep_folds = math.pow((float(1) + poisson_rep_factor(rep_f_exp)), rep_times)
     return int(init_seq_nums * rep_folds)

def poisson_rep_factor(exp):
    ranmd = int(exp*10)
    rep_f = float(np.random.poisson(ranmd)/10)
    if rep_f > 1:
        rep_f = 1
    return rep_f



# 2020/06/27
def sta_value_hash(hs,  file, ladder=10):
    max = max_value_hash(hs)
    hs_sta = {}

    i = 0
    while (i <= int(max / ladder + 1)):
        hs_sta[i * ladder] = 0
        i = i + 1

    for aKey in hs:
        range_value = int(hs[aKey] / ladder) * ladder
        hs_sta[range_value] += 1

    hashToFile(hs_sta, file)


def test_deG(file, deG, byte_len =60):
    seqinfo = file_to_array(file)

    fd_ids = []
    found_num = 0
    corrupt_num = 0
    multi_num = 0
    wrong_num = 0
    i = 0
    for seq in seqinfo:
        i = i + 1
        arr = seq.split('\t')
        # print('.', end='')
        deG.find_droplet_DNA(int(arr[0]), byte_len)
        if len(deG.crc_paths) > 1:
            multi_num = multi_num + 1
        else:
            if len(deG.crc_paths) > 0:
                if deG.crc_paths[0] in arr[1]:
                    found_num = found_num + 1
                    print(arr[0])
                    fd_ids.append(arr[0])
                else:
                    wrong_num = wrong_num + 1

            else:
                corrupt_num = corrupt_num + 1

        # if i >= 1000:
        #     print(found_num,end='\t')
        #     print(multi_num, end= '\t')
        #     print(corrupt_num)
        #     i = 0
    return fd_ids
    # return found_num



def test_deG1(dnas, deG, byte_len =35):
    # seqinfo = dnas
    found_num = 0
    corrupt_num = 0
    multi_num = 0
    wrong_num = 0
    i = 0
    for anum in dnas:
        i = i + 1
        deG.find_droplet_DNA(anum, byte_len)
        if len(deG.crc_paths) > 1:
            multi_num = multi_num + 1
        else:
            if len(deG.crc_paths) > 0:
                if deG.crc_paths[0] in dnas[anum]:
                    found_num = found_num + 1
                else:
                    wrong_num = wrong_num + 1

            else:
                corrupt_num = corrupt_num + 1

        if i >= 1000:
            print(found_num,end='\t')
            print(multi_num, end= '\t')
            print(corrupt_num)
            i = 0
    return found_num


def test_dump_file(dump_file,seqinfo_file,kmer_len = 21, byte_len=17):
    deG = DeBruijnGraph()
    deG.kmer_len = kmer_len
    deG.veri_kmers = False
    deG.open_dump(dump_file)
    return test_deG(seqinfo_file,deG, byte_len)


def test_strand_dropouts(file, seq_file,long_kmer_len=81):
    # 2020-05-16
    deG = DeBruijnGraph()
    deG.kmer_len = long_kmer_len
    deG.veri_kmers = False
    deG.open_fasta(seq_file)

    seqinfo = file_to_array(file)
    found_num = 0
    for seq in seqinfo:
        arr = seq.split('\t')
        dnaseq = arr[2]
        key_kmers = key_kmers_of_strand(dnaseq,long_kmer_len)
        if kmers_in_dict(key_kmers, deG.kmers):
            found_num = found_num + 1
    print(found_num)
    return deG


def key_kmers_of_strand(dnastr, kmer_len, overlap=20):
    kmers = []
    str_len=len(dnastr)
    kmers.append(kmers_of_position(dnastr, kmer_len, 0))
    kmers.append(kmers_of_position(dnastr,kmer_len,overlap))
    kmers.append(kmers_of_position(dnastr, kmer_len, str_len-overlap-kmer_len))
    kmers.append(kmers_of_position(dnastr, kmer_len, str_len - kmer_len))
    return kmers


def hashToFile(hs, file):
# 2020/06/30
    out = open(file,'tw')
    for a in hs:
        out.write(str(a))
        out.write("\t")
        out.write(str(hs[a]))
        out.write("\n")
    out.close()


def hashToFile2(hs, file):
# 2020/06/30
    out = open(file,'tw')
    for a in hs:
        for b in hs[a]:
            out.write(str(a))
            out.write("\t")
            out.write(str(b))
            out.write("\t")
            out.write(str(hs[a][b]))
            out.write("\n")
    out.close()

def max_value_hash(hs):
    max = 0
    for aKey in hs:
        if hs[aKey] > max:
            max = hs[aKey]
    return max


def sta_value_hash(hs,  file, ladder=10):
    # 2020/06/27
    max = max_value_hash(hs)
    hs_sta = {}

    i = 0
    while(i <= int(max/ladder + 1)):
        hs_sta[i*ladder] = 0
        i = i + 1

    for aKey in hs:
        range_value = int(hs[aKey]/ladder) * ladder
        hs_sta[range_value] +=1

    hashToFile(hs_sta,file)
    return hs_sta


def hasInitCodes(dnastr):
    if "ATG" in dnastr:
        return True
    if "CAT" in dnastr:
        return True
    if "GTG" in dnastr:
        return True
    if "CAC" in dnastr:
        return True
    return False


def strand_recover_rate(droplets, DNAs, kmer_length=21, cov_cutoff=1):
    # 20200829 First Implementation
    # 20200909 Adding theoritical decoding num
    # droplet_all_copy = copy.deepcopy(droplets)
    droplet_all = droplets
    data_block_length = len(droplets[0].data)

    print("Analyzing droplets...")
    droplet_IDs = []
    droplet_ID_DNA = {}
    maxID = 1

    for dps in droplet_all:
        droplet_ID_DNA[dps.head_index] = dps.to_DNA_CRC()
        droplet_IDs.append(dps.head_index)
        if dps.head_index > maxID:
            maxID = dps.head_index
        # deG.addSeq(dps.to_DNA_CRC())

    deG = DeBruijnGraph()
    deG.kmer_len = kmer_length

    print("Adding DNA seqs to deG object")
    deG.add_seqs(DNAs)


    # i = min_index
    theory_num = 0
    find_path_num = 0
    fail_path_num = 0
    multi_pathAB_num = 0
    multi_path_crc_num = 0
    multi_path_IDs = []
    fail_droplet_IDs = []

    deGnoCut = copy.deepcopy(deG)


    deG.remove_low_cov_kmers(cov_cutoff)
    print("Recovering paths...")
    for id in droplet_IDs:
        print("Finding", end="\t")
        print(id, end="\t")
        print(data_block_length)
        deG.find_droplet_DNA(id, data_block_length)
        # print("Path nums:", end="\t")
        # print(len(deG.crc_paths), end="\t")
        # print(len(deG.pathAB), end="\t")
        if (len(deG.crc_paths) > 0):
            if (len(deG.crc_paths) > 1):
                multi_path_crc_num = multi_path_crc_num + 1
                multi_pathAB_num = multi_pathAB_num + 1
                multi_path_IDs.append(id)
            else:
                if (deG.crc_paths[0] == droplet_ID_DNA[id]):
                    find_path_num = find_path_num + 1
                    if len(deG.pathAB) > 1:
                        multi_pathAB_num = multi_pathAB_num + 1
                else:
                    fail_path_num = fail_path_num + 1
                    fail_droplet_IDs.append(id)
        else:
            fail_path_num = fail_path_num + 1
            fail_droplet_IDs.append(id)

        kmers = kmers_of_str(droplet_ID_DNA[id], kmer_length)
        if kmers_in_dict(kmers, deGnoCut.kmers):
            theory_num += 1

    return [theory_num, find_path_num, multi_pathAB_num, multi_path_crc_num]


def strand_theory_recover_rate(droplets, DNAs, kmer_length=21, cov_cutoff=0):
    # droplet_all_copy = copy.deepcopy(droplets)
    droplet_all = droplets
    data_block_length = len(droplets[0].data)

    #print("Analyzing droplets...")
    droplet_IDs = []
    droplet_ID_DNA = {}
    maxID = 1

    for dps in droplet_all:
        droplet_ID_DNA[dps.head_index] = dps.to_DNA_CRC()
        droplet_IDs.append(dps.head_index)
        if dps.head_index > maxID:
            maxID = dps.head_index
        # deG.addSeq(dps.to_DNA_CRC())

    deG = DeBruijnGraph()
    deG.kmer_len = kmer_length

    #print("Adding DNA seqs to deG object")
    deG.add_seqs(DNAs)


    # i = min_index
    theory_num = 0
    find_path_num = 0
    fail_path_num = 0
    multi_pathAB_num = 0
    multi_path_crc_num = 0
    multi_path_IDs = []
    fail_droplet_IDs = []

    deGnoCut = copy.deepcopy(deG)


    # deG.remove_low_cov_kmers(cov_cutoff)
    #print("Testing integrity of paths...")
    for id in droplet_IDs:
        #print(r'.', end='')
        # print("Testing integrity of ", end="\t")
        # print(id, end="\t")
        # print(data_block_length)

        kmers = kmers_of_str(droplet_ID_DNA[id], kmer_length)
        if kmers_in_dict(kmers, deGnoCut.kmers):
            theory_num += 1

    return theory_num


def strand_in_graph(strand, deG):
    kmer_length = deG.kmer_len
    kmers = kmers_of_str(strand, kmer_length)
    if kmers_in_dict(kmers, deG.kmers):
        return True
    else:
        return False


def strand_num_in_graph(strands, deG):
    a = 0
    for strd in strands:
        if strand_in_graph(strd, deG):
            a = a + 1
    return a



###-----------------------------------------------
###Some Old Functions, maybe removed later--------------------------
#An old function 20200829


def recover_num(droplets, error_type="sub", error_rate=0.01, seq_num=10, seq_rand=False, kmer_length=21):

    droplet_all_copy = copy.deepcopy(droplets)
    data_block_length = len(droplets[0].data)
    dna_handler = DNAHandler()

    primerF = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
    primerE = 'CTGACACTGATGCATCCG'  # complement seq of P2

    if seq_rand:
        min_DNA_copies = 0
        max_DNA_copies = seq_num * 2
    else:
        min_DNA_copies = seq_num
        max_DNA_copies = seq_num

    deG = DeBruijnGraph()
    deG.kmer_len = kmer_length

    all_DNAs = []
    droplet_IDs = []
    droplet_ID_DNA = {}
    maxID = 1
    print("Copying DNAs with Errors...")
    for dps in droplet_all_copy:
        deG.add_seqs(
            dna_handler.copyDNAsWithSubErr([primerF + dps.to_DNA_CRC() + primerE], min_DNA_copies, max_DNA_copies,
                                           error_rate))
        # deG.addSeq(primerF + dps.to_DNA_CRC() + primerE)
        droplet_ID_DNA[dps.head_index] = dps.to_DNA_CRC()
        all_DNAs.append(dps.to_DNA_CRC())
        droplet_IDs.append(dps.head_index)
        if (dps.head_index > maxID):
            maxID = dps.head_index
        # deG.addSeq(dps.to_DNA_CRC())

    print("Recovering paths...")
    # i = min_index
    find_path_num = 0
    fail_path_num = 0
    find_multi_path_num = 0
    multi_path_IDs = []

    fail_droplet_IDs = []
    for id in droplet_IDs:
        deG.find_droplet_DNA(id, data_block_length)
        if (len(deG.crc_paths) > 0):
            if (len(deG.crc_paths) > 1):
                find_multi_path_num = find_multi_path_num + 1
                multi_path_IDs.append(id)
            else:
                if (deG.crc_paths[0] == droplet_ID_DNA[id]):
                    find_path_num = find_path_num + 1
                else:
                    fail_path_num = fail_path_num + 1
                    fail_droplet_IDs.append(id)
        else:
            fail_path_num = fail_path_num + 1
            fail_droplet_IDs.append(id)

    return find_path_num


def test_deGraphsize(data_bytes, fountain_seed, min_DNA_copies, max_DNA_copies, err_rates, number_of_droplets=10000, fountain_init_index=1, data_block_length=30, kmer_length=21):
    # test DeGraphSize
    primerF = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
    primerE = 'CTGACACTGATGCATCCG'  # complement seq of P2

    fdna1 = DNAFountain(data_bytes, data_block_length, fountain_init_index, fountain_seed)
    print("Generating droplets!")
    droplet_all = get_droplets_CheckKmer(number_of_droplets, fdna1)
    all_DNAs = []
    for dps in droplet_all:
        all_DNAs.append(primerF + dps.to_DNA_CRC() + primerE)

    deG = DeBruijnGraph()
    deG.kmer_len = kmer_length
    deG.add_seqs(all_DNAs)

    deGE = DeBruijnGraph()

    results = {}
    print("Testing SubErrors....")
    results['sub'] = {}
    results['ins'] = {}
    results['del'] = {}

    dna_handler = DNAHandler()
    for err in err_rates:
        deGE = DeBruijnGraph()
        deGE.kmer_len = kmer_length
        deGE.add_seqs(dna_handler.copy_randnum_with_sub(all_DNAs, min_DNA_copies, max_DNA_copies, err))

        results['sub'][str(err)] = {}
        results['sub'][str(err)]['Ori'] = len(deG.kmers)
        results['sub'][str(err)]['Err'] = len(deGE.kmers)
        # results['sub'][str(err)]['Veri'] = len(deGE.kmerList)


def kmer_num(seqs, kmer_length=21):
    deG = DeBruijnGraph()
    deG.kmer_len = kmer_length
    deG.add_seqs(seqs)
    return len(deG.kmers)

def kdkn(seq_seqs, ori_seqs, kmer_length=21):
    deGori = DeBruijnGraph()
    deGori.kmer_len = kmer_length
    deGori.add_seqs(ori_seqs)

    deGseq = DeBruijnGraph()
    deGseq.kmer_len = kmer_length
    deGseq.add_seqs(seq_seqs)

    corr_kmers = 0
    for km in deGori.kmers:
        if km in deGseq.kmers:
            corr_kmers = corr_kmers + 1

    kd = (len(deGori.kmers) - corr_kmers) / len(deGori.kmers)
    kn = (len(deGseq.kmers) - corr_kmers) / corr_kmers

    sr = 0
    for seq in ori_seqs:
        kmers = kmers_of_str(seq, kmer_length)
        if kmers_in_dict(kmers, deGseq.kmers):
            sr += 1

    return [kn, kd, sr]

def knkd_analysis(rs, cut_off=0):

    for cov in [1, 2, 4,  6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 64, 128, 256]:
        print(str(cov) + '\t', end='')
        for item in ['kd', 'kn', 'found']:
            # print(item + '\t')
            aa = item_std_dec(rs, cov, item, cut_off)
            print(str(aa[0]) + '\t', end='')
            print(str(aa[1]) + '\t', end='')
            print(str(aa[2]) + '\t', end='')

        print('\n')



def item_values(rs, cov, it, cut_off=0):
    vals = []
    for a in rs:
        vals.append(a[cov][1][cut_off][it])
    return vals

def item_std_dec(rs, cov, it, cut_off=0):
    vals = []
    cut = []
    for a in rs:
        if cut_off == -1:
            best_cut = len(a[cov][1]) - 2
            vals.append(a[cov][1][best_cut][it])
        else:
            best_cut = cut_off
            vals.append(a[cov][1][cut_off][it])
        cut.append(best_cut)

    return [np.mean(vals), np.std(vals,ddof=1), np.mean(cut)]

#2020-12-09
def kmer_dis_simulation(f_seed, data, number_of_droplets = 10000,exp_seq_num=50,  err_rate = 0.04, kmer_length=21):
    # file1 = open('2.pdf', 'rb')
    fountain_seed = f_seed
    run_name = 'seed' + str(fountain_seed) + '_'
    fountain_init_index = 1

    fdna1 = DNAFountain(data, 30, 1, fountain_seed)
    dna_handler = DNAHandler()
    degree_table = get_degrees(fdna1.num_of_chunks, int(fdna1.num_of_chunks * 5), fountain_seed)

    print("Generating droplets!")
    # dropout_rate = dropout_rate(number_of_droplets,fdna1)
    droplet_all = get_droplets_check_repeat_kmer(number_of_droplets, fdna1)

    deG = DeBruijnGraph()
    deG.kmer_len = kmer_length

    all_DNAs = []
    droplet_IDs = []
    droplet_ID_DNA = {}
    maxID = 1
    print("Degradating DNAs...")
    for dps in droplet_all:
        droplet_ID_DNA[dps.head_index] = dps.to_DNA_CRC()
        all_DNAs.append(dps.to_DNA_CRC())
        droplet_IDs.append(dps.head_index)
        if (dps.head_index > maxID):
            maxID = dps.head_index
        # deG.addSeq(dps.to_DNA_CRC())

    deG.add_seqs(all_DNAs)

    dna_hd = DNAHandler()
    droplet_DNAs = all_DNAs

    deGori = DeBruijnGraph()
    deGori.kmer_len = kmer_length
    deGori.kmers = deG.kmers

    print('Generating 50 copies error-rich seqs')
    eDNAs = dna_hd.copy_seq_poisson(droplet_DNAs, exp_seq_num)
    eDNAs = dna_hd.add_rand_sub_new(eDNAs, err_rate / 2)
    eDNAs = dna_hd.add_rand_ins_new(eDNAs, err_rate / 4)
    eDNAs = dna_hd.add_rand_del_new(eDNAs, err_rate / 4)

    deGE = DeBruijnGraph()
    deGE.kmer_len = kmer_length

    print('Adding 10 copies error-rich seqs to deG')
    deGE.add_seqs(eDNAs)

    for a in deGori.kmers:
        if a in deGE.kmers:
            deGori.kmers[a] = deGE.kmers[a]
        else:
            deGori.kmers[a] = 0
    res = []
    res.append(sta_value_hash(deGE.kmers, run_name + '.10', 1))
    res.append(sta_value_hash(deGori.kmers, run_name + '.10.ori', 1))
    return res


def majority_merge(reads, weight = 0.4):
    # Function used in grass's study
    # assume reads have the same length
    res = ""
    reads_upper = []
    for rd in reads:
        reads_upper.append(rd.upper())

    for i in range(len(reads_upper[0])):
        #aseq = str.upper()
        counts = {'A':0,'C':0,'G':0,'T':0,'-':0,'N':0}
        for j in range(len(reads_upper)):
            counts[reads_upper[j][i]] +=1
        counts['-'] *= weight
        mv = max(counts.items(), key=operator.itemgetter(1))[0]
        if mv != '-':
            res += mv
    return res

    # 20201218 read aln file
def read_fasta(file):
    f = open(file, "r")
    seqs = []
    matchLineA = re.compile('^(>)')

    line = f.readline()

    seq = ''
    while line.strip():
        if not matchLineA.match(line):
            seq = seq + line.strip()
        else:
            if seq != '':
                seqs.append(seq)
                seq = ''
        line = f.readline()
    f.close()

    if seq != '':
        seqs.append(seq)

    return seqs


def sim_file_seqs(file):

    f = open(file)
    arr = []
    line = f.readline()
    while line.strip():
        a = line.strip()
        arr.append(a.split('\t')[1])
        line = f.readline()

    return arr

def read_sim(file):

    f = open(file)
    hs = {}
    line = f.readline()
    while line.strip():
        a = line.strip()
        hs[int(a.split('\t')[0])] = a.split('\t')[1]
        line = f.readline()
    return hs

def hash_to_sim(data, file):
    f = open(file, 'tw')
    for id in data:
        f.write(str(id))
        f.write('\t')
        f.write(data[id])
        f.write("\n")
    f.close()



def read_starcode_clusters(file, min_seq_num):

    f = open(file)
    clu_seqs = {}
    line = f.readline()
    clu_num = 1
    while line.strip():
        a = line.strip()
        arr = a.split('\t')

        arr2 = arr[2].split(',')
        if len(arr2) >= min_seq_num:
            clu_seqs[clu_num] = arr2
        line = f.readline()
        clu_num = clu_num + 1
    return clu_seqs

def seq_in_seqs(seq, seqs):
    for s in seqs:
        if seq in s:
            return True
    return False

def cut_num(cov, p1=2, p2=25):
    if cov <= p1:
        return 0
    if cov <= p2:
        return 1
    cut = math.log2(cov/p2) + 2
    return int(cut)


def cut_num_old(cov, p1=0.1, p2=2.5):
    cut = cov * p1 / p2 ** (math.log10(cov) - 1)
    return int(cut)

def hash_keys_values(hs):
    for a in hs:
        print(a,end="\t")
        print(hs[a])


def sub_seqs(seqs, ff=0, ee=-1):
    if ee < 0:
        ee = len(seqs[0]) -1
    arr = []
    for a in seqs:
        arr.append(a[ff:ee])
    return arr

def kmer_editing_dist(k1, k2):
    assert len(k1) == len(k2), "length of kmer 1 and 2 are not same"
    dis = 0
    for i in range(0,len(k1)):
        if k1[i] != k2[i]:
            dis = dis + 1
    return dis

def kmer_editing_dist_seq_deG(seq, deG):
    k_len = len(deG.kmer_len)
    min_dist = k_len
    for km1 in kmers_of_str(seq, k_len):
        for km2 in deG.kmers:
            dist = editing_dist_kmer(km1, km2)
            if dist < min_dist:
                min_dist = dist
    return min_dist


def chunks_data(ft, chunk_nums):
    data = None
    for num in chunk_nums:
        if data is None:
            data = ft.chunk(num)
        else:
            data = xor(data, ft.chunk(num))
    return data


def filter_con_seqs(con_seqs, dps_seqs, p1 = 50, p2 = 150):
    new_con_seqs = []
    deG = DeBruijnGraph()
    for dp in dps_seqs:
        deG.add_seq(dps_seqs[dp])

    for seq in con_seqs:
        if strand_in_graph(seq[p1:p2], deG):
            new_con_seqs.append(seq)

    return new_con_seqs

def check_cons(con_seqs, dps_seqs):
    kmer_clu_seqs = {}

    print('Building kmer clu_con_seq rels')
    for i, seq in enumerate(con_seqs):
        kms = kmers_of_str(seq, 16)
        for km in kms:
            if km in kmer_clu_seqs:
                kmer_clu_seqs[km].append(i)
            else:
                kmer_clu_seqs[km] = [i]

    ss = 0

    print('Checking dps one by one')
    for dp in dps_seqs:
        #print('.', end='')
        id_km = dps_seqs[dp][0:16]
        if id_km in kmer_clu_seqs:
            for dgs in kmer_clu_seqs[id_km]:
                if len(con_seqs[dgs]) == 200:
                    if dps_seqs[dp] in con_seqs[dgs]:
                        ss = ss + 1
                        #print('\^')
                        break
    return ss

def seqs_to_fasta(seqs, file):
    f = open(file, 'tw')
    i = 1
    for seq in seqs:
        f.write(">" + str(i) + "\n")
        #f.write('\t')
        f.write(seq)
        f.write("\n")
        i = i + 1
    f.close()