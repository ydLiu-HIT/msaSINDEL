#!/usr/bin/env python
# coding=utf-8

import os
import sys
import gc
import logging
import subprocess
import shlex
import argparse
from intervaltree import IntervalTree
from multiprocessing import Pool
from collections import Counter, defaultdict
from time import time
from utils import *
from check import *
from getHCVariant import get_hc_variantList
from subCall import mergeCoVar, most_likely_genotype, getExactMatch, makeCandidate
from poa_module import POA
from ksw_module import ksw_aligner

GAP = 50000
LARGECUT = 50 #20
MAXR = 250

#majorContigs = {"chr"+str(a) for a in list(range(0,23))+["X", "Y"]}.union({str(a) for a in list(range(0,23))+["X", "Y"]})
majorContigs = {"chr"+str(a) for a in list(range(0,23))}.union({str(a) for a in list(range(0,23))})

def setup_logging():
    """
    Default logger
    """
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(level=logging.INFO, format=log_format)

def getHCV(fp, threads, perror_for_indel, perror_for_snp):
    logging.info("Getting HC variants")
    return get_hc_variantList(fp, threads, perror_for_indel, perror_for_snp, OUT=False)

def can_combine(svT1, svT2, S1, S2):
    if svT1 == 'S' and svT2 == 'S':
        return True
    if svT1 != 'S' and svT2 != 'S':
        return True
    if S2-S1 < 20:
        return True
    return False

def getCandidateFromFile(args):
    AllRen = {}
    fin = open(args.canfix+".lc.info", 'r')
    for line in fin:
        line = line.split()
        ctg = line[0]
        #if (args.chrName != None and ctg != args.chrName) or ctg not in majorContigs:
            #continue
        if ctg not in majorContigs: continue
        if ctg not in AllRen:
            AllRen[ctg] = []
        pos = int(line[1])
        if args.chrStart != None and pos < args.chrStart: continue
        if args.chrEnd != None and pos > args.chrEnd: continue
        #cal ave indel length rC(readCount)
        rC = int(line[2]); refB = line[3]; refC = int(line[4])
        l = 5
        while l < len(line):
            svT = line[l]; cnt = int(line[l+1]); l += 2
            if svT == "I" or svT == "D":
                mostAllele = Counter(line[l:l+cnt]).most_common(1)
                ave_l = sum([len(x) for x in line[l:l+cnt]]) / cnt
                #ave_l = sum([int(x) for x in line[l:l+cnt]]) / cnt
                AllRen[ctg].append((pos, svT, int(ave_l), rC-cnt+1, rC+1, mostAllele[0][0]))
                l += cnt
            else:
                AllRen[ctg].append((pos, 'S', 0, refC+1, rC+1, svT)) 

    fin.close()
    return AllRen

def mergeCandidate(AllRen, slideWindow):
    CanRen = {}
    TotalSite = 0
    for ctg in AllRen:
        CanRen[ctg] = list()
        Candi = AllRen[ctg]
        Start, svT = Candi[0][0:2]
        i = 1; j = 0; preS = Start
        while i < len(Candi):
            if Candi[i][0] - Start < slideWindow and can_combine(svT, Candi[i][1], preS, Candi[i][0]):
                preS = Candi[i][0]
                i += 1
            else:
               L = [m for m in Candi[j:i]]
               CanRen[ctg].append(L)
               Start, svT = Candi[i][0:2]; preS = Start
               TotalSite += 1
               j = i; i += 1

        # last
        L = [m for m in Candi[j:i]]
        CanRen[ctg].append(L)
        TotalSite += 1

    return CanRen, TotalSite

def get_reference_sequence(ref_path, chrName, chrStart, chrEnd):
    refernce_sequences = []

    samtools_faidx_process = subprocess.Popen(shlex.split("samtools faidx %s %s:%d-%d" % (ref_path, chrName, chrStart, chrEnd) ), stdout=subprocess.PIPE, bufsize=8388608, universal_newlines=True)  # 2308

    while True:
        row = samtools_faidx_process.stdout.readline()
        is_finish_reading_output = row == '' and samtools_faidx_process.poll() is not None
        if is_finish_reading_output:
            break
        if row:
            refernce_sequences.append(row.rstrip())

    # first line is reference name ">xxxx", need to be ignored
    reference_sequence = "".join(refernce_sequences[1:])

    # uppercase for masked sequences
    reference_sequence = reference_sequence.upper()

    samtools_faidx_process.stdout.close()
    samtools_faidx_process.wait()
    if samtools_faidx_process.returncode != 0:
        return None

    return reference_sequence

def main_ctrl(args):
    CheckFileExist(args.fin_bam)
    CheckFileExist(args.fin_ref)
    CheckFileExist(args.canfix+".lc.info")
    
    AllRen = getCandidateFromFile(args)
    CanRen, TotalSite = mergeCandidate(AllRen, args.slideWindow) 
    del AllRen
    gc.collect()

    logging.info("Total %d region-sites need to be processed for indel detection" %(TotalSite))
    # get Ref
    logging.info("Loading reference genome")
    global REF_dict
    REF_dict = LoadRef(args.fin_ref) 
    if args.threads > 1:
        analysis_pools = Pool(processes=int(args.threads))
    logging.info("Run the script with %d threads" %(args.threads))
    semiResult = list(); finalResult = list()
    finalResult.extend(getHCV(args.canfix+".hc.info", args.threads, args.perror_for_indel, args.perror_for_snp))
    prefixDict = defaultdict(str)
    for ctg in CanRen:
        siteStart = 1; end = 0; start = 0; cnt = 0
        LL = len(CanRen[ctg])
        PE = CanRen[ctg][0][-1][0]
        while siteStart < LL:
            if CanRen[ctg][siteStart][0][0] - PE < GAP and cnt < args.ChunkSite:
                end += 1; cnt += 1
            else:
                if end > LL: end = LL - 1
                fp_prefix = "%s_%d_%d_" %(ctg, start, end)
                prefixDict[fp_prefix] = end - start
                start = siteStart; end = siteStart; cnt = 0

            PE = CanRen[ctg][siteStart][-1][0]
            siteStart += 1
        fp_prefix = "%s_%d_%d_" %(ctg, start, end)
        prefixDict[fp_prefix] = end - start
 
    fplist = sorted(prefixDict.items(), key = lambda x: x[1], reverse=True)
    logging.info("Getting LC variants")
    for item in fplist:
        fp_prefix = item[0]; ctg, start, end = fp_prefix.split('_')[0:3]
        #logging.info("current fp_prefix: %s" %(fp_prefix))
        if args.threads > 1:
            param = [(fp_prefix, CanRen[ctg][int(start):int(end) + 1], ctg, args)]
            semiResult.append(analysis_pools.map_async(run_parallel, param))
        else:
            param = (fp_prefix, CanRen[ctg][int(start):int(end) + 1], ctg, args)
            finalResult.extend(run_parallel(param))
    if args.threads > 1:
        analysis_pools.close()
        analysis_pools.join()
        for res in semiResult:
            try:
                finalResult.extend(res.get()[0])
            except:
                pass

    logging.info("Writing variant to disk")
    # writing header
    vcfOut = open(args.fout_vcf, 'w')
    PrintVCFHeader(vcfOut, args.fin_ref, "HG002")

    finalResult = sorted(finalResult, key=lambda var: (str(var.chr), int(var.pos)))
    Generate_Output(vcfOut, finalResult)

    vcfOut.close()

def getFlanking(candi):
    indL = 0
    for ci in candi:
        indL = ci[2] if indL < ci[2] else indL
    
    # 新方法,flanking的影响很小，但是会增加MSA和call的速度
    flanking = 100 
    if indL == 0: #snp
        flanking = 50
    else:
        flanking = 100
    #elif indL < 10:
    #    flanking = 100 ##
    #elif indL >= 10 and indL < 20:
    #    flanking = 150
    #elif indL >= 20 and indL < 30:
    #    flanking = 200
    #else:
    #    flanking = 300

    return flanking
    
def run_parallel(args):
    return run_core(*args)

def getLocalSeq(ref, rS, queries, candi, Start, End):
    query = list(); qual = list(); target = ''
    cnt = 0
    for l in queries:
        POS = l.begin
        SEQ, CIGAR, QUAL = l.data

        refStart = POS; readStart = 0
        ExtractStart = readStart; ExtractEnd = 0 
        flagS = False; flagE = True; inner = True
        if refStart < Start:
            flagS = True
        elif refStart-Start > LARGECUT:
            inner = False

        op_l = 0; lastop = 0; lastopl = 0
        for m in str(CIGAR):
            if m.isdigit():
                op_l = op_l * 10 + int(m)
                continue

            if m == "S":
                readStart += op_l
            elif m == "M" or m == "=" or m == "X": # M, X, =
                refStart += op_l
                readStart += op_l
            elif m == "I": # I
                readStart += op_l
            elif m == "D": # D
                refStart += op_l
            lastop = m; lastopl = op_l
            op_l = 0
            
            # cal start
            if flagS and refStart >= Start:
                if m == "M" or m == "=" or m == "X":
                    ExtractStart = readStart - (refStart - Start)
                else:
                    ExtractStart = readStart
                flagS = False
            # cal end
            if refStart >= End:
                if m == "M" or m == "=" or m == "X":
                    ExtractEnd = readStart - (refStart - End)
                else:
                    ExtractEnd = readStart
                flagE = False
                break
        if flagE:
            if End - refStart > LARGECUT:
                inner = False
            else:
                if lastop == "S" or lastop == "I":
                    ExtractEnd = readStart - lastopl
                else:
                    ExtractEnd = readStart

        if End - refStart > LARGECUT:
            ExtractEnd = readStart
            inner = False

        Lread = ExtractEnd - ExtractStart
        Lref = End - Start

        #print(readStart, ExtractStart, ExtractEnd, Lread, Lref)
        if Lread > Lref/2 and Lread < 3*Lref/2 and inner and cnt < MAXR:
            query.append(bytes(SEQ[ExtractStart: ExtractEnd], 'utf-8'))
            qual.append(QUAL[ExtractStart: ExtractEnd])
            cnt += 1
    
    name = "RefSeq_" + "id_" +str(Start) + "_" + str(End) + "_" + "_".join(["%d|%s|%d|%d|%s" %(x[0],x[1],x[3],x[4],x[5]) for x in candi])
    target = str(ref[Start+1-rS:End+1-rS])

    return query, qual, target, name

def ExtractSeqAroundVariant(chrStart, chrEnd, chrName, args):
    p2 = subprocess.Popen(shlex.split("samtools view -F 2316 -q %d %s %s:%d-%d" % (args.minMQ, args.fin_bam, chrName, chrStart, chrEnd) ), stdout=subprocess.PIPE, bufsize=8388608, universal_newlines=True)  # 2308
    ReadTree = IntervalTree()
    for l in p2.stdout:
        l = l.strip().split()
        if l[0][0] == "@":
            continue
        CIGAR = l[5]
        refStart = int(l[3]) - 1
        refEnd = refStart
        advance = 0
        for m in str(CIGAR):
            if m.isdigit():
                advance = advance * 10 + int(m)
                continue
            if m == "M" or m == "=" or m == "X" or m == "D":
                refEnd += advance
            advance = 0
        ReadTree.addi(refStart, refEnd, (l[9].upper(), CIGAR, l[10]))
    p2.stdout.close()
    p2.wait()

    return ReadTree

def run_core(fp_prefix, CanRen, chrName, args): 
    chrStart = CanRen[0][0][0] - 1000
    chrEnd = CanRen[-1][-1][0] + 1000
    ##be care! ref seq by samtools and REF_dict are different. samtools startPos = 0, REF_dict startPos = 1
    #reference_sequence = get_reference_sequence(args.fin_ref, chrName, chrStart, chrEnd)
    reference_sequence = str(REF_dict[chrName].seq[chrStart-1:chrEnd])
    ReadTree = ExtractSeqAroundVariant(chrStart, chrEnd, chrName, args)
    finalList = list()
    sid = 0
    for candi in CanRen:
        flanking = getFlanking(candi)
        Start = candi[0][0] - flanking
        End = candi[-1][0] + flanking
        IsID = "S" if flanking == 50 else "ID"

        queries = sorted(ReadTree.overlap(Start, End))
        if len(queries) < args.minCNT: continue
        
        query, qual, target, name = getLocalSeq(reference_sequence, chrStart, queries, candi, Start, End)
        if len(query) == 0: continue 

        # debug
        #if name != "RefSeq_id_21520822_21520922_21520872|S|1|6|G": continue
        #with open("test.fa", 'w') as fw:
        #    for t in range(len(query)):
        #        fw.write('>read_%s\n' %(str(t)))
        #        fw.write(bytes.decode(query[t]))
        #        fw.write('\n')

        consensus = POA(query, 2, 0.2 if IsID == "S" else 0.3)
        canList = ksw_core(consensus,target, chrName, name, reference_sequence, chrStart, args.shift)
        if len(canList) == 0: continue

        canList = mergeCoVar(canList, target, sid)
        
        most_likely_genotype(canList, query, qual)

        finalList.extend(canList)

        del query; del qual; del queries; del canList 

    del reference_sequence 
    del ReadTree
    #del queries
    #del query
    #del qual
    #del canList

    gc.collect()
    logging.info("[localMSA.py] Finish prefix: %s" %(fp_prefix))
    return finalList

def ksw_core(consensus, target, chrName, RNAME, reference_sequence, rS, shift):
    canList = []
    for i in range(len(consensus)):
        msa = bytes.decode(consensus[i])
        alignment = ksw_aligner(msa, target)
        alignment[1] = getExactMatch(alignment[1], target, msa, 0)
        if alignment[0] < 0: continue
        canList.extend(makeCandidate(chrName, alignment, RNAME.split('_'), reference_sequence, rS, msa, shift, "LC"))
    
    return canList


def run():
    parser = argparse.ArgumentParser(description="MSA examples using abPOA")
    parser.add_argument("--fin_bam", type=str, default="input.bam",
            help="Sorted bam file, default: %(default)s")
    parser.add_argument("--fin_ref", type=str, default="ref.fa",
            help="Reference fasta input")
    parser.add_argument("--fout_vcf", type=str, default="detect.vcf",
            help="vcf path for detected variants, default: %(default)s")
    parser.add_argument("--canfix", type=str, default="candidate",
            help="Candidate file prefix used in ExtractVariant.py, default: %(default)s")
    parser.add_argument("--chrName", type=str, default=None,
            help="The name of reference to be processed")
    parser.add_argument("--chrStart", type=int, default=None,
            help="The 1-based starting positions of the reference to be processed")
    parser.add_argument("--chrEnd", type=int, default=None,
            help="The inclusive ending positions of the reference to be processed")
    parser.add_argument("--perror_for_snp", type=float, default="0.1",
            help="P-error is the probability of observing a heterozygote SNP, default:%(default)f")
    parser.add_argument("--perror_for_indel", type=float, default="0.1",
            help="P-error is the probability of observing a heterozygote indel, default:%(default)f")
    parser.add_argument("--minMQ", type=int, default=10,
            help="Minimum Mapping Quality. Mapping quality lower than the setting will be filtered, default:%(default)d")
    parser.add_argument("--minCNT", type=int, default=3,
            help="Minimum read counts required to call a variant, default:%(default)d")
    parser.add_argument("--ChunkSite", type=int, default=50000,
            help="Divide job with smaller candidate site count for parallelism, default:%(default)d")
    parser.add_argument("--threads", type=int, default=1,
            help="Number of threads to use, default:%(default)d")
    parser.add_argument("--shift", type=int, default=5,
            help="The distance bewteen the detect variant and the activate region, default:%(default)d")
    parser.add_argument("--flanking", type=int, default=50,
            help="Flanking base pairs around variant site, default: %(default)d")
    parser.add_argument("--slideWindow", type=int, default=100,
            help="The window lenth for two candidate can be merge as one to speed up, default: %(default)d")

    args = parser.parse_args()    
    
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    setup_logging()
    main_ctrl(args)
    
if __name__ == "__main__":
    t_s = time()
    run()
    t_e = time()
    logging.info("Finish the script in %f seconds", t_e - t_s)
