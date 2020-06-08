#!/usr/bin/env python
# coding=utf-8

import os
import sys
import gc
import logging
import subprocess
import shlex
import re
import argparse
from intervaltree import IntervalTree
from multiprocessing import Pool
from collections import Counter
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from utils import LoadRef
from check import *

GAP = 50000

cigarRe = r"(\d+)([MIDNSHP=X])"
majorContigs = {"chr"+str(a) for a in list(range(0,23))+["X", "Y"]}.union({str(a) for a in list(range(0,23))+["X", "Y"]})

def setup_logging():
    """
    Default logger
    """
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(level=logging.INFO, format=log_format)

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
    fin = open(args.fin_can, 'r')
    for line in fin:
        line = line.split()
        ctg = line[0]
        if ctg not in majorContigs:
            continue
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

def main_ctrl(args):
    if not os.path.exists(args.workDir):
        os.makedirs(args.workDir)

    CheckFileExist(args.fin_bam)
    CheckFileExist(args.fin_ref)
    CheckFileExist(args.fin_can)
    
    AllRen = getCandidateFromFile(args)
    CanRen, TotalSite = mergeCandidate(AllRen, args.slideWindow) 
    
    del AllRen
    gc.collect()

    logging.info("Total %d site need to be processed for indel detection" %(TotalSite))
    # get Ref
    REF_dict = LoadRef(args.fin_ref)

    analysis_pools = Pool(processes=int(args.threads))
    logging.info("Run the script with %d threads" %(args.threads))

    prefixList = []
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
                prefixList.append(fp_prefix)
                param = [(fp_prefix, CanRen[ctg][start:end + 1], ctg, REF_dict[ctg].seq, start, args)]
                analysis_pools.map_async(run_parallel, param)
                #param = (fp_prefix, CanRen[ctg][start:end + 1], ctg, REF_dict[ctg].seq, start, args)
                #run_parallel(param)
                start = siteStart; end = siteStart; cnt = 0

            PE = CanRen[ctg][siteStart][-1][0]
            siteStart += 1

        fp_prefix = "%s_%d_%d_" %(ctg, start, end)
        prefixList.append(fp_prefix)
        param = [(fp_prefix, CanRen[ctg][start:end + 1], ctg, REF_dict[ctg].seq, start, args)]
        analysis_pools.map_async(run_parallel, param)
        #param = (fp_prefix, CanRen[ctg][start:end + 1], ctg, REF_dict[ctg].seq, start, args)
        #run_parallel(param)

    analysis_pools.close()
    analysis_pools.join()

    # write prefixList
    with open(args.workDir + "prefixList", 'w') as fw:
        for v in prefixList:
            fw.write(v + '\n')
    
    logging.info("Finishing extracting regional reads for MSA and BWA-mem")

def getFlanking(candi):
    indL = 0
    for ci in candi:
        indL = ci[2] if indL < ci[2] else indL
    
    # 新方法,flanking的影响很小，但是会增加MSA和call的速度
    flanking = 100 
    if indL == 0: #snp
        flanking = 50
    elif indL < 10:
        flanking = 100 ##
    elif indL >= 10 and indL < 20:
        flanking = 150
    elif indL >= 20 and indL < 30:
        flanking = 200
    else:
        flanking = 300

    return flanking
    
def run_parallel(args):
    return ExtractSeqAroundVariant(*args)

def ExtractSeqAroundVariant(fp_prefix, CanRen, chrName, ref, candiS, args):
    ###
    #   should adjust flanking bases according to indel length for speeding
    ###
    chrStart = CanRen[0][0][0] - 1000
    chrEnd = CanRen[-1][-1][0] + 1000
    p2 = subprocess.Popen(shlex.split("samtools view -F 2316 -q %d %s %s:%d-%d" % (args.minMQ, args.fin_bam, chrName, chrStart, chrEnd) ), stdout=subprocess.PIPE, bufsize=8388608, universal_newlines=True)  # 2308
    ReadTree = IntervalTree()
    for l in p2.stdout:
        l = l.strip().split()
        if l[0][0] == "@":
            continue
        CIGAR = l[5]
        refStart = int(l[3]) - 1
        refEnd = refStart
        skipBase = 0
        totalAlnPos = 0
        for m in re.finditer(cigarRe, CIGAR):
            advance = int(m.group(1))
            totalAlnPos += advance
            if m.group(2)  == "S":
                skipBase += advance
            elif m.group(2) in ("M", "=", "X", "D"):
                refEnd += advance
        #if 1.0 - float(skipBase) / (totalAlnPos + 1) < 0.55:
        #    continue
        ReadTree.addi(refStart, refEnd, (l[9].upper(), CIGAR, l[10]))
    p2.stdout.close()
    p2.wait()
    
    ExtractRead = {"ID":[], "S": []}
    RefSeq = {"ID":[], "S": []}
    ReadCount = {"ID":[0], "S": [0]}
    SumCount = {"ID":0, "S":0}
    candidx = {"ID":candiS, "S":candiS}

    #ExtractRead_ID = list(); ExtractRead_S = list()
    #RefSeq_ID = list(); RefSeq_S = list()
    #ReadCount_ID = [0]; ReadCount_S = [0]
    #SumCount_ID = 0; SumCount_S = 0
    #candidx = candiS
    for candi in CanRen:
        flanking = getFlanking(candi)
        Start = candi[0][0] - flanking
        End = candi[-1][0] + flanking
        IsID = "S" if flanking == 50 else "ID"

        queries = sorted(ReadTree.overlap(Start, End))
        if len(queries) < args.minCNT:
            continue

        CurrentList = []; readidx = 0
        for l in queries:
            POS = l.begin
            SEQ, CIGAR, QUAL = l.data

            refStart = POS
            readStart = 0; ExtractStart = readStart; ExtractEnd = 0 
            flagS = 0
            inner = False
            if refStart < Start:
                flagS = 1
            elif refStart-Start > 20:
                inner = True

            for m in re.finditer(cigarRe, CIGAR):
                op_l = int(m.group(1))
                if m.group(2) == "S":
                    readStart += op_l
                    continue
                if m.group(2) in ("M", "=", "X"): # M, X, =
                    refStart += op_l
                    readStart += op_l
                elif m.group(2) == "I": # I
                    readStart += op_l
                elif m.group(2) == "D": # D
                    refStart += op_l
                
                # cal start
                if flagS == 1 and refStart >= Start:
                    if m.group(2) in ("M", "=", "X"):
                        ExtractStart = readStart - (refStart - Start)
                    else:
                        ExtractStart = readStart
                    flagS = 0
                # cal end
                if refStart >= End:
                    if m.group(2) in ("M", "=", "X"):
                        ExtractEnd = readStart - (refStart - End)
                    else:
                        ExtractEnd = readStart
                    break
            if End - refStart > 20:
                ExtractEnd = readStart
                inner = True

            Lread = ExtractEnd - ExtractStart
            Lref = End - Start

            #print readStart, ExtractStart, ExtractEnd, Lread, Lref
            if Lread > Lref/2 and Lread < 3*Lref/2:
                name = "Read_"+str(candidx[IsID])+"_" + str(ExtractStart) + "_" + str(ExtractEnd) + "_" + str(inner) + "_" +str(readidx)
                ql = [ord(i)-33 for i in QUAL[ExtractStart: ExtractEnd]] 
                Record = SeqIO.SeqRecord(seq=Seq(SEQ[ExtractStart:ExtractEnd]),name=name,id=name,description=str(inner),letter_annotations={"phred_quality":ql})
                #ExtractRead[IsID].append(Record)
                #SumCount[IsID] += 1
                CurrentList.append(Record)
                readidx += 1
        
        if readidx < 8: #remvoe some reads located inner the region
                CurrentList, readidx = remove_Inner_reads(CurrentList)
        if readidx > 1: #readidx > args.minCNT 
            SumCount[IsID] += readidx
            ReadCount[IsID].append(SumCount[IsID])
            ExtractRead[IsID].extend(CurrentList)
            name = "RefSeq_"+str(candidx[IsID])+"_" + str(Start) + "_" + str(End) + "_" + "_".join(["%d|%s|%d|%d|%s" %(x[0],x[1],x[3],x[4],x[5]) for x in candi])
            Record = SeqIO.SeqRecord(seq=ref[Start:End], name=name, id=name, description=name)
            RefSeq[IsID].append(Record)
            candidx[IsID] += 1
        del CurrentList

    filelist = {"ID":[], "S":[]}
    # ID
    for ty in ["ID", "S"]:
        for i in range(1, len(ReadCount[ty])):
            filename = fp_prefix + ty + "_" + "Region_" + str(candiS+i-1)+".fa"
            path = args.workDir + filename
            filelist[ty].append(path)
            SeqIO.write(ExtractRead[ty][ReadCount[ty][i-1]:ReadCount[ty][i]], path, "fastq")
    for ty in ["ID", "S"]:
        with open(args.workDir + fp_prefix + ty + "_" + "fileList", 'w') as fw:
            fw.write('\n'.join(filelist[ty]))
            fw.write('\n')
    for ty in ["ID", "S"]:
        SeqIO.write(RefSeq[ty], args.workDir + fp_prefix + ty + "_" + "RefSeq.fa", "fasta") 

    del ReadTree
    del ExtractRead
    del RefSeq
    del ReadCount
    del queries

    gc.collect()

    logging.info("[ExtractRead.py] Finish prefix: %s" %(fp_prefix))

def remove_Inner_reads(Record):
    CurrentList = []; Cnt = 0
    for item in Record:
        if item.description == "True":
            continue
        CurrentList.append(item)
        Cnt += 1
    return CurrentList, Cnt

def run():
    parser = argparse.ArgumentParser(description="MSA examples using abPOA")
    parser.add_argument("--fin_bam", type=str, default="input.bam",
            help="Sorted bam file, default: %(default)s")
    parser.add_argument("--fin_can", type=str, default="candi.info",
            help="Candidates of potential variants genenrated by ExtractVariantCandidate.py")
    parser.add_argument("--fin_ref", type=str, default="ref.fa",
            help="Reference fasta input")
    parser.add_argument("--varType", type=str, default="all",
            help="Extract candidates of SNP, indel or all of small variant,[snp, indel, all], default: %(default)s")
    parser.add_argument("--chrName", type=str, default="22",
            help="The name of reference to be processed, default: %(default)s")
    parser.add_argument("--chrStart", type=int, default=None,
            help="The 1-based starting positions of the reference to be processed")
    parser.add_argument("--chrEnd", type=int, default=None,
            help="The inclusive ending positions of the reference to be processed")
    parser.add_argument("--minMQ", type=int, default=10,
            help="Minimum Mapping Quality. Mapping quality lower than the setting will be filtered, default:%(default)d")
    parser.add_argument("--minCNT", type=int, default=3,
            help="Minimum read counts required to call a variant, default:%(default)d")
    parser.add_argument("--ChunkSite", type=int, default=50000,
            help="Divide job with smaller candidate site count for parallelism, default:%(default)d")
    parser.add_argument("--workDir", type=str, default="tempDir",
            help="Temporary working path for Multiple sequence alignment, default:%(default)s")
    parser.add_argument("--threads", type=int, default=1,
            help="Number of threads to use, default:%(default)d")
    parser.add_argument("--flanking", type=int, default=50,
            help="Flanking base pairs around variant site, default: %(default)d")
    parser.add_argument("--slideWindow", type=int, default=100,
            help="The window lenth for two candidate can be merge as one to speed up, default: %(default)d")

    args = parser.parse_args()    
    
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    setup_logging()
    cwd = os.getcwd()
    args.workDir = cwd + '/' + args.workDir + '/'

    main_ctrl(args)

if __name__ == "__main__":
    t_s = time()
    run()
    t_e = time()
    logging.info("Finish the script in %f seconds", t_e - t_s)





