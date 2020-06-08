#!/usr/bin/env python
# coding=utf-8
import os
import sys
import argparse
import logging
import intervaltree
import param
import re
import shlex
import subprocess
from collections import Counter
import gc
from Bio import SeqIO
from time import time
from collections import defaultdict
from multiprocessing import Pool
from utils import LoadRef
from check import CheckFileExist

Multi = True

cigarRe = r"(\d+)([MIDNSHP=X])"
#majorContigs = {"chr"+str(a) for a in list(range(0,23))+["X", "Y"]}.union({str(a) for a in list(range(0,23))+["X", "Y"]})
majorContigs = {"chr"+str(a) for a in list(range(0,23))}.union({str(a) for a in list(range(0,23))})
def setup_logging():
    """
    Default logger
    """
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(level=logging.INFO, format=log_format)

def split_reference(REF_dict, args):
    totalLen = 0
    refProcessList = list()
    order = 0
    if args.chrStart != None and args.chrEnd != None:
        s = args.chrStart
        while s < args.chrEnd:
            e = s + args.bandWidth
            if e > args.chrEnd:
                e = args.chrEnd
            refProcessList.append((args.chrName, s, e, order))
            s = e
            order += 1
    else:
        for key in REF_dict:
            if key not in majorContigs:
                continue
            L = len(REF_dict[key].seq)
            s = 0
            while s < L:
                e = s + args.bandWidth
                if e > L:
                    e = L
                refProcessList.append((key, s, e, order))
                s = e
                order += 1

    refProcessList = sorted(refProcessList, key = lambda x: (x[2]-x[1]), reverse = True)
    #for item in refProcessList:
    #    print item
    return refProcessList

def main_ctrl(args):
    CheckFileExist(args.fin_bam + ".bai")
    REF_dict = LoadRef(args.fin_ref)
    refProcessList = split_reference(REF_dict, args)

    can_region = {}
    if args.fin_bed != None:
        with open(args.fin_bed, 'r') as f:
            for line in f:
                line = line.strip().split()
                contig = line[0]
                if contig not in can_region:
                    can_region[contig] = intervaltree.Intervaltree()
                begin = int(line[1])
                end = int(line[2])
                can_region[contig].addi(begin, end)
            if args.chrName not in can_region:
                logging.error("chrName %s is not in the bed file, are you using the correct bed file (%s)" %(args.chrName, args.fin_bed))
                sys.exit(1)
    ## out put
    can_hc_fp = open(args.fout_can+".hc.info", 'w')
    can_lc_fp = open(args.fout_can+".lc.info", 'w')

    ResultList = list(); Order = []
    semi_result = list()

    if Multi:
        analysis_pools = Pool(processes=int(args.threads))
        logging.info("Run the script with %d threads" %(args.threads))

        for batch in refProcessList:
            contig = batch[0]
            Order.append(batch[3])
            param = [(REF_dict[contig], contig, int(batch[1]), int(batch[2]), can_region, args)]
            ResultList.append(analysis_pools.map_async(run_parallel, param))

        analysis_pools.close()
        analysis_pools.join()

        for res in ResultList:
            try:
                semi_result.append(res.get()[0])
            except:
                pass
    else:
        for batch in refProcessList:
            contig = batch[0]
            Order.append(batch[3])
            param = (REF_dict[contig], contig, int(batch[1]), int(batch[2]), can_region, args)
            semi_result.append(run_parallel(param))

    reverseOrder = sorted(list(range(len(Order))), key = lambda k:Order[k])
    
    logging.info("Writing into disk.")
    for o in reverseOrder:
        for item in semi_result[o]:
            which = divide_candidate(item, args.minCov_for_snp, args.minCov_for_indel)
            if which == "HC":
                can_hc_fp.write(' '.join([str(x) for x in item]))
                can_hc_fp.write('\n')
            elif which == "LC":
                can_lc_fp.write(' '.join([str(x) for x in item]))
                can_lc_fp.write('\n')

    can_hc_fp.close()
    can_lc_fp.close()

SNP2S = {"A":"S", "C":"S", "G":"S", "T":"S", "I":"I", "D":"D", "N":"S"}
def divide_candidate(item, minCov_for_snp, minCov_for_indel):
    chrN, pos, readC, refB, refC, typ, alleleC = item[0:7]
    if readC < 1:
        return "Filter"

    typ = SNP2S[typ] 
    if typ == "S":
        if alleleC < minCov_for_snp: return "Filter"
        if len(item) > 7: return "LC"  # multiallel
        ratio = float(alleleC) / readC
        if readC >= 12:
            return "HC" if ratio > 0.3 else "LC"   #0.4, 0.3(best)
        elif readC > 5:
            return "HC" if ratio > 0.8 else "LC"
        else:
            return "LC"
    else:
        if alleleC == 1 or (alleleC < minCov_for_indel and float(alleleC)/readC < 0.75): return "Filter"
        if len(item) > alleleC + 7: return "LC"  # multiallel
        mostAllele = Counter(item[7:]).most_common(2)
        if len(mostAllele) == 1:
            ratio = float(mostAllele[0][1]) / readC
            if readC >= 12:
                return "HC" if ratio > 0.4 else "LC"  #0.4  F1:99.52
            elif readC > 5:
                return "HC" if ratio > 0.8 else "LC"
            else:
                return "LC"
        else:
            most, second = mostAllele
            if float(second[1]) / most[1] < 0.2 and float(most[1]) / readC > 0.4:
                return "HC"
            else:
                return "LC"

def run_parallel(args):
    return MakeCandidates(*args)

def MakeCandidates(REF, chrName, chrStart, chrEnd, can_region, args):
    chrStart += 1
    refStart = chrStart; refEnd = chrEnd
    refStart -= param.expandReferenceRegion
    refStart = 1 if refStart < 1 else refStart
    refEnd += param.expandReferenceRegion
    refEnd = len(REF.seq) if refEnd > len(REF.seq) else refEnd
    RefSeq = str(REF.seq[refStart:refEnd]).upper()

    Result_list = list() 
    #pileup = {}
    #pileup_l = {}
    pileup = defaultdict(lambda: {"A":0,"C":0,"D":0,"G":0,"I":0,"N":0,"T":0})
    pileup_l = defaultdict(lambda: {"I":[], "D":[]})
    sweep = 0

    p2 = subprocess.Popen(shlex.split("samtools view -F 2316 -q %d %s %s:%d-%d" % (args.minMQ, args.fin_bam, chrName, chrStart, chrEnd) ), stdout=subprocess.PIPE, bufsize=8388608, universal_newlines=True)   #2308

    for l in p2.stdout:
        l = l.strip().split()
        if l[0][0] == "@":
            continue

        QNAME = l[0]
        RNAME = l[2]

        if RNAME != chrName:
            continue
        
        
        FLAG = int(l[1])
        POS = int(l[3]) - 1 # switch from 1-base to 0-base to match sequence index
        MQ = int(l[4])
        CIGAR = l[5]
        SEQ = l[9].upper()
        refPos = POS
        queryPos = 0

        #if MQ < args.minMQ:
        #    continue
        #skipBase = 0
        #totalAlnPos = 0
        #for m in re.finditer(cigarRe, CIGAR):
        #    advance = int(m.group(1))
        #    totalAlnPos += advance
        #    if m.group(2)  == "S":
        #        skipBase += advance

        #if 1.0 - float(skipBase) / (totalAlnPos + 1) < 0.55: # skip a read less than 55% aligned
        #    continue

        for m in re.finditer(cigarRe, CIGAR):
            advance = int(m.group(1))
            if m.group(2) == "S":
                queryPos += advance
                continue
            if m.group(2) in ("M", "=", "X"):
                for _ in range(advance):
                    pileup[refPos][SEQ[queryPos]] += 1
                    refPos += 1
                    queryPos += 1
            elif m.group(2) == "I":
                if advance < 53:
                    pileup[refPos-1]["I"] += 1
                    pileup_l[refPos-1]["I"].append(SEQ[queryPos:queryPos+advance])
                queryPos += advance
            elif m.group(2) == "D":
                if advance < 53:
                    pileup[refPos-1]["D"] += 1
                    pileup_l[refPos-1]["D"].append(RefSeq[refPos - refStart: refPos - refStart + advance])
                refPos += advance

        while sweep < POS:
            flag = pileup.get(sweep)
            if flag is None:
                sweep += 1; continue
            baseCount = list(pileup[sweep].items())
            baseCount_l = pileup_l[sweep]
            refBase = RefSeq[sweep - refStart]
            outFlag = 0
            out = None
            if sweep >= chrStart and sweep <= chrEnd:
                if args.fin_bed != None:
                    if chrName in can_region and len(can_region[chrName].search(sweep)) != 0:
                        outFlag = 1
                else:
                    outFlag = 1
            
            if outFlag == 1:
                out = CheckCandaidate(chrName, sweep, baseCount, baseCount_l, refBase, args.minRatio, args.variantType)
            if out != None:
                totalCount, outline = out
                Result_list.append(outline)
            del pileup[sweep]
            del pileup_l[sweep]
            sweep += 1;
    p2.stdout.close()
    p2.wait()

    # check remaining bases
    remainder = list(pileup.keys())
    remainder.sort()
    for sweep in remainder:
        baseCount = list(pileup[sweep].items())
        baseCount_l = pileup_l[sweep]
        refBase = RefSeq[sweep - refStart]
        outFlag = 0
        out = None
        if sweep >= chrStart and sweep <= chrEnd:
            if args.fin_bed != None:
                if chrName in can_region and len(can_region[chrName].search(sweep)) != 0:
                    outFlag = 1
            else: 
                outFlag = 1
        
        if outFlag == 1:
            out = CheckCandaidate(chrName, sweep, baseCount, baseCount_l, refBase, args.minRatio, args.variantType)
        if out != None:
            totalCount, outline = out
            Result_list.append(outline)
 
    del pileup
    del pileup_l
    gc.collect()

    logging.info("[ExtractVariant] Fininsh region (%s:%d-%d)" %(chrName, chrStart, chrEnd))
    return Result_list 

def CheckCandaidate(chrName, pos, baseCount, baseCount_l, refBase, threshold, variantType):
    if refBase == "N":
        return None
    totalCount = 0; readCount = 0; indelCount = 0; refCount = 0
    for x in baseCount:
        totalCount += x[1]
        if x[0] not in ["I", "D"]: 
            readCount += x[1]
            if x[0] == refBase:
                refCount += x[1]
        else:
            indelCount += x[1]

    denominator = totalCount
    if denominator == 0:
        denominator = 1
    baseCount.sort(key = lambda x:-x[1])  # sort baseCount descendingly
    
    #print(pos+1, totalCount, refBase, readCount, refCount, indelCount, baseCount, baseCount_l)

    p0 = float(baseCount[0][1]) / denominator
    
    sig = False
    if baseCount[0][0] == refBase and p0 > 1.0 - threshold:
        return None
    outC = [chrName, pos+1, readCount, refBase, refCount]
    for baseC in baseCount[:3]:
        p1 = float(baseC[1]) / denominator
        if p1 >= threshold and baseC[0] != refBase:
            typ = baseC[0]
            if typ in ["I", "D"]:  # indel candidate
                if variantType == "indel" or variantType == "all":
                    outC.extend([typ, len(baseCount_l[typ])])
                    outC.extend(baseCount_l[typ])
                    sig = True
            else: # SNP candidate
                if variantType == "snp" or variantType == "all":
                    outC.extend([typ, baseC[1]])
                    sig = True

    if sig: return totalCount, outC

    return None

def main():
    parser = argparse.ArgumentParser(description="Generate variant candidates using alignments")
    parser.add_argument("--fin_bam", type=str, default="input.bam",
            help="Sorted bam file, default: %(default)s")
    parser.add_argument("--fin_bed", type=str, default=None,
            help="Call variants only in these regions")
    parser.add_argument("--fin_ref", type=str, default="ref.fa",
            help="Reference fasta input, default: %(default)s")
    parser.add_argument("--fout_can", type=str, default="candidate",
            help="Variant candidate output prefix, default: %(default)s")
    parser.add_argument("--minMQ", type=int, default=10,
            help="Minimum Mapping Quality. Mapping quality lower than the setting will be filtered, default:%(default)d")
    parser.add_argument("--minCov_for_snp", type=int, default=2,
            help="Minimum read counts required to call a snp, default:%(default)d")
    parser.add_argument("--minCov_for_indel", type=int, default=4,
            help="Minimum read counts required to call a indel, default:%(default)d")
    parser.add_argument("--minRatio", type=float, default=0.15,   ## previouse use 0.125 as default
            help="Minimum variant supported read count ratio, default:%(default)f")
    parser.add_argument("--chrName", type=str, default=None,
            help="The name of reference to be processed, default:%(default)s")
    parser.add_argument("--chrStart", type=int, default=None,
            help="The 1-based starting positions of the reference to be processed")
    parser.add_argument("--chrEnd", type=int, default=None,
            help="The inclusive ending positions of the reference to be processed")
    parser.add_argument("--variantType", type=str, default="all",
            help="Extract candidates of SNP, indel or all of small variant,[snp, indel, all], default: %(default)s")
    parser.add_argument("--threads", type=int, default=1,
            help="Number of threads to use, default:%(default)d")
    parser.add_argument("--bandWidth", type=int, default=10000000,
            help="Reference length to detect candidate in one loop, default:%(default)d")
    parser.add_argument("--tools", type=str, default="samtools",
            help="Method to process sorted bam file, samtools or pysam, default:%(default)s")
    parser.add_argument("--DEBUG", type=param.str2bool, nargs='?', const=True, default=False,
            help="Output some results for debugging, defualt:%(default)s")
    args = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    setup_logging()
    #MakeCandidates(args)
    main_ctrl(args)
    

if __name__ == "__main__":
    t_s = time()
    main()
    t_e = time()
    logging.info("Finish the script in %f seconds", t_e - t_s)
