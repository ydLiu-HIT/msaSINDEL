#!/usr/bin/env python
# coding=utf-8

###
#  get genotype likelihood by realignment mapping score
###

import os
import sys
import gc
import logging
import re
import math
import numpy as np
import argparse
from multiprocessing import Pool
from time import time
from utils import *
from check import *
import tempfile
from ksw_module import ksw_aligner
from math_func import *
from Canclass import *
from reMSA import reCall
from getHCVariant import get_hc_variantList

cigarRe = r"(\d+)([MIDNSHP=X])"
majorContigs = {"chr"+str(a) for a in list(range(0,23))+["X", "Y"]}.union({str(a) for a in list(range(0,23))+["X", "Y"]})
_MIN_MQ = 1e-10
_MAX_LIKELIHOODSCORE = 1e-250
_MAX_PROB = 1.0

#penalty
# optimal: GAPO=6, GAPE=2, MATCH=2
# optimal: GAPO=6, GAPE=3, MATCH=3, MQ = 5/10  基本没差别，个别例子
# optimal: GAPO=6, GAPE=3, MATCH=3, flank = 5/10, MQ = 5/10  基本没差别，个别例子
MATCH = 3
MISMATCH = -6
GAPO = -6
GAPE = -3
DELP = -7.8 #-7.8
INSP = -8.8
EXTFLANK = 5
AMPLIFIER = 10.0

CDEBUG = False

def setup_logging():
    """
    Default logger
    """
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(level=logging.INFO, format=log_format)

def main_ctrl(args):
    if not os.path.exists(args.workDir):
        sys.exit("Error: work folder %s not found" %(args.workDir))

    CheckFileExist(args.fin_ref)
    fpList = args.workDir+"prefixList"
    CheckFileExist(fpList)
    
    prefixList = list()
    with open(fpList, 'r') as fin:
        for line in fin:
            prefixList.append(line.strip())

    # get Ref
    REF_dict = LoadRef(args.fin_ref)

    ResultList = list(); semi_result = list()
    logging.info("Run the script with %d threads" %(args.threads))
    ## get HC variant from HC file candidates
    logging.info("Getting HC variants.")
    if not DEBUG: ResultList.extend(get_hc_variantList(args.canfix+".hc.info", args.perror_for_indel, args.perror_for_snp))

    # use multiple threads
    analysis_pools = Pool(processes=int(args.threads))
    for v in prefixList:
        ctg = v.split('_')[0]
        if ctg not in majorContigs:
            continue
        candiS = int(v.split('_')[1])
        param = [(v, ctg, args.workDir, REF_dict[ctg].seq, args.shift)]
        semi_result.append(analysis_pools.map_async(run_parallel, param))

    analysis_pools.close()
    analysis_pools.join()

    for res in semi_result:
        try:
            ResultList.extend(res.get()[0])
        except:
            pass

    # without multiple thread
    #for v in prefixList:
    #    ctg = v.split('_')[0]
    #    if ctg not in majorContigs:
    #        continue
    #    candiS = int(v.split('_')[1])
    #    param = (v, ctg, args.workDir, REF_dict[ctg].seq, args.shift)
    #    ResultList.extend(run_parallel(param)) 

    logging.info("Writing variant to disk")
    # writing header
    vcfOut = open(args.fout_vcf, 'w')
    PrintVCFHeader(vcfOut, args.fin_ref, "HG002")

    ResultList = sorted(ResultList, key=lambda var: int(var.pos))
    Generate_Output(vcfOut, ResultList)

    vcfOut.close()
    
def run_parallel(args):
    return ParseMSAandREF(*args)


def haveMatch(pos, typ, IndelCan, idx, shift):
    mindis = shift; minid = 0; sig = False
    for i in range(len(IndelCan)):
        r = int(IndelCan[i][0])
        t = IndelCan[i][1]
        if abs(pos - r) < mindis and typ == t:
            mindis = abs(pos-r); minid = i; sig = True

    return sig, minid

def get_reverse_complememt(seq):
    RC = {"A":"T", "T":"A", "G":"C", "C":"G"}
    return ''.join([RC[b] for b in seq[::-1]])

def getRefList(fpath):
    fg = open(fpath, 'r')
    RefSeqDict = dict()
    RefidDict = dict()
    rN = ""; Cnt = 0
    for line in fg:
        if line.startswith('>'):
            rN = line.strip(); Cnt += 1
            RefSeqDict[rN] = ""
            RefidDict[rN.split('_')[1]] = rN
        else:
            RefSeqDict[rN] += line.strip()
    fg.close()
    return RefSeqDict, RefidDict, Cnt
    
def getMSAList(fpath):
    fm = open(fpath, 'r')
    MsaSeqDict = dict()
    for line in fm:
        if line.startswith('>'):
            mid = line.strip().split('_')[2]
            MsaSeqDict.setdefault(mid, [])
        else:
            MsaSeqDict[mid].append(line.strip())
    fm.close()
    return MsaSeqDict

def compare2cigarsignal(alt, cigarsignal):
    if alt[1:] == cigarsignal:
        return alt
    else:
        return alt[0] + cigarsignal

def makeCandidate(chrName, alignment, RNAME, REF, SEQ, shift, source):
    quePos = 0; TS = int(RNAME[2]); refPos = TS
    IndelCan = []; Candidates = []
    for i in RNAME[4:]:
        IndelCan.append(i.split('|'))
    LL = len(IndelCan); idx = 0
    CIGAR = alignment[1]
    for m in re.finditer(cigarRe, CIGAR):
        op_l = int(m.group(1))
        if m.group(2) == "S":
            quePos += op_l
            continue
        if m.group(2) in ("M", "="): # M, X, =
            refPos += op_l
            quePos += op_l
        elif m.group(2) == "X":
            for _ in range(op_l):
                refPos += 1; quePos += 1;
                sig, idx = haveMatch(refPos, "S", IndelCan, idx, 2)
                if sig == True:
                    rseq = REF[refPos-1]; qseq = SEQ[quePos-1]
                    cand = [chrName, refPos, "S", str(rseq), qseq, source, alignment[0], 0, int(IndelCan[idx][2])-1, int(IndelCan[idx][3])-1, quePos, refPos-TS, SEQ]
                    Candidates.append(cand)
                    del IndelCan[idx]; HaveVar = True
        elif m.group(2) == "I":
            sig, idx = haveMatch(refPos, "I", IndelCan, idx, shift)
            if sig == True:
                rseq = REF[refPos - 1]; qseq = SEQ[quePos-1:quePos+op_l]
                #qseq = compare2cigarsignal(qseq, IndelCan[idx][4])
                cand = [chrName, refPos, "INS", str(rseq), qseq, source, alignment[0], op_l, int(IndelCan[idx][2])-1, int(IndelCan[idx][3])-1, quePos, refPos-TS, SEQ]
                Candidates.append(cand)
                del IndelCan[idx]; HaveVar = True
            quePos += op_l
        elif m.group(2) == "D":
            sig, idx = haveMatch(refPos, "D", IndelCan, idx, shift)
            if sig == True:
                rseq = REF[refPos-1:refPos+op_l]; qseq = SEQ[quePos-1]
                #rseq = compare2cigarsignal(rseq, IndelCan[idx][4])
                cand = [chrName, refPos, "DEL", str(rseq), qseq, source, alignment[0], op_l, int(IndelCan[idx][2])-1, int(IndelCan[idx][3])-1, quePos, refPos-TS, SEQ]
                Candidates.append(cand)
                del IndelCan[idx]; HaveVar = True
            refPos += op_l
    del IndelCan

    return Candidates 

def ParseMSAandREF(fp_prefix, chrName, workDir, REF, shift):
    CandidateList = list()
    s = int(fp_prefix.split('_')[1])
    outliers = list()
    for ty in["ID", "S"]:
        RefSeqDict, RefidDict, Cnt = getRefList(workDir + fp_prefix + ty + "_RefSeq.fa") 
        MsaSeqDict = getMSAList(workDir + fp_prefix + ty + "_MSA.fa")

        CandidateDict = dict()
        e = s + Cnt
        for i in range(s, e):
            i = str(i); rN = RefidDict[i]; CandidateDict[i] = []
            #if rN != ">RefSeq_1986_19934333_19934534_19934433|D|11|22|T_19934434|S|1|11|G": continue
            #print(rN)
            RefCons = RefSeqDict[rN]
            MSACons = MsaSeqDict[i]
            RNAME = rN.split('_')
            haveIns = False
            for msa in MSACons:
                alignment = ksw_aligner(msa, RefCons)
                alignment[1] = getExactMatch(alignment[1], RefCons, msa, 0)
                if DEBUG:
                    print("ref:", RefCons)
                    print("msa:", msa)
                    print(alignment)

                if alignment[0] < 0: continue
                instance = makeCandidate(chrName, alignment, RNAME, REF, msa, shift, "LC")
                if len(instance) > 0: haveIns = True
                CandidateDict[i].extend(instance)
            if not haveIns and ty == "ID": outliers.append([ty, i])
        #if ty == "ID":
        #    logging.info("[call.py] Re-MSA for prefix: %s" %(fp_prefix))
        #    fp_outliers = reCall(workDir, fp_prefix, outliers, 0.2)
        #    MsaSeqDict = getMSAList(fp_outliers)
        #    for i in range(s, e):
        #        i = str(i); flag = MsaSeqDict.get(i)
        #        if flag is None: continue
        #        rN = RefidDict[i]; CandidateDict[i] = []   
        #        RefCons = RefSeqDict[rN]
        #        MSACons = MsaSeqDict[i]
        #        RNAME = rN.split('_')
        #        for msa in MSACons:
        #            alignment = ksw_aligner(msa, RefCons)
        #            alignment[1] = getExactMatch(alignment[1], RefCons, msa, 0)
        #            if DEBUG:
        #                print "ref:", RefCons
        #                print "msa:", msa
        #                print alignment

        #            if alignment[0] < 0: continue
        #            CandidateDict[i].extend(makeCandidate(chrName, alignment, RNAME, REF, msa, shift, "reMSA"))

        CanList = mergeCoVar(CandidateDict, s, e, RefidDict, RefSeqDict)
        
        #sort
        CanList = sorted(CanList, key=lambda x:x.pos) 

        #cal genotype quality
        most_likely_genotype(CanList, workDir, fp_prefix, ty) 
        CandidateList.extend(CanList)
        
        del CanList
        del CandidateDict
        del RefSeqDict
        del RefidDict
        gc.collect()

    if len(CandidateList) == 0:
        logging.info("[call.py] No variant founded.")
    logging.info("[call.py] Fininsh prefix: %s" %(fp_prefix))
    
    return CandidateList

def re_calc_MQ_locally(cigars, qual, ts, vl, extflank):
    start = ts - extflank
    end = ts + extflank + vl + 1
    rpos = 0; qpos = 0; prob = 0
    if CDEBUG: print("cigar:", cigars)
    for m in re.finditer(cigarRe, cigars):
        advance = int(m.group(1))
        op = m.group(2)
        if op in ("=", "X", "D"):
            tp = advance + rpos - start
            rpos += advance; qpos += 0 if op == "D" else advance
            if tp > 0:
                tp = tp if rpos < end else tp+end-rpos
                qpos = qpos if rpos < end else qpos if op == "D" else qpos+end-rpos
                if op == "=":
                    if CDEBUG: print("=", qpos, tp, end=' ')
                    for i in range(qpos-tp, qpos):
                        Q = ord(qual[i])-33; 
                        if CDEBUG: print(Q, end=' ')
                        prob += math.log10(1 - math.pow(10,-Q/10.0))
                    if CDEBUG: print('\n')
                elif op == "X":
                    if CDEBUG: print("X", qpos, tp, end=' ')
                    for i in range(qpos-tp, qpos):
                        Q = ord(qual[i])-33; 
                        if CDEBUG: print(Q,)
                        prob += (-Q/10.0)
                    if CDEBUG: print('\n')
                elif op == "D":
                    if CDEBUG: print("D", qpos, tp, end=' ')
                    Q_max = 0
                    s = qpos-2 if qpos-2>0 else 0
                    e = qpos+3 if qpos+3<len(qual) else len(qual)
                    for i in range(s, e):
                        Q = ord(qual[i])-33
                        Q_max = Q if Q > Q_max else Q_max
                    prob += (-Q_max / 10.0) * tp
                start = rpos
            if rpos >= end:
                break
        elif op == "I":
            if rpos < start:
                qpos += advance; continue
            if CDEBUG: print("I", qpos, advance,)
            for i in range(qpos, qpos + advance):
                Q = ord(qual[i])-33; 
                if CDEBUG: print(Q, end=' ')
                prob += (-Q/10.0)
            qpos += advance
            if CDEBUG: print("\n")
    if CDEBUG: print("prob=", prob) 
    return prob

def re_calc_AS_locally(cigars, ts, vl, extflank, match=MATCH, mismatch=MISMATCH, gapo=GAPO, gape=GAPE):
    start = ts - extflank
    end = ts + extflank + vl
    rpos = 0
    score_A = 0; score_p = 0; score_log = 0; misLen = 0
    for m in re.finditer(cigarRe, cigars):
        advance = int(m.group(1))
        op = m.group(2)
        if op in ("=", "X","M", "D"):
            tp = advance + rpos - start
            rpos += advance
            if tp > 0:
                score_A += (tp*match if op in ("M", "=") else gapo+tp*gape if op == "D" else tp*mismatch)
                if op == "D":
                    score_p += gapo+tp*gape
                    misLen += tp
                start = rpos
            if rpos > end:
                tp = rpos - end
                score_A -= (tp*match if op in ("M", "=") else gapo+tp*gape if op == "D" else tp*mismatch)
                if op == "D":
                    score_p -= gapo+tp*gape
                    misLen -= tp
                break
        elif op == "I":
            if rpos < start:
                continue
            score_A += (gapo+advance*gape)
            score_p += (gapo+advance*gape)
            misLen += advance

    return score_A, score_p, misLen

def getExactMatch(Cigars, target, query, start):
    Ts = start; Qs = 0
    newCigar = ''
    for m in re.finditer(cigarRe, Cigars):
        advance = int(m.group(1))
        op = m.group(2)
        if op == "D":
            Ts += advance
            newCigar += (str(advance)+op)
        elif op == "I" or op =="S":
            Qs += advance
            newCigar += (str(advance)+op)
        elif op == "M":
            tM = 0; tX = 0
            for i in range(advance):
                if target[Ts+i] == query[Qs+i]:
                    if tX > 0: newCigar += (str(tX)+"X"); tX = 0
                    tM += 1
                else:
                    if tM > 0: newCigar += (str(tM)+"="); tM = 0
                    tX += 1

            newCigar += (str(tM)+"=") if tM > 0 else (str(tX)+"X")
            
            Ts += advance; Qs += advance
    return newCigar

def realignWithoutQuality(fileList, Haplot, var, Param):
    likeliScore = []
    with open(fileList[str(var.SeqID)], 'r') as flocal:
        querySeq = ''
        inner = ''
        pre_inner = ''
        for line in flocal:
            if line.startswith('>'):
                inner = line.split('_')[4]
                if querySeq == '' or pre_inner == "True":
                    pre_inner = inner
                    querySeq = ''
                    continue
                
                ## sw
                rsc = []; od = 0
                for hap in Haplot:
                    TS, TL, LH = Param[od]; od += 1;
                    alignment = ksw_aligner(querySeq, hap)
                    newCigar = getExactMatch(alignment[1], hap, querySeq, 0); alignment[1] = newCigar
                    alignment[0], score_p, misLen = re_calc_AS_locally(alignment[1], TS, TL, LH)
                    err = MATCH*LH - alignment[0]
                    #print LH, alignment[0], err, score_p, alignment[1] 
                    alignment[0] = min(max(math.pow(10, -err/10.0), _MIN_MQ), _MAX_PROB) # method1
                    #alignment[0] = max(math.pow(10, score_p/10.0), _MIN_MQ) # method2
                    #alignment[0] = math.pow(detal, misLen) # method3
                    rsc.append(alignment)
                else:
                    likeliScore.append(rsc)
                querySeq = ''
                pre_inner = inner
            else:
                querySeq += line.strip()
        else: # last read
            if pre_inner == "False":
                rsc = []; od = 0
                for hap in Haplot:
                    TS, TL, LH = Param[od]; od += 1
                    alignment = ksw_aligner(querySeq, hap)
                    newCigar = getExactMatch(alignment[1], hap, querySeq, 0); alignment[1] = newCigar
                    alignment[0], score_p, misLen = re_calc_AS_locally(alignment[1], TS, TL, LH)
                    err = MATCH*LH - alignment[0]
                    alignment[0] = min(max(math.pow(10, -err/10.0), _MIN_MQ), _MAX_PROB) # method1
                    #alignment[0] = max(math.pow(10, score_p/10.0), _MIN_MQ) # method2
                    #alignment[0] = math.pow(detal, misLen) # method3
                    rsc.append(alignment)
                else:
                    likeliScore.append(rsc) 

    return likeliScore

def realignWithQuality(fileList, Haplot, var, Param, amplifier, extflank):
    if DEBUG:
        print("amplifier = ", amplifier)
    likeliScore = []
    with open(fileList[str(var.SeqID)], 'r') as flocal:
        querySeq = ''; qual = ''
        inner = "True"
        idx = 0
        name = ''
        for line in flocal:
            idx += 1
            linum = idx % 4
            if linum == 3: continue
            if linum == 2: # query read
                querySeq = line.strip()
                continue
            if linum == 0: # qual
                qual = line.strip()
                continue
            if linum == 1: # read name
                if inner == "True":
                    inner = line.split('_')[4]
                    continue
                inner = line.split('_')[4]
                name = line.strip()
            # ksw
            rsc = []; od = 0
            for hap in Haplot:
                TS, TL, LH = Param[od]; od += 1;
                alignment = ksw_aligner(querySeq, hap)
                newCigar = getExactMatch(alignment[1], hap, querySeq, 0); alignment[1] = newCigar
                prob = re_calc_MQ_locally(alignment[1], qual, TS, TL, extflank); alignment[0] = prob
                alignment[0] = min(max(math.pow(10, prob/amplifier), _MIN_MQ), _MAX_PROB) # method1
                rsc.append(alignment)
            else:
                likeliScore.append(rsc)
        else: # last read
            if inner == "False":
                rsc = []; od = 0
                for hap in Haplot:
                    TS, TL, LH = Param[od]; od += 1
                    alignment = ksw_aligner(querySeq, hap)
                    newCigar = getExactMatch(alignment[1], hap, querySeq, 0); alignment[1] = newCigar
                    prob = re_calc_MQ_locally(alignment[1], qual, TS, TL, extflank); alignment[0] = prob
                    alignment[0] = min(max(math.pow(10, prob/amplifier), _MIN_MQ), _MAX_PROB) # method1
                    rsc.append(alignment)
                else:
                    likeliScore.append(rsc) 

    return likeliScore

def realignToHaplotype(fileList, Haplot, var): 
    Param = []
    VLen = str(var.varlen).split(',')
    QPos = str(var.quePos).split(',')
    Vtype = var.type.split(',')
    extflank = 1 if var.type == "S" else EXTFLANK
    TS = var.refPos; TL = max([int(i) for i in VLen]) if "DEL" in Vtype else 0; LH = extflank*2 + TL
    Param.append((TS, TL, LH))
    for l in range(len(VLen)):
        TS_alt = int(QPos[l]); TL_alt = int(VLen[l]) if Vtype[l] == "INS" else 0; LH_alt = extflank*2 + TL_alt
        Param.append((TS_alt, TL_alt, LH_alt))
    #get local reads
    if DEBUG:
        print(Param)
    
    amplifier = 3.5 if var.type == 'S' else AMPLIFIER
    #likeliScore = realignWithoutQuality(fileList, Haplot, var, Param)
    likeliScore = realignWithQuality(fileList, Haplot, var, Param, amplifier, extflank)

    return likeliScore

def most_likely_genotype(CandidateList, workDir, fp_prefix, ty):
    fileList = {}

    with open(workDir+fp_prefix + ty +"_fileList", 'r') as flist:
        for line in flist:
            fileList.setdefault(line.split('_')[-1].split('.')[0], line.strip())

    for var in CandidateList:
        if DEBUG:
            print(var.pos, var.type, var.ALT, var.varlen, var.quePos, var.refPos)
        #get halplotypes
        Haplot = []
        Haplot.append(var.ref)  # ref
        Haplot.extend(var.alt)

        if DEBUG:
            print("ref:", Haplot[0])
            print("alt:", Haplot[1:])
        
        #######    method 1  ###################
        #       Based on alignment score #
        likeliScore = realignToHaplotype(fileList, Haplot, var)
        if DEBUG:
            print(np.array(likeliScore))
        #
        GenotypeLikeliHood = getGenotypeLikeliHood(likeliScore, len(Haplot))
        log10_probs = normalize_log10_probs([math.log10(p) for p in GenotypeLikeliHood])
        real_probs = toRealSpace(log10_probs)
        predictions = round_gls(real_probs)
        gq, qual, gt = computer_quals(predictions, len(Haplot))
        
        if DEBUG:
            print("genotype likelihood:", GenotypeLikeliHood)
            print("log10_prob", log10_probs)
            print("real_prob:", real_probs)
            print("rounded probs:", predictions)
            print(gq, qual, gt)
        
        var.qual = qual
        var.GQ = gq
        var.GT = str(gt[0]) + '/' + str(gt[1])

    del fileList

def confidence_read(prob):
    pV = [p[0] for p in prob]
    Midx = pV.index(max(pV))
    midx = pV.index(min(pV))
    #if pV[Midx] - pV[midx] < 0.5:
    #    return None
    if pV[Midx] < pV[midx] * 3 and pV[Midx] > 0.99 and pV[midx] < 0.7:  #0.999 get max 
        prob[midx][0] = 0.1

    return prob

def getGenotypeLikeliHood(likeliScore, H):
    HL = sum(x+1 for x in range(H))
    GenotypeLikeHood = [1.0 for i in range(HL)]
    breaksig = False

    for prob in likeliScore:
        prob = confidence_read(prob)
        if prob == None:
            continue

        gt = 0
        for i in range(H):
            for j in range(0, i + 1):
                p = (prob[i][0] + prob[j][0])/2
                if i != j: GenotypeLikeHood[gt] *= p if p < 0.6 else 0.6
                else: GenotypeLikeHood[gt] *= p
                gt += 1

        for t in GenotypeLikeHood:
            if t < _MAX_LIKELIHOODSCORE:
                breaksig = True; break
        if breaksig:
            break
    
    return GenotypeLikeHood

def createVariantCandidate(value, ref, sid, *alt):
    var = variantsCan(value)
    
    var.set_seqID(sid)
    for a in alt:
        var.set_alt(a)

    qpos = str(var.quePos).split(',')[0]
    var.set_ref(get_ref_haplotype(alt[0], int(qpos), var.REF, var.ALT.split(',')[0]))
    var.refPos = int(qpos)
    
    return var

def get_ref_haplotype(alt_hap, pos, REF, ALT):
    ref_hap = alt_hap[0:pos-1]
    ref_hap += REF
    ref_hap += alt_hap[pos+len(ALT)-1:]
    
    return ref_hap

def mergeCoVar(CandidateDict, s, e, RefidDict, RefSeqDict):
    CandidateList = list()
    #merge variant in the same position
    for i in range(int(s), int(e)):
        i = str(i)
        Lcan = len(CandidateDict[i])
        if Lcan == 0:
            continue
        elif Lcan == 1:
            var = createVariantCandidate(CandidateDict[i][0][:-1], RefSeqDict[RefidDict[str(i)]], i, CandidateDict[i][0][-1])
            CandidateList.append(var)
        else: # maybe multiploid
            CandidateDict[i] = sorted(CandidateDict[i], key=lambda x:x[1])
            ins1 = CandidateDict[i][0]
            t = 1
            while t < Lcan:
                ins2 = CandidateDict[i][t]
                if ins2[1] != ins1[1]:
                    var = createVariantCandidate(ins1[:-1], RefSeqDict[RefidDict[str(i)]], i, ins1[-1])
                    CandidateList.append(var)
                    ins1 = CandidateDict[i][t]
                    t += 1
                elif ins2[2] != ins1[2]: #the same pos, ins and del
                    ins1[3] = ins1[3] if len(ins1[3]) > len(ins2[3]) else ins2[3]
                    dl = ins1[3][1:]
                    dA = ins1[4]+dl if len(ins1[4]) > 1 else ins1[4]
                    dA += ','
                    dA += ins2[4]+dl if len(ins2[4]) > 1 else ins2[4]
                    ins1[2] = ins1[2] + ',' + ins2[2]
                    ins1[4] = dA
                    ins1[7] = str(ins1[7]) + ',' + str(ins2[7])
                    ins1[10] = str(ins1[10]) + ',' + str(ins2[10])
                    
                    var = createVariantCandidate(ins1[:-1], RefSeqDict[RefidDict[str(i)]], i, ins1[-1], ins2[-1])
                    CandidateList.append(var)

                    if t < Lcan - 1:
                        ins1 = CandidateDict[i][t+1]
                    t += 2
                else: # the same pos, MNP
                    if ins2[4] != ins1[4]:
                        ins1[2] = "INS,INS"
                        ins1[4] = ins1[4] + ',' + ins2[4]
                        ins1[7] = str(ins1[7]) + ',' + str(ins2[7])
                        ins1[10] = str(ins1[10]) + ',' + str(ins2[10])
                        var = createVariantCandidate(ins1[:-1], RefSeqDict[RefidDict[str(i)]], i, ins1[-1], ins2[-1])
                    elif ins1[2] == "DEL" and ins1[3] != ins2[3]:
                        ins1[2] = "DEL,DEL"
                        ins1[3] = ins1[3] if len(ins1[3]) > len(ins2[3]) else ins2[3]
                        ins1[4] = ins1[3][0:len(ins1[3])-ins1[7]] + ',' + ins1[3][0:len(ins1[3])-ins2[7]]
                        ins1[7] = str(ins1[7]) + ',' + str(ins2[7])
                        ins1[10] = str(ins1[10]) + ',' + str(ins2[10])
                        var = createVariantCandidate(ins1[:-1], RefSeqDict[RefidDict[str(i)]], i, ins1[-1], ins2[-1])
                    else: # the same variant record  
                        var = createVariantCandidate(ins1[:-1], RefSeqDict[RefidDict[str(i)]], i, ins1[-1])

                    CandidateList.append(var)

                    if t < Lcan - 1:
                        ins1 = CandidateDict[i][t+1]
                    t += 2
            if t == Lcan:
                var = createVariantCandidate(ins1[:-1], RefSeqDict[RefidDict[str(i)]], i, ins1[-1])
                CandidateList.append(var)

    return CandidateList

def run():
    parser = argparse.ArgumentParser(description="MSA examples using abPOA")
    parser.add_argument("--fin_ref", type=str, default="ref.fa",
            help="Reference fasta input")
    parser.add_argument("--fout_vcf", type=str, default="detect.vcf",
            help="vcf path for detected variants, default: %(default)s")
    parser.add_argument("--workDir", type=str, default="tempDir",
            help="Temporary working path for Multiple sequence alignment, default:%(default)s")
    parser.add_argument("--canfix", type=str, default="candidate",
            help="Candidate file prefix used in ExtractVariant.py, default: %(default)s")
    parser.add_argument("--shift", type=int, default=5,
            help="The distance bewteen the detect variant and the activate region, default:%(default)d")
    parser.add_argument("--perror_for_snp", type=float, default="0.1",
            help="P-error is the probability of observing a heterozygote SNP, default:%(default)f")
    parser.add_argument("--perror_for_indel", type=float, default="0.1",
            help="P-error is the probability of observing a heterozygote indel, default:%(default)f")
    parser.add_argument("--threads", type=int, default=1,
            help="Number of threads to use, default:%(default)d")
    parser.add_argument("--debug", type=bool, default=False,
            help="debugs, default:%(default)s")
    args = parser.parse_args()    
    
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)
    setup_logging()
    cwd = os.getcwd()
    args.workDir = cwd + '/' + args.workDir + '/'

    global DEBUG
    DEBUG = True if args.debug else False

    main_ctrl(args)

if __name__ == "__main__":
    t_s = time()
    run()
    t_e = time()
    logging.info("Finish the script in %f seconds", t_e - t_s)





