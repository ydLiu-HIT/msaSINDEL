#!/usr/bin/env python
# coding=utf-8

import os
import sys
import subprocess
import logging
#import shlex
import argparse
from multiprocessing import Pool
from time import time
from check import *

def setup_logging():
    """
    Default logger
    """
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(level=logging.INFO, format=log_format)

def main_ctrl(args):
    if not os.path.exists(args.workDir):
        sys.exit("Error: work folder %s not found" %(args.workDir))

    fpList = args.workDir+"prefixList"
    CheckFileExist(fpList)
    
    prefixList = list()
    with open(fpList, 'r') as fin:
        for line in fin:
            prefixList.append(line.strip())

    analysis_pools = Pool(processes=int(args.threads))
    logging.info("Run the script with %d threads" %(args.threads))

    for v in prefixList:
        candiS = int(v.split('_')[1])
        param = [(v, args.workDir, candiS, args.freq_for_indel, args.freq_for_SNP)]
        analysis_pools.map_async(run_parallel, param)

    analysis_pools.close()
    analysis_pools.join()

    #for v in prefixList:
    #    candiS = int(v.split('_')[1])
    #    param = (v, args.workDir, candiS, args.freq_for_indel, args.freq_for_SNP)
    #    run_parallel(param)
    
    logging.info("Finishing MSA and BWA-mem")
    
def run_parallel(args):
    return call_MSA(*args)

def call_MSA(fp_prefix, workDir, candiS, freq_for_indel, freq_for_SNP): 
    for ty in ["ID", "S"]:
        path = workDir + fp_prefix + ty + "_MSA.fa"
        fw = open(path, 'w')
        freq = freq_for_indel if ty == "ID" else freq_for_SNP
        cmd = "abPOA -c -m 2 -f %f -l %s" %(freq, workDir+fp_prefix+ty+"_fileList")
        f = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
        idx = 0; orde = 0
        for row in f.stdout:
            if idx % 2 == 0:
                if len(row.split('_')) > 2:
                    if int(row.split('_')[2]) == 1:
                        row = row.replace(row, ">Consensus_sequence_" + str(orde + candiS) + "_" + "1" + "\n")
                    else:
                        row = row.replace(row, ">Consensus_sequence_" + str(orde + candiS) + "_" + "2" + "\n")
                        orde += 1
                else:
                    row = row.replace(row, ">Consensus_sequence_"+str(orde + candiS) + "\n")
                    orde += 1

            fw.write(row)
            idx += 1
        f.stdout.close()
        f.wait()
        fw.close()

    #logging.info("Run BWA for aligning consensus to reference sequence...")
    #SAMout = workDir + fp_prefix + "Alignment.sam"

    #cmd = "bwa index %s 2>%s && bwa mem -B 6 %s %s > %s 2>%s" %(workDir+fp_prefix+"RefSeq.fa", workDir+"BWAinde.log", workDir+fp_prefix+"RefSeq.fa", path, SAMout, workDir+"BWAmem.log")   # bwa -t speed up 

    #subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    logging.info("[runMSA.py] Finish prefix: %s" %(fp_prefix))
    #os.system(cmd)

def run():
    parser = argparse.ArgumentParser(description="MSA examples using abPOA")
    parser.add_argument("--workDir", type=str, default="tempDir",
            help="Temporary working path for Multiple sequence alignment, default:%(default)s")

    parser.add_argument("--freq_for_indel", type=float, default=0.3,
            help="minimum frequency of each haploid, default:%(default)f")

    parser.add_argument("--freq_for_SNP", type=float, default=0.2,
            help="minimum frequency of each haploid, default:%(default)f")

    parser.add_argument("--threads", type=int, default=1,
            help="Number of threads to use, default:%(default)d")

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





