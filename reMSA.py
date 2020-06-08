#!/usr/bin/env python
# coding=utf-8

import subprocess

def reCall(workDir, prefix, outliers, freq):
    outlier_fileList = workDir+prefix+"ID_outliers_fileList"
    pfix = workDir + prefix + "ID_Region_"
    with open(outlier_fileList, 'w') as fw:
        for out in outliers:
            f = pfix + str(out[1]) + ".fa"
            fw.write(f+"\n")

    #re MSA
    fw = open(workDir + prefix + "ID_outliers_MSA.fa", 'w')
    cmd = "/home/ydliu/abPOA/bin/abPOA -c -m 2 -f %f -l %s" %(freq, outlier_fileList)
    f = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    idx = 0; orde = 0
    for row in f.stdout:
        if idx % 2 == 0:
            if len(row.split('_')) > 2:
                if int(row.split('_')[2]) == 1:
                    row = row.replace(row, ">Consensus_sequence_" + str(outliers[orde][1]) + "_" + "1" + "\n")
                else:
                    row = row.replace(row, ">Consensus_sequence_" + str(outliers[orde][1]) + "_" + "2" + "\n")
                    orde += 1
            else:
                row = row.replace(row, ">Consensus_sequence_"+str(outliers[orde][1]) + "\n")
                orde += 1

        fw.write(row)
        idx += 1
    f.stdout.close()
    f.wait()   
    fw.close()

    return workDir + prefix + "ID_outliers_MSA.fa"


