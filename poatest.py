#!/usr/bin/env python
# coding=utf-8
from poa_module import POA
import sys

fp = sys.argv[1]
Rlist = list()
with open(fp, 'r') as f:
    for line in f:
        if line.startswith('>'):
            continue
        Rlist.append(bytes(line.strip(), 'gbk'))

consensus = POA(Rlist, 2, 0.2)
print(consensus)


