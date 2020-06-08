#!/usr/bin/env python
# coding=utf-8

class variantsCan:
    def __init__(self, value):
        self.SeqID = ''
        self.chr = value[0]
        self.pos = value[1]   # variant pos
        self.type = value[2] # variant type
        self.REF = value[3]  # ref base
        self.ALT = value[4]  # Alt base
        self.qname = value[5]  # variant from 
        self.MQ = value[6] # mapping quality
        self.varlen = value[7] # variant length
        self.supportR = value[8] # reads supported reference
        self.supportT = value[9] # reads supported allel
        self.quePos = value[10]   # variant pos in reads
        self.refPos = value[11]   # variant pos in ref

        self.qual = 0.0
        self.GT = ""
        self.ref = ''
        self.alt = []

    def set_ref(self, ref):
        self.ref = ref

    def set_alt(self, alt):
        self.alt.append(alt)

    def set_seqID(self, sid):
        self.SeqID = sid

    def set_GT(self, gt):
    	self.GT = gt

    def set_qual(self, qual):
        self.qual = qual

class hapSeq:
    def __init__(self, nid):
        self.nid = nid
        self.refHap = None
        self.alleleHap = []
    def set_ref_hap(self, refhap):
        self.refHap = refhap
    def set_allele_hap(self, allelehap):
        self.alleleHap.append(allelehap)
