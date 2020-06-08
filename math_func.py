#!/usr/bin/env python
# coding=utf-8

import math
import numpy as np

_MAX_CONFIDENCE = 1.0 - 1e-16

# rounded to, for numerical stability.
_GL_PRECISION = 10
_QUAL_PRECISION = 1
LOG_10 = math.log(10.0)


def ptrue_to_bounded_phred(pref, ptotal, max_prob=_MAX_CONFIDENCE):
    ratio = float(pref) / ptotal
    prob = math.pow(ratio, pref)
    if not 0 < prob <= 1:
        raise ValueError('ptrue must be between zero and one: {}'.format(ptrue))
    return int(np.around(perror_to_phred(min(prob, max_prob))))

def perror_to_phred(perror):
    return -10 * math.log10(perror)

def ptrue_to_bounded_phred_dv(ptrue, max_prob = _MAX_CONFIDENCE):
    if not 0 <= ptrue <= 1:
        raise ValueError("ptrue must be between zero and one: {}".format(ptrue))

    return perror_to_phred(1.0 - min(ptrue, max_prob))

def computer_quals(predictions, n_alleles):
    """Computes GQ and QUAL values from a set of prediction probabilities."""
    ## GQ is prob(GQ) / prob(all genotypes)
    ## GQ is rounded to the nearest integer to comply with the VCF spec
    index, gt = get_max_index(predictions, n_alleles)
    gq = int(np.around(ptrue_to_bounded_phred_dv(predictions[index])))
    qual = ptrue_to_bounded_phred_dv(min(sum(predictions[1:]), 1.0))
    rounded_qual = round(qual, _QUAL_PRECISION)

    return gq, rounded_qual, gt

def get_max_index(predictions, n_alleles = 2):
    index_of_max = np.argmax(predictions)

    index = 0
    for h1 in range(0, n_alleles):
        for h2 in range(0, h1 + 1):
            if index == index_of_max:
                return index, [h2, h1]
            index += 1
    raise ValueError("No corresponding genotype for predictions", predictions)


def round_gls(gls, precision = _GL_PRECISION):
    """Returns genotype likelihoods rounded to the desired precision level.
    
    Returns genotype likelihoods rounded to the desired precision level.
    
    Args:
        gls: A list of floats. The input genotype likelihoods at any precision.
        precision: Positive int. The number of places past the decimal point to
                round to. If None, no rounding is performed.
    """
    if abs(sum(gls) - 1) > 1e-6:
        raise ValueError('Invalid genotype likelihoods do not sum to one: sum({}) = {}'.format(gls, sum(gls)))

    min_ix = 0
    min_gl = gls[0]
    for ix, gl in enumerate(gls):
        if gl < min_gl:
            min_gl = gl
            min_ix = ix

    rounded_gls = [round(gl, precision) for gl in gls]

    rounded_gls[min_ix] = max(
            0.0,
            round(1 - sum(rounded_gls[:min_ix] + rounded_gls[min_ix + 1:]), 
            precision))

    return rounded_gls

def toRealSpace(log10_probs):
    """
    sum(10^log10_probs) close to one
    """
    return [math.pow(10.0, x) for x in log10_probs]

def log10sumexp(log10_probs):
    """Returns log10(sum(10^log10_probs)) computed in a numerically-stable way"""

    m = max(log10_probs)

    return m + math.log10(sum(pow(10.0, x-m) for x in log10_probs))

def normalize_log10_probs(log10_probs):
    """
    Approximately normalizes log10 probabilities
    """
    #log10_probs = [math.log10(p) for p in probs]
    log10_probs = np.array(log10_probs)
    if np.max(log10_probs) > 0.0:
        raise ValueError("log10_probs all must be <=0", log10_probs)

    lse = log10sumexp(log10_probs)

    return np.minimum(log10_probs - lse, 0.0)

def calc_reference_confidence(n_ref, n_total, p_error):
    if n_ref < 0:
        raise ValueError("n_ref={} must be >=0".format(n_ref))
    if n_total < n_ref:
        raise ValueError("n_total={} must be >= n_ref={}".format(n_total, n_ref))
    n_alts = n_total - n_ref
    logp = math.log(p_error) / LOG_10
    log1p = math.log1p(-p_error) / LOG_10

    p_ref = n_ref * log1p + n_alts * logp
    p_het = -n_total * math.log(2) / LOG_10
    p_hom_alt= n_ref *logp + n_alts * log1p

    #print [p_ref, p_het, p_hom_alt]
    #return [p_ref, p_het, p_hom_alt]

    return normalize_log10_probs([p_ref, p_het, p_hom_alt])

def rescale_read_counts(n_ref, n_total, max_allowed_reads=100):
    """Ensures that n_total <= max_allowed_reads, rescaling if necessary."""
    if n_total > max_allowed_reads:
        ratio = n_ref / (1.0 * n_total)
        n_ref = int(math.ceil(ratio * max_allowed_reads))
        n_total = max_allowed_reads

    return n_ref, n_total
