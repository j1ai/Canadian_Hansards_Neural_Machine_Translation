# Copyright 2020 University of Toronto, all rights reserved

'''Calculate BLEU score for one reference and one hypothesis

You do not need to import anything more than what is here
'''

from math import exp  # exp(x) gives e^x


def grouper(seq, n):
    '''Extract all n-grams from a sequence

    An n-gram is a contiguous sub-sequence within `seq` of length `n`. This
    function extracts them (in order) from `seq`.

    Parameters
    ----------
    seq : sequence
        A sequence of words or token ids representing a transcription.
    n : int
        The size of sub-sequence to extract.

    Returns
    -------
    ngrams : list
    '''
    ngrams = []
    for i in range(len(seq)):
        if i + (n-1) < len(seq):
            ngrams.append(seq[i:i+(n-1)])
    return ngrams


def n_gram_precision(reference, candidate, n):
    '''Calculate the precision for a given order of n-gram

    Parameters
    ----------
    reference : sequence
        The reference transcription. A sequence of words or token ids.
    candidate : sequence
        The candidate transcription. A sequence of words or token ids
        (whichever is used by `reference`)
    n : int
        The order of n-gram precision to calculate

    Returns
    -------
    p_n : float
        The n-gram precision. In the case that the candidate has length 0,
        `p_n` is 0.
    '''
    if len(candidate) == 0:
        return 0.0
    ref_ngrams = grouper(reference, n)
    cand_ngrams = grouper(candidate, n)
    #No capping from tutorial slides
    #Only 1 reference and 1 candidate at a time based on tutorial slides
    count_cand_in_ref = 0
    for cand_ngram in cand_ngrams:
        if cand_ngram in ref_ngrams:
            count_cand_in_ref += 1
    return float(count_cand_in_ref) / len(cand_ngrams)
    
    


def brevity_penalty(reference, candidate):
    '''Calculate the brevity penalty between a reference and candidate

    Parameters
    ----------
    reference : sequence
        The reference transcription. A sequence of words or token ids.
    candidate : sequence
        The candidate transcription. A sequence of words or token ids
        (whichever is used by `reference`)

    Returns
    -------
    BP : float
        The brevity penalty. In the case that the candidate transcription is
        of 0 length, `BP` is 0.
    '''
    if len(candidate) == 0:
        return 0.0
    ref_length = len(reference)
    cand_length = len(candidate)
    brevity = float(ref_length) / cand_length
    if cand_length <= ref_length:
        return exp(1 - brevity)
    else:
        return 1.0


def BLEU_score(reference, hypothesis, n):
    '''Calculate the BLEU score

    Parameters
    ----------
    reference : sequence
        The reference transcription. A sequence of words or token ids.
    candidate : sequence
        The candidate transcription. A sequence of words or token ids
        (whichever is used by `reference`)
    n : int
        The maximum order of n-gram precision to use in the calculations,
        inclusive. For example, ``n = 2`` implies both unigram and bigram
        precision will be accounted for, but not trigram.

    Returns
    -------
    bleu : float
        The BLEU score
    '''
    brevity = brevity_penalty(reference, hypothesis)
    n_gram_precisions = 1
    for num in range(1,n+1):
        n_gram_precisions = n_gram_precisions * n_gram_precision(reference, hypothesis, num)
    return brevity * (n_gram_precisions ** (1/n))
        
    
