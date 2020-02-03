# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Nabih Estefan

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
import doctest


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # TODO: implement this

    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'G':
        return 'C'

    pass

#Passed
#doctest.run_docstring_examples(get_complement, globals(), verbose=True)


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO: implement this
    complement = ''
    i = len(dna) - 1
    while i >= 0:
        complement += get_complement(dna[i])
        i -= 1

    return complement
    pass

#Passed
#doctest.run_docstring_examples(get_reverse_complement, globals(), verbose=True)


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    # TODO: implement this
    i = 0
    l = len(dna)
    rest = ''
    while i < l-3:
        if(dna[i:i+3] == 'TAG' or dna[i:i+3] == 'TAA' or dna[i:i+3] == 'TGA'):
            break
        rest += dna[i:i+3]
        i += 3
    if(i >= l-3):
        if(dna[i:i+3] != 'TAG' and dna[i:i+3] != 'TAA' and dna[i:i+3] != 'TGA'):
            rest += dna[i:l]

    return rest
    pass

#Passed
#doctest.run_docstring_examples(rest_of_ORF, globals(), verbose=True)


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    # TODO: implement this
    i = 0
    l = len(dna)
    ORFs = []
    singleORF = ''

    while i < l-3:
        if(dna[i:i+3] == 'ATG'):
            singleORF = rest_of_ORF(dna[i:])
            i += len(singleORF) + 3
            ORFs.append(singleORF)
        else:
            i += 3

    return ORFs

    pass
#Passed
#doctest.run_docstring_examples(find_all_ORFs_oneframe, globals(), verbose=True)


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # TODO: implement this
    list = []
    tempo = []
    i = 0

    tempo = find_all_ORFs_oneframe(dna)
    while i < len(tempo):
        list.append(tempo[i])
        i += 1
    i = 0

    tempo = find_all_ORFs_oneframe(dna[1:])
    while i < len(tempo):
        list.append(tempo[i])
        i += 1
    i = 0

    tempo = find_all_ORFs_oneframe(dna[2:])
    while i < len(tempo):
        list.append(tempo[i])
        i += 1
    i = 0

    """
    list.append(find_all_ORFs_oneframe(dna))
    list.append(find_all_ORFs_oneframe(dna[1:]))
    list.append(find_all_ORFs_oneframe(dna[2:]))
    """

    return list

    pass

#Passed
#doctest.run_docstring_examples(find_all_ORFs, globals(), verbose=True)


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """

    list = []
    tempo = []
    i = 0

    revdna = get_reverse_complement(dna)

    tempo = find_all_ORFs(dna)
    while i < len(tempo):
        list.append(tempo[i])
        i += 1
    i = 0


    tempo = find_all_ORFs(revdna)
    while i < len(tempo):
        list.append(tempo[i])
        i += 1
    i = 0
    # TODO: implement this

    return list
    pass

#Passed
#doctest.run_docstring_examples(find_all_ORFs_both_strands, globals(), verbose=True)


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass


"""
if __name__ == "__main__":
    import doctest
    doctest.testmod()
"""
