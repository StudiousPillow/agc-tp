#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    amplicon_file = isfile(amplicon_file)
    with gzip.open(amplicon_file, "rb") as handle:
        lines = iter(handle.readlines())
        seq = ''
        for line in lines:
            line = line.decode('utf-8').strip()
            if line.startswith('>'):
                if len(seq)>=minseqlen:
                    yield seq
                seq = ''
            else:
                seq = seq + str(line)
        if len(seq)>=minseqlen:
            yield seq


def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    sequences = (read_fasta(amplicon_file, minseqlen))
    count = Counter(sequences)
    count = count.most_common()
    for item in count: ## item is [sequence, count]
        if item[1]>=mincount:
            yield list(item)
        

def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    nb_identique = 0
    len_align = len(alignment_list[0])
    for base_idx in range(len_align):
        if alignment_list[0][base_idx]==alignment_list[1][base_idx]:
            nb_identique += 1
    id = nb_identique/len_align*100
    return(id)

def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int, kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    sequences = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    OTU = []
    n = 0
    for idxseq in range(len(sequences)):
        threshold = sequences[idxseq][1]
        # print("t", threshold)
        for idxseq2 in range(len(sequences[idxseq:])):
            if sequences[idxseq2][1]>threshold: ## si plus abondante
                # print("abondante")
                alignement = nw.global_align(sequences[idxseq][0],
                                             sequences[idxseq2][0]) ## returns a set with the two aligned sequences : (A-GT,ACG-)
                if not get_identity(list(alignement))>97:
                    
                    OTU.append(sequences[idxseq])
                    n+=1
                    print(n, end="\r")
                    
    return OTU
                    

def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    n=0
    with open(output_file, 'w') as handle:
        for otu in OTU_list:
            n+=1
            handle.write(f'>OTU_{n} occurrence:{otu[1]}\n')
            handle.write(textwrap.fill(otu[0],width=80))
            if n != len(OTU_list):
                handle.write('\n')

#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    # args = get_arguments()
    # Votre programme ici
    # print(read_fasta("data/amplicon.fasta.gz",200))
    
    # seq1 = "ACTACGGGGCGCAGCAGTAGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGTATGAAGAAGGTTTTCGGATCGTAAAGTACTGTTGTTAGAGAAGAACAAGGATAAGAGTAACTGCTTGTCCCTTGACGGTATCTAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAGTTAGTGGGCGTAAAGCGCGCGCAGGCGGTCTTTTAAGTCTGATGTCAAAGCCCCCGGCTTAACCGGGGAGGGTCATTGGAAACTGGAAGACTGGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGATATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAA"
    # seq2 = "TAGGGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTCTTCGGATCGTAAAACTCTGTTATTAGGGAAGAACATATGTGTAAGTAACTGTGCACATCTTGACGGTACCTAATCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTACAGCGCG"
    # print(seq1)
    # print("--")
    # OTU_list = abundance_greedy_clustering('tests/test_sequences.fasta.gz', 200, 3, 50, 8)
    # print(seq1 in OTU_list)
    
    otu = [("TCAGCGAT", 8), ("TCAGCGAA", 8), ("ACAGCGAT", 8), ("ACAGCGAA", 8)]
    write_OTU(otu, 'results/test.fasta')
    import hashlib
    with open('results/test.fasta', 'rb') as f:
        print(hashlib.md5(f.read()).hexdigest())
        print('0a7caf3d43ba5f0c68bc05cb74782dbb')
    with open('results/test.fasta', 'r') as f:
        print(f.read())
    with open('results/test.fasta', 'r') as f:
        print(repr(f.read()))
        # print(f.read())
    # write_OTU(OTU_list, "results/OTU_found.fasta")
    

    
if __name__ == '__main__':
    main()
