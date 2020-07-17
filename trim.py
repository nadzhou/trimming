
import numpy as np
from Bio import AlignIO
from Bio import SeqIO

from identical_sequence_parser import IdenticalSequencesParser
from emboss import emboss_water



def find_best_hits(aligned_file): 
    water_file = AlignIO.parse(aligned_file, "msf")
    result_record = []

    for aln in water_file: 
        parsed = IdenticalSequencesParser(aln[1], aln[0])

        result = parsed.highly_identical_seqs()

        if result: 
            result_record.append(result)
    return result_record



def write_best_hit_list(result_record): 
    return SeqIO.write(result_record, "best_core_seqs.fasta", "fasta")


def main(): 
    seq_file_1 = input("Input sequence file 1: ")
    seq_file_2 = input("Enter sequence file 2: ")


    print("\nInitiating local alignment...")
    emboss_water(seq_file_1, seq_file_2, "water.fasta")
    print("\n Local alignment: done.")

    print("Intiating trimming...")
    result = find_best_hits("water.fasta")

    print("\nHits now being written to file...")
    write_best_hit_list(result)
    print("Fasta file writing: done.")



if __name__ == '__main__': 
    main()