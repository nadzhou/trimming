from Bio.Emboss.Applications import WaterCommandline
import subprocess


def emboss_water(seq_a_file: str, seq_b_file: str, out_file: str):
    """ Do a global pairwise alignment using EMBOSS

        Args: 
            seq_a_file: First sequence
            seq_b_file: second sequence
            out_file: Output file

        Returns: 
            r [subprocess object]: Execute the commandline command for EMBOSS
        
    """
    needle_cline = WaterCommandline(asequence=seq_a_file,
                                        bsequence=seq_b_file,
                                        outfile=out_file,
                                        verbose=True,
                                        gapextend=1,
                                        gapopen=10)

    cmd = str(needle_cline)
    cmd = cmd.split(" ")
    cmd.append("-aformat=msf")

    r = subprocess.Popen(cmd)
    if r.communicate():
        print("Global alignment done.")


def main():
    seq_a = "polished_core_1.fasta"
    seq_b = "polished_core_2.fasta"
    out_file = "water.fasta"

    emboss_water(seq_a, seq_b, out_file)


if __name__ == '__main__':
    main()
