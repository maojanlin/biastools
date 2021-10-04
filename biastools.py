import argparse
import sys
from subprocess import call

def main(fn_vcf, fn_fas, fn_id, fn_n, fn_c):


    for command in (#"bcftools view -i \'GT~\"1|.\"\' -o hapA.vcf " + fn_vcf,
                    #"bcftools view -i \'GT~\".|1\"\' -o hapB.vcf " + fn_vcf,
                    "mason_simulator -iv hapA.vcf -ir " + fn_fas + " -o hapA1.fq -or hapA2.fq -oa hapA.sam -n " + fn_n,
                    #"mason_simulator -iv hapB.vcf -ir " + fn_fas + " -o hapB1.fq -or hapB2.fq -oa hapB.sam  -n " + fn_n,
                    #"samtools sort hapA.sam -o sorted_hapA.sam",
                    #"python VCF_Processing.py -v hapA.vcf -o hapA_het.vcf",
                    "python3 ref_bi.py -v hapA_het.vcf -s sorted_hapA.sam -f " + fn_fas + " -o hapA_ref_bi.txt"):#,
                    #"ipython create_ref_bi_graph.py"):
        toDisplay = "I am about to run " + command
        #print(toDisplay, file = sys.stdout)
        call(command, shell=True)

    #message = "For each haplotype, two fastq files with simulated reads were produced and a SAM with alignments to reference genome."
    #print(message, file=sys.stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help = 'vcf file')
    parser.add_argument('-f', '--fas', help='reference genome fasta file')
    parser.add_argument('-i', '--id', help='Individual id')
    parser.add_argument('-n', '--numReads', help='number of reads to generate')
    parser.add_argument('-c', '--coverage', help='average coverage for simulated reads')

    # What is the best way to implement the "stretch goal" - allowing the user to specify additional parameters
    # that are passed through to mason_simulator

    args = parser.parse_args()
    fn_vcf = args.vcf
    fn_fas = args.fas
    fn_id = args.id
    fn_n = args.numReads
    fn_c = args.coverage

    print("fn_vcf: ", fn_vcf)
    print("fn_fas: ", fn_fas)
    print("fn_id: ", fn_id)
    print("fn_n: ", fn_n)
    print("fn_c: ", fn_c)

main(fn_vcf, fn_fas, fn_id, fn_n, fn_c)
