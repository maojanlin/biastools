import argparse
import pysam


def main():
    parser = argparse.ArgumentParser(description="Generate a BED file from a VCF file.")
    parser.add_argument('-v', '--vcf', help='list of vcf files for input', required=True)
    parser.add_argument('-o', '--out', help='output bed file', required=True)
    args = parser.parse_args()

    vcf_path = args.vcf
    vcf = pysam.VariantFile(vcf_path)
    fo = open(args.out, 'w')
    for record in vcf:
        chrom = record.chrom
        start = record.start
        end = record.stop
        fo.write(f"{chrom}\t{start}\t{end}\n")
    fo.close()


if __name__ == "__main__":
    main()
