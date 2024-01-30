#program to find HET sites from VCF
import argparse
import pysam

def parse_het_site(fn_vcf, fn_output):
    in_vcf_file  = pysam.VariantFile(fn_vcf, 'r')
    out_vcf_file = pysam.VariantFile(fn_output, 'w', header=in_vcf_file.header)
    for segment in in_vcf_file:
        #hap_info = str(segment).split()[9].split('|') # "0|0", "1|0", "0|1" tag
        #if hap_info[0] != hap_info[1]:
        phase_info = segment.samples[0]['GT']
        if len(phase_info) != 2:
            print("WARNING! non diploid haplotype.")
            continue
        hap_0, hap_1 = phase_info
        if hap_0 == None or hap_1 == None:
            print("WARNING! one haplotype information is missing.")
            continue
        if hap_0 + hap_1 != 0:
            out_vcf_file.write(segment)
    in_vcf_file.close()
    out_vcf_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf/vcf.gz file for chromosomes')
    parser.add_argument('-o', '--out', help='output vcf.gz file with HET sites')
    args = parser.parse_args()
    
    fn_vcf = args.vcf
    fn_output = args.out
    
    parse_het_site(fn_vcf, fn_output)
