#program to find HET sites from VCF
import argparse

def main(fn_vcf, fn_output):
    #t = time.time()
    file = open(fn_vcf, 'r')
    f = open(fn_output, 'a+')

    index = -1
    toAdd = []
    should = True
    for line in file:
        if should and line.startswith("##"):        # "##" is the header
            should = should
            f.write(line)
        else:
            if line.startswith("#"):
                categories = line.split()
                for i in range(len(categories)):
                    if categories[i] == 'NA12878':  # 'NA12878' is the last item but seems to be not universal
                        index = i
                toAdd = [categories[i] for i in range(8)]  # leaves only 9 columns
                toAdd.append(categories[index])
                f.write(line)
            else:
                spl = line.split()

                if len(set(spl[index].split("|")) )>1:  # there are differences in the 0|0 entry
                    toAdd = [spl[i] for i in range(8)]
                    toAdd.append(spl[index])
                    f.write("\t".join(toAdd))
                    f.write("\n")
        #print("time: ", time.time() - t)
    f.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf file for chromosomes')
    parser.add_argument('-o', '--out', help='output vcf file with HET sites')
    args = parser.parse_args()
    fn_vcf = args.vcf
    fn_output = args.out
    print('vcf', fn_vcf)
    print('output', fn_output)
    main(fn_vcf, fn_output)
