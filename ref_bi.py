import argparse
import re
import pickle
import os.path
from os import path


def main(fn_vcf, fn_sam, fn_fasta, fn_output):
    ref_file = open(fn_fasta, 'r')
    reference = ''
    for line in ref_file:
        if not line.startswith('>'):
            reference += line.strip()

    file = open(fn_vcf, 'r')
    count_line = 0
    index_ref = 0
    index_pos = 0
    index_alt = 0
    index_cat = 0

    vcf_file_name = fn_output+'.chr_vcf.pickle'
    sam_file_name = fn_output+'.chr_sam.pickle'

    f = open(fn_output, 'w')
    
    # read the fn_vcf file, build the pickle file it is not exist
    count_het = 0
    if not path.exists(vcf_file_name):
        chr_vcf = {}
        for line in file:
            #if not line.startswith("##") and line.startswith("#"): # the "#" category line
            if line.startswith("#"): # the "#" category line
                if line.startswith("##"):
                    pass
                else:
                    categories = line.split()
                    index_cat = count_line
                    for i in range(len(categories)):
                        if categories[i] == 'REF':
                            index_ref = i
                        elif categories[i] == 'ALT':
                            index_alt = i
                        elif categories[i] == 'POS':
                            index_pos = i
            elif count_line > index_cat and index_cat != 0:
                now = line.split()
                #print("line: ", line)
                if line == '\n':
                    if not chr_vcf:
                        chr_vcf[chr] = [[],[]]#this works since we are doing individual chromosomes
                    else:
                        chr_vcf[chr][0].append('N/A')
                        chr_vcf[chr][1].append(['N/A', 'N/A'])
                    continue
                # chr = int(now[0])
                chr = now[0]
                if not chr_vcf:
                    chr_vcf[chr] = [[], []]
                elif chr not in chr_vcf.keys():
                    chr_vcf[chr] = [[], []]

                #print(line)
                if (len(now[index_alt]) == 1 and len(now[index_ref]) == 1):
                    chr_vcf[chr][0].append(int(now[index_pos]) - 1)
                    chr_vcf[chr][1].append([now[index_ref], now[index_alt]])
                elif ',' in now[index_alt]:
                    possib = now[index_alt].split(',')
                    chr_vcf[chr][0].append(int(now[index_pos]) - 1)
                    chr_vcf[chr][1].append([now[index_ref], ''.join(possib)])
            count_line += 1
        #pickle.dump(chr_vcf, open(vcf_file_name, 'wb'))
        print ('Dump to {}'.format(vcf_file_name))
    else:
        #print ('Load from {}'.format(vcf_file_name))
        chr_vcf = pickle.load(open(vcf_file_name, "rb"))

    file = open(fn_sam, 'r')
    f.write("CHR\tHET_SITE\tREFERENCE_BIAS\tREF_COUNT\tALT_COUNT\tGAP_COUNT\tOTHER_COUNT\tNUM_READS\tSUM_MAPQ\tREAD_DISTRIBUTION")#modified
    f.write("\n")
    count_line = 0
    chr_sam = {}

    if not path.exists(sam_file_name):
        for line in file:
            if not line.startswith('@'):
                spl = line.split() 
                tag = int(spl[1])
                if (tag & 4):
                    continue
                chr = spl[2]
                #: ignores haplotype suffixes
                if chr.endswith('A') or chr.endswith('B'):
                    chr = chr[:-1]
                # chr = spl[2][0:2]
                # chr = int(spl[2][0:2])
                cigar = spl[5]
                mapq = int(spl[4])
                start_pos = int(spl[3]) - 1 #TODO
                sequence = spl[9]
                rg_tag = spl[-1]# added this okay
                #print("rg tag is ")
                #print(rg_tag)
                #print("\n")
                mod_sequence = ''
                #: ignores unmapped reads
                #: tag 4: unaligned
                if tag & 4:
                    continue
                if not cigar == (str(len(sequence))+'M') and not (tag & 4):
                    change = 0
                    start = start_pos
                    count_del = 0
                    count_ins = 0
                    for num1, idm in re.findall('(\d+)([IDMS])', cigar):
                        # print(start_pos)
                        if idm == 'M':
                            mod_sequence += sequence[change:change + int(num1)]
                        elif idm == 'D':
                            count_del += int(num1)
                            for i in range(int(num1)):
                                mod_sequence += '-'
                        elif idm == 'I':
                            count_ins += int(num1)
                        elif idm == 'S':
                            count_ins += int(num1)
                        else:
                            print ('error: unexpected cigar letter', cigar)
                            exit ()

                        if idm != 'D':
                            change += int(num1)
                    #print("tag: {0}, cigar: {1}".format(tag, cigar))
                    ref = reference[start_pos:start_pos + len(sequence) + count_del - count_ins]
                    try:
                        assert len(ref) == len(mod_sequence)
                    except:
                        if len(ref) == 0:
                            continue
                        else:
                            print("read name: ", spl[0])
                            print ("ref      seq", len(ref))
                            print ("modified seq", mod_sequence)
                            print ("start_pos", start_pos)
                            print ("count_ins", count_ins)
                            print ("count_del", count_del)
                            print ("cigar", cigar)
                            print ("flag", tag) 
                            exit()
                else:
                    mod_sequence = sequence
                if not chr_sam:
                    chr_sam[chr] = [[], [], [], [], [], []]
                elif chr not in chr_sam.keys():
                    chr_sam[chr] = [[], [], [], [], [], []]
                chr_sam[chr][0].append(start_pos) #position
                chr_sam[chr][1].append(mod_sequence) #sequence
                chr_sam[chr][2].append(tag)
                chr_sam[chr][3].append(cigar)
                chr_sam[chr][4].append(mapq)
                chr_sam[chr][5].append(rg_tag)#added this 
            count_line += 1
        #pickle.dump(chr_sam, open(sam_file_name, 'wb'))
        print ('Dump to {}'.format(sam_file_name))
    else:
        #print ('Load from {}'.format(sam_file_name))
        chr_sam = pickle.load(open(sam_file_name, 'rb'))

    chr_list = list(chr_vcf.keys())
    chr_list.sort()
    #count = 0
    for chr in chr_list:
        het_site_list = chr_vcf[chr][0]
        options = chr_vcf[chr][1]
        sam_pos = chr_sam[chr][0]
        sam_reads = chr_sam[chr][1]
        list_mapq = chr_sam[chr][4]
        rg_tag = chr_sam[chr][5]#added this
        have_started = False
        starting_point = 0
        count_pos = 0

        total_ref_count = 0
        total_alt_count = 0
        total_gap_count = 0
        total_other_count = 0
        hap_map = ref_hap_map(fn_vcf)

        for pos in het_site_list:
            #count+=1
            #if count > 5:
            #    break
            ref_count = 0
            alt_count = 0
            gap_count = 0
            other_count = 0
            reads_at_het = 0
            sum_mapq = 0
            # reference_hap = find_ref_hap(pos, fn_vcf)#added this
            reference_hap = hap_map[pos] # replace reading every time to a map
            count_a = 0.0#added this
            count_b = 0.0#added this
            read_dis = 0.0
            if pos == 'N/A':
                f.write('\n')
                count_pos += 1
                continue
            for i in range(starting_point, len(sam_pos)):
                align = sam_pos[i]
                ran = range(align, align + len(sam_reads[i]))
                
                #print(rg_tag[i])
#                if 'hapA' in rg_tag[i]:#added this                                                                                                       
#                    count_a += 1.0
#                elif 'hapB' in rg_tag[i]:       
#                    count_b += 1.0#up to here

                if pos in ran:
                    #print("in here. pos is: ", pos, "   and align is: ", align)
                    #print("sam[align]: ", sam_reads[i])
                    #print("rg tag is: ", rg_tag[i])
                    #print("here, it overlaps")
                    if not have_started:
                        have_started = True
                        starting_point = i
                    reads_at_het += 1
                    sum_mapq += list_mapq[i]
                    if 'hapA' in rg_tag[i]:#added this                                                                                                                                     
                        count_a += 1.0
                    elif 'hapB' in rg_tag[i]:
                        count_b += 1.0#up to here 
                    try:
                        allele = sam_reads[i][pos - align]
                        #print("allele: ", allele)
                        if allele in options[count_pos][0]:
                            ref_count += 1
                            total_ref_count += 1
                        elif allele in options[count_pos][1]:
                            alt_count += 1
                            total_alt_count += 1
                        elif allele == '-':
                            gap_count += 1
                            total_gap_count += 1
                        else:
                            other_count += 1
                            total_other_count += 1
                    except:
                        print("in here. pos is: ", pos, "   and align is: ", align)
                        print("sam[align]: ", sam_reads[i])
                else:
                    if align > pos:
                        break
            have_started = False
            f.write(str(chr))
            f.write("\t")
            f.write(str(pos + 1)) #: outputs in 1-based format
            # f.write(str(pos))
            f.write("\t")
            if ref_count + alt_count == 0:
                f.write("N/A")
            else:
                f.write(format(ref_count / float(ref_count + alt_count), '.5f')) 
            f.write(f'\t{ref_count}\t{alt_count}\t{gap_count}\t{other_count}\t{reads_at_het}\t{sum_mapq}\t')#modified this \t vs \n
            #f.write("\t" + str(ref_count) + "\t" + str(alt_count) + "\t" + str(gap_count) + "\t" + (other_count) + "\t" + str(reads_at_het) + "\t" + str(sum_mapq) + "\t")
            if count_a + count_b == 0:#added this up to 
                f.write("N/A")
            else:
                if reference_hap =='hapA':
                    f.write(format((count_a)/(count_a+count_b), '.5f'))
                else:
                    f.write(format((count_b)/(count_a+count_b), '.5f'))
            f.write("\n")#up to here
            count_pos += 1

    f.close()


def find_ref_hap(het_site, fn_vcf):
    file_in = open(fn_vcf, 'r')
    for line in file_in:
        if line.startswith("#"):
            continue
        else:
            spl = line.split()
            if int(spl[1]) == het_site+1:#+1 because 0-base to 1-base                                                                                                            
                hap = spl[9].split("|")
                if hap[0] == '0':
                    return 'hapA'
                else:
                    return 'hapB'
    print("error, het site not found")


def ref_hap_map(fn_vcf):
    """
    dict_hap_map{}
    - key: het_site
    - value: hapA or hapB
    """
    dict_hap_map = {}
    file_in = open(fn_vcf, 'r')
    for line in file_in:
        if line.startswith("#"):
            continue
        else:
            spl = line.split()
            site = int(spl[1]) - 1 # transform the 1-base into 0-base
            if dict_hap_map.get(site):
                continue
            if spl[9][0] == '0':
                dict_hap_map[site] = 'hapA'
            else:
                dict_hap_map[site] = 'hapB'
    file_in.close()
    return dict_hap_map



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf file')
    parser.add_argument('-s', '--sam', help='sam file')
    parser.add_argument('-f', '--fasta', help='reference fasta file')
    parser.add_argument('-o', '--out', help='output file')
    args = parser.parse_args()
    fn_vcf = args.vcf
    fn_sam = args.sam
    fn_fasta = args.fasta
    fn_output = args.out
    print('Input vcf:', fn_vcf)
    print('Input sam:', fn_sam)
    print('Input fasta:', fn_fasta)
    print('Output file:', fn_output)
    main(fn_vcf, fn_sam, fn_fasta, fn_output)
