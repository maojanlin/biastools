import argparse
import numpy as np


def read_list_bed(list_bed):
    dict_chr_bias = {}
    for fn_bed in list_bed:
        f = open(fn_bed)
        f.readline()
        for line in f:
            fields = line.split()
            contig = fields[0]
            start = int(fields[1])
            stop  = int(fields[2])
            if dict_chr_bias.get(contig):
                dict_chr_bias[contig].append((start, stop))
            else:
                dict_chr_bias[contig] = [(start, stop)]
    return dict_chr_bias


def compare_bias_regions(dict_target, dict_improve, dict_lowRd, out_file):
    assert sorted(dict_target.keys()) == sorted(dict_improve.keys()), "discrepancy on the reference of the two lists"
    f_o = open(out_file, 'w')
    f_o.write('#chrom\tchromStart\tchromEnd\tname(%;initial;improve;lowRd)\n')

    total_100 = []
    total_75 = []
    total_50 = []
    total_25 = []
    total_under_25 = []
    for contig in sorted(dict_target.keys()):
        local_100 = []
        local_75 = []
        local_50 = []
        local_25 = []
        local_under_25 = []

        region_target  = dict_target [contig]
        region_improve = dict_improve[contig]
        region_lowRd   = dict_lowRd  [contig]
        idx_2 = 0
        idx_3 = 0
        for region in region_target:
            start_1, stop_1 = region
            #if stop_1 - start_1 < 1000:
            #    continue
            contain_region_2 = []
            contain_lowRd    = []
            for idx in range(idx_2, len(region_improve)):
                start_2, stop_2 = region_improve[idx]
                if stop_2 < start_1:
                    continue
                elif start_2 < stop_1:
                    contain_region_2.append(region_improve[idx])
                else:
                    idx_2 = idx-1
                    break
            for idx in range(idx_3, len(region_lowRd)):
                start_3, stop_3 = region_lowRd[idx]
                if stop_3 < start_1:
                    continue
                elif start_3 < stop_1:
                    contain_lowRd.append(region_lowRd[idx])
                else:
                    idx_3 = idx-1
                    break

            len_region_1 = stop_1 - start_1
            len_region_2 = sum([ele[1]-ele[0] for ele in contain_region_2])
            len_region_3 = sum([ele[1]-ele[0] for ele in contain_lowRd])
            improve_len = len_region_1 - len_region_2 - len_region_3
            if improve_len == len_region_1:
                local_100.append(region)
                f_o.write(contig + '\t' + str(start_1) + '\t' + str(stop_1) + '\t' + '100;' + str(len_region_1) + ';' + str(len_region_2) + ';' + str(len_region_3) + '\n')
            elif improve_len >= len_region_1*0.75:
                local_75.append(region)
                f_o.write(contig + '\t' + str(start_1) + '\t' + str(stop_1) + '\t' + '75;' + str(len_region_1) + ';' + str(len_region_2) + ';' + str(len_region_3) + '\n')
            elif improve_len >= len_region_1*0.5:
                local_50.append(region)
                f_o.write(contig + '\t' + str(start_1) + '\t' + str(stop_1) + '\t' + '50;' + str(len_region_1) + ';' + str(len_region_2) + ';' + str(len_region_3) + '\n')
            elif improve_len >= len_region_1*0.25:
                local_25.append(region)
                f_o.write(contig + '\t' + str(start_1) + '\t' + str(stop_1) + '\t' + '25;' + str(len_region_1) + ';' + str(len_region_2) + ';' + str(len_region_3) + '\n')
            else:
                local_under_25.append(region)
        total_100 += local_100
        total_75  += local_75
        total_50  += local_50
        total_25  += local_25
        total_under_25 += local_under_25
        print(contig, len(local_100), len(local_75), len(local_50), len(local_25), len(local_under_25))
    f_o.close()
    len_total = sum([len(total_100), len(total_75), len(total_50), len(local_25), len(total_under_25)])
    print(len(total_100), round(len(total_100)/len_total,3), \
          len(total_75),  round(len(total_75)/len_total,3),  \
          len(total_50),  round(len(total_50)/len_total,3),  \
          len(total_25),  round(len(total_25)/len_total,3),  \
          len(total_under_25), round(len(total_under_25)/len_total,3))







if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-lt',  '--list_target',  nargs='+', required=True, help='the first list of scanning bias bed report')
    parser.add_argument('-li',  '--list_improve', nargs='+', required=True, help='the second list of scanning bias bed report, the region should contain in list 1')
    parser.add_argument('-lrd', '--list_lowRd',   nargs='+', required=True, help='the second list of scanning bias bed report, the region should contain in list 1')
    parser.add_argument('-out', '--output_improve', help="output the improve regions")
    args = parser.parse_args()

    list_target  = args.list_target
    list_improve = args.list_improve
    list_lowRd   = args.list_lowRd
    out_file   = args.output_improve

    dict_target  = read_list_bed(list_target)
    dict_improve = read_list_bed(list_improve)
    dict_lowRd   = read_list_bed(list_lowRd)

    compare_bias_regions(dict_target, dict_improve, dict_lowRd, out_file)

