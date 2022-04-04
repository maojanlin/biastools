import argparse


def merge_and_report(
    fn_bias     :str,
    fn_golden   :str,
    fn_output   :str
    ) -> None:
    """
    Comment
    """
    f_b = open(fn_bias, "r")
    f_g = open(fn_golden, "r")
    f_o = open(fn_output, "w")
    headline = f_b.readline().strip()
    goldline = f_g.readline().split()
    gold_idx = 2
    ref_idx  = 3
    alt_idx  = 4
    headline += '\t' + goldline[gold_idx] + '\tGOLDEN_' + goldline[ref_idx] + '\tGOLDEN_' + goldline[alt_idx] + '\n'
    f_o.write(headline)
    for line in f_b:
        fields = line.split()
        goldline = f_g.readline().split()
        assert(fields[0] == goldline[0])
        assert(fields[1] == goldline[1])
        f_o.write(line.strip() + '\t' + goldline[gold_idx] + '\t' + goldline[ref_idx] + '\t' + goldline[alt_idx] + '\n')
    f_b.close()
    f_g.close()
    f_o.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bias_report',   help='bias report file')
    parser.add_argument('-g', '--golden_report', help='golden report file')
    parser.add_argument('-o', '--out', help='merged report file')
    args = parser.parse_args()
    
    fn_bias   = args.bias_report
    fn_golden = args.golden_report
    fn_output = args.out

    merge_and_report(fn_bias, fn_golden, fn_output)
