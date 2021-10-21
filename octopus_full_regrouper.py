import argparse
import pysam


def regroup_sam(fn_primary, fn_secondary, fn_output):
    print("Start recording the sam file", fn_primary)
    in_sam_file = pysam.AlignmentFile(fn_primary, "r")
    out_sam_file = pysam.AlignmentFile(fn_output, "w", header=in_sam_file.header)
    set_primary = set()
    for segment in in_sam_file:
        seq_name  = segment.query_name
        rg_tag    = segment.get_tag("RG")
        set_primary.add((seq_name, rg_tag))
        out_sam_file.write(segment)
    in_sam_file.close()
    
    in_sam_file = pysam.AlignmentFile(fn_secondary, "r")
    for segment in in_sam_file:
        seq_name  = segment.query_name
        rg_tag    = segment.get_tag("RG")
        identity = (seq_name, rg_tag)
        if identity in set_primary:
            pass
        else:
            out_sam_file.write(segment)
    in_sam_file.close()

    out_sam_file.close()
    return set_primary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--primary',   help='primary sam file')
    parser.add_argument('-s', '--secondary', help='secondary sam file')
    parser.add_argument('-o', '--output',    help='output merged sam file')
    args = parser.parse_args()
    
    fn_primary = args.primary
    fn_secondary = args.secondary
    fn_output = args.output

    regroup_sam(fn_primary, fn_secondary, fn_output)
