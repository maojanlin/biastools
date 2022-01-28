#program to separate true positive alignment from the whole sam file
import argparse
import pysam


def separate_mapping(
        hap_map      :list,
        f_target_sam :pysam.AlignmentFile,
        threshold    :int,
        output_path  :str
        ) -> None:
    """
    compare if each read is within the threshold of golden position
    """
    f_true  = pysam.AlignmentFile(output_path + '.true.sam',  'w', header=f_target_sam.header)
    f_false = pysam.AlignmentFile(output_path + '.false.sam', 'w', header=f_target_sam.header)
    f_miss  = pysam.AlignmentFile(output_path + '.missed.sam', 'w', header=f_target_sam.header)
    for segment in f_target_sam:
        if (segment.flag & 4): # bitwise AND 4, segment unmapped
            continue
        ref_name  = segment.reference_name
        seq_name  = segment.query_name
        pos_start = segment.reference_start
        rg_tag    = segment.get_tag("RG")
        if rg_tag == "hapA":
            hap_idx = 0
        elif rg_tag == "hapB":
            hap_idx = 1
        else:
            print("Error in haplotype information!")
        if (segment.flag & 64):
            pair_idx = 0
        elif (segment.flag & 128):
            pair_idx = 1
        else:
            print("Error in pair-end information!")

        try:
            golden_start = hap_map[hap_idx][ref_name][seq_name][pair_idx]
            if abs(golden_start - pos_start) > threshold:
                f_false.write(segment)
                #print(ref_name, pos_start, golden_start - pos_start)
            else:
                f_true.write(segment)
        except:
            f_miss.write(segment)
    f_true.close()
    f_false.close()
    f_miss.close()


def golden_map(
        f_sam0 :pysam.AlignmentFile,
        f_sam1 :pysam.AlignmentFile
        ) -> dict:
    """
    hap_map is the map recording all reads' golden start position
    """
    hap_map = [{},{}]
    for contig in f_sam0.references:
        hap_map[0][contig] = {}
        hap_map[1][contig] = {}

    for segment in f_sam0:
        if (segment.flag & 4): # bitwise AND 4, segment unmapped
            continue
        ref_name  = segment.reference_name
        seq_name  = segment.query_name
        pos_start = segment.reference_start
        if hap_map[0][ref_name].get(seq_name):
            pass
        else:
            hap_map[0][ref_name][seq_name] = [0,0]
        if (segment.flag & 64): # bitwise AND 64, first pair
            hap_map[0][ref_name][seq_name][0] = pos_start
        elif (segment.flag & 128):
            hap_map[0][ref_name][seq_name][1] = pos_start
        else:
            print("Error! Pair-end information error!!")
    for segment in f_sam1:
        if (segment.flag & 4): # bitwise AND 4, segment unmapped
            continue
        ref_name  = segment.reference_name
        seq_name  = segment.query_name
        pos_start = segment.reference_start
        if hap_map[1][ref_name].get(seq_name):
            pass
        else:   
            hap_map[1][ref_name][seq_name] = [0,0]
        if (segment.flag & 64): # bitwise AND 64, first pair
            hap_map[1][ref_name][seq_name][0] = pos_start
        elif (segment.flag & 128):
            hap_map[1][ref_name][seq_name][1] = pos_start
        else:
            print("Error! Pair-end information error!!")
    return hap_map



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s0', '--hap0_sam', help='golden sam file of haplotype 0 information')
    parser.add_argument('-s1', '--hap1_sam', help='golden sam file of haplotype 1 information')
    parser.add_argument('-st', '--target_sam', help='target sam file for separation')
    parser.add_argument('-th', '--threshold', type=int, default=15, help='the toleratant distance for correct mapping')
    parser.add_argument('-o',  '--out', help='output path')
    args = parser.parse_args()
    
    fn_sam0 = args.hap0_sam
    fn_sam1 = args.hap1_sam
    fn_target_sam = args.target_sam
    threshold = args.threshold
    fn_output = args.out
    
    f_sam0 = pysam.AlignmentFile(fn_sam0)
    f_sam1 = pysam.AlignmentFile(fn_sam1)
    print("Building the golden mapping positions...")
    hap_map = golden_map(
        f_sam0=f_sam0,
        f_sam1=f_sam1
        )
    print("Separate the true and false alignment...")
    f_target_sam = pysam.AlignmentFile(fn_target_sam)
    separate_mapping(
        f_target_sam=f_target_sam,
        hap_map=hap_map,
        threshold=threshold,
        output_path=fn_output
        )

