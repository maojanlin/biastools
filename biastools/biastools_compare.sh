path_out=$1
sample_id=$2
run_id=$3
target_bed=$4
improve_bed=$5
improve_lowRd=$6
path_module=$7
prefix=${path_out}/${sample_id}

bedtools subtract  -a ${improve_bed} -b ${improve_lowRd} > ${prefix}.improve.goodRd.bias.bed
bedtools intersect -a ${target_bed}  -b ${improve_lowRd} > ${prefix}.improve.skipped.bias.bed

python3 ${path_module}compare_bias_with_RD.py -lt ${target_bed} -li ${prefix}.improve.goodRd.bias.bed -lrd ${prefix}.improve.skipped.bias.bed -out ${prefix}.${run_id}.improve.bias.bed
#python3 check_inside_centromere.py -lr1 centromere_extend.bed -lr2 ${prefix}.${run_id}.improve.bias.bed
