mkdir lift_sam/
mkdir mapped_sam/
python3 ../levioSAM/scripts/chain_invert.py -c ref2hapA.chain -o hapA2ref.chain
python3 ../levioSAM/scripts/chain_invert.py -c ref2hapB.chain -o hapB2ref.chain
leviosam index -t 4 -c hapA2ref.chain -F GRCh38_chr21.hapA.fa.fai -p ./lift_sam/golden_hapA_to_chr21
leviosam index -t 4 -c hapB2ref.chain -F GRCh38_chr21.hapB.fa.fai -p ./lift_sam/golden_hapB_to_chr21
leviosam lift -C ./lift_sam/golden_hapA_to_chr21.clft -a hapA.sorted.bam -p ./lift_sam/golden_hapA_to_chr21 -O bam
leviosam lift -C ./lift_sam/golden_hapB_to_chr21.clft -a hapB.sorted.bam -p ./lift_sam/golden_hapB_to_chr21 -O bam

#samtools sort  ./lift_sam/golden_hapA_to_chr21.bam > ./lift_sam/golden_hapA_to_chr21.sorted.bam
#samtools index ./lift_sam/golden_hapA_to_chr21.sorted.bam
#samtools sort  ./lift_sam/golden_hapB_to_chr21.bam > ./lift_sam/golden_hapB_to_chr21.sorted.bam
#samtools index ./lift_sam/golden_hapB_to_chr21.sorted.bam
#python3 separate_sam.py -s0 lift_sam/golden_hapA_to_chr21.sorted.bam -s1 lift_sam/golden_hapB_to_chr21.sorted.bam -st bt2.sorted.het.bam -th 20 -o mapped_sam/bt2.het 
python3 separate_sam.py -s0 lift_sam/golden_hapA_to_chr21.bam -s1 lift_sam/golden_hapB_to_chr21.bam -st bt2.sorted.het.bam -th 20 -o mapped_sam/bt2.het 
