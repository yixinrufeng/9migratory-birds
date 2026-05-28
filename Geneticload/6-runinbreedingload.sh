#!/bin/bash
#Select scaffolds > 1000000 bp
awk -v minlen=1000000 'BEGIN{OFS="\t"} $2 > minlen {print $1, 0, $2}' ../../02reference/GCF_001431845.1_ASM143184v1_genomic.fna.fai > ./scaffolds.gt2Mb.original.bed
awk -v minlen=1000000 '$2 > minlen {print $1}' ../../02reference/GCF_001431845.1_ASM143184v1_genomic.fna.fai > ./scaffolds.gt2Mb.original.list
#plink can calculate 95 chrom at same time, we need split all scaffolds to A.list/bed B.list/bed C.list/bed
bcftools norm -m -any /path/all.snpeff.ann.vcf -Oz -o all.norm.vcf.gz
bcftools index all.norm.vcf.gz
for i in A B C
do
    # 1. 创建染色体重命名映射文件
    awk 'BEGIN{OFS="\t"} {print $1, NR}' ${i}.list > ${i}.scaffold.original_to_number.map
    awk 'BEGIN{OFS="\t"} {print NR, $1}' ${i}.list > ${i}.scaffold.number_to_original.map

    # 2. 提取指定区域的变异
    bcftools view all.norm.vcf.gz \
        --regions-file ${i}.bed \
        -Oz -o ${i}.subset.vcf.gz
    bcftools index ${i}.subset.vcf.gz

    # 3. 重命名染色体（将原始scaffold名改为数字）
    bcftools annotate ${i}.subset.vcf.gz \
        --rename-chrs ${i}.scaffold.original_to_number.map \
        -Oz -o ${i}.renamed.vcf.gz
    bcftools index ${i}.renamed.vcf.gz

    # 4. 排序（可选，重命名后可能需要）
    bcftools sort ${i}.renamed.vcf.gz \
        -Oz -o ${i}.renamed.sorted.vcf.gz
    bcftools index ${i}.renamed.sorted.vcf.gz
#5 plink convert binary
plink --vcf ${i}.renamed.sorted.vcf.gz --double-id --keep-allele-order --allow-extra-chr --chr-set 95 --make-bed --out $i.plink.rename123
plink --bfile $i.plink.rename123 \
        --allow-extra-chr \
        --chr-set 95 \
  --homozyg-window-snp 20 \
  --homozyg-window-het 1 \
  --homozyg-snp 10 \
  --homozyg-kb 100 \
  --homozyg-density 10 \
  --homozyg-window-missing 5 \
  --homozyg-het 1 \
  --homozyg-window-threshold 0.05 \
  --homozyg-group \
  --out $i
done
#manualy select the threshold for long and short roh
#mkdir long 
#cd long
#for i in A B C
#do
#awk '
#NR==1 {
 # for (i=1; i<=NF; i++) {
  #  if ($i=="KB") kb=i
  #}
  #print
  #next
#}
#$kb >= 500
#' ../$i.hom > $i.hom
# if run the >=500, then the non roh regin also include the roh less than 500, we should care about this
STRICT_REF_ANCESTRAL=1 ./inbreedingload.sh ./$i.hom ../$i.renamed.sorted.vcf.gz ../outgroup.txt ../all $i
done
