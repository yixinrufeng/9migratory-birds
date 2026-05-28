#!/bin/bash
#extract vcf that include outgroup
for i in all
do
        vcftools --vcf /path/all.snpeff.ann.vcf --keep ./$i --recode --recode-INFO-all --out $i
        vcftools --vcf ./$i.recode.vcf --recode --recode-INFO-all --max-missing 1.0 --out $i.nomissing
        bcftools norm -m -any $i.nomissing.recode.vcf -Oz -o xx.split.vcf.gz
        bcftools index xx.split.vcf.gz
./homratio.py \
  --vcf xx.split.vcf.gz \
  --target-list Terek.txt \
  --outgroup-list outgroup.txt \
  --prefix lof \
  --nrep 10000 \
  --bins 20 \
  --seed 12345 \
  --outgroup-max-missing 0.2 \
  --target-max-missing 0.5
##cat
rm $i.recode.vcf
rm $i.nomissing.recode.vcf
rm xx.split.vcf.gz
rm xx.split.vcf.gz.csi
done
