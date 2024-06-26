## Following are the shell commands for running HATCHet
## Here's an example of codes to run HATCHet on patient sample A021

# Run freebayes on germline sample to generate heterozygous SNPs
freebayes \
	-f /mnt/projects/huangcy/huangcy/Ref/GRCh38p7.fa \
	-b /mnt/projects/huangcy/huangcy/ITH2/A021/CHL668.bam \
	-C 4 -m 20 --min-coverage 20 > CHL668.freebayes.vcf 

bcftools view -Oz -o /mnt/projects/huangcy/huangcy/ITH2/A021/CHL668.freebayes.vcf.gz /mnt/projects/huangcy/huangcy/ITH2/A021/CHL668.freebayes.filtered.vcf
tabix -p vcf /mnt/projects/huangcy/huangcy/ITH2/A021/CHL668.freebayes.vcf.gz

# Run HATCHet with a config file
export HATCHET_COMPUTE_CN_SOLVER=cbc

hatchet run ./hatchet.A021.ini

hatchet plot-cn results.diploid.n2.bbc.ucn.tsv -n A021 -sG "(16,4)"
