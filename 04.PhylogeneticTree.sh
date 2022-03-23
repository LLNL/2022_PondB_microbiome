# 04. Phylogenetic Tree
# Ref for SINA: Pruesse, E., Peplies, J. & Glöckner, F. O. SINA: Accurate high-throughput multiple sequence alignment of ribosomal RNA genes. Bioinformatics 28, 1823–1829 (2012).
# Ref for trimal: Capella-Gutiérrez, S., Silla-Martínez, J. M. & Gabaldón, T. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics 25, 1972–1973 (2009).
# Ref for FastTree2: Price, M. N., Dehal, P. S. & Arkin, A. P. FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments. PLOS ONE 5, e9490 (2010).
# Water and sediment sequences were analyzed separately

wd=<path to working dir>
cd ${wd}

#### Make fasta file with sequences of OTUs that are being kept ####
input=${wd}/tax_table_for_tree_ps_filt5.txt
output=${wd}/sequences_for_tree_ps_filt5.fa

python3 <path to tab2fasta.py> ${input} 2 1 > ${output}
#python tab2fasta.py <tab-file> <sequence column> <header column 1> <header column 2> <header column n>  > <outfile>
sed -i -e '1,2d' ${output} #delete first two lines from file
sed 's/_NA//g' ${output} > ${output}2 #delete all the NAs
rm ${output}-e ${output}
mv ${output}2 ${output}

#### SINA alignment (v1.6.0) ####
OTU_input=${wd}/sequences_for_tree_ps_filt5.fa
concate_align=${ASV_input%.*}.clean_ps_filt5.fa
ref_db=<path to SILVA_138_SSURef_NR99_05_01_20_opt.arb>

# Clean the names
tr "[ -%,;\(\):=\.\\\[]\"\']" "_" < ${OTU_input} > ${concate_align}

# sina alignment (v1.7.2)
sina -i ${concate_align} -o ${concate_align%.*}_SINAaligned.fasta -r ${ref_db} \
-v --log-file ${concate_align%.*}_SINAaligned_log.txt --fasta-write-dna \
--search --meta-fmt=csv \
--search-max-result 1 --lca-fields=tax_slv

#### remove gaps in columns containing only gaps ####
trimal -in ${concate_align%.*}_SINAaligned.fasta -out ${concate_align%.*}_SINAaligned_removedgaps.fasta -htmlout ${concate_align%.*}_SINAaligned_removedgaps.html -fasta -noallgaps

sed '/^>/ s/ .*//' ${concate_align%.*}_SINAaligned_removedgaps.fasta > ${concate_align%.*}_SINAaligned_removedgaps2.fasta
mv ${concate_align%.*}_SINAaligned_removedgaps2.fasta ${concate_align%.*}_SINAaligned_removedgaps.fasta

#### Fast tree 2 ####
fasttree -nt -gtr -log ${concate_align%.*}_SINAaligned_removedgaps_tree.log ${concate_align%.*}_SINAaligned_removedgaps.fasta > ${concate_align%.*}_SINAaligned_removedgaps_tree.tre

# slow option
fasttree -nt -gtr -slow -log ${concate_align%.*}_SINAaligned_removedgaps_tree_slow.log ${concate_align%.*}_SINAaligned_removedgaps.fasta > ${concate_align%.*}_SINAaligned_removedgaps_tree_slow.tre

# Check using FigTree; save as a newick file
