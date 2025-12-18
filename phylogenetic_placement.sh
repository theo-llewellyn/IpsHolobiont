#Phylogenetic placement of Ips ASVs


###############
# STEP 1: Prepare input files
cd EPA
mkdir Microascales Ophiostomatales
grep 'Microascales' OTUstaxonomy_minusNC.txt | cut -f 1 > Microascales/Microascales_OTUs.txt
grep 'Ophiostomatales' OTUstaxonomy_minusNC.txt | cut -f 1 > Ophiostomatales/Ophiostomatales_OTUs.txt

#filter for OTUs identified as Ophiostomatoids

#pull Ophiostomatoid taxa from ASV file
cat ../OTUs/OTUssequences.fasta | seqkit grep -f Microascales/Microascales_OTUs.txt > MicroascalesMicroascales_OTUs.fasta
cat ../OTUs/OTUssequences.fasta | seqkit grep -f Ophiostomatales/Ophiostomatales_OTUs.txt > Ophiostomatales/Ophiostomatales_OTUs.fasta

#Alignment - on MAFFT web page
#align asvs to ITS alignment on the mafft website use the --add sequences to existing alignment. And choose the fragment option so it doesnt try to align across the whole length
#same for Ophiostomatales
#split the alignment into the OTUs (Microascales_OTUs_msa.fa) and references (ITS_Microascales_trimmed.fas)

#drop tips from reference tree that dont have ITS in R see droptipforEPA.R
#make a file with a list of the taxa in the MSA to keep for the tree
grep '>' ../../MAFFT/Microascales_TBAS_noASVs_nosp/Microascales_7genes_TBAS_noASVs_nosp_msa_p3_ITS-out.fas| sed 's/)/_/g;s/(/_/g;s/</_/g;s/\[/_/g;s/\]/_/g;s/>//g' > ITS_tips_2_keep.txt

#edit MSA taxa to match those in the tree
gsed -i 's/)/_/g;s/(/_/g;s/</_/g;s/\[/_/g;s/\]/_/g' Microascales/ITS_Microascales_trimmed.fas
gsed -i 's/)/_/g;s/(/_/g;s/</_/g' Ophiostomatales/ITS_Ophiostomatales_trimmed.fas

#rename OTUs with the dynamic clustering IDs
awk -F'\t' 'NR > 1 {gsub(/ /, "_", $9); print $1, $1"_"$9}' OFS='\t' ../OTUs/OTUstaxonomy.txt > rename.map
awk 'BEGIN { while ((getline < "rename.map") > 0) map[$1]=$2 }
     /^>/ { id=substr($0, 2); print ">" (id in map ? map[id] : id); next }
     { print }' Microascales/Microascales_OTUs_msa.fa > Microascales/Microascales_OTUs_msa.renamed.fa
awk 'BEGIN { while ((getline < "rename.map") > 0) map[$1]=$2 }
     /^>/ { id=substr($0, 2); print ">" (id in map ? map[id] : id); next }
     { print }' Ophiostomatales/Ophiostomatales_OTUs_msa.fa > Ophiostomatales/Ophiostomatales_OTUs_msa.renamed.fa

### EPA-ng
#best model can check iqtree output for its gene tree, need to add 1.0 to the end as iqtree format only lists the 5 substitution rates and assumed G-T is 1.0

conda activate epa-ng-env

epa-ng \
 --ref-msa Microascales/ITS_Microascales_trimmed.fas \
 --tree Microascales/Microascales_7genes_TBAS_noASVs_nosp_414T.treefile \
 --query Microascales/Microascales_OTUs_msa.renamed.fa \
 --model 'GTR{1.36567,2.57501,2.17568,1.17365,3.71831,1.0}+FU{0.230917,0.269643,0.243777,0.255663}' \
 --outdir Microascales -T 8 \
 --redo
 
#same for Ophiostomatales
epa-ng \
 --ref-msa Ophiostomatales/ITS_Ophiostomatales_trimmed.fas \
 --tree Ophiostomatales/Ophiostomatales_7genes_MagnaporthalesOG_nosp_496T.treefile \
 --query Ophiostomatales/Ophiostomatales_OTUs_msa.renamed.fa \
 --model 'GTR{1.53225,2.19523,1.63412,1.55091,4.06576,1.0}+FU{0.191974,0.335201,0.276926,0.195898}' \
 --outdir Ophiostomatales \
 -T 8 --redo

### GAPPA 

#summarise and place seqs as new branches on tree
conda activate gappa-env
cd Microascales
gappa examine assign --jplace-path epa_result.jplace --taxon-file Gappa_taxonomy_formatted_treenames.txt --allow-file-overwriting --per-query-results --best-hit
gappa examine graft --jplace-path epa_result.jplace --fully-resolve --allow-file-overwriting

cd ../Ophiostomatales
gappa examine assign --jplace-path epa_result.jplace --taxon-file Gappa_taxonomy_formatted_Ophiostomatales.txt --allow-file-overwriting --per-query-results --best-hit
gappa examine graft --jplace-path epa_result.jplace --allow-file-overwriting

#LWR is likelihood weight ratio of that top placement, aLWR incorporates accumulated LWR of all hits so if the placement is split across two taxa in same genus aLWR will combine their values for the placement within that genus. fract is the fraction of placements to that taxon and same with afract
