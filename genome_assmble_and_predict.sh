#Kmer
mkdir survey
nohup /opt/biosoft/jellyfish-2.3.0/bin/jellyfish count -m 21 -t 32 -s 6G -C -o ls_k21.jf <(zcat ../L-6_L2_139A39.R1.fastq.gz) <(zcat ../L-6_L2_139A39.R2.fastq.gz) &
/opt/biosoft/jellyfish-2.3.0/bin/jellyfish histo ls_k21.jf > kmer_21_hist.tsv
nohup /opt/biosoft/jellyfish-2.3.0/bin/jellyfish count -m 19 -t 32 -s 6G -C -o ls_k19.jf <(zcat ../L-6_L2_139A39.R1.fastq.gz) <(zcat ../L-6_L2_139A39.R2.fastq.gz) &
/opt/biosoft/jellyfish-2.3.0/bin/jellyfish histo ls_k19.jf > kmer_19_hist.tsv 
nohup /opt/biosoft/jellyfish-2.3.0/bin/jellyfish count -m 17 -t 32 -s 6G -C -o ls_k17.jf <(zcat ../L-6_L2_139A39.R1.fastq.gz) <(zcat ../L-6_L2_139A39.R2.fastq.gz) &
/opt/biosoft/jellyfish-2.3.0/bin/jellyfish histo ls_k19.jf > kmer_17_hist.tsv 

/opt/biosoft/genomescope2.0-1.0.0/genomescope.R -i kmer_21_hist.tsv -o ./ -k 21 -p 2 > genomescope_21.out
/opt/biosoft/genomescope2.0-1.0.0/genomescope.R -i kmer_19_hist.tsv -o ./ -k 19 -p 2 > genomescope_19.out
/opt/biosoft/genomescope2.0-1.0.0/genomescope.R -i kmer_19_hist.tsv -o ./ -k 17 -p 2 > genomescope_17.out

#genome assembly

cd L_HIFI/CCS
conda activate canu
ccs /nfs_genome1/yanyujie/L_HIFI/P05TYD23817003-1_r64116_20230310_065804_1_A02.subreads.bam ls_A02.ccs.bam &
ccs /nfs_genome1/yanyujie/L_HIFI/P05TYD23817003-1_r64118_20230307_074800_1_F01.subreads.bam ls_F01.ccs.bam &
nohup bam2fasta -o ls_A02.hifi_reads P05TYD23817003-1_r64116_20230310_065804_1_A02.hifi_reads.bam &
nohup bam2fasta -o ls_F01.hifi_reads P05TYD23817003-1_r64118_20230307_074800_1_F01.hifi_reads.bam &
hifiasm -o yyj -t 20 -k 45 -r 2 -a 2 -m 2000000 -p 20000 -n 3 -x 0.8 -y 0.2 -s 0.75 -u --primary -l 0 worm.ccs.fasta 2> hifiasm.log

###genome annotation#####
#
mkdir 00-RepeatMask
RepeatMasker -pa 20 -species Polychaeta -html -gff -xsmall -dir 00-RepeatMask ./purge_yyj_70.fasta

BuildDatabase -name lin -engine ncbi purge_yyj_70.fasta
RepeatModeler -database lin -pa 10 -LTRStruct 

RepeatMasker -lib lin-families.fa -s -gff -parallel 10 -xsmall -alignments -dir RepeatMask ./purge_yyj_70.fasta 

#transcript
####hisat-rnaseq-####
cd /nfs_genome1/yanyujie/L_HIFI/annoatation
mkdir 02-rnaseq
cd 02-rnaseq
hisat2-build /nfs_genome1/yanyujie/L_HIFI/annotation/01-RepeatModeler/RepeatMask/purge_yyj_70.fasta.masked lin
for id in $(ls *_1.clean.fq.gz | sed "s/_1.clean.fq.gz//"); 
do base=`basename $id`; 
echo "hisat2 -x lin -p 10 -1 /nfs_genome1/yanyujie/rnaseq/"$id"_1.clean.fq.gz -2 /nfs_genome1/yanyujie/rnaseq/"$id"_2.clean.fq.gz | samtools sort -@ 10 -O BAM -o $id.bam"; 
done >hisat2.sh
sh hisat2.sh
samtools merge AB.bam AB-1.bam AB-2.bam AB-3.bam AB-4.bam AB-5.bam AB-6.bam 
samtools view -h AB.bam

#de-novo
Trinity --seqType fq --samples_file AB.txt --CPU 20 --max_memory 1024G --output trinity_AB
Trinity --seqType fq --samples_file CE.txt --CPU 20 --max_memory 1024G --output trinity_CE
Trinity --seqType fq --samples_file DB.txt --CPU 20 --max_memory 1024G --output trinity_DB
Trinity --seqType fq --samples_file PB.txt --CPU 20 --max_memory 1024G --output trinity_PB
Trinity --seqType fq --samples_file VB.txt --CPU 20 --max_memory 1024G --output trinity_VB

#transcript_predict
nohup stringtie AB.bam -l AB -o AB.gtf -p 32 >AB.out &
nohup stringtie PB.bam -l PB -o PB.gtf -p 32 >PB.out &
nohup stringtie VB.bam -l VB -o VB.gtf -p 32 >VB.out &
nohup stringtie Ce.bam -l Ce -o Ce.gtf -p 16 >Ce.out &
nohup stringtie DB.bam -l DB -o DB.gtf -p 16 >DB.out &

stringtie --merge AB.gtf PB.gtf VB.gtf Ce.gtf DB.gtf -o rna_merge.gtf
util/gtf_genome_to_cdna_fasta.pl rna_merge.gtf genome.fasta > transcripts.fasta
util/gtf_to_alignment_gff3.pl transcripts.gtf > transcripts.gff3
TransDecoder.Predict -t transcripts.fasta
cdna_alignment_orf_to_genome_orf.pl \
     transcripts.fasta.transdecoder.gff3 \
     transcripts.gff3 \
     transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3
gff3_to_sequences.pl --out_prefix transdecoder.genome.gff3 --only_gene_sequences --only_coding_gene_sequences --only_first_isoform purge_yyj_70.fasta.masked transcripts.fasta.transdecoder.genome.gff3



#braker

nohup braker.pl --species=Lindaspio --genome=purge_yyj_70.fasta.masked --softmasking --bam=AB.bam,DB.bam,Ce.bam,PB.bam,VB.bam --threads 16 --GENEMARK_PATH=/nfs_genome/home/yanyujie01/software/gmes_linux_64_4 --workingdir=/nfs_genome1/yanyujie/L_HIFI/annotation/03-braker >braker.out & 
cd /nfs_genome1/yanyujie/L_HIFI/annotation/03-braker
perl -e 'while (<>) { if (m/\tCDS\t/) { print; s/\tCDS\t/\texon\t/; print; } elsif (m/\ttranscript\t/) { next; } elsif (m/^#/) { next;} elsif (m/\tgene\t/)	{ next; } elsif (m/\tintron\t/) { next; } else { print; } }' braker.gtf > braker_1.gtf
gtfToGff3.pl braker.gtf > braker3.gff3
GFF3Clear --genome genome.fasta --no_attr_add --GFF3_source BRAKER --gene_prefix braker braker3.gff3 > out; mv out braker3.gff3
gff3_to_sequences.pl --out_prefix braker3 --only_gene_sequences --only_coding_gene_sequences --only_first_isoform genome.fasta braker3.gff3

###maker

cd /nfs_genome1/yanyujie/L_HIFI/annotation/04-maker
cp /nfs_genome1/yanyujie/rnaseq/trinity_AB/Trinity.fasta /nfs_genome1/yanyujie/L_HIFI/annotation/04-maker/AB_trinity.fasta 
cp /nfs_genome1/yanyujie/rnaseq/trinity_CE/Trinity.fasta /nfs_genome1/yanyujie/L_HIFI/annotation/04-maker/CE_trinity.fasta
cp /nfs_genome1/yanyiujie/rnaseq/trinity_DB/Trinity.fasta /nfs_genome1/yanyujie/L_HIFI/annotation/04-maker/DB_trinity.fasta
cp /nfs_genome1/yanyujie/rnaseq/trinity_PB/Trinity.fasta /nfs_genome1/yanyujie/L_HIFI/annotation/04-maker/PB_trinity.fasta
cp /nfs_genome1/yanyujie/rnaseq/trinity_VB/Trinity.fasta /nfs_genome1/yanyujie/L_HIFI/annotation/04-maker/VB_trinity.fasta
cat *_trinity.fasta > merged_trinity.fasta
nohup cd-hit-est -i merged_trinity.fasta -o merged_trinity_rmdup_0.9.fasta -c 0.9 -n 10 -d 0 -M 64000 -T 8 &

#protein homology evidence
perl -p -i -e 'if (m/^>/) { s/\./_/g; }' Riftia_pachyptila.protein.fa
perl -p -e 'if (m/^>/) { s/[\|-]/_/g; }' Lamellibrachia_luymsi.proteins.fa >Lamellibrachia_luymsi.proteins_reform.fa
perl -p -i -e 'if (m/^>/) { s/\s+.*//; s/\./_/g; }' Dimorphilus_gyrociliatus.protein.faa
perl -p -i -e 'if (m/^>/) { s/\s+.*//; s/\./_/g; }'  Helobdella_robusta.protein.faa
perl -p -i -e 'if (m/^>/) { s/\s+.*//; s/\./_/g; }' Capitella.teleta_protein.faa
perl -p -i -e 'if (m/^>/) { s/\s+.*//; s/\./_/g; }' Eisenia_andrei.Protein.faa
cat Dimorphilus_gyrociliatus.protein.faa Helobdella_robusta.protein.faa Lamellibrachia_luymsi.proteins_reform.fa Capitella.teleta_protein.faa Eisenia_andrei.Protein.faa Riftia_pachyptila.protein.fa >proteins_new.fasta

nohup mpiexec -n 50 maker > maker_01.out &
fasta_merge -d purge_yyj_70_master_datastore_index.log
gff3_merge -d purge_yyj_70_master_datastore_index.log

#evm
# 对小份数据进行EVM并行运算
~/software/EVidenceModeler-v2.1.0/EvmUtils/write_EVM_commands.pl --genome genome.fasta --gene_predictions gene_predictions.gff3 --transcript_alignments transcript_alignments.gff3 --repeats genome.repeat.gff3  --weights `pwd`/weights.txt --partitions partitions_list.out --output_file_name evm_03.out > commands.write_EVM_commands.list
ParaFly -c commands.write_EVM_commands.list -CPU 64
~/software/EVidenceModeler-v2.1.0/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm_03.out
~/software/EVidenceModeler-v2.1.0/EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output_file_name evm_03.out --genome genome.fasta

cat gene_prediction/*/evm_03.out.gff3 > evm.out.gff3
GFF3Clear --gene_prefix evm --genome genome.fasta --GFF3_source EVM evm.out.gff3 > evm.gff3
gff3_to_sequences.pl --out_prefix evm.filter --only_gene_sequences --only_coding_gene_sequences --only_first_isoform genome.fasta evm.filter.gff3
evm_genes_filtering.pl evm.gff3 evidence.gff3 0.3 genome.fasta 1e-6 0.3 8 > evm.filter.gff3
busco -i ../protein_file/evm_02.protein.fasta -c 64 -o evm_02.protein.fasta -m proteins -l /nfs_genome1/yanyujie/Database/metazoa_odb10 --offline
gff3_to_sequences.pl --out_prefix evm.filter.adjust --only_gene_sequences --only_coding_gene_sequences --only_first_isoform genome.fasta evm.filter.adjust.gff3
busco -i ../protein_file/evm.filter.adjust.protein.fasta -c 64 -o evm.filter.adjust.protein.fasta -m proteins -l /nfs_genome1/yanyujie/Database/metazoa_odb10 --offline

#function_annotation
mkdir /nfs_genome1/yanyujie/L_HIFI/annotation/07-functional_annotation
cd /nfs_genome1/yanyujie/L_HIFI/annotation/07-functional_annotation
perl -p -e 's/\*$//' /nfs_genome1/yanyujie/L_HIFI/annotation/06-evm/02/adjust/evm.filter.adjust.protein.fasta > proteins.fasta
sed 's/^\(>[^ ]*\).*/\1/' proteins.fasta > pro.fasta 
rm proteins.fasta
mv pro.fasta proteins.fasta
#nr
mkdir 01.Nr
cd 01.Nr
ln -sf /nfs_genome_old/BioDB/nr.dmnd ./
diamond blastp --db nr --query ../stben.proteins.fasta --out Nr.xml --outfmt 5 --sensitive --max-target-seqs 20 --evalue 1e-5 --id 10 --index-chunks 1
/nfs_genome1/yanyujie/train_scripts/bin/parsing_blast_result.pl --out-hit-confidence --suject-annotation Nr.xml > Nr.tab
/nfs_genome1/yanyujie/train_scripts/bin/nr_species_distribution.pl Nr.tab > Nr_species_distribution.txt
/nfs_genome1/yanyujie/train_scripts/bin/gene_annotation_from_Nr.pl Nr.tab > Nr.txt

#Swiss-Prot
mkdir -p /home/train/13.functional_annotation/02.Swiss-Prot
cd 02.Swiss-Prot
ln -sf  /nfs_genome1_old/yanyujie/Database/Swiss-prot/uniprot_sprot.dmnd  ./
ln -sf /disks/node1_software_R5_16T/opt/opt/biosoft/bioinfomatics_databases/Swiss-Prot/uniprot_sprot.dmnd ./
diamond blastp --db uniprot_sprot --query ../proteins.fasta --out uniprot_sprot.xml --outfmt 5 --sensitive --max-target-seqs 20 --evalue 1e-5 --id 10 --index-chunks 1
parsing_blast_result.pl --out-hit-confidence --suject-annotation uniprot_sprot.xml > uniprot_sprot.tab
gene_annotation_from_SwissProt.pl uniprot_sprot.tab > SwissProt.txt
cp uniprot_sprot.tab ../functional_annotation.Swiss-Prot.tab
cp SwissProt.txt ../functional_annotation.Swiss-Prot.txt
cd ..

#COG / KOG
mkdir 03.COG
cd 03.COG
ln -sf /disks/node1_software_R5_16T/opt/opt/biosoft/bioinfomatics_databases/COG/kog.dmnd ./
vim diamond.sh
diamond blastp --db kog --query ../proteins.fasta --out kog.xml --outfmt 5 --sensitive --max-target-seqs 200 --evalue 1e-5 --id 10 --tmpdir /dev/shm --index-chunks 1
ssub.pl -n 50 diamond.sh node7
cog_from_xml.pl --coverage 0.2 --evalue 1e-5  --db-fasta /disks/node1_software_R5_16T/opt/opt/biosoft/bioinfomatics_databases/COG/kog.fasta --db-class /disks/node1_software_R5_16T/opt/opt/biosoft/bioinfomatics_databases/COG/kog --fun-txt /disks/node1_software_R5_16T/opt/opt/biosoft/bioinfomatics_databases/COG/fun.kog.txt kog.xml
cut -f 1,3,4 out.annot | gene_annotation_from_table.pl - > KOG.txt
cog_R.pl --title "KOG Function Classification of Whole Genome Genes of Lindaspio polybranchiata" --y-name "Number of Genes" out.class
cp out.annot ../functional_annotation.KOG.tab
cp KOG.txt ../functional_annotation.KOG.txt
cp out.pdf ../functional_annotation.KOG.pdf
cp out.png ../functional_annotation.KOG.png
cd ..

#eggNOG 
mkdir 04.eggNOG
cd 04.eggNOG
ln -sf out.emapper.annotations eggNOG.annot
cp eggNOG.annot ../functional_annotation.eggNOG.tab
grep -v -P "^#" out.emapper.annotations | cut -f 1,8 > eggNOG.txt
cp eggNOG.txt ../functional_annotation.eggNOG.txt
cd ..
#Pfam 注释
mkdir 06.Pfam
cd 06.Pfam
vim pfam.sh
para_hmmscan.pl --outformat --cpu 50 --hmm_db /disks/node1_software_R5_16T/opt/opt/biosoft/bioinfomatics_databases/Pfam/Pfam-A.hmm ../proteins.fasta > Pfam.tab
cut -f 1,2,7 Pfam.tab | perl -e '<>; print <>' | gene_annotation_from_table.pl - > Pfam.txt
cp Pfam.tab ../functional_annotation.Pfam.tab
cp Pfam.txt ../functional_annotation.Pfam.txt
cd ..
# GO 注释
mkdir 07.GO
cd 07.GO
go_from_eggNOG_and_interpro.pl ../04.eggNOG/out.emapper.annotations > go.annot
go_reducing_go_number_para.pl /disks/node1_software_R5_16T/opt/opt/biosoft/go_class/bin/go-basic.obo go.annot 8 > go_reduced.annot
sort go_reduced.annot > go.annot; rm go_reduced.annot
gene_annotation_from_table.pl go.annot > GO.txt
annot2wego.pl go.annot > go.wego
get_Genes_From_GO.pl /disks/node1_software_R5_16T/opt/opt/biosoft/go_class/bin/go-basic.obo go.wego > go_class.tab
go_svg.pl --outdir ./ --name out --color "green" --mark "Whole Genome Genes" --note "GO Class of whole genome genes of Lindaspio polybranchiata" go.wego
perl -p -i -e 's/MovePer:0.125/MovePer:0.5/; s/FontSize:.*/FontSize:24/;' out.lst
/disks/node1_software_R5_16T/opt/opt/biosoft/go_class/svg/distributing_svg.pl out.lst out.svg
/disks/node1_software_R5_16T/opt/opt/biosoft/go_class/bin/changsvgsize.pl out.svg 150 -100
convert out.svg out.png
cp go.annot ../functional_annotation.GO.tab
cp GO.txt ../functional_annotation.GO.txt
cp out.svg ../functional_annotation.GO_class.svg
cp out.png ../functional_annotation.GO_class.png
cp go_class.tab ../functional_annotation.GO_class.tab
cd ..

#KAAS 注释
mkdir 08.KAAS
cd 08.KAAS
# http://www.genome.jp/kaas-bin/kaas_main
gene_annotation_from_kaas.pl query.ko > KEGG.txt
cp KEGG.txt ../functional_annotation.KEGG.txt
cd ..
# integrate
gene_annotation_from_all.pl --all-protein proteins.fasta --html --table-width 3500 --column-width-geneID 150 01.Nr/Nr.txt 02.Swiss-Prot/SwissProt.txt 03.COG/KOG.txt 04.eggNOG/eggNOG.txt 06.Pfam/Pfam.txt 07.GO/GO.txt 08.KAAS/KEGG.txt > functional_annotation.All.html
#Nr.txt	20051	72.26%
#SwissProt.txt	14425	51.98%
#KOG.txt	9391	33.84%
#eggNOG.txt	13937	50.23%
#Pfam.txt	18180	65.52%
#GO.txt	9364	33.74%
#KEGG.txt	8685	31.3%
#All	21462	77.35%
