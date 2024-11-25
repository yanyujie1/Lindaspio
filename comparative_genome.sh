#mkdir -p a.preparing_data/reform_date
for i in $(ls *.genome.fasta | sed "s/.genome.fasta//"); 
do base=`basename $i`; 
echo "GFF3Clear --coverage 0.3  --gene_code_length 5 --gene_prefix $i --genome $i.genome.fasta $i.genome.gff > $i.geneModels.gff3 2> $i.GFF3Clear.log; gff3ToGtf.pl $i.genome.fasta $i.geneModels.gff3 > $i.geneModels.gtf 2> $i.gff3ToGtf.log; eukaryotic_gene_model_statistics.pl $i.geneModels.gtf $i.genome.fasta $i > $i.statistics"
done > command.geneModels.list
   command.geneModels.list -CPU 30
   
mkdir b.OrthoFinder
mkdir compliantFasta
cd compliantFasta
for i in `cut -f 1 ../../a.preparing_data/reform_date/source.txt`
do
    echo "perl -p -e 's/>/>$i\|/' ../../a.preparing_data/reform_date/$i.protein.fasta > $i.fasta"
done | sh
cd ../

mkdir compliantFasta_CDS
cd compliantFasta_CDS
for i in `cut -f 1 ../../a.preparing_data/reform_date/source.txt`
do
    echo "perl -p -e 's/>/>$i\|/' ../../a.preparing_data/reform_date/$i.CDS.fasta > $i.fasta"
done | sh

for i in $(ls *.fasta | sed "s/.fasta//"); 
do base=`basename $i`; 
echo "grep '>' $i.fasta | wc -l";
done > gene_number.sh

orthofinder -f compliantFasta -o OrthoFinder -op > command.diamond.list

perl -p -i -e 'if (m/diamond blastp/ ) { s/$/ --outfmt 5 --max-target-seqs 500 --id 10/; s/--compress 1//; s/.txt/.xml/; s/ -e \S+/ --evalue 1e-5/; s/ -p 1 / -p 4 /; } else { s/^.*\n$//; }' command.diamond.list
ParaFly -c command.diamond.list -CPU 30
perl -e 'while (<>) { print "parsing_blast_result.pl --no-header --max-hit-num 500 --evalue 1e-9 --CIP 0.3 --query-coverage 0.5 --db-query 0.5 $1.xml | gzip -c - > $1.txt.gz\n" if m/(\S+).xml/; }' command.diamond.list > command.parsing_blast_result.list
ParaFly -c command.parsing_blast_result.list -CPU 30

OrthoFinderWorkingDir=`head -n 1 command.diamond.list | perl -ne 'print $1 if m/-d (\S+)\//'`
/nfs_genome_old/home/yanyujie01/software/OrthoFinder/orthofinder -b $OrthoFinderWorkingDir -og
orthomcl_mcl2Groups.pl OrthoFinder/*/WorkingDirectory/OrthoFinder/*/Orthogroups/Orthogroups.txt
ln -sf groups.OG.txt groups.txt

cat groups.OG.txt groups.PG.txt groups.filtered.txt > groups.All.txt
orthomcl_genes_number_stats.pl groups.All.txt compliantFasta > genes_number_stats.txt
cd ..

#c. constract tree
mkdir c.species_tree
cd c.species_tree

orthomcl_extract_ortholog_seqs.pl --out_directory orthologGroups_CDS --species_ratio 1 --single_copy_species_ratio 1 --copy_num 10 --max_seq_length 6000 ../b.OrthoFinder/groups.txt ../b.OrthoFinder/compliantFasta_CDS
orthomcl_extract_ortholog_seqs.pl --out_directory orthologGroups_Protein --species_ratio 1 --single_copy_species_ratio 1 --copy_num 10 --max_seq_length 6000 --compliantFasta_dir_for_calculating_seq_length  ../b.OrthoFinder/compliantFasta_CDS ../b.OrthoFinder/groups.txt ../b.OrthoFinder/compliantFasta
perl -p -i -e 's/\*$//; s/\*/X/g' orthologGroups_Protein/*.fasta
ls orthologGroups_Protein/ | perl -pe 's/.fasta//' > orthologGroups.txt

for i in `cat orthologGroups.txt`
do
    echo "linsi orthologGroups_Protein/$i.fasta > orthologGroups_Protein/$i.fasta.align"
done > command.mafft.list
ParaFly -c command.mafft.list -CPU 32

for i in `cat orthologGroups.txt`
do
    echo "proteinAlignment2CDSAlignemnt.pl orthologGroups_Protein/$i.fasta.align orthologGroups_CDS/$i.fasta > orthologGroups_CDS/$i.fasta.align"
done > command.proteinAlignment2CDSAlignemnt.list
ParaFly -c command.proteinAlignment2CDSAlignemnt.list -CPU 32

for i in `cat orthologGroups.txt`
do
    echo "/opt/biosoft/Gblocks_0.91b/Gblocks orthologGroups_CDS/$i.fasta.align -t=c; if [ -f orthologGroups_CDS/$i.fasta.align-gb ]; then echo $i completed; fi"
    echo "/opt/biosoft/Gblocks_0.91b/Gblocks orthologGroups_Protein/$i.fasta.align -t=p; if [ -f orthologGroups_Protein/$i.fasta.align-gb ]; then echo $i completed; fi"
done > command.Gblocks.list
ParaFly -c command.Gblocks.list -CPU 32

for i in `cat orthologGroups.txt`
do
    perl -p -e 's/(^>\w+).*/$1/; s/ //g' orthologGroups_CDS/$i.fasta.align-gb > orthologGroups_CDS/$i.fasta.align-gb.fasta
done

for i in `cat orthologGroups.txt`
do
    echo "/opt/biosoft/RAxML-8.2.12/usefulScripts/convertFasta2Phylip.sh orthologGroups_CDS/$i.fasta.align-gb.fasta > orthologGroups_CDS/$i.phy"
done > command.convertFasta2Phylip.list
ParaFly -c command.convertFasta2Phylip.list -CPU 32

perl -e 'while (<orthologGroups_CDS/*.align-gb>) { open IN, $_ or die $!; while (<IN>) { if (m/^>(\w+)/) { $seq_id = $1; } else { s/\s+//g; $seq{$seq_id} .= $_; } } } foreach (sort keys %seq) { print ">$_\n$seq{$_}\n"; }' > allSingleCopyOrthologsAlign.Codon.fasta

mkdir RAxML
cd RAxML
/opt/biosoft/RAxML-8.2.12/usefulScripts/convertFasta2Phylip.sh ../allSingleCopyOrthologsAlign.Codon.fasta > allSingleCopyOrthologsAlign.Codon.phy
mpirun -np 20 /opt/biosoft/RAxML-8.2.12/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -# 100 -m PROTGAMMAGTRX -s merge_protein.phy -n ex -T 8 &

cp RAxML_bipartitions.out_codon tree_abbr.RAxML
cp tree_abbr.RAxML tree_fullName.RAxML
cut -f 1,2 ../../a.preparing_data/reform_date/source.txt | perl -p -e 's/\t/\t\"/; s/$/\"/;' | perl -p -e 's#\s+#/#; s#^#perl -p -i -e "s/#; s#\n$#/" tree_fullName.RAxML\n#;' | perl -p -e 's#/\"#/\\\"#; s#\"/#\\\"/#;' | sh
cp tree* ../

# d. mcmctree
mkdir d.divergence_time
cd d.divergence_time
/opt/biosoft/RAxML-8.2.12/usefulScripts/convertFasta2Phylip.sh ../c.species_tree/allSingleCopyOrthologsAlign.Codon.fasta > input.txt
perl -p -i -e 's/\s+/  /' input.txt

perl -p -e 's/mtCDNApri123/input/; s/mtCDNApri/input/; s/<1.0/<8.0/; s/^(\s*ndata =) .*/$1 1/; s/usedata = .*/usedata = 3/; s/model = .*/model = 7/; s/alpha = .*/alpha = 0.5/; s/ncatG = .*/ncatG = 5/; s/cleandata = .*/cleandata = 0/; s/burnin = .*/burnin = 10000000/; s/sampfreq = .*/sampfreq = 100/; s/nsample = .*/nsample = 500000/;' /opt/biosoft/paml4.9i/examples/DatingSoftBound/mcmctree.ctl > mcmctree.ctl
mcmctree mcmctree.ctl
cp out.BV in.BV

perl -p -i -e 's/usedata = .*/usedata = 2    \* 0:/;' mcmctree.ctl
mkdir run01 run02
cp input.txt input.trees mcmctree.ctl in.BV run01
cp input.txt input.trees mcmctree.ctl in.BV run02
echo 'cd run01; mcmctree mcmctree.ctl &> mcmctree.log
cd run02; mcmctree mcmctree.ctl &> mcmctree.log' > command.mcmctree.list
ParaFly -c command.mcmctree.list -CPU 2
perl -n -e 'my $out; while (s/(.*?(\d\.\d+))//) { my $info = $1; my $value = $2; my $new = $value * 100; $info =~ s/$value/$new/;  $out .= $info; } $out .= $_; print $out' run01/FigTree.tre > tree01.nex
perl -n -e 'my $out; while (s/(.*?(\d\.\d+))//) { my $info = $1; my $value = $2; my $new = $value * 100; $info =~ s/$value/$new/;  $out .= $info; } $out .= $_; print $out' run02/FigTree.tre > tree02.nex
perl -e 'while (<>) { if (s/\s*UTREE.*?=\s*//) { s/\s*\[.*?\]//g; print; } }' tree01.nex > tree01.txt
perl -e 'while (<>) { if (s/\s*UTREE.*?=\s*//) { s/\s*\[.*?\]//g; print; } }' tree02.nex > tree02.txt
calculating_branchLength_bias_percentage_of_two_trees.pl --no_normalization_of_total_branch_length tree01.txt tree02.txt > bias_of_2runs.txt



## e.cafe
mkdir e.cafe
cd e.cafe
orthomcl_extract_ortholog_seqs.pl --out_directory orthologGroups_CDS --out_tab_for_CAFE orthomcl2cafe.tab --species_ratio 0.4 --single_copy_species_ratio 0.0 --copy_num 1000 --max_seq_length 100000 ../b.OrthoFinder/groups.txt ../b.OrthoFinder/compliantFasta_CDS
chmod 755 cafe_command
/caferror.py -i cafe_command -e 0.1
parsing_cafeOut.pl caferror_1/cafe_final_report.cafe
/nfs_genome1/yanyujie/train_scripts/bin/parsing_cafeOut.pl caferror_1/cafe_final_report.cafe --out_prefix cafe_out_0.05 --p_value 0.05
python /nfs_genome/BioInformatics/perls/cafetutorial_report_analysis.py -i ./caferror_1/cafe_final_report.cafe -o summary_report
cd ..

## f. PAML
mkdir f.PSG_analysis
cd f.PSG_analysis
orthomcl_extract_ortholog_seqs.pl --out_directory orthologGroups_CDS --species_ratio 0.6 --single_copy_species_ratio 0.0 --copy_num 1000 --max_seq_length 100000 ../b.OrthoFinder/groups.txt ../b.OrthoFinder/compliantFasta_CDS
orthomcl_extract_ortholog_seqs.pl --out_directory orthologGroups_Protein --species_ratio 0.6 --single_copy_species_ratio 0.0 --copy_num 1000 --max_seq_length 100000 --compliantFasta_dir_for_calculating_seq_length ../b.OrthoFinder/compliantFasta_CDS/ ../b.OrthoFinder/groups.txt ../b.OrthoFinder/compliantFasta

perl -p -i -e 's/\*$//; s/\*/X/g' orthologGroups_Protein/*.fasta
ls orthologGroups_Protein/ | perl -pe 's/.fasta//' > orthologGroups.txt
for i in `cat orthologGroups.txt`
do
    echo "/nfs_genome/anaconda/envs/metaphlan/bin/linsi orthologGroups_Protein/$i.fasta > orthologGroups_Protein/$i.fasta.align"
done > command.mafft.list

ssub.pl command.mafft.list

for i in `cat orthologGroups.txt`
do
    echo "proteinAlignment2CDSAlignemnt.pl orthologGroups_Protein/$i.fasta.align orthologGroups_CDS/$i.fasta > orthologGroups_CDS/$i.fasta.align"
done > command.proteinAlignment2CDSAlignemnt.list
ssub.pl command.proteinAlignment2CDSAlignemnt.list node7

perl -p -i -e 's/^(>\w+).*/$1/' orthologGroups_CDS/*.fasta.align

for i in `cat orthologGroups.txt`
do
    perl -p -e 's/(^>\w+).*/$1/; s/ //g' orthologGroups_CDS/$i.fasta.align > orthologGroups_CDS/$i.fasta.align.fasta
done

for i in `cat orthologGroups.txt`
do
    echo "/disks/node1_software_R5_16T/opt/opt/biosoft/RAxML-8.2.12/usefulScripts/convertFasta2Phylip.sh orthologGroups_CDS/$i.fasta.align.fasta > orthologGroups_CDS/$i.phy"
done > command.convertFasta2Phylip.list

ParaFly -c command.convertFasta2Phylip.list -CPU 32
rm -rf orthologGroups_CDS/*.fasta*

# YN00 dn/ds
mkdir dnds_yn00
for i in `cat orthologGroups.txt`
do
    echo "calculating_omega_by_yn00.pl --omega_for_PSG 1.0 --input_tree tree.txt orthologGroups_CDS/$i.phy > dnds_yn00/$i.txt 2> dnds_yn00/$i.stats"
done > command.calculating_omega_by_yn00.list
ParaFly -c command.calculating_omega_by_yn00.list -CPU 32

rename PSGyes PSGyes.dnds_yn00 orthologGroups_CDS/*.phy.PSGyes
ls orthologGroups_CDS/*.phy.PSGyes.dnds_yn00 | perl -ne 'print "$1\n" if m/(OG\d+)/' > PSG_yn00.list
 
mkdir dnds_ML
for i in `cat orthologGroups.txt`
do
    echo "calculating_omega_by_codeml.pl --omega_for_PSG 1.0 orthologGroups_CDS/$i.phy tree.txt > dnds_ML/$i.txt 2> dnds_ML/$i.stats"
done > command.calculating_omega_by_ML.list

ParaFly -c command.calculating_omega_by_ML.list -CPU 32

rename PSGyes PSGyes.dnds_ML orthologGroups_CDS/*.phy.PSGyes
ls orthologGroups_CDS/*.phy.PSGyes.dnds_ML | perl -ne 'print "$1\n" if m/(OG\d+)/' > PSG_ML.list

cat PSG_yn00.list PSG_ML.list | sort | uniq > candidate_PSG.list
rm -rf */dnds_* dnds_ML dnds_yn00


mkdir -p BS_lipol/
for i in `cat candidate_PSG.list`
    do
        echo "paml_branch-site_model_analysis.pl --target_branch_species lipol orthologGroups_CDS/$i.phy tree.txt BS_lipol/$i > BS_lipol/$i.p"
    done > command.paml_branch-site_model.list

ParaFly -c command.paml_branch-site_model.list -CPU 32


#Chemoreceptor and photoreceptor gene annotation
makeblastdb -in lipol.genome.fasta -dbtype nucl -title lipol.genome -parse_seqids -out lipol.genome -logfile lipol.genome.log
cd lipol; mkdir GRL; cd GRL;
tblastn -query ../../PLCB.pep -db ../lipol.genome  -evalue 1e-8 -out lipol.out -outfmt 7 -num_threads 16;
seqkit subseq --bed id.bed -u 10000 -d 10000 -o lipol.PLCB.can.fa ../lipol.genome.fasta;
cp ~/script/genewise.sh ./;
sh genewise_command_list.txt;
para_hmmscan.pl --outformat --cpu 4 --hmm_db /opt/biosoft/bioinfomatics_databases/Pfam/Pfam-A.hmm ../PLCB.protein.fasta > Pfam.tab;

for i in `ls ./*.fasta`
do 
	x=${i/.fasta/}
	x=${x/*\//}
	echo "makeblastdb -in $x.fasta -dbtype nucl -title $x -parse_seqids -out $x.gene -logfile $x.fasta.log"
done > make_db.sh

for i in `ls ./*.fasta`
do 
	x=${i/.fasta/}
	x=${x/*\//}
	echo "tblastn -query ./pep/gene_loss.pep -db $x.gene  -evalue 1e-8 -out $x.out -outfmt '7 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads 16"
done > tblastn_db.sh






















