#!/bin/bash
#PipeRiboMNase.sh
#This is a pipeline for MNase treated Riboseq data analysis
#This script requires the lncRNA $7=$8 in ${genome}.RefSeq.reduced.bed12
#Inputs: Trimmed (e.g use PipeSETrimmer.sh) Ribo-seq single-end raw reads in FASTQ format
#This pipeline will remove rRNA and miRNA reads, and use rRNA/miRNA free reads for genome mapping and downstream analysis
#This pipeline also predicts P site offset based on the reads mapped to STOP codon, using random forest method
#The R scripts involved in P site learning and prediction require several packages
#Outputs: 
#Example: PipeRiboseq.sh -i Embryo.RPF.fastq.gz -g hg38 -noqc -normCDS
#Version1 basic write up: 2021-11, Li Xie and Yu Sun

array=( "$@" )

#Usage
if [ ! -n "$1" ]
then
  echo "********************************************************************************************"
  echo "*               PipeRiboMNase: pipeline for MNase treated Ribo-seq analysis.               *"
  echo "*                                Version 1, 2021-11-15                                     *"
  echo "* Usage: `basename $0`                                                                  *"
  echo "*        Required (single end data):                                                       *"
  echo "*                  -i [Data.fastq or Data.fastq.gz]                                        *"
  echo "*                  -g [mm10/hg38/rn6/GRCm38/dm6]                                           *"
  echo "*        Optional: -p [Number of CPUs, default=1]                                          *"
  echo "*                  -m [Genome mapping mismatch, default=2]                                 *"
  echo "*                  -s [Salmon strand type, default=A, automatically detect]                *"
  echo "*                  -pre Sequences.fa [Run pre-mapping to Sequences.fa]                     *"
  echo "*                  -cufflinksrun [Run cufflinks]                                           *"
  echo "*                  -t [Cufflinks library type, default=fr-firststrand]                     *"
  echo "*                  -noqc [Suppress fastqc]                                                 *"
  echo "*                  -plotRNA FILE [A list of mRNAs]                                         *"
  echo "*                  -normCDS [Normalize by CDS counts]                                      *"
  echo "*                  -normFactor [Num] [User defined normalization factor, to divide]        *"
  echo "*                             -normCDS and -normFactor cannot run at the same time         *"
  echo "*                  -PsiteModel FILE [Rdata model file for random forest]                   *"
  echo "* Inputs: The fastq file needs to be trimmed, and has .fastq (or.fastq.gz) as suffix       *"
  echo "* Run: Default to run fastqc, rRNA & genomic mapping, RiboQC, featureCounts, salmon        *"
  echo "*      Figures will be generated in /plots folder, and bigWig files in /tracks folder      *"
  echo "* Outputs: All output files will be generated in the same folder as the pipeline submitted *"
  echo "* This pipeline requires additional scripts in PipelineHomeDir/bin folder                  *"
  echo "* If PsiteModel provided, it will use it to predict P site, or it will train using the data*"
  echo "* The pipeline sets default P site offset to be 14 in case the prediction step fails       *"
  echo "********************************************************************************************"
  exit 1
fi

echo "*****************************************************************************"
echo "*              PipeRiboseq: pipeline for Ribo-seq analysis.                 *"
echo "*                         Version 1, 2021-11-15                             *"
echo "*****************************************************************************"
echo "0. Loading softwares:"
#Get current time
St=`date`

#Load softwares: This part may need to be changed when using different cluster
#Main tools: STAR, bowtie, bowtie2, cufflinks, salmon, featureCounts (from Subread package)
#Other tools: samtools, bedtools
#R packages required: dplyr, mlr, mlr3, ranger, ggplot2, forcats, parallelMap, pheatmap

#Module load for bsub cluster:
#if [ -f /etc/profile.d/modules.sh ]; then
#source /etc/profile.d/modules.sh 
#fi
#module load STAR-2.5.2a 
#module load cufflinks-2.2.1
#module load samtools-1.1
#module load bowtie2-2.1.0
#module load FastQC-0.11.2

#Module load for Bluehive:
module load bowtie
module load bowtie2
module load rnastar
module load subread
#module load cufflinks      #Optional. Only required when -cufflinksrun specified.
#module load salmon
module load bedtools
module load samtools
module load fastqc         #Optional when -noqc specified.
#module load kentutils     #The two tools used in this pipeline have been put in /bin folder
module load r/3.5.0/b1     #A set of R packages must be installed

#Get pipeline directory
HomeDir=$(dirname `readlink -f $0`)
#For mac local check only:
#HomeDir="/PipelineHomeDir"
echo "   Home Directory:"
echo "   "$HomeDir

echo "1. Resolving inputs:"

#Get parameters
Data="unassigned"
datastats=0
genome="unassigned"
genomestats=0
CPU=1
GenomeMM=2               #Default Genome mapping MM
LibraryType="fr-firststrand"
StrandType="unassigned"
runfastqc=1              #You can change this default to allow/suppress fastqc
runcufflinks=0
premap=0
premapData="unassigned"
harrdata=0
normCDS=0
PsiteModel="unassigned"
plotRNA="unassigned"
normFactor="unassigned"

for arg in "$@"
do
 if [[ $arg == "-i" ]]
  then
    Data=${array[$counter+1]}
    echo '   Single end data: '$Data
 elif [[ $arg == "-g" ]]
  then
    genome=${array[$counter+1]}
    echo '   Genome: '$genome
 elif [[ $arg == "-p" ]]
  then
    CPU=${array[$counter+1]}
    echo '   CPU: '$CPU
 elif [[ $arg == "-m" ]]
  then
    GenomeMM=${array[$counter+1]}
    echo '   Genome mapping mismatch: '$GenomeMM
 elif [[ $arg == "-t" ]]
  then
    LibraryType=${array[$counter+1]}
    echo '   Library type: '$LibraryType
 elif [[ $arg == "-s" ]]
  then
    StrandType=${array[$counter+1]}
    echo '   Salmon strand type: '$StrandType
 elif [[ $arg == "-pre" ]]
  then
    premapData=${array[$counter+1]}
    premap=1
    echo '   Pre mapping file: '$premapData
 elif [[ $arg == "-cufflinksrun" ]]
  then
    runcufflinks=1
    echo '   Run cufflinks'
 elif [[ $arg == "-noqc" ]]
  then
    runfastqc=0
    echo '   Suppress fastqc quality control'
 elif [[ $arg == "-plotRNA" ]]
  then
    plotRNA=${array[$counter+1]}
    echo '   Plot a set of user-provided mRNAs'
 elif [[ $arg == "-normCDS" ]]
  then
    normCDS=1
    echo '   Use CDS mapping reads for normalization'
 elif [[ $arg == "-normFactor" ]]
  then
    normFactor=${array[$counter+1]}
    echo '   Use user defined norm factor: '$normFactor
 elif [[ $arg == "-PsiteModel" ]]
  then
    PsiteModel=${array[$counter+1]}
    echo '   Use model file: '$PsiteModel
 fi
  let counter=$counter+1
done

calculating(){ awk "BEGIN { print "$*" }"; }

CountFASTALength() {
  Data=$1
  awk 'NR%2==1' $Data | sed 's/^>//' > $Data.header
  awk 'NR%2==0' $Data > $Data.seq
  awk '{print length($0)}' $Data.seq > $Data.seq.len
  paste $Data.header $Data.seq.len > $Data.len
  rm -rf $Data.header && rm -rf $Data.seq && rm -rf $Data.seq.len
}

#Get current directory and create folders or files
[ $runfastqc == "1" ] && FastqcDir=fastqc && mkdir -p $FastqcDir
[ $premap == "1" ] && PremapDir=pre_mapping && mkdir -p $PremapDir
GenomeMappingDir=genome_mapping && mkdir -p $GenomeMappingDir
[ $runcufflinks == "1" ] && CufflinksDir=cufflinks_results && mkdir -p $CufflinksDir
SalmonOutputDir=salmon_results && mkdir -p $SalmonOutputDir
TxMappingDir=tx_mapping && mkdir -p $TxMappingDir
FeatureDir=feature_counts && mkdir -p $FeatureDir
FiguresDir=plots && mkdir -p $FiguresDir
TracksDir=tracks && mkdir -p $TracksDir

#Check data
if [ $Data == "unassigned" ];then
  echo "     >>> [Error]: Please input data!"
else
  datastats=1
fi

#Getting Output Suffix and default strand types
OutputSuffix=$(echo $Data|sed 's/.fastq.*//g')
TABLE=${OutputSuffix}.summary
if [ $StrandType == "unassigned" ];then
  StrandType="A"
fi

#Check genome
if [ $genome == "unassigned" ];then
  echo "     >>> [Error]: Please assign genome file!"
else
  if [ -d $HomeDir/$genome ]; then
    echo "     >>> This genome is supported."
    genomestats=1
  else
    echo "     >>> [Error]: Genome unsupported!"
  fi
fi

#Export checking status
if [ $datastats == "1" -a $genomestats ==  "1" ];then
    echo "   Output Suffix: "$OutputSuffix
  echo "   Status check: pass"
else
  echo "   Status check: failed, stop the pipeline"
  exit 1
fi

echo "*****************************************************************************"
echo "2. Checking dependencies:"
echo "Annotations:"
#Genome files and index check
#This check includes folders: /Annotation, /Sequence, /Index
#                    files: ${genome}.RefSeq.gtf, ${genome}.RefSeq.bed12
#                    files: $genome.fa, $genome.RefSeq.fa
dependenciescount=0

#Annotation and FASTA files
if [ -s $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12 ];then
  echo "   BED12 Annotation: "${genome}.RefSeq.reduced.bed12
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] BED12 Annotation lost "
fi
if [ -s $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.mRNA.bed12 ];then
  echo "   BED12 mRNA Annotation: "${genome}.RefSeq.reduced.mRNA.bed12
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] BED12 mRNA Annotation lost, generating..."
  awk '$7!=$8' $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12 > \
      $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.mRNA.bed12 && echo "   BED12 mRNA Annotation generated" && \
      let dependenciescount=$dependenciescount+1
fi
if [ -s $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12.geneid.gtf ];then
  echo "   GTF Annotation: "${genome}.RefSeq.reduced.bed12.geneid.gtf
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] GTF Annotation lost "
fi
if [ -s $HomeDir/$genome/Annotation/${genome}.uniqMatching.txt ];then
  echo "   RefSeq to refFlat: "${genome}.uniqMatching.txt
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] RefSeq to refFlat Matching lost "
fi
if [ -s $HomeDir/$genome/Sequence/${genome}.fa ];then
  echo "   Genome: "${genome}.fa
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] Genome lost "
fi
if [ -s $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ];then
  echo "   ChromInfo: "${genome}.ChromInfo.txt
  let dependenciescount=$dependenciescount+1
else
  $HomeDir/bin/faSize -tab -detailed $HomeDir/$genome/Sequence/${genome}.fa > $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt && \
  echo "   [Error] ChromInfo lost, generate a new one.." && let dependenciescount=$dependenciescount+1
fi
if [ -s $HomeDir/$genome/Sequence/${genome}.RefSeq.reduced.bed12.fa ];then
  echo "   Transcriptome: "${genome}.RefSeq.reduced.bed12.fa
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] Transcriptome lost"
fi
if [ -s $HomeDir/$genome/Sequence/${genome}.RefSeq.reduced.bed12.fa.len ];then
  echo "   Transcriptome length: "${genome}.RefSeq.reduced.bed12.fa.len
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] Transcriptome length file lost, generating..."
  $HomeDir/bin/faSize -tab -detailed $HomeDir/$genome/Sequence/${genome}.RefSeq.reduced.bed12.fa > $HomeDir/$genome/Sequence/${genome}.RefSeq.reduced.bed12.fa.len
  let dependenciescount=$dependenciescount+1
fi
if [ -s $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.mRNA.cds.bed12.tx.STOP ];then
  echo "   STOP codon annotation: "${genome}.RefSeq.reduced.mRNA.cds.bed12.tx.STOP
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] STOP codon annotation lost, generating..."
  if [ ! -s ${genome}.RefSeq.reduced.mRNA.cds.bed12 ];then
    bash $HomeDir/bin/BED12Extractor.sh -a cds -i $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.mRNA.bed12 -o $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.mRNA.cds.bed12
  fi
  if [ ! -s ${genome}.RefSeq.reduced.mRNA.cds.bed12.tx ];then
    bash $HomeDir/bin/BED12TranslateGenomePos2Tx.sh -i $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.mRNA.bed12 -t $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.mRNA.cds.bed12 -o $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.mRNA.cds.bed12.tx
  fi
  awk '{OFS="\t";print $1,$3-3,$3}' $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.mRNA.cds.bed12.tx > $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.mRNA.cds.bed12.tx.STOP
  let dependenciescount=$dependenciescount+1
fi

###Index folders
echo "Index files:"
if [ -d $HomeDir/$genome/Index/STARIndex ];then
  echo "   STARIndex exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] STARIndex lost"
fi
if [ -d $HomeDir/$genome/Index/SalmonIndex ];then
  echo "   SalmonIndex exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] SalmonIndex lost "
fi
if [ -d $HomeDir/$genome/Index/rRNAIndex ];then
  echo "   rRNAIndex exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] rRNAIndex lost "
fi
if [ -d $HomeDir/$genome/Index/miRNAIndex ];then
  echo "   miRNAIndex exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] miRNAIndex lost "
fi
if [ -d $HomeDir/$genome/Index/TxIndex ];then
  echo "   TxIndex exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] TxIndex lost "
fi

#echo $dependenciescount
if [ $dependenciescount == "14" ];then
  echo "   Dependencies check: pass"
else
  echo "   [Error] Some dependencies lost "
fi

echo "*****************************************************************************"
echo "3. Run fastqc quality control:"
if [ $runfastqc == "1" ];then
echo "   Running fastqc single-end mode"
fastqc \
    -f fastq \
    -o $FastqcDir \
    $Data \
    2> $FastqcDir/${OutputSuffix}.fastqc.log && \
    echo "   Done fastqc"
elif [ $runfastqc == "0" ];then
echo "   Skipping fastqc..."
fi

echo "*****************************************************************************"
echo "4. rRNA mapping and removal for further analysis:"
echo "   Running bowtie2 single-end mode"
bowtie2 \
    -x $HomeDir/$genome/Index/rRNAIndex/rRNAIndex \
    -U $Data \
    -q \
    --very-sensitive-local \
    -k 1 \
    -p $CPU \
    -S /dev/null \
    --un-gz $GenomeMappingDir/${OutputSuffix}.No_rRNA.fastq.gz \
    2> $GenomeMappingDir/${OutputSuffix}.rRNA.log && \
    echo "   Done rRNA mapping"
	rRNAReads=`head -4 $GenomeMappingDir/${OutputSuffix}.rRNA.log | tail -1 | awk '{print $1}'`

echo "*****************************************************************************"
echo "5. miRNA mapping and removal for further analysis:"
echo "   Running bowtie2 single-end mode"
bowtie2 \
    -x $HomeDir/$genome/Index/miRNAIndex/miRNAIndex \
    -U $GenomeMappingDir/${OutputSuffix}.No_rRNA.fastq.gz \
    -q \
    --very-sensitive-local \
    -k 1 \
    -p $CPU \
    -S /dev/null \
    --un-gz $GenomeMappingDir/${OutputSuffix}.No_rRNA.No_miRNA.fastq.gz \
    2> $GenomeMappingDir/${OutputSuffix}.miRNA.log && \
    echo "   Done miRNA mapping"
	miRNAReads=`head -4 $GenomeMappingDir/${OutputSuffix}.miRNA.log | tail -1 | awk '{print $1}'`

if [ $premap == "1" ];then
echo "5.a Running bowtie2 single-end mode for extra sequences"
echo "   Building index for input sequence using bowtie2"
echo "   Those mapped reads won't be removed for further analysis"
bowtie2-build $premapData $PremapDir/ExtraSeq 2> $PremapDir/${OutputSuffix}.index.log 1> $PremapDir/temp && rm -rf $PremapDir/temp
bowtie2 \
    -x $PremapDir/ExtraSeq \
    -U $Data \
    -q \
    -a \
    --very-sensitive-local \
    -p $CPU \
    -S $PremapDir/${OutputSuffix}.premap.sam \
    2> $PremapDir/${OutputSuffix}.premap.log && \
    echo "   Done extra sequences mapping"
	premapReads=`head -4 $PremapDir/${OutputSuffix}.premap.log | tail -1 | awk '{print $1}'`

echo "   Summarizing extra sequences mapping results"
samtools view -F 4 $PremapDir/${OutputSuffix}.premap.sam > $PremapDir/${OutputSuffix}.premap.sam.hits
samtools view -bS $PremapDir/${OutputSuffix}.premap.sam > $PremapDir/${OutputSuffix}.premap.bam && rm -rf $PremapDir/${OutputSuffix}.premap.sam
rm -rf $PremapDir/${OutputSuffix}.premap.summary.txt && echo -e "Gene\tReads" > $PremapDir/${OutputSuffix}.premap.summary.txt
grep '>' $premapData | sed 's/>//' > $PremapDir/FastaList
awk '{print $3}' $PremapDir/${OutputSuffix}.premap.sam.hits |sort|uniq -c|awk '{OFS="\t";print $2,$1}' > $PremapDir/${OutputSuffix}.premap.summary.temp
awk '{print $1}' $PremapDir/${OutputSuffix}.premap.summary.temp > $PremapDir/${OutputSuffix}.premap.summary.tempList
cat $PremapDir/FastaList $PremapDir/${OutputSuffix}.premap.summary.tempList | sort|uniq -c|awk '$1==1'|awk '{OFS="\t";print $2,0}' > $PremapDir/${OutputSuffix}.premap.summary.tempfill
cat $PremapDir/${OutputSuffix}.premap.summary.temp $PremapDir/${OutputSuffix}.premap.summary.tempfill > $PremapDir/${OutputSuffix}.premap.summary.tempall
for name in `cat $PremapDir/FastaList`;do awk -v t=$name '{OFS="\t";if ($1==t) print $0}' $PremapDir/${OutputSuffix}.premap.summary.tempall;done >> $PremapDir/${OutputSuffix}.premap.summary.txt
rm -rf $PremapDir/${OutputSuffix}.premap.summary.temp* $PremapDir/FastaList 
fi

echo "*****************************************************************************"
echo "6. Genome mapping using STAR:"
echo "   Running STAR with "${GenomeMM}" mismatch(es) allowed"
STAR \
    --runMode alignReads \
    --genomeDir $HomeDir/$genome/Index/STARIndex \
    --readFilesIn $GenomeMappingDir/${OutputSuffix}.No_rRNA.No_miRNA.fastq.gz \
    --readFilesCommand zcat \
    --runThreadN $CPU \
    --outSAMattributes All \
    --outFilterMismatchNmax $GenomeMM \
    --alignEndsType EndToEnd \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --outMultimapperOrder Random \
    --outFileNamePrefix $GenomeMappingDir/${OutputSuffix}.${genome}. \
    --outSAMtype BAM Unsorted 2>&1 1> $GenomeMappingDir/${OutputSuffix}.STAR.log && \
    echo "   Done STAR mapping"

echo "*****************************************************************************"
echo "7. Post-mapping processing:"
echo "   Sorting the bam file..."
samtools sort -T $GenomeMappingDir/aln.sorted $GenomeMappingDir/${OutputSuffix}.${genome}.Aligned.out.bam -o $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam
echo "   Indexing..."
samtools index $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam
samtools view -h $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam
echo "   Extracting unique mapping RPF reads, to get *RPF.unique.bam"
samtools view -H $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.header
samtools view -q 10 $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.uMAPQ10
#Expand the RPF length to 40nt
awk '{if (length($10) >25 && length($10) <41) print $0}' $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.uMAPQ10 > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.uMAPQ10.MNase
cat $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.header $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.uMAPQ10.MNase > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.unique.sam
samtools view -b $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.unique.sam > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.unique.bam
samtools index $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.unique.bam
rm -rf $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.header $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.uMAPQ10
rm -rf $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.uMAPQ10.MNase $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.unique.sam
rm -rf $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam

echo "   Getting statistics:"
rm -rf $TABLE
InputReads_rRNA_miRNA=`grep 'Number of input reads' $GenomeMappingDir/$OutputSuffix.${genome}.Log.final.out | awk '{print $NF}'`
InputReads=$((InputReads_rRNA_miRNA+rRNAReads+miRNAReads))
UniquReads=`grep 'Uniquely mapped reads number' $GenomeMappingDir/$OutputSuffix.${genome}.Log.final.out | awk '{print $NF}'`
MultiReads=`grep 'Number of reads mapped to multiple loci' $GenomeMappingDir/$OutputSuffix.${genome}.Log.final.out | awk '{print $NF}'`
AllMapReads=$((UniquReads+MultiReads))
UnMapReads=$((InputReads_rRNA_miRNA-UniquReads-MultiReads))
GMappingRate=`awk -v mapping=$AllMapReads -v total=$InputReads 'BEGIN{print (mapping/total*100)}'`
ContaminationRate=`awk -v mapping=$rRNAReads -v total=$InputReads 'BEGIN{print (mapping/total*100)}'`
ContaminationMIRate=`awk -v mapping=$miRNAReads -v total=$InputReads 'BEGIN{print (mapping/total*100)}'`
MappingRate=`awk -v m1=$AllMapReads -v m2=$rRNAReads -v m3=$miRNAReads -v total=$InputReads 'BEGIN{print ((m1+m2+m3)/total*100)}'`

echo -e "   total input reads:\t${InputReads}"
echo -e "   rRNA_reads:\t${rRNAReads}"
echo -e "   miRNA_reads:\t${miRNAReads}"
if [ $premap == "1" ];then
  echo -e "   premap_reads:\t${premapReads}"
fi
echo -e "   genomic_mapped_reads:\t${AllMapReads}"
echo -e "   genomic_unique_mapped_reads:\t${UniquReads}"
echo -e "   genomic_multiple_mapped_reads:\t${MultiReads}"
echo -e "   genomic_unmappable_reads:\t${UnMapReads}"
echo -e "   genomic mapping rate(%):\t${GMappingRate}"
echo -e "   rRNA mapping rate(%):\t${ContaminationRate}"
echo -e "   Mapping rate(%):\t${MappingRate}"
echo -e "total input reads:\t${InputReads}" >> $TABLE
echo -e "rRNA_reads:\t${rRNAReads}" >> $TABLE
echo -e "miRNA_reads:\t${miRNAReads}" >> $TABLE
if [ $premap == "1" ];then
  echo -e "premap_reads:\t${premapReads}" >> $TABLE
fi
echo -e "genomic_mapped_reads:\t${AllMapReads}" >> $TABLE
echo -e "genomic_unique_mapped_reads:\t${UniquReads}" >> $TABLE
echo -e "genomic_multiple_mapped_reads:\t${MultiReads}" >> $TABLE
echo -e "genomic_unmappable_reads:\t${UnMapReads}" >> $TABLE
echo -e "genomic mapping rate(%):\t${GMappingRate}" >> $TABLE
echo -e "rRNA mapping rate(%):\t${ContaminationRate}" >> $TABLE
echo -e "Mapping rate(%):\t${MappingRate}" >> $TABLE

echo "*****************************************************************************"
echo "8. Transcriptome assembly and abundance calculation using cufflinks:"
if [ $runcufflinks == "1" ];then
echo "   Running cufflinks single-end mode"
echo "   Using library type: "$LibraryType
cufflinks \
    -o $CufflinksDir \
    -p $CPU \
    -G $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12.geneid.gtf \
    -b $HomeDir/$genome/Sequence/${genome}.fa \
    -u \
    --library-type $LibraryType \
    --compatible-hits-norm \
    --no-update-check \
    $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam \
                2> $CufflinksDir/${OutputSuffix}.${genome}.cufflinks.log
mv $CufflinksDir/genes.fpkm_tracking $CufflinksDir/${OutputSuffix}.${genome}.genes.fpkm_tracking
mv $CufflinksDir/isoforms.fpkm_tracking $CufflinksDir/${OutputSuffix}.${genome}.isoforms.fpkm_tracking
mv $CufflinksDir/transcripts.gtf $CufflinksDir/${OutputSuffix}.${genome}.transcripts.gtf
echo "   Done running cufflinks"
else
  echo "   Skipping cufflinks..."
fi

echo "*****************************************************************************"
echo "9. Direct transcriptome mapping using Salmon:"
echo "   Salmon will automatically detect library type by default (-l A)"
echo "   Using strand type: "$StrandType
module load salmon
salmon quant \
    -i $HomeDir/$genome/Index/SalmonIndex \
    -p $CPU \
    -l $StrandType \
    -r $Data \
    -o $SalmonOutputDir \
                2> $SalmonOutputDir/${OutputSuffix}.${genome}.salmon.log
echo "   Done Salmon quantification."
echo "   You can look at the strand info in lib_format_counts.json file"
mv $SalmonOutputDir/quant.sf $SalmonOutputDir/${OutputSuffix}.${genome}.quant.sf
mkdir $SalmonOutputDir/${OutputSuffix}.salmon
mv $SalmonOutputDir/aux_info $SalmonOutputDir/${OutputSuffix}.salmon/ && mv $SalmonOutputDir/cmd_info.json $SalmonOutputDir/${OutputSuffix}.salmon/
mv $SalmonOutputDir/lib_format_counts.json $SalmonOutputDir/${OutputSuffix}.salmon/ && mv $SalmonOutputDir/libParams $SalmonOutputDir/${OutputSuffix}.salmon/
mv $SalmonOutputDir/logs $SalmonOutputDir/${OutputSuffix}.salmon/
module unload salmon

echo "*****************************************************************************"
echo "10. Quantify gene abundance using featureCounts from Subread package:"
echo "   Running featureCounts, single-end mode, for all reads and unique MNase RPF reads separately"
echo "   featureCounts does not count reads overlapping with more than one feature"
echo "   feature_counts/*gene.txt file contains tidy gene counts"
featureCounts \
    -T $CPU \
    -t exon \
    -g gene_id \
    -a $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12.geneid.gtf \
    -o $FeatureDir/${OutputSuffix}.${genome}.featureCounts.FullTable.txt \
    $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam \
    2> $FeatureDir/${OutputSuffix}.${genome}.featureCounts.log
    
featureCounts \
    -T $CPU \
    -t exon \
    -g gene_id \
    -a $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12.geneid.gtf \
    -o $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.FullTable.txt \
    $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.unique.bam \
    2> $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.log

featureCounts \
    -T $CPU \
    -t CDS \
    -g gene_id \
    -a $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12.geneid.gtf \
    -o $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.mRNA.CDS.FullTable.txt \
    $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.unique.bam \
    2> $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.mRNA.CDS.log

awk 'NR>1' $FeatureDir/${OutputSuffix}.${genome}.featureCounts.FullTable.txt | awk 'BEGIN{OFS="\t"}{print $1,$6,$7}' > \
    $FeatureDir/${OutputSuffix}.${genome}.featureCounts.gene.txt
awk 'NR>1' $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.FullTable.txt | awk 'BEGIN{OFS="\t"}{print $1,$6,$7}' > \
    $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt
awk 'NR>1' $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.mRNA.CDS.FullTable.txt | awk 'BEGIN{OFS="\t"}{print $1,$6,$7}' > \
    $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.mRNA.CDS.gene.txt

awk '{print $1}' $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt > $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt.namesall
awk 'NR>1' $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt|awk '{print $1}' > $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt.names
awk 'NR>1' $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.mRNA.CDS.gene.txt|awk '{print $1}' > $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.mRNA.CDS.gene.txt.names
cat $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt.names $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.mRNA.CDS.gene.txt.names|sort|uniq -c|awk '$1==1'|awk '{print $2}' > \
	$FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt.names.lncRNA
python $HomeDir/bin/ReorderByList.py $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt.names.lncRNA \
  $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt.lncRNA
cat $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.mRNA.CDS.gene.txt $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt.lncRNA > \
	$FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.All.CDS.gene.txt.unorder
python $HomeDir/bin/ReorderByList.py $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.All.CDS.gene.txt.unorder $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt.namesall \
  $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.All.CDS.gene.txt
rm -rf $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt.namesall $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt.names $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt.names.lncRNA
rm -rf $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.mRNA.CDS.gene.txt.names
mv $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.gene.txt.lncRNA $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.lncRNA.gene.txt
rm -rf $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.MNase.All.CDS.gene.txt.unorder

echo "   Done running featureCounts"

echo "*****************************************************************************"
echo "11. Summarize Ribo-seq mapping results and generate bed13"
#Convert bam to bed13, without size selection
echo "   Processing bam to bed12"
bedtools bamtobed -bed12 -ed -i $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bed12
samtools view $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam
awk '{print $10}' $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.seq
#Length distribution
awk '{print $1,$10}' $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.col10
sort $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.col10|uniq|awk '{print length($2)}' | sort|uniq -c |awk 'BEGIN{OFS="\t"}{if ($2>17) print $2,$1}'|awk '$1<51' > $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis && \
    rm -rf $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.col10
awk '{print $1}' $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis > $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis.len
for i in {18..50};do echo $i;done > $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis.ref
cat $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis.len $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis.ref|sort|uniq -c |awk '$1==1'|awk '{print $2}'|awk 'BEGIN{OFS="\t"}{print $1,0}' > \
    $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis.fil
cat $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis.fil|sort -k1,1 -n > $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis.c
rm -rf $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis.len $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis.ref $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis.fil \
    && mv $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis.c $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis
echo "   Generating bed13"
paste $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bed12 $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.seq > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bed13t
python $HomeDir/bin/BED13Finalizer.py $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bed13t $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bed13 && \
    rm -rf $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bed13t $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bed12

echo "*****************************************************************************"
echo "12. Filter MNase RPF [26, 40] nt and run direct transcriptome mapping"
echo "   Filtering MNase RPF length"
awk '{if (length($13)>25 && length($13)<41) print $0}' $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bed13 > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13
cd $TxMappingDir && ln -s ../$GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13 && cd ../
python $HomeDir/bin/ConvertBED13ToUniqFasta.py $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13 $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa
echo "   Running bowtie direct transcriptome mapping, allowing 2 mismatches"
bowtie \
  -S \
  -k 1 \
  -v 2 \
  -p $CPU \
  --best \
  --no-unal \
    $HomeDir/$genome/Index/TxIndex/TxIndex \
    -f $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa \
    --un $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.unaligned \
    1> $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.sam \
    2> $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.log
echo "   Post-mapping processing"
samtools view -b $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.sam > $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bam
bedtools bamtobed -i $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bam > $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed6
samtools view $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bam | awk '{print $10}' > $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.seq
#The col7 sequence has - strand mapping reads, with reversed sequence. But we will discard them.
paste $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed6 $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.seq > $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7
awk '$6=="+"' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7 > $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus

echo "*****************************************************************************"
echo "13. Prepare data for P site training and prediction using machine learning"
if [ $PsiteModel == "unassigned" ];then
  echo "   No model file provided. Prepare features around 5end and 3end of each read"
  bedtools intersect -a $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus -b $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.mRNA.cds.bed12.tx.STOP -wa -wb -F 1 > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads
  awk '$1==$8' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads > $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res
  
  echo "   Retrieving transcript length"
  python $HomeDir/bin/RetrieveLegnthInfo.py $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res 1 $HomeDir/$genome/Sequence/${genome}.RefSeq.reduced.bed12.fa.len $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len
  awk '{OFS="\t";print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$2-8,$11-$3-8,$9-$2,$3-$10}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext
  awk '$12>=0 && $13>=0 && $14>=12 && $15>=12' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext > $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil

  echo "   Retrieving 8nt features"
  awk '{OFS="\t";print $1,$2-8,$2}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.5minus
  awk '{OFS="\t";print $1,$3,$3+8}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.3plus
  bedtools getfasta -s -name -fi $HomeDir/$genome/Sequence/${genome}.RefSeq.reduced.bed12.fa -bed $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.5minus | grep -v ">"|awk '{print toupper($0)}' > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.5minus.seq
  bedtools getfasta -s -name -fi $HomeDir/$genome/Sequence/${genome}.RefSeq.reduced.bed12.fa -bed $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.3plus | grep -v ">"|awk '{print toupper($0)}' > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.3plus.seq
  echo "   Extracting 8nt before and after each read, for training reads"
  awk '{print $7}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil |cut -c1 > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.5end
  awk '{print $7}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil |rev|cut -c1 > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.3end
  awk '{print $7}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil |cut -c2-9 > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.5plus
  awk '{print $7}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil |rev|cut -c2-9|rev > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.3minus
  paste $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.5minus.seq \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.5end $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.5plus \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.3minus $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.3end \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.3plus.seq > $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.merged
  awk '{split($16,a,"");for (i=1;i<=length(a);i++){printf a[i]" "};printf $17" ";split($18,b,"");for (i=1;i<=length(b);i++){printf b[i]" "};split($19,c,"");for (i=1;i<=length(c);i++){printf c[i]" "};printf $20" ";split($21,d,"");for (i=1;i<=length(d);i++){printf d[i]" "};print length($7)" "$14}' \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.merged > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.merged.ml
fi

echo "   Extracting 8nt before and after each read, for prediction reads"
python $HomeDir/bin/RetrieveLegnthInfo.py $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus 1 $HomeDir/$genome/Sequence/${genome}.RefSeq.reduced.bed12.fa.len $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len
awk '{OFS="\t";print $1,$2,$3,$4,$5,$6,$7,$8,$2-8,$8-$3-8}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext
awk '$9>=0 && $10>=0' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil
awk '{OFS="\t";print $1,$2-8,$2}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil > $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.5minus
awk '{OFS="\t";print $1,$3,$3+8}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil > $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.3plus
bedtools getfasta -s -name -fi $HomeDir/$genome/Sequence/${genome}.RefSeq.reduced.bed12.fa -bed $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.5minus | grep -v ">"|awk '{print toupper($0)}' > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.5minus.seq
bedtools getfasta -s -name -fi $HomeDir/$genome/Sequence/${genome}.RefSeq.reduced.bed12.fa -bed $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.3plus | grep -v ">"|awk '{print toupper($0)}' >  \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.3plus.seq
awk '{print $7}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil |cut -c1 > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.5end
awk '{print $7}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil |rev|cut -c1 > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.3end
awk '{print $7}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil |cut -c2-9 > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.5plus
awk '{print $7}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil |rev|cut -c2-9|rev > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.3minus
paste $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.5minus.seq \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.5end $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.5plus \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.3minus $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.3end \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.3plus.seq > $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged
awk '{split($11,a,"");for (i=1;i<=length(a);i++){printf a[i]" "};printf $12" ";split($13,b,"");for (i=1;i<=length(b);i++){printf b[i]" "};split($14,c,"");for (i=1;i<=length(c);i++){printf c[i]" "};printf $15" ";split($16,d,"");for (i=1;i<=length(d);i++){printf d[i]" "};print length($7)}' \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged > $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged.ml

echo "14. Learn and predict MNase treated Riboseq P sites"
if [ $PsiteModel == "unassigned" ];then
    echo "   No model file provided. Learn from the current data"
    echo "   Running random forest learning"
    Rscript $HomeDir/bin/RiboMNaseML_learn.R $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.merged.ml $TxMappingDir/${OutputSuffix}.${genome}.STOPTrainedRF 2> $TxMappingDir/${OutputSuffix}.${genome}.STOPTrainedRF.log
    if [ -s $TxMappingDir/${OutputSuffix}.${genome}.STOPTrainedRF.RFModel.35features.heatmap.pdf ];then
      mv $TxMappingDir/${OutputSuffix}.${genome}.STOPTrainedRF.RFModel*pdf $FiguresDir/
    fi  
    echo "   Model file generated: ./tx_mapping/"${OutputSuffix}.${genome}.STOPTrainedRF.TunedModel.rds
    echo "   Running prediction"
    Rscript $HomeDir/bin/RiboMNaseML_predict.R $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged.ml $TxMappingDir/${OutputSuffix}.${genome}.STOPTrainedRF.TunedModel.rds $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged.ml.predict.txt
else
    echo "   Use provided model file: "$PsiteModel
    echo "   Running prediction"
    Rscript $HomeDir/bin/RiboMNaseML_predict.R $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged.ml $PsiteModel $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged.ml.predict.txt
fi
if [ ! -s $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged.ml.predict.txt ];then
    echo "   P site prediction file is empty. Use 14nt as default P offset"
    TotalReadNum=`wc -l $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged.ml|awk '{print $1}'`
    for i in $(eval echo {1..$TotalReadNum});do echo "14";done > $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged.ml.predict.txt
fi

echo "   Getting offset nucleotide values"
awk '{print $4}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged > $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged.readid
paste $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged.readid $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged.ml.predict.txt > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged.readid.matching
python $HomeDir/bin/RetrieveLegnthInfo.py $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13 4 $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged.readid.matching $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.Poffset

echo "   Correcting P offset"
python $HomeDir/bin/Shift5endOffset_BED14RPF_MNase.py $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.Poffset $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.Poffset.XShifted
awk 'BEGIN{OFS="\t"}{if ($6=="+") print $1,$2,$2+1,$4,$5,$6,$2,$2+1,$9,"1","1","0",$13;else print $1,$3-1,$3,$4,$5,$6,$3-1,$3,$9,"1","1","0",$13}' $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.Poffset.XShifted > \
    $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.Poffset.XShifted.5end

#rm -rf $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.Poffset
mv $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.Poffset.XShifted $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.PShifted.bed13
mv $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.Poffset.XShifted.5end $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.PShifted.bed13.5end


echo "*****************************************************************************"
echo "15. Prepare files for visualization"
echo "   Fetching mRNAs annotation"
if [ $plotRNA == "unassigned" ];then
  echo "   No user provided list, plot top 20 abundant mRNAs"
  sort -rk5,5 -n $SalmonOutputDir/${OutputSuffix}.${genome}.quant.sf |awk '{print $1}' > $FiguresDir/${OutputSuffix}.${genome}.AbSortedlist
  touch $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.bed12
  for name in `cat $FiguresDir/${OutputSuffix}.${genome}.AbSortedlist`;do
    awk -v name=$name 'BEGIN{OFS="\t"}{if ($4==name) print $0}' $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.mRNA.bed12 >> \
        $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.bed12
    AggLines=`wc -l $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.bed12|awk '{print $1}'`
    if [ $AggLines -gt "19" ];then
      break
    fi
  done
else
  touch $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.bed12
  for name in `cat $plotRNA`;do
    awk -v name=$name 'BEGIN{OFS="\t"}{if ($4==name) print $0}' $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.mRNA.bed12 >> \
        $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.bed12
    AggLines=`wc -l $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.bed12|awk '{print $1}'`
  done
  ToPlot=`wc -l $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.bed12 |awk '{print $1}'`
  echo "   Plot user provided mRNA list: "$ToPlot" mRNAs"
fi
awk '{print $4}' $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.bed12 > $FiguresDir/${OutputSuffix}.${genome}.mRNAlist

echo "   Preparing annotation around AUG and STOP, default extension: 60 nt"
ExtendedLength=60
RealExtLength=$((ExtendedLength+5))
Frag=`calculating 2*$ExtendedLength+1`
python $HomeDir/bin/BED12BorderExtender.py $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.bed12 $RealExtLength $RealExtLength $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12 > /dev/null
bash $HomeDir/bin/BED12Extractor.sh -a cds -i $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12 -o $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds > /dev/null
bash $HomeDir/bin/BED12TranslateGenomePos2Tx.sh -i $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12 -t $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds -o \
    $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds.tx > /dev/null
Left=$ExtendedLength
Right=$ExtendedLength
awk -v L=$Left -v R=$Right -v Ext=$RealExtLength 'BEGIN{OFS="\t"}{print $1,$2-L,$2+L+1}' $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds.tx > \
    $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds.tx.extaug
awk -v L=$Left -v R=$Right -v Ext=$RealExtLength 'BEGIN{OFS="\t"}{print $1,$3-L-3,$3+L+1-3}' $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds.tx > \
    $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds.tx.extstop
awk -v L=$Left -v R=$Right -v Ext=$RealExtLength 'BEGIN{OFS="\t"}{print $1,$2-L,$3+L}' $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds.tx > \
    $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds.tx.extboth
bash $HomeDir/bin/BED12TranslateTx2GenomePos.sh -i $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12 -t $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds.tx.extaug -o \
    $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.AUG60ext.bed12 > /dev/null
bash $HomeDir/bin/BED12TranslateTx2GenomePos.sh -i $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12 -t $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds.tx.extstop -o \
    $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.STOP60ext.bed12 > /dev/null
bash $HomeDir/bin/BED12TranslateTx2GenomePos.sh -i $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12 -t $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds.tx.extboth -o \
    $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.Both60ext.bed12 > /dev/null
rm -rf $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12 $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds.tx
rm -rf $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds.tx.ext

echo "*****************************************************************************"
echo "16. Plotting figures..."
LineNum=`wc -l $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.AUG60ext.bed12|awk '{print $1}'`
#Plot four figures for each transcript: Full with intron, Tx without intron, Around AUG, Around STOP.
OutputPrefix_Full=$FiguresDir/${OutputSuffix}.mRNA
OutputPrefix_Both=$FiguresDir/${OutputSuffix}.mRNA.Both
OutputPrefix_AUG=$FiguresDir/${OutputSuffix}.mRNA.AUG
OutputPrefix_STOP=$FiguresDir/${OutputSuffix}.mRNA.STOP
Data5endRPF=$GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.PShifted.bed13.5end
DataFullRPF=$GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.PShifted.bed13
Norm=1
if [ $normCDS == "1" ];then
  echo "   Use CDS reads as normalization factor"
  bash $HomeDir/bin/BED12Extractor.sh -a cds -i $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.mRNA.bed12 -o $GenomeMappingDir/${genome}.RefSeq.reduced.mRNA.cds.bed12 > /dev/null
  bedtools intersect -s -u -wa -split -b $GenomeMappingDir/${genome}.RefSeq.reduced.mRNA.cds.bed12 -a $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.PShifted.bed13 > \
      $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.PShifted.bed13.cds
  rm -rf $GenomeMappingDir/${genome}.RefSeq.reduced.mRNA.cds.bed12
  awk '{print $4}' $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.PShifted.bed13.cds|sort|uniq > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.PShifted.bed13.cds.id
  Reads=`wc -l $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.PShifted.bed13.cds.id|awk '{print $1}'`
  echo "   Reads mapped to CDS (multi-mapping counts once): "$Reads
  Norm=`calculating $Reads/1000000`
  if [ $Reads -lt "100" ];then
    echo "   Too few reads (<100)... Cancel normalization step."
    Norm=1
  fi
  rm -rf $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.PShifted.bed13.cds.id
else
  echo "   Not normalized to CDS"
fi
if [ $normFactor != "unassigned" ];then
  echo "   User defined normalization factor exists"
  Norm=$normFactor
fi
echo "   Normalization factor for further analysis: "$Norm
echo -e "Normalization Factor:\t${Norm}" >> $TABLE

echo "   Plotting length distribution of mapped Ribo-seq reads"
Rscript $HomeDir/bin/RiboPlot.R lendis $GenomeMappingDir/${OutputSuffix}.${genome}.Lendis $Norm $FiguresDir/${OutputSuffix}.${genome}.Lendis.pdf
ScalingPlus=`calculating 1/$Norm`
ScalingMinus=`calculating -1/$Norm`
rm -rf ${OutputPrefix_AUG}.bedGraph.summary ${OutputPrefix_AUG}.bedGraph.summary.t && touch ${OutputPrefix_AUG}.bedGraph.summary.t
rm -rf ${OutputPrefix_STOP}.bedGraph.summary ${OutputPrefix_STOP}.bedGraph.summary.t && touch ${OutputPrefix_STOP}.bedGraph.summary.t
#Generate bigWig format for the data, raw reads:
echo "   Generating track files"
bedtools genomecov -bga -split -strand + -i ${Data5endRPF} -scale $ScalingPlus -g $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${Data5endRPF}.plus.bedGraph
bedtools genomecov -bga -split -strand - -i ${Data5endRPF} -scale $ScalingMinus -g $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${Data5endRPF}.minus.bedGraph
bedtools genomecov -bga -split -strand + -i ${DataFullRPF} -scale $ScalingPlus -g $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${DataFullRPF}.plus.bedGraph
bedtools genomecov -bga -split -strand - -i ${DataFullRPF} -scale $ScalingMinus -g $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${DataFullRPF}.minus.bedGraph
$HomeDir/bin/bedGraphToBigWig ${Data5endRPF}.plus.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${Data5endRPF}.plus.bedGraph.bw
$HomeDir/bin/bedGraphToBigWig ${Data5endRPF}.minus.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${Data5endRPF}.minus.bedGraph.bw
$HomeDir/bin/bedGraphToBigWig ${DataFullRPF}.plus.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${DataFullRPF}.plus.bedGraph.bw
$HomeDir/bin/bedGraphToBigWig ${DataFullRPF}.minus.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${DataFullRPF}.minus.bedGraph.bw

touch $FiguresDir/${OutputSuffix}.${genome}.mRNAlistgene
echo "   Plotting quadruple figures for individual transcripts"
for name in `cat $FiguresDir/${OutputSuffix}.${genome}.mRNAlist`
do
    #Get gene name
    GeneName=`awk -v TxName=$name '{OFS="\t";if ($1==TxName) print $2}' $HomeDir/$genome/Annotation/${genome}.uniqMatching.txt`
    echo -e $name"\t"$GeneName
    echo -e $name"\t"$GeneName >> $FiguresDir/${OutputSuffix}.${genome}.mRNAlistgene
    Curr=`awk -v GeneName=$name '{OFS="\t";if ($4==GeneName) print $0}' $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.bed12`
    mRNAStart=`echo $Curr|awk '{print $2}'`
    mRNAEnd=`echo $Curr|awk '{print $3}'`
    awk -v GeneName=$name 'BEGIN{OFS="\t"}{if ($4==GeneName) print $0}' $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.Both60ext.bed12 > ${OutputPrefix_Both}.${name}.bed12
    awk -v GeneName=$name 'BEGIN{OFS="\t"}{if ($4==GeneName) print $0}' $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.AUG60ext.bed12 > ${OutputPrefix_AUG}.${name}.bed12
    bedtools getfasta -s -split -fi $HomeDir/$genome/Sequence/${genome}.fa -bed ${OutputPrefix_AUG}.${name}.bed12 -fo ${OutputPrefix_AUG}.${name}.bed12.fa
    awk -v GeneName=$name 'BEGIN{OFS="\t"}{if ($4==GeneName) print $0}' $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.STOP60ext.bed12 > ${OutputPrefix_STOP}.${name}.bed12
    bedtools getfasta -s -split -fi $HomeDir/$genome/Sequence/${genome}.fa -bed ${OutputPrefix_STOP}.${name}.bed12 -fo ${OutputPrefix_STOP}.${name}.bed12.fa

    #Check chromosome, and filter out range                                                                                                                                                                                     
    TempChr=`awk '{print $1}' ${OutputPrefix_AUG}.${name}.bed12`
    awk -v chr=$TempChr 'BEGIN{OFS="\t"}{if ($1==chr) print $0}' $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${OutputPrefix_Both}.${name}.chr
    awk -v chr=$TempChr 'BEGIN{OFS="\t"}{if ($1==chr) print $0}' $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${OutputPrefix_AUG}.${name}.chr
    awk -v chr=$TempChr 'BEGIN{OFS="\t"}{if ($1==chr) print $0}' $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${OutputPrefix_STOP}.${name}.chr

    #Intersect reads to get only overlapped reads                                                                                                                                                                               
    bedtools intersect -wa -a ${Data5endRPF} -b ${OutputPrefix_Both}.${name}.bed12 > ${OutputPrefix_Both}.${name}.RPF
    bedtools intersect -wa -a ${Data5endRPF} -b ${OutputPrefix_AUG}.${name}.bed12 > ${OutputPrefix_AUG}.${name}.RPF
    bedtools intersect -wa -a ${Data5endRPF} -b ${OutputPrefix_STOP}.${name}.bed12 > ${OutputPrefix_STOP}.${name}.RPF

    #Get bedGraph from bw files
    $HomeDir/bin/bigWigToBedGraph ${Data5endRPF}.plus.bedGraph.bw -chrom=$TempChr -start=$mRNAStart -end=$mRNAEnd ${OutputPrefix_Full}.${name}.plus.bedGraph
    $HomeDir/bin/bigWigToBedGraph ${Data5endRPF}.minus.bedGraph.bw -chrom=$TempChr -start=$mRNAStart -end=$mRNAEnd ${OutputPrefix_Full}.${name}.minus.bedGraph
    #Get bedGraph for each transcript                                                                                                                                                                                           
    bedtools genomecov -bga -split -strand + -i ${OutputPrefix_Both}.${name}.RPF -scale $ScalingPlus -g ${OutputPrefix_Both}.${name}.chr > ${OutputPrefix_Both}.${name}.plus.bedGraph
    bedtools genomecov -bga -split -strand - -i ${OutputPrefix_Both}.${name}.RPF -scale $ScalingMinus -g ${OutputPrefix_Both}.${name}.chr > ${OutputPrefix_Both}.${name}.minus.bedGraph
    bedtools genomecov -bga -split -strand + -i ${OutputPrefix_AUG}.${name}.RPF -scale $ScalingPlus -g ${OutputPrefix_AUG}.${name}.chr > ${OutputPrefix_AUG}.${name}.plus.bedGraph
    bedtools genomecov -bga -split -strand - -i ${OutputPrefix_AUG}.${name}.RPF -scale $ScalingMinus -g ${OutputPrefix_AUG}.${name}.chr > ${OutputPrefix_AUG}.${name}.minus.bedGraph
    bedtools genomecov -bga -split -strand + -i ${OutputPrefix_STOP}.${name}.RPF -scale $ScalingPlus -g ${OutputPrefix_STOP}.${name}.chr > ${OutputPrefix_STOP}.${name}.plus.bedGraph
    bedtools genomecov -bga -split -strand - -i ${OutputPrefix_STOP}.${name}.RPF -scale $ScalingMinus -g ${OutputPrefix_STOP}.${name}.chr > ${OutputPrefix_STOP}.${name}.minus.bedGraph

    #The worst case, if there are no reads in the given region, fill 0:
    if [ ! -s ${OutputPrefix_Both}.${name}.plus.bedGraph ];then
        awk '{OFS="\t";print $1,"0",$2,"0"}' ${OutputPrefix_Both}.${name}.chr > ${OutputPrefix_Both}.${name}.plus.bedGraph
    fi
    if [ ! -s ${OutputPrefix_Both}.${name}.minus.bedGraph ];then
        awk '{OFS="\t";print $1,"0",$2,"-0"}' ${OutputPrefix_Both}.${name}.chr > ${OutputPrefix_Both}.${name}.minus.bedGraph
    fi

    if [ ! -s ${OutputPrefix_AUG}.${name}.plus.bedGraph ];then
        awk '{OFS="\t";print $1,"0",$2,"0"}' ${OutputPrefix_AUG}.${name}.chr > ${OutputPrefix_AUG}.${name}.plus.bedGraph
    fi
    if [ ! -s ${OutputPrefix_AUG}.${name}.minus.bedGraph ];then
        awk '{OFS="\t";print $1,"0",$2,"-0"}' ${OutputPrefix_AUG}.${name}.chr > ${OutputPrefix_AUG}.${name}.minus.bedGraph
    fi
  
    if [ ! -s ${OutputPrefix_STOP}.${name}.plus.bedGraph ];then
        awk '{OFS="\t";print $1,"0",$2,"0"}' ${OutputPrefix_STOP}.${name}.chr > ${OutputPrefix_STOP}.${name}.plus.bedGraph
    fi
    if [ ! -s ${OutputPrefix_STOP}.${name}.minus.bedGraph ];then
        awk '{OFS="\t";print $1,"0",$2,"-0"}' ${OutputPrefix_STOP}.${name}.chr > ${OutputPrefix_STOP}.${name}.minus.bedGraph
    fi

    #Calculate single base coverage
    python $HomeDir/bin/BedGraphRegionFiller.py ${OutputPrefix_Full}.${name}.plus.bedGraph ${OutputPrefix_Full}.${name}.plus.bedGraph.bb
    python $HomeDir/bin/BedGraphRegionFiller.py ${OutputPrefix_Full}.${name}.minus.bedGraph ${OutputPrefix_Full}.${name}.minus.bedGraph.bb
    python $HomeDir/bin/BedGraphDirectTxRegionExtractor.py ${OutputPrefix_Both}.${name}.bed12 ${OutputPrefix_Both}.${name}.plus.bedGraph ${OutputPrefix_Both}.${name}.plus.bedGraph.bb
    python $HomeDir/bin/BedGraphDirectTxRegionExtractor.py ${OutputPrefix_Both}.${name}.bed12 ${OutputPrefix_Both}.${name}.minus.bedGraph ${OutputPrefix_Both}.${name}.minus.bedGraph.bb
    python $HomeDir/bin/BedGraphDirectTxRegionExtractor.py ${OutputPrefix_AUG}.${name}.bed12 ${OutputPrefix_AUG}.${name}.plus.bedGraph ${OutputPrefix_AUG}.${name}.plus.bedGraph.bb
    python $HomeDir/bin/BedGraphDirectTxRegionExtractor.py ${OutputPrefix_AUG}.${name}.bed12 ${OutputPrefix_AUG}.${name}.minus.bedGraph ${OutputPrefix_AUG}.${name}.minus.bedGraph.bb
    python $HomeDir/bin/BedGraphDirectTxRegionExtractor.py ${OutputPrefix_STOP}.${name}.bed12 ${OutputPrefix_STOP}.${name}.plus.bedGraph ${OutputPrefix_STOP}.${name}.plus.bedGraph.bb
    python $HomeDir/bin/BedGraphDirectTxRegionExtractor.py ${OutputPrefix_STOP}.${name}.bed12 ${OutputPrefix_STOP}.${name}.minus.bedGraph ${OutputPrefix_STOP}.${name}.minus.bedGraph.bb

    #Get sense coverage (no matter strand direction, always report 5'-end to 3'-end coverage).                                                         
    python $HomeDir/bin/BedGraphDirectTxRegionSenseCoverage.py ${OutputPrefix_Both}.${name}.bed12 ${OutputPrefix_Both}.${name}.plus.bedGraph.bb ${OutputPrefix_Both}.${name}.minus.bedGraph.bb ${OutputPrefix_Both}.${name}.sense.bedGraph
    python $HomeDir/bin/BedGraphDirectTxRegionSenseCoverage.py ${OutputPrefix_AUG}.${name}.bed12 ${OutputPrefix_AUG}.${name}.plus.bedGraph.bb ${OutputPrefix_AUG}.${name}.minus.bedGraph.bb ${OutputPrefix_AUG}.${name}.sense.bedGraph
    python $HomeDir/bin/BedGraphDirectTxRegionSenseCoverage.py ${OutputPrefix_STOP}.${name}.bed12 ${OutputPrefix_STOP}.${name}.plus.bedGraph.bb ${OutputPrefix_STOP}.${name}.minus.bedGraph.bb ${OutputPrefix_STOP}.${name}.sense.bedGraph

    #Cleaning up, you can optimize the items to be deleted here
    rm -rf ${OutputPrefix_Both}.${name}.chr ${OutputPrefix_Both}.${name}.RPF
    rm -rf ${OutputPrefix_Both}.${name}.plus.bedGraph ${OutputPrefix_Both}.${name}.minus.bedGraph
    rm -rf ${OutputPrefix_AUG}.${name}.chr ${OutputPrefix_AUG}.${name}.RPF
    rm -rf ${OutputPrefix_AUG}.${name}.plus.bedGraph ${OutputPrefix_AUG}.${name}.minus.bedGraph
    rm -rf ${OutputPrefix_STOP}.${name}.chr ${OutputPrefix_STOP}.${name}.RPF
    rm -rf ${OutputPrefix_STOP}.${name}.plus.bedGraph ${OutputPrefix_STOP}.${name}.minus.bedGraph

    #Check the output file, if empty, add one column of Nt
    if [ ! -s ${OutputPrefix_AUG}.bedGraph.summary.t ];then
        #Get the Ns
        TotalLines=`wc -l ${OutputPrefix_AUG}.${name}.sense.bedGraph|awk '{print $1}'`
        Extend=`calculating $TotalLines/2-0.5`
        Right=`calculating $Extend-3+1`
        echo "Nucleotide" > ${OutputPrefix_AUG}.bedGraph.summary.t
        for i in $(eval echo {1..$Extend});do echo "N";done >> ${OutputPrefix_AUG}.bedGraph.summary.t
        echo "A" >> ${OutputPrefix_AUG}.bedGraph.summary.t
        echo "T" >> ${OutputPrefix_AUG}.bedGraph.summary.t
        echo "G" >> ${OutputPrefix_AUG}.bedGraph.summary.t
        for i in $(eval echo {1..$Right});do echo "N";done >> ${OutputPrefix_AUG}.bedGraph.summary.t
    fi
    if [ ! -s ${OutputPrefix_STOP}.bedGraph.summary.t ];then
        #Get the Ns                                                                                                                                                                                                         
        TotalLines=`wc -l ${OutputPrefix_STOP}.${name}.sense.bedGraph|awk '{print $1}'`
        Extend=`calculating $TotalLines/2-0.5`
        Right=`calculating $Extend-3+1`
        echo "Nucleotide" > ${OutputPrefix_STOP}.bedGraph.summary.t
        for i in $(eval echo {1..$Extend});do echo "N";done >> ${OutputPrefix_STOP}.bedGraph.summary.t
        echo "E" >> ${OutputPrefix_STOP}.bedGraph.summary.t
        echo "N" >> ${OutputPrefix_STOP}.bedGraph.summary.t
        echo "D" >> ${OutputPrefix_STOP}.bedGraph.summary.t
        for i in $(eval echo {1..$Right});do echo "N";done >> ${OutputPrefix_STOP}.bedGraph.summary.t
    fi

    #Summarize ATG
    awk '{print $4}' ${OutputPrefix_AUG}.${name}.sense.bedGraph > ${OutputPrefix_AUG}.${name}.sense.bedGraph.temp
    echo ${name} > ${OutputPrefix_AUG}.${name}.sense.bedGraph.name
    cat ${OutputPrefix_AUG}.${name}.sense.bedGraph.name ${OutputPrefix_AUG}.${name}.sense.bedGraph.temp > ${OutputPrefix_AUG}.${name}.sense.bedGraph.withname
    paste ${OutputPrefix_AUG}.bedGraph.summary.t ${OutputPrefix_AUG}.${name}.sense.bedGraph.withname > ${OutputPrefix_AUG}.bedGraph.summary
    rm -rf ${OutputPrefix_AUG}.bedGraph.summary.t
    mv ${OutputPrefix_AUG}.bedGraph.summary ${OutputPrefix_AUG}.bedGraph.summary.t
    rm -rf ${OutputPrefix_AUG}.${name}.sense.bedGraph.temp ${OutputPrefix_AUG}.${name}.sense.bedGraph.withname ${OutputPrefix_AUG}.${name}.sense.bedGraph.name
    #clean up
    rm -rf ${OutputPrefix_AUG}.${name}.bed12 ${OutputPrefix_AUG}.${name}.plus.bedGraph.bb ${OutputPrefix_AUG}.${name}.minus.bedGraph.bb
    
    #Summarize STOP
    awk '{print $4}' ${OutputPrefix_STOP}.${name}.sense.bedGraph > ${OutputPrefix_STOP}.${name}.sense.bedGraph.temp
    echo ${name} > ${OutputPrefix_STOP}.${name}.sense.bedGraph.name
    cat ${OutputPrefix_STOP}.${name}.sense.bedGraph.name ${OutputPrefix_STOP}.${name}.sense.bedGraph.temp > ${OutputPrefix_STOP}.${name}.sense.bedGraph.withname
    paste ${OutputPrefix_STOP}.bedGraph.summary.t ${OutputPrefix_STOP}.${name}.sense.bedGraph.withname > ${OutputPrefix_STOP}.bedGraph.summary
    rm -rf ${OutputPrefix_STOP}.bedGraph.summary.t
    mv ${OutputPrefix_STOP}.bedGraph.summary ${OutputPrefix_STOP}.bedGraph.summary.t
    rm -rf ${OutputPrefix_STOP}.${name}.sense.bedGraph.temp ${OutputPrefix_STOP}.${name}.sense.bedGraph.withname ${OutputPrefix_STOP}.${name}.sense.bedGraph.name
    #clean up
    rm -rf ${OutputPrefix_STOP}.${name}.bed12 ${OutputPrefix_STOP}.${name}.plus.bedGraph.bb ${OutputPrefix_STOP}.${name}.minus.bedGraph.bb
    
    #Clean up Both
    rm -rf ${OutputPrefix_Both}.${name}.plus.bedGraph.bb ${OutputPrefix_Both}.${name}.minus.bedGraph.bb ${OutputPrefix_Both}.${name}.bed12 
    
    #Clean up Full
    rm -rf ${OutputPrefix_Full}.${name}.plus.bedGraph ${OutputPrefix_Full}.${name}.minus.bedGraph
    
    #Plotting figures
    Rscript --vanilla $HomeDir/bin/RiboPlot.R quadplot ${OutputPrefix_Full} ${name} $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.bed12 ${GeneName} >> $FiguresDir/${OutputSuffix}.${genome}.mRNAlistgene $Norm
    
    #Last clean up, for debugging use
    rm -rf ${OutputPrefix_AUG}.${name}.bed12.fa ${OutputPrefix_STOP}.${name}.bed12.fa 
    rm -rf ${OutputPrefix_Both}.${name}.sense.bedGraph ${OutputPrefix_AUG}.${name}.sense.bedGraph ${OutputPrefix_STOP}.${name}.sense.bedGraph
    rm -rf ${OutputPrefix_Full}.${name}.plus.bedGraph.bb ${OutputPrefix_Full}.${name}.minus.bedGraph.bb
done

if [ $normCDS == "1" ];then
  #Rename the files that have been normalized:
  mv ${Data5endRPF}.plus.bedGraph ${Data5endRPF}.plus.NormCDS.bedGraph
  mv ${Data5endRPF}.minus.bedGraph ${Data5endRPF}.minus.NormCDS.bedGraph
  mv ${DataFullRPF}.plus.bedGraph ${DataFullRPF}.plus.NormCDS.bedGraph
  mv ${DataFullRPF}.minus.bedGraph ${DataFullRPF}.minus.NormCDS.bedGraph
  mv ${Data5endRPF}.plus.bedGraph.bw ${Data5endRPF}.plus.NormCDS.bedGraph.bw
  mv ${Data5endRPF}.minus.bedGraph.bw ${Data5endRPF}.minus.NormCDS.bedGraph.bw
  mv ${DataFullRPF}.plus.bedGraph.bw ${DataFullRPF}.plus.NormCDS.bedGraph.bw
  mv ${DataFullRPF}.minus.bedGraph.bw ${DataFullRPF}.minus.NormCDS.bedGraph.bw
  echo "   Generating track files using raw data"
  bedtools genomecov -bga -split -strand + -i ${Data5endRPF} -scale 1 -g $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${Data5endRPF}.plus.bedGraph
  bedtools genomecov -bga -split -strand - -i ${Data5endRPF} -scale -1 -g $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${Data5endRPF}.minus.bedGraph
  bedtools genomecov -bga -split -strand + -i ${DataFullRPF} -scale 1 -g $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${DataFullRPF}.plus.bedGraph
  bedtools genomecov -bga -split -strand - -i ${DataFullRPF} -scale -1 -g $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${DataFullRPF}.minus.bedGraph
  $HomeDir/bin/bedGraphToBigWig ${Data5endRPF}.plus.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${Data5endRPF}.plus.bedGraph.bw
  $HomeDir/bin/bedGraphToBigWig ${Data5endRPF}.minus.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${Data5endRPF}.minus.bedGraph.bw
  $HomeDir/bin/bedGraphToBigWig ${DataFullRPF}.plus.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${DataFullRPF}.plus.bedGraph.bw
  $HomeDir/bin/bedGraphToBigWig ${DataFullRPF}.minus.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${DataFullRPF}.minus.bedGraph.bw
  rm -rf ${Data5endRPF}.plus.bedGraph ${Data5endRPF}.minus.bedGraph ${DataFullRPF}.plus.bedGraph ${DataFullRPF}.minus.bedGraph
fi
if [ $normFactor != "unassigned" ];then
  #Rename the files that have been normalized:
  mv ${Data5endRPF}.plus.bedGraph ${Data5endRPF}.plus.Norm.bedGraph
  mv ${Data5endRPF}.minus.bedGraph ${Data5endRPF}.minus.Norm.bedGraph
  mv ${DataFullRPF}.plus.bedGraph ${DataFullRPF}.plus.Norm.bedGraph
  mv ${DataFullRPF}.minus.bedGraph ${DataFullRPF}.minus.Norm.bedGraph
  mv ${Data5endRPF}.plus.bedGraph.bw ${Data5endRPF}.plus.Norm.bedGraph.bw
  mv ${Data5endRPF}.minus.bedGraph.bw ${Data5endRPF}.minus.Norm.bedGraph.bw
  mv ${DataFullRPF}.plus.bedGraph.bw ${DataFullRPF}.plus.Norm.bedGraph.bw
  mv ${DataFullRPF}.minus.bedGraph.bw ${DataFullRPF}.minus.Norm.bedGraph.bw
  echo "   Generating track files using raw data"
  bedtools genomecov -bga -split -strand + -i ${Data5endRPF} -scale 1 -g $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${Data5endRPF}.plus.bedGraph
  bedtools genomecov -bga -split -strand - -i ${Data5endRPF} -scale -1 -g $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${Data5endRPF}.minus.bedGraph
  bedtools genomecov -bga -split -strand + -i ${DataFullRPF} -scale 1 -g $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${DataFullRPF}.plus.bedGraph
  bedtools genomecov -bga -split -strand - -i ${DataFullRPF} -scale -1 -g $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${DataFullRPF}.minus.bedGraph
  $HomeDir/bin/bedGraphToBigWig ${Data5endRPF}.plus.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${Data5endRPF}.plus.bedGraph.bw
  $HomeDir/bin/bedGraphToBigWig ${Data5endRPF}.minus.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${Data5endRPF}.minus.bedGraph.bw
  $HomeDir/bin/bedGraphToBigWig ${DataFullRPF}.plus.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${DataFullRPF}.plus.bedGraph.bw
  $HomeDir/bin/bedGraphToBigWig ${DataFullRPF}.minus.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${DataFullRPF}.minus.bedGraph.bw
  rm -rf ${Data5endRPF}.plus.bedGraph ${Data5endRPF}.minus.bedGraph ${DataFullRPF}.plus.bedGraph ${DataFullRPF}.minus.bedGraph
fi

echo "*****************************************************************************"
echo "17. Generating meta-figures"
mv ${OutputPrefix_AUG}.bedGraph.summary.t ${OutputPrefix_AUG}.bedGraph.summary.txt
mv ${OutputPrefix_STOP}.bedGraph.summary.t ${OutputPrefix_STOP}.bedGraph.summary.txt
Rscript --vanilla $HomeDir/bin/RiboPlot.R metaAUG ${OutputPrefix_AUG}
Rscript --vanilla $HomeDir/bin/RiboPlot.R metaSTOP ${OutputPrefix_STOP}

echo "*****************************************************************************"
echo "18. Finishing:"
#clean up fastq files:
rm -rf $GenomeMappingDir/${OutputSuffix}.No_rRNA.fastq.gz $GenomeMappingDir/${OutputSuffix}.No_rRNA.No_miRNA.fastq.gz
#Remove mapping files:
rm -rf $GenomeMappingDir/${OutputSuffix}.${genome}.Aligned.out.bam $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.seq
rm -rf $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.sam $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bam
#Remove Tx direct mapping related files:
rm -rf $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed6 $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7
rm -rf $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext
rm -rf $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.[35]*
rm -rf $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.len.ext.fil.merged.readid*
rm -rf $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res
rm -rf $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext
rm -rf $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil
rm -rf $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.bed7.plus.STOPreads.res.len.ext.fil.[35]* $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.tx.seq
rm -rf $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa.unaligned $TxMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13.fa
#If you want, you can keep the following files:
#rm -rf $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bed13 $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.RPF.bed12 $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.MNase.bed13
mv $GenomeMappingDir/${OutputSuffix}.*bedGraph* $TracksDir/
#You can also keep the following annotation files for plotting:
rm -rf $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds.tx.extaug $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds.tx.extboth $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.ext.bed12.cds.tx.extstop
rm -rf $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.AUG60ext.bed12 $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.STOP60ext.bed12 $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.Both60ext.bed12
rm -rf $FiguresDir/${OutputSuffix}.${genome}.mRNAlist.bed12 $FiguresDir/${OutputSuffix}.${genome}.mRNAlist $FiguresDir/${OutputSuffix}.${genome}.AbSortedlist
rm -rf Rplots.pdf

echo "Time Started: "$St
Ed=`date`
echo "Time Ended:   "$Ed
echo "*                           End of the pipeline                             *"
echo "*****************************************************************************"
