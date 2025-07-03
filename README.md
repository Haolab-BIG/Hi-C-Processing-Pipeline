# Hi-C-Processing-Pipeline
This is the pipeline of downstream analysis for HiC data, including quality contrl (QC) and TAD, Loop, Compartment identification. 
# Part I Introduction
## i. Workflow
Here stands an throughout workflow of Hi-C data analysis.
![hicprocedure](https://github.com/user-attachments/assets/08f722df-ef5d-4cd8-9f2a-fd163913ecee)

As illustrated in the figure,
(i) Solid-line boxes represent the files we input/output;
(ii) Dashed-line boxes represent the QC results after processing at each step.

We will proceed by structuring our workflow according to (ii).

## File Structure
```
Hi-C
├─ 1.rawdata
│    ├─ *rep1_R1_001.fastq.gz
│    ├─ *rep1_R2_001.fastq.gz
│    ├─ *rep2_R1_001.fastq.gz
│    └─ *rep2_R2_001.fastq.gz
├─ 2.trim
│    ├─ *rep1_R1_trimed.fastq
│    ├─ *rep1_R2_trimed.fastq
│    ├─ *rep2_R1_trimed.fastq
│    └─ *rep2_R2_trimed.fastq
├─ 3.bowtie 
│    ├─ *rep1_R1.hg38.sorted.bam
│    ├─ *rep1_R2.hg38.sorted.bam
│    ├─ *rep2_R1.hg38.sorted.bam
│    ├─ *rep2_R2.hg38.sorted.bam
│    ├─ *mergedR1.bam
│    └─ *mergedR2.bam
├─ 4.homertag
│    ├─ *rep1
│    ├─ *rep2
│    └─ *HicMerge
├─ 5.transform
│    ├─ multiBigWigSummary
│    ├─ PCA 
│    └─ correlation
├─ 6.OnTAD
│    ├─ multiBigWigSummary
│    ├─ PCA 
│    └─ correlation
└─ 7.loop
       ├─ computeMatrix 
       └─ heatmap
```

## iii. Software Required
```
homer
fastqc
bowtie2
straw
bedtools
3DChromatin_ReplicateQC
samtools
deeptools
OnTAD
fancplot
fithic2
fanc
hicexplorer
```
# Part II Codes
## Step 1: Preprocessing
In this section, you will convert the raw FASTQ files into .hic files, which can be used for subsequent analysis.
```
workdir=/2.trim
botiemap=/3.bowtie
homertag=/4.homertag
rawdatadir=/1.rawdata
bowtie2index=/bowtie2_index   ## Prepare in advance
cd ${rawdatadir}
for fileName in `ls -1 *.fastq.gz | awk -F_ '{NF-=2; print $0}' OFS=_ | sort | uniq`;do
	echo ${fileName}
	echo "fastq trimming and reads alignment"
	cd ${workdir}
	homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 ${rawdatadir}/${fileName}_R1_001.fastq.gz > ${workdir}/${fileName}_R1_trim.txt
	homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 ${rawdatadir}/${fileName}_R2_001.fastq.gz > ${workdir}/${fileName}_R2_trim.txt
	mv ${rawdatadir}/${fileName}_R1_001.fastq.gz.trimmed ${workdir}/${fileName}_R1_trimed.fastq
	mv ${rawdatadir}/${fileName}_R2_001.fastq.gz.trimmed ${workdir}/${fileName}_R2_trimed.fastq
	fastqc -t 20 ${workdir}/${fileName}_R1_trimed.fastq -o ${workdir}
	fastqc -t 20 ${workdir}/${fileName}_R2_trimed.fastq -o ${workdir}
	filenameshort=$(echo ${fileName} | rev | cut -d"_" -f3- | rev)
	bowtie2 -p 20 -x ${bowtie2index} -U ${workdir}/${fileName}_R1_trimed.fastq | samtools view -b - | samtools sort - -@ 20 -o ${botiemap}/${filenameshort}_R1.hg38.sorted.bam
	bowtie2 -p 20 -x ${bowtie2index} -U ${workdir}/${fileName}_R2_trimed.fastq | samtools view -b - | samtools sort - -@ 20 -o ${botiemap}/${filenameshort}_R2.hg38.sorted.bam
	samtools index ${botiemap}/${filenameshort}_R1.hg38.sorted.bam
	samtools index ${botiemap}/${filenameshort}_R2.hg38.sorted.bam
	echo "Create Hi-C Tag Directory with HOMER"
 	makeTagDirectory ${homertag}/HicExp1TagDir_${fileName}/ ${botiemap}/${filenameshort}_R1.hg38.sorted.bam,${botiemap}/${filenameshort}_R2.hg38.sorted.bam -tbp 1
 	echo "Check visulization files for juicer"
 	tagDir2hicFile.pl ${homertag}/HicExp1TagDir_${fileName}/ -juicer auto -genome hg38 -p 50
done
```

## Step 2: Combine deplications
Check duplication correlation
```
#### contact matrix
matrixdir=/5.transform
homerdir=/4.homertag
for fileName in XXX ;do
	mkdir -p ${matrixdir}/${fileName}
	for chrA in {1..22} X;do
	  if [[ "$chrA" == "X" ]]; then
	    chrB_list="X Y"
	  elif [[ "$chrA" == "Y" ]]; then
	    chrB_list="Y"
	  else
	    chrB_list="$(seq $chrA 22) X Y"
	  fi
	  for chrB in $chrB_list;do
	    straw NONE ${homerdir}/${fileName}/${fileName}.hic ${chrA} ${chrB} BP 100000 | awk -v chrA="$chrA" -v chrB="$chrB" '{print chrA "\t" $1 "\t" chrB "\t" $2 "\t" $3}' > ${matrixdir}/${fileName}/${fileName}.chr.${chrA}.${chrB}.100kb.txt
	  done
	done
	cat ${matrixdir}/${fileName}/${fileName}.chr.* > ${matrixdir}/${fileName}/${fileName}.merged.100kb.txt
	gzip ${matrixdir}/${fileName}/${fileName}.merged.100kb.txt
	rm ${matrixdir}/${fileName}/${fileName}.chr.*
done

#### run 3DChromatin_ReplicateQC
homerdir=/4.homertag
transformdir=/5.transform
genomesize=XX    ## Prepare in advance
bedtools makewindows -g ${genomesize} -w 100000 | awk '{print $1, $2, $3, $2}' OFS="\t" - > ${transformdir}/Bins.w100000.bed
gzip ${transformdir}/Bins.w100000.bed
for filename in XXX;
do
	3DChromatin_ReplicateQC run_all --metadata_samples ${transformdir}/${filename}.metadata.samples --metadata_pairs ${transformdir}/${filename}.metadata.pairs --bins ${transformdir}/Bins.w100000.bed.gz --outdir ${transformdir}/${filename}.ReplicateQC
done
```
Combine deplications
```
botiemap=/3.bowtie
homertag=/4.homertag
for filenameshort in XX ;do
	echo ${filenameshort3}
	samtools merge -@ 20 ${botiemap}/${filenameshort}_mergedR1.bam ${botiemap}/${filenameshort}_rep1_R1.hg38.sorted.bam ${botiemap}/${filenameshort}_rep2_R1.hg38.sorted.bam
	samtools merge -@ 20 ${botiemap}/${filenameshort}_mergedR2.bam ${botiemap}/${filenameshort}_rep1_R2.hg38.sorted.bam ${botiemap}/${filenameshort}_rep2_R2.hg38.sorted.bam
	cd ${findtad}
	echo "Create Hi-C Tag Directory with HOMER"
	makeTagDirectory ${homertag}/HicMerge_${filenameshort}/ ${botiemap}/${filenameshort}_mergedR1.bam,${botiemap}/${filenameshort}_mergedR2.bam -tbp 1
	echo "Check visulization files for juicer"
	tagDir2hicFile.pl ${homertag}/HicMerge_${filenameshort}/ -juicer auto -genome hg38 -p 50
done
```
## Step 3: Identifying TAD and loop
Identifying TAD
```
onTADdir=XX    ### path to onTAD software
ontadout=/6.OnTAD
homerdir=/4.homertag
genomesize=XX    ## Prepare in advance
for filename in XXX ;do
	echo "Using OnTAD to do TAD identification"
	for i in  $(seq 1 22) "X" "Y";do
		temp="chr"${i}
		num=$(grep -w "${temp}" ${genomesize} |cut -f 2)
		echo ${i}
		echo ${temp}
		echo ${num}
		${onTADdir}/OnTAD ${homerdir}/${filename}/${filename}.hic -bedout ${i} ${num} 25000 -o ${ontadout}/${filename}_chr${i} -maxsz 400 -penalty 0.3
		echo 'obtain the largest regions'
		grep '56,108,176' ${ontadout}/${filename}_chr${i}.bed |awk 'BEGIN{FS=OFS="\t"}{if(NR%2==1){print $1,$2,$3,".",1,"+"}else{print $1,$2,$3,".",2,"+"}}' > ${ontadout}/${filename}_chr${i}use.bed
	done
	cat ${ontadout}/${filename}_chr*use.bed > ${ontadout}/${filename}_use.bed
	echo "draw several examples"
	fancplot --width 12  -o ${ontadout}/${filename}.chr4:50mb-60mb_25kb.pdf 4:50mb-60mb -p triangular ${homerdir}/${filename}/${filename}.hic@25kb -vmin 0 -vmax 20 -p bar -y 0 2.5 ${ontadout}/${filename}_use.bed
	fancplot --width 12  -o ${ontadout}/${filename}.chr4:60mb-90mb_25kb.pdf 4:60mb-90mb -p triangular ${homerdir}/${filename}/${filename}.hic@25kb -vmin 0 -vmax 20 -p bar -y 0 2.5 ${ontadout}/${filename}_use.bed
	fancplot --width 12  -o ${ontadout}/${filename}.chr4:140mb-170mb_25kb.pdf 4:140mb-170mb -p triangular ${homerdir}/${filename}/${filename}.hic@25kb -vmin 0 -vmax 20 -p bar -y 0 2.5 ${ontadout}/${filename}_use.bed
	fancplot --width 12  -o ${ontadout}/${filename}.chr7:50mb-68mb_25kb.pdf 7:50mb-68mb -p triangular ${homerdir}/${filename}/${filename}.hic@25kb -vmin 0 -vmax 20 -p bar -y 0 2.5 ${ontadout}/${filename}_use.bed
	fancplot --width 12  -o ${ontadout}/${filename}.chr7:148mb-158mb_25kb.pdf 7:148mb-158mb -p triangular ${homerdir}/${filename}/${filename}.hic@25kb -vmin 0 -vmax 20 -p bar -y 0 2.5 ${ontadout}/${filename}_use.bed
	fancplot --width 12  -o ${ontadout}/${filename}.chr20:6mb-17mb_25kb.pdf 20:6mb-17mb -p triangular ${homerdir}/${filename}/${filename}.hic@25kb -vmin 0 -vmax 20 -p bar -y 0 2.5 ${ontadout}/${filename}_use.bed
	fancplot --width 12  -o ${ontadout}/${filename}.chr16:36mb-50mb_25kb.pdf 16:36mb-50mb -p triangular ${homerdir}/${filename}/${filename}.hic@25kb -vmin 0 -vmax 20 -p bar -y 0 2.5 ${ontadout}/${filename}_use.bed
	fancplot --width 12  -o ${ontadout}/${filename}.chr2:2mb-20mb_25kb.pdf 2:2mb-20mb -p triangular ${homerdir}/${filename}/${filename}.hic@25kb -vmin 0 -vmax 20 -p bar -y 0 2.5 ${ontadout}/${filename}_use.bed
done
```
QC for TAD
```
ontadout=/6.OnTAD
homerdir=/4.homertag
for filename in XXX ;do
	echo "check the signal distribution along TADs"
	awk 'BEGIN{FS=OFS="\t"}{print "chr"$0}' ${ontadout}/${filename}_use.bed > ${ontadout}/${filename}_use_chr.bed
	echo "draw CTCF distribution"
	bwoutdir=CTCF.bw  ## path to CTCF.bw
	computeMatrix scale-regions -S ${bwoutdir} -R ${ontadout}/${filename}_use_chr.bed --beforeRegionStartLength 500000 --regionBodyLength 500000 --afterRegionStartLength 500000 --skipZeros -o ${ontadout}/${filename}_CTCF_TADs.mat.gz -p 20
	plotHeatmap -m ${ontadout}/${filename}_CTCF_TADs.mat.gz -out ${ontadout}/${filename}_CTCF_TADs.pdf --colorMap inferno --heatmapHeight 20 --perGroup --zMax 1.0
	echo "draw conhensin distribution"
	bwoutdir=Rad21.bw  ## path to Rad21.bw
	computeMatrix scale-regions -S ${bwoutdir} -R ${ontadout}/${filename}_use_chr.bed --beforeRegionStartLength 500000 --regionBodyLength 500000 --afterRegionStartLength 500000 --skipZeros -o ${ontadout}/${filename}_conhensin_TADs.mat.gz -p 20
	plotHeatmap -m ${ontadout}/${filename}_conhensin_TADs.mat.gz -out ${ontadout}/${filename}_conhensin_TADs.pdf --colorMap inferno --heatmapHeight 20 --perGroup --
	echo "Check the TAD strength"
	source activate /mnt/share/software/FAN-C/env
	fanc aggregate ${homerdir}/${filename}/${filename}.hic@50kb ${ontadout}/${filename}_use.bed ${ontadout}/${filename}_use_50kb.agg -p ${ontadout}/${filename}_use_50kb.png --tads  --tad-strength ${ontadout}/${filename}.tadstrength_50kb.txt
done
```
Identifying Loop
```
###### generate the fragment mappability file, for a fixed-size dataset
outputdir=/7.loop
fithic2dir=XX   ## path to fithic2
genomesize=XX    ## Prepare in advance
python3 ${fithic2dir}/fithic/utils/createFitHiCFragments-fixedsize.py --chrLens ${genomesize} --resolution 10000 --outFile ${outputdir}/fragments.10kb.txt.gz
zcat ${outputdir}/fragments.10kb.txt.gz | sed 's/^chr//g' - | awk -F'\t' 'BEGIN{OFS="\t"} {tmp=$2; $2=$3; $3=tmp; print $0}' - | gzip - > ${outputdir}/fragments.10kb.mod.txt.gz

###### bias
homerdir=/4.homertag
insdir=/5.transform/matrix
for filename in `ls ${homerdir} | head -n 9 | tail -n 1`;
do
	echo ${filename}
	mkdir -p ${outputdir}/${filename}
	python3 ${fithic2dir}/fithic/utils/HiCKRy.py -i ${insdir}/${filename}/${filename}_matrix_10000_combine.txt.gz -f ${outputdir}/fragments.10kb.mod.txt.gz -o ${outputdir}/${filename}/${filename}_fragments.10kb.mod.bias.0.1.txt.gz -x 0.1
	echo "run intra"
	python3 ${fithic2dir}/fithic/fithic.py -i ${insdir}/${filename}/${filename}_matrix_10000_combine.txt.gz -f ${outputdir}/fragments.10kb.mod.txt.gz -t ${outputdir}/${filename}/${filename}_fragments.10kb.mod.bias.0.1.txt.gz -r 10000 -L 30000 -v -x intraOnly -o ${outputdir}/${filename}/intra -l ${filename} -U 100000000 -p 2 > ${outputdir}/${filename}/${filename}_intra.log
	echo "run inter"
	python3 ${fithic2dir}/fithic/fithic.py -i ${insdir}/${filename}/${filename}_matrix_10000_combine.txt.gz -f ${outputdir}/fragments.10kb.mod.txt.gz -t ${outputdir}/${filename}/${filename}_fragments.10kb.mod.bias.0.1.txt.gz -r 10000 -L 30000 -v -x interOnly -o ${outputdir}/${filename}/inter -l ${filename} -U 100000000 -p 2 > ${outputdir}/${filename}/${filename}_inter.log
	echo "obtain significantly HiC loops for Intra"
	cd ${fithic2dir}/fithic/utils
	bash ${fithic2dir}/fithic/utils/merge-filter.sh ${outputdir}/${filename}/intra/${filename}.spline_pass2.res10000.significances.txt.gz 10000 ${outputdir}/${filename}/intra/${filename}.spline_pass2.res10000.significances.mergedfilter.txt.gz 0.01 > ${outputdir}/${filename}/intra/mergedfilter.log
	echo "output bedpe for visulization"
	echo "Do filteration:FDR<0.01 & minumCotact>5 & minumdist>100kb"
	zcat ${outputdir}/${filename}/intra/${filename}.spline_pass2.res10000.significances.txt.gz|awk 'BEGIN{FS=OFS="\t"}{if($4>5 && $7<0.01){if(($2>$4 && ($2-$4)>100000) || ($2<$4 && ($4-$2)>100000))print}}' > ${outputdir}/${filename}/intra/${filename}_loops_filter1
	awk 'BEGIN{FS=OFS="\t"}{print $1,$2-500,$2+500,$3,$4-500,$4+500}' ${outputdir}/${filename}/intra/${filename}_loops_filter1 > ${outputdir}/${filename}/intra/${filename}_loops_filter1.bedpe
	echo "We will remove the singlet-loops"
	zcat ${outputdir}/${filename}/intra/${filename}.spline_pass2.res10000.significances.mergedfilter.txt.gz |awk 'BEGIN{FS=OFS="\t"}{if(($9-$8)>10000 || ($11-$10)>10000){print $1,$8,$9,$1,$10,$11}}'|sort|uniq|awk 'BEGIN{FS=OFS="\t"}{print $0,"Loop_"NR}' > ${outputdir}/${filename}/intra/${filename}_temp.bedpe
	cut -f 1-3 ${outputdir}/${filename}/intra/${filename}_loops_filter1.bedpe|sort -k1,1 -k2,2n|uniq > ${outputdir}/${filename}/intra/${filename}_temp1.bed
	cut -f 4-6 ${outputdir}/${filename}/intra/${filename}_loops_filter1.bedpe|sort -k1,1 -k2,2n|uniq > ${outputdir}/${filename}/intra/${filename}_temp2.bed
	intersectBed -wo -a ${outputdir}/${filename}/intra/${filename}_temp1.bed -b ${outputdir}/${filename}/intra/${filename}_temp.bedpe > ${outputdir}/${filename}/intra/${filename}_temp1Overlap.bed
	cut -f 4- ${outputdir}/${filename}/intra/${filename}_temp.bedpe |intersectBed -wo -a ${outputdir}/${filename}/intra/${filename}_temp2.bed -b - > ${outputdir}/${filename}/intra/${filename}_temp2Overlap.bed
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if(A[$1$2$3]!=""){A[$1$2$3]=A[$1$2$3]";"$10}else{A[$1$2$3]=$10}}NR>FNR{if(A[$1$2$3]!=""){split(A[$1$2$3],a,";");for(i in a){print $0,a[i]}}}' ${outputdir}/${filename}/intra/${filename}_temp1Overlap.bed ${outputdir}/${filename}/intra/${filename}_loops_filter1.bedpe |awk 'BEGIN{FS=OFS="\t"}NR==FNR{if(A[$1$2$3]!=""){A[$1$2$3]=A[$1$2$3]";"$7}else{A[$1$2$3]=$7}}NR>FNR{if(A[$4$5$6]!=""){split(A[$4$5$6],a,";");for(i in a){print $0,a[i]}}}' ${outputdir}/${filename}/intra/${filename}_temp2Overlap.bed -|awk 'BEGIN{FS=OFS="\t"}{if($7==$8){print}}' > ${outputdir}/${filename}/intra/${filename}_loops_filter2.bedpe
	rm *temp*
	echo "Do visilization in QC_TADAndLoop.sh"
	echo "Do filteration analysis for inter"
	bash ${fithic2dir}/fithic/utils/merge-filter.sh ${outputdir}/${filename}/inter/${filename}.spline_pass1.res10000.significances.txt.gz 10000 ${outputdir}/${filename}/inter/${filename}.spline_pass1.res10000.significances.mergedfilter.txt.gz 0.01 > ${outputdir}/${filename}/inter//mergedfilter.log
done
```
QC for Loop
```
homerdir=/4.homertag
insdir=/5.transform/matrix
outputdir=/7.loop
for filename in `ls ${homerdir} | tail -n 4`;
do
	echo ${filename}
	echo "Check loops"
	pyGenomeTracks --tracks ${outputdir}/${filename}.pyGenome.ini --region chr1:39000000-46000000 --out ${outputdir}/${filename}_Chr1_39-46M_loops_tracks.pdf
	pyGenomeTracks --tracks ${outputdir}/${filename}.pyGenome.ini --region chr1:40000000-43000000 --out ${outputdir}/${filename}_Chr1_40-43M_loops_tracks.pdf
	pyGenomeTracks --tracks ${outputdir}/${filename}.pyGenome.ini --region chr20:4000000-18000000 --out ${outputdir}/${filename}_Chr20_4-18M_loops_tracks.pdf
	#
	fanc aggregate ${homerdir}/${filename}/${filename}.hic@25kb ${outputdir}/${filename}/intra/${filename}_loops_filter2.bedpe ${outputdir}/${filename}/intra/${filename}_loops_filter2.agg -p ${outputdir}/${filename}_loops_filter2.pdf --loops --pixels 30
	echo "Check CTCF and Cohensin binding"
	awk 'BEGIN{FS=OFS="\t"}{print "chr"$1,$2+500,$3+500;print "chr"$4,$5+500,$6+500}' ${outputdir}/${filename}/intra/${filename}_loops_filter2.bedpe |sort|uniq > ${outputdir}/${filename}/intra/${filename}_loop_anchor.bed
	bwoutdir=CTCF.bw  ## path to CTCF.bw
	computeMatrix reference-point -S ${bwoutdir} -R ${outputdir}/${filename}/intra/${filename}_loop_anchor.bed --beforeRegionStartLength 500000 --afterRegionStartLength 500000 --skipZeros -o ${outputdir}/${filename}_CTCF_LoopAnchor.mat.gz -p 20
	plotHeatmap -m ${outputdir}/${filename}_CTCF_LoopAnchor.mat.gz -out ${outputdir}/${filename}_CTCF_LoopAnchor.pdf --colorMap inferno --heatmapHeight 20
	echo "draw conhensin distribution"
	bwoutdir=Rad21.bw  ## path to Rad21.bw
	computeMatrix reference-point -S ${bwoutdir} -R ${outputdir}/${filename}/intra/${filename}_loop_anchor.bed --beforeRegionStartLength 500000 --afterRegionStartLength 500000 --skipZeros -o ${outputdir}/${filename}_conhensin_LoopAnchor.mat.gz -p 20
	plotHeatmap -m ${outputdir}/${filename}_conhensin_LoopAnchor.mat.gz -out ${outputdir}/${filename}_conhensin_LoopAnchor.pdf --colorMap inferno --heatmapHeight 20
	#echo "draw APA plot for inter-chromatin loops"
done
```
## Step 4: Identifying compartment
```
###### generate the PC1
homerdir=/4.homertag
for filename in XX ;
do
	runHiCpca.pl auto ${homerdir}/${filename} -res 50000 -window 50000 -genome hg38 -cpu 10
	awk 'BEGIN{FS=OFS="\t"} NR > 1 && ($2 ~ /^chr[1-9]$|^chr1[0-9]$|^chr2[0-2]$|^chrX$|^chrY$/) {if ($6 > 0) $5 = "A"; else if ($6 < 0) $5 = "B"; print $2, $3, $4, $5, $6}' ${homerdir}/${filename}/${filename}.50x50kb.PC1.txt | sort -V -k1,1 -k2,2 - > ${homerdir}/${filename}/${filename}.50x50kb.PC1_AB.txt
done
```





