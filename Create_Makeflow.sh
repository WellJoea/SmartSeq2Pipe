INpath=$1
SID=$2
Module=$3

Lanemkflow()
{
    INp=$1
    OUp=$2
    cat $INp >> $OUp
    cat>>$OUp<<EOF

CATEGORY="fastqclink.${AID}"
CORES=7
MEMORY=20000
DISK=10
${RawFQpre}_r1.fq.gz ${RawFQpre}_r2.fq.gz : ${R1} ${R2}
    sh $WORKFLOW_DIR/FastqLink $INp

CATEGORY="cutadapt.trim.${AID}"
CORES=12
MEMORY=20000
DISK=10
${CUTADAPTpre}_trimmed_1.fq.gz  ${CUTADAPTpre}_trimmed_2.fq.gz : ${RawFQpre}_r1.fq.gz  ${RawFQpre}_r2.fq.gz
    sh $WORKFLOW_DIR/CUTAdapt $INp

CATEGORY="fastqc.${AID}"
CORES=7
MEMORY=20000
DISK=10
${FASTQCpre}_r1_fastqc.html ${FASTQCpre}_trimmed_2_fastqc.html : ${RawFQpre}_r1.fq.gz  ${RawFQpre}_r2.fq.gz ${CUTADAPTpre}_trimmed_1.fq.gz ${CUTADAPTpre}_trimmed_2.fq.gz
    sh $WORKFLOW_DIR/FastQC $INp

CATEGORY="trim_galore.${AID}"
CORES=7
MEMORY=20000
DISK=10
${TRIMGALpre}_val_1.fq.gz  ${TRIMGALpre}_val_2.fq.gz : ${RawFQpre}_r1.fq.gz  ${RawFQpre}_r2.fq.gz
    sh $WORKFLOW_DIR/TRIMGalore $INp 

CATEGORY="star.${AID}"
CORES=20
MEMORY=30000
DISK=10
${STARpre}.Aligned.toTranscriptome.out.bam ${STARpre}.Aligned.sortedByCoord.out.bam : ${CUTADAPTpre}_trimmed_1.fq.gz  ${CUTADAPTpre}_trimmed_2.fq.gz
    sh $WORKFLOW_DIR/STAr $INp  ${CUTADAPTpre}_trimmed_1.fq.gz ${CUTADAPTpre}_trimmed_2.fq.gz

CATEGORY="hisat2.${AID}"
CORES=20
MEMORY=30000
DISK=10
${HISAT2pre}.hisat2.bam : ${CUTADAPTpre}_trimmed_1.fq.gz  ${CUTADAPTpre}_trimmed_2.fq.gz
    sh $WORKFLOW_DIR/HISAT2  $INp  ${CUTADAPTpre}_trimmed_1.fq.gz ${CUTADAPTpre}_trimmed_2.fq.gz




EOF
}

Repmkflow()
{
    INp=$1
    OUp=$2
    cat $INp >> $OUp
    cat>>$OUp<<EOF


CATEGORY="bwa.trans.star.${TID}"
CORES=10
MEMORY=20000
DISK=10
${BAMTranSM} ${BAMTranSMN} ${BAMCoordSM} ${BAMCoordSMN} : ${BAMTran[@]}  ${BAMCoord[@]} 
    sh $WORKFLOW_DIR/BAMTranStar $INp

CATEGORY="bwa.trans.hisat.${TID}"
CORES=10
MEMORY=20000
DISK=10
${BAMHisSM} ${BAMHisSMN} : ${BAMHis[@]}
    sh $WORKFLOW_DIR/BAMTranHisat $INp

CATEGORY="rsemce.${TID}"
CORES=10
MEMORY=20000
DISK=10
${RSEMcepre}.rsemce.isoforms.results ${RSEMcepre}.rsemce.genes.results : $BAMTranSMN
    sh $WORKFLOW_DIR/RSEMce $INp  $BAMTranSMN

CATEGORY="featureCounts.${TID}"
CORES=10
MEMORY=20000
DISK=10
${FTcountspre}.fCounts.txt : $BAMCoordSM
    sh $WORKFLOW_DIR/FeatureCounts $INp $BAMCoordSM

CATEGORY="asgal.splicing.${TID}"
CORES=10
MEMORY=20000
DISK=10
$ASGALpre/ASGAL.csv : ${FQEachRep[@]}
    sh $WORKFLOW_DIR/ASGAL $INp

CATEGORY="spladder.splicing.${TID}"
CORES=10
MEMORY=20000
DISK=10
$SPLADDERout/$(basename ${BAMCoordSM/%.bam/})_exon_skip_C3.txt.gz : ${BAMCoordSM}
    sh $WORKFLOW_DIR/SplAdder $INp ${BAMCoordSM}

CATEGORY="cufflinks.${TID}"
CORES=10
MEMORY=20000
DISK=10
${CUFFLinkspre}.combined.gtf : ${BAMCoordSM}
    sh $WORKFLOW_DIR/CuffLinks $INp ${BAMCoordSM}



EOF
}


INSf=( $INpath/${SID}__D*${Module}.input )
. ${INSf[0]}
WorkSH=$AUT/WorkShell
SampMK=$WorkSH/${SID}_${Module}.mf

if [ -f $SampMK ];then
    rm $SampMK
fi
mkdir -p $WorkSH

#########each lane
for i in ${INSf[@]};do
    . $i
    . $WORKFLOW_DIR/env.Variable
    cp  $i $WorkSH/${AID}_${Module}.input
    Lanemkflow $WorkSH/${AID}_${Module}.input $SampMK
done

#########each rep
for j in ${TIDs[@]};do
    Tidf=( $WorkSH/${j}*${Module}.input )
    . ${Tidf[0]}
    . $WORKFLOW_DIR/env.Variable
    Repmkflow ${Tidf[0]} $SampMK
done

