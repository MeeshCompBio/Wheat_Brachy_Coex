#BrachiGenome STAR index
STAR \
    --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir BdistachyonBd21-3v1.1_STAR \
    --genomeFastaFiles BdistachyonBd21-3v1.1.fasta \
    --sjdbGTFfile BdistachyonBd21-3v1.1.gene.gff3 \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbGTFfeatureExon CDS \
    --sjdbOverhang 125

#Same for Wheat
STAR \
    --runThreadN 24 \
    --runMode genomeGenerate \
    --genomeDir 161010_Chinese_Spring_v1.0_STAR \
    --genomeFastaFiles 161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta \
    --sjdbGTFfile iwgsc_refseqv1.0_HighConf_2017Mar13-transposed.gff3 \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbOverhang 125

