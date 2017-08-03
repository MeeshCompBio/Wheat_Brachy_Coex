#Please note that there are hard coded file paths for GNU parallel usage
#Multi-threading was also used in multiple scripts. Please check that your
####machine has the available resourses
#Please modify the scripts accordingly
#Wheat genome needs 164GB of RAM to run with STAR(regardless of threads)

#Index both genomes using STAR
bash STARIndex_Wheat_Brachy.sh

#Run cutadapt to check for Illumina adapters and trim the end of reads with low quality
parallel --jobs 8 < CutCommands_Brachy.txt
parallel --jobs 8 < CutCommands_Wheat.txt

#Align reads to the genome using STAR
parallel --jobs 3 < STAR_Commands_Brachy.sh
#only run one at a time for wheat due to resources
bash STARCommands_Wheat.sh

#STAR outputs SAM files so we will convert to bam and filter 
####for unique reads
parallel --jobs 8 < SamCommands_Brachy.txt
parallel --jobs 8 < SamCommands_Wheat.txt

#Generate FPKM using cufflinks
parallel --jobs 2 < FPKMCommands_Brachy.txt
parallel --jobs 2 < FPKMCommands_Wheat.txt
