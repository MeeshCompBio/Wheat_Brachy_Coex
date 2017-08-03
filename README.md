# Wheat_Brachy_Coex
Run the following command to go from fastq files to generating FPKM values from uniquely mapped reads. 

```
bash Run_Pipeline.sh
```
## Please note
* There are hard coded file paths for GNU parallel usage to fastq files and outputs.
* Multi-threading was also used in multiple scripts. Please check that your machine has the available resources.
* Wheat genome needs 164GB of RAM to run with STAR(regardless of threads).
* Please modify the scripts accordingly.

## Software versions
```
GNU parallel 20160822
FastQC 0.11.5
Cutadapt 1.8.1
STAR 2.5.3a
Samtools 1.2 (using htslib 1.2.1)
```

