# Wheat_Brachy_Coex
Run the following command to go from fastq files to generating FPKM values from uniquely mapped reads. 

```
bash Run_Pipeline.sh
```
## Please note
* There are hard coded file paths for GNU parallel usage to fastqfiles and outputs.
* Multi-threading was also used in multiple scripts. Please check that your machine has the available resources.
* Wheat genome needs 164GB of RAM to run with STAR(regardless of threads).
* Please modify the scripts accordingly.

