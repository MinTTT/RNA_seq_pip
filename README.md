# RNA_seq_pip

These codes are written for analyzing the next generation data, not only for mRNA seq data.

We have verified the data types for analyzing: 
* RNA seq
* Deep seq

### RNA seq
the script `RNA_seq_pip_CP` was written for mRNA-seq analysis





### RNA seq



### trimmed adapter



```shell
cutadapt -q 20,20 --minimum-length 100:100 --max-n 3 --pair-filter=any -o "${sample}"_QF_R1.fastq.gz -p "${sample}"_QF_R2.fastq.gz "${sample}"_R1.fastq.gz "${sample}"_R2.fastq.gz > "${sample}"_logs/QF.log

# Demultiplexing /////////////////////////////////////////////////////////////
# when the library quality is high, say, each adaptor can be found strictly anchored at the 5' of R1
cutadapt --pair-filter=any --no-indels --minimum-length 20 --times 1 -e 0.15 --overlap 7 -g file:adaptor_R1.fasta -A file:adaptor_R2.fasta -o end5-{name}.R1.fastq.gz -p end5-{name}.R2.fastq.gz "${sample}"_QF_R1.fastq.gz "${sample}"_QF_R2.fastq.gz > "${sample}"_logs/demultiplexing.log

```



#### alignment

* bowtie 2:
  * `-X`: The maximum fragment length for valid paired-end alignments. set to `1000`
  * `-I `: The minimum fragment length for valid paired-end alignments. set to `18`
  * `--no-mixed`: 
  * `--no-discordant`:  
  * 
