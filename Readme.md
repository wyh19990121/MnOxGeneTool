
MnOxGeneTool is a bioinformatics tool comprising two modules (MnOxGeneTool-hmm and MnOxGeneTool-abundance) based on a Mn(II)-oxidizing protein database containing homologous proteins of 15 non-redundant reported Mn(II)-oxidizing genes.
MnOxGeneTool-hmm supports the input form of long DNA sequences (complete genomic sequences, contigs) or complete protein sequences for identifying gene sequences encoding Mn(II)-oxidizing proteins or Mn(II)-oxidizing proteins.
MnOxGeneTool-abundance supports the input form of single-end or either end of paired-end metagenomic short reads data for calculating the abundance of Mn(II)-oxidizing genes in metagenomic datasets.

MnOxGeneTool is a product of the paper titled `"MnOxGeneTool: A Comprehensive Tool for Identifying and Quantifying Mn(II)-oxidizing Genes, Revealing Phylogenetic Diversity and Environmental Drivers of Mn(II)-oxidizers"` (currently under review). You can find more detailed information of the tool in the paper.

If you have any questions, please create an issue, or contact wyh (1365298466@qq.com).





### **Prerequisites**
A basic Linux shell environment is required to run the software. Additionally, users need to install the following software in advance:
- `Diamond>=2.0.15` [Diamond](https://github.com/bbuchfink/diamond)
- `Blast>=2.15` [Blast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- `Samtools>=1.15` [Samtools](https://github.com/samtools/samtools)
- `Prodigal>=2.6.3` [Prodigal](https://github.com/hyattpd/Prodigal)
- `Hmmer>=3.4` [Prodigal](https://github.com/EddyRivasLab/hmmer)
- `Bwa-mem2` [Bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
- `Seqkit`[Seqkit](https://github.com/shenwei356/seqkit)
- `FragGeneScanRs` [FragGeneScanRs](https://github.com/unipept/FragGeneScanRs)
- `Uproc` [Uproc](https://github.com/gobics/uproc)

Once the required software is installed, you can clone the repository using:
```bash
git clone https://github.com/wyh19990121/MnOxGeneTool.git
```
Then, navigate into the directory and run the following command:
```bash
bash initial.sh
```
This command will create launch files under `/usr/local/bin`. Once the process completes, you can execute the corresponding functions anywhere in the system using `MnOxGeneTool-hmm` or `MnOxGeneTool-abundance`. 

If you do not wish to create launch files under `/usr/local/bin`, you can use the `-i` option to store the launch files in any preferred folder. Afterward, you just need to add the path of that folder to your environment variable:
```bash
bash initial.sh -i /Path/You/Prefer
```

---

### **Usage**

#### **MnOxGeneTool-hmm**
- `-i`: Input file in FASTA format
- `-o`: Output directory where the results will be saved
- `-f`: Input type, either `fa` (for genomic/contig sequences) or `faa` (for protein sequences). If the input is a genome or contig sequence, choose `fa`; if it is protein sequence, choose `faa`.
- `-r`: If you do not want to keep intermediate files, select `-r` to delete them after execution. The process files are saved in the `tmp` folder by default.
- `-t`: Number of threads to use

**Example Usage:**
```bash
MnOxGeneTool-hmm -i ./testdata/sequence.fa -o ./testdata -f fa -t 180
```
In this example, the input file contains the full genomes of two microorganisms, `NC_010322.1` and `NZ_CP019304.1`. As it is genome data, the `-f` option is set to `fa`. The output is specified to be saved in the `testdata` folder, with the output file named automatically based on the input file name (e.g., `sequence_HMM_result`). The process will use 180 threads and retain intermediate files.

The output file `sequence_HMM_result` will look like this:
```
NC_010322.1 katG 496.7
NC_010322.1 mcoA 626.2
NC_010322.1 mnxG_P 1031.3
NC_010322.1 mopA_P 1281.2
NZ_CP019304.1 boxA 437.0
NZ_CP019304.1 boxA 441.8
NZ_CP019304.1 boxA 446.7
NZ_CP019304.1 katG 481.0
```
Here, the first column represents the sequence information from the input FASTA file, the second column shows the identified gene type, and the third column is the corresponding score. In the genome `NZ_CP019304.1`, `boxA` appears three times, indicating three copies of the gene in the genome.

---

#### **MnOxGeneTool-abundance**
- `-i`: Input file in FASTA format
- `-o`: Output directory where the results will be saved
- `-s`: Sample size for extraction. By default, no extraction is performed. If you want to extract a specific number of sequences for analysis, use `-s`. For example, `-s 5000000` will extract 5 million reads from the original sample for analysis. This can help reduce processing time, but it is not recommended to set it below 5 million.
- `-a`: Select either `-a 16s` to calculate abundance based on 16S or `-a cell` to calculate abundance based on cell counts. By default, both will be calculated, and you can exclude one to reduce unnecessary calculations.
- `-r`: If you do not want to keep intermediate files, select `-r` to delete them after execution. The process files are saved in the `tmp` folder by default.
- `-t`: Number of threads to use (default is 1)

**Example Usage:**
```bash
MnOxGeneTool-abundance -i ./SRR5739200_200w.fa -o ./testdata -t 180
```
In this example, the input file `SRR5739200_200w.fa` contains 2 million reads from environmental metagenomic data. The output will be saved in the `testdata` folder. As the `-a` option is not specified, both abundance calculations for 16S and cell count will be performed, and two files (`SRR5739200_200w_abundance_16s` and `SRR5739200_200w_abundance_cell`) will be generated.

The `SRR5739200_200w_abundance` output file will look like this:
```
boxA    0.00229495
katG    0.274935
mcoA    0.0147713
mnxG_P  0.0163207
mokA    0.00842726
mopA_A  0.000297223
mopA_P  0.0146553
mopA_R  0.000986554
moxA    0.0302352
```
In this file, the first column shows the Mn(II)-oxidizing gene types detected in the input file, and the second column represents their corresponding abundance.

---

Feel free to reach out if you encounter any issues or need further assistance!
