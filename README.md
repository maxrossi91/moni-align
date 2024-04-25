<!-- [![Release](https://img.shields.io/github/release/maxrossi91/moni.svg)](https://github.com/maxrossi91/moni/releases)
[![Downloads](https://img.shields.io/github/downloads/maxrossi91/moni/total?logo=github)](https://github.com/maxrossi91/moni/archive/master.zip) -->

# MONI-Align
```console
___  ___            _         ___  _ _             
|  \/  |           (_)       / _ \| (_)            
| .  . | ___  _ __  _ ______/ /_\ \ |_  __ _ _ __  
| |\/| |/ _ \| '_ \| |______|  _  | | |/ _` | '_ \ 
| |  | | (_) | | | | |      | | | | | | (_| | | | |
\_|  |_/\___/|_| |_|_|      \_| |_/_|_|\__, |_| |_|
                                        __/ |      
                                       |___/       
                                            ver 0.1.0
```
A Read Aligner with Multi-Genome References.

MONI index uses the prefix-free parsing of the text [2][3] to build the Burrows-Wheeler Transform (BWT) of the reference genomes, the suffix array (SA) samples at the beginning and at the end of each run of the BWT, and the threshold positions of [1]. The MONI index can be built from (a) one or more FASTA file(s) or (b) a reference FASTA and one or more VCF files.  

Maximal Exact Matches (MEMs) are extracted between the reads and the multi-genome reference using the MONI index. The reads are aligned to the multi-genome reference by using the MEMs as seeds to guide read alignment, a similar strategy employed by other popular read aligners [4]. The alignments are reported in reference coordinates by lifting over the multi-genome alignments to the reference genome [5]. A SAM file is produced that can be used in downstream analysis software such as variant callers. 

# Install

There are two ways to be able to run `moni-align`. The easiest way is to run the tool from the `Docker` or `Singularity/Apptainer` image. The other way is to build it from source (Linux only).  

## Docker

1. Pull the image from DockerHub
```
docker pull rvarki/moni-align:davide
```
2. Run the `moni` help command. 
```
docker run rvarki/moni-align moni -h
```

> [!NOTE] 
> Make sure to include the -v option in the docker run command before attempting to run any `moni` commands on your data. The -v option will allow you to [mount](https://docs.docker.com/storage/bind-mounts/) your host directory into the Docker container so that `moni` can interact with your data. Users are encouraged to mount into `\mnt` directory in the Docker containter. 

## Singularity/Apptainer

If working on a high performance computing cluster, you might not have access to `Docker` due to it needing root permissions. Instead, you can use the `Singularity/Apptainer` version of the tool.

1. Build the SIF image of `moni-align`
```
singularity pull moni.sif docker://rvarki/moni-align:davide
```
or
```
apptainer build moni.sif docker://rvarki/moni-align:davide
```
2. Run the `moni` help command

```
./moni.sif moni -h
```

## Building from Source
Currently, building from source only works on Linux machines. 

`moni-align` has the following dependencies:

- `gcc` tool (version 9.3.0/9.4.0)
- `cmake`tool (version 3.15 or greater) 

> [!WARNING]
> To build from source, it is absolutely necessary to have one of those versions of gcc installed.

After installing the dependencies, the following steps can be run to build the project.

1. Run the build steps 

```
git clone https://github.com/maxrossi91/moni-align.git
cd moni-align
mkdir build
cd build
cmake ..
cmake ..
make
```

> [!WARNING]
> Make sure to run the cmake command twice as specified.

2. Set the LD_LIBRARY_PATH environment variable
```
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$PWD/thirdparty/lib"
```

> [!IMPORTANT]
> Make sure to run the export LD_LIBRARY_PATH command prior to running any `moni` commands if building the project from source.

3. Run the `moni` help command
```
./moni -h
```

# Example

We provide test data located in the `data/mouse` directory to test the commands. We use this test data to showcase how to perform some common operations.

## Building the Moni index

Prior to running the alignment step, the `moni` index of the references must first be created. The references can stored either in concatenated FASTA format or reference FASTA and VCF file. The only major difference between building with either option is that the alignment liftover step is currently only possible when the index is built with a reference FASTA and VCF file. If the index with a concatenated FASTA, the downstream alignments will be reported to the haplotype it aligned to instead. Here we show building the index with both type of input.

1. Building the index with reference FASTA file + VCF file.
```
./moni build -r ../data/mouse/ref/mouse.chr19.fa.gz -v ../data/mouse/vcf/mouse.chr19.subset.vcf.gz -S ../data/mouse/vcf/mouse_samples.txt -H12 -o ../data/mouse/index/mouse 
```

Here we provide a brief summary of the options used in the above command.

- -r : The location of the reference FASTA file. The reference FASTA is expected to have the .gz extension.

- -v : The location of the corresponding VCF file. The file should be structured such that the reference alleles at all positions match the reference FASTA file.

- -S : A file containing the names of the samples in the VCF to build the index from. This is optional to provide, however `moni` assumes to use all samples in the VCF if not provided.

- -H : The haplotype information from the VCF file to include. Write -H1 to include only the information from haplotype 1, -H2 to include only the information from haplotype 2, or -H12 to include information from both haplotypes.

- -o : The location of the output directory and the prefix of the output files. For the example, the files will be written to this location: `../data/mouse/index/` with the files having the `mouse` prefix. 


2. Building the index with only a FASTA file.

```
./moni build -f ../data/mouse/ref/mouse.chr19.fa.gz -o ../data/mouse/index/mouse
```

Here we provide a brief summary of the options used in the above command.

- -f : The location of the reference FASTA file. The reference FASTA is expected to have the .gz extension.

- -o : The location of the output directory and the prefix of the output files. For example, in this example, the files will be written to this location: `../data/mouse/index/` with the files having the `mouse` prefix.

Regardless of which method used to build the index, the following files should be produced in the `data/mouse/index` directory: `mouse.ldx`, `mouse.lidx`, `mouse.moni.log`, `mouse.plain.slp`, `mouse.slcp`, and `mouse.thrbv.full.lcp.ms`

> [!NOTE]
> To see the other build options available, type ./moni build -h 

## Aligning the reads

After building the `moni` index, reads can be aligned to the pangenome. The align function accepts both single and paired-end reads in either FASTA or FASTQ format. Here we show how to run the alignment step with both paired-end and single-end reads.

1. Aligning with paired-end reads.
```
./moni align -i ../data/mouse/index/mouse -1 ../data/mouse/reads/mouse.chr19.R1.fastq -2 ../data/mouse/reads/mouse.chr19.R2.fastq -o ../data/mouse/output/mouse.sam
```

Here we provide a brief summary of the options used in the above command.

- -i : The path to where the previously created index files are located and the prefix of these files. For the example, the index files are located in this location: `../data/mouse/index/` with the `mouse` prefix.

- -1 : The location of the mate 1 file.
- -2 : The location of the mate 2 file.

- -o : The location of where to write the output SAM file and the name of the SAM file. For the example, the SAM file will be written to  this location: `../data/mouse/index/` with the name: `mouse.sam`.

2. Aligning with single-end reads.

```
./moni align -i ../data/mouse/index/mouse -p ../data/mouse/reads/mouse.chr19.R1.fastq -o ../data/mouse/output/mouse.sam
```
Here we provide a brief summary of the options used in the above command.

- -i : The path to where the previously created index files are located and the prefix of these files. For the example, the index files are located in this location: `../data/mouse/index/` with the `mouse` prefix.

- -p : The location of the single-end read file.

- -o : The location of where to write the output SAM file and the name of the SAM file. For the example, the SAM file will be written to  this location: `../data/mouse/index/` with the name: `mouse.sam`.


Either command should produce the `mouse.sam` file in the `data/mouse/output` directory.

> [!NOTE]
> To see the other build options available, type ./moni align -h

> [!IMPORTANT]
> By default, there are two types of filters applied: A MEM occurance filter (-S) and MEM orientation (-D) filter. These filters were designed to speed up the alignment step, however both can potentially lead to a loss of sensitivity or specificity in the read alignments. To disable the MEM occurance filter, add the -s option. To disable the MEM orientation filter, add the -d option.

## Writing MEMs to SAM

If only the MEM occurances in the pangenome are desired, we provide an option to write the MEMs of the read to the SAM file instead of the full read alignments. This can be done with both single and paired-end reads in either FASTA or FASTQ format. Here we show how to do it with paired-end reads, but the command can easily be modified for single-end reads.  

1. Writing the MEMs of length 25 or greater to the SAM file for all reads. 
```
./moni align -i ../data/mouse/index/mouse -1 ../data/mouse/reads/mouse.chr19.R1.fastq -2 ../data/mouse/reads/mouse.chr19.R2.fastq -o ../data/mouse/output/mouse.sam --report_mems -d -s -l 25
```

Here we provide a brief summary of the options used in the above command.

- -i : The path to where the previously created index files are located and the prefix of these files. For the example, the index files are located in this location: `../data/mouse/index/` with the `mouse` prefix.

- -1 : The location of the mate 1 file.
- -2 : The location of the mate 2 file.

- -o : The location of where to write the output SAM file and the name of the SAM file. For the example, the SAM file will be written to  this location: `../data/mouse/index/` with the name: `mouse.sam`.


- --report_mems: This option will cause `moni` write the MEMs of mate1 and mate2 found in the forward and reverse direction to the SAM file rather than the read alignments. 

- -s : Disable the MEM occurance filter.

- -d : Disable the MEM orientation filter.

- -l : Set the MEM length cutoff filter. 

Running this command will write all MEMs of length 25 or greater to the SAM file for all reads. This will likely cause the SAM file created to be much larger in size than if reporting the read alignments. This is because each read can have multiple entries in the SAM file depending on how many MEMs are found in the pangenome. The same MEM from the same read can have multiple entires if it is found in multiple of the genomes. 


# External resources

* [Big-BWT](https://github.com/alshai/Big-BWT.git)
    * [gSACA-K](https://github.com/felipelouza/gsa-is.git)
    * [malloc_count](https://github.com/bingmann/malloc_count)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)
    * [Divsufsort](https://github.com/simongog/libdivsufsort.git)
* [klib](https://github.com/attractivechaos/klib)
* [r-index](https://github.com/maxrossi91/r-index.git)
* [pfp-thresholds](https://github.com/maxrossi91/pfp-thresholds.git)
* [bigrepair](https://gitlab.com/manzai/bigrepair.git)
* [shaped_slp](https://github.com/maxrossi91/ShapedSlp.git)
* [htslib](https://github.com/samtools/htslib)
* [levioSAM](https://github.com/alshai/levioSAM)
<!-- * [Google Benchmark](https://github.com/google/benchmark.git)
    * [Google Test](https://github.com/google/googletest) -->

<!-- # Citation 

Please, if you use this tool in an academic setting cite the following paper:

    @article{BoucherCGHMNR20,
    author    = {Christina Boucher and
                Ondřej Cvacho and
                Travis Gagie and
                Jan Holub and
                Giovanni Manzini and
                Gonzalo Navarro and
                Massimiliano Rossi},
    title     = {PFP Data Structures},
    journal   = {CoRR},
    volume    = {abs/xxxx.xxxxx},
    year      = {2020},
    url       = {https://arxiv.org/abs/xxxx.xxxxx},
    archivePrefix = {arXiv},
    eprint    = {xxxx.xxxxx},
    } -->


# Authors

### Theoretical results:

* Christina Boucher
* Ben Langmead
* [Massimiliano Rossi](https://github.com/maxrossi91)

### Implementation:

* [Massimiliano Rossi](https://github.com/maxrossi91)

### Experiments

* [Rahul Varki](https://github.com/rvarki)
* [Eddie Ferro](https://github.com/EddieFerro)
* [Marco Oliva](https://github.com/marco-oliva)


# Why "MONI"?

**Moni** is the Finnish word for *multi*.

# References

[1] Hideo Bannai, Travis Gagie, and Tomohiro I, *"Refining ther-index"*, Theoretical Computer Science, 812 (2020), pp. 96–108

[2] Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini, *"Prefix-Free Parsing for Building Big BWTs"*, In Proc. of the 18th International Workshop on Algorithms in Bioinformatics (WABI 2018).

[3] Christina Boucher, Travis Gagie, Alan Kuhnle, Ben Langmead, Giovanni Manzini, and Taher Mun. *"Prefix-free parsing for building big BWTs."*, Algorithms for Molecular Biology 14, no. 1 (2019): 13.

[4] Heng Li. *"Minimap2: pairwise alignment for nucleotide sequences."* Bioinformatics, 34(18):3094–3100,
2018.

[5] Taher Mun, Nae-Chyun Chen, and Ben Langmead. *"LevioSAM: fast lift-over of variant-aware reference alignments."* Bioinformatics, 37(22):4243–4245, 2021.
