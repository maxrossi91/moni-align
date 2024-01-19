[![Release](https://img.shields.io/github/release/maxrossi91/moni.svg)](https://github.com/maxrossi91/moni/releases)
[![Downloads](https://img.shields.io/github/downloads/maxrossi91/moni/total?logo=github)](https://github.com/maxrossi91/moni/archive/master.zip)

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

There are two ways to be able to run Moni-align. The easiest way is to run the tool from the Docker or Singularity/Apptainer image. The other way is to build it from source (Linux only).  

## Docker

1. Pull the image from DockerHub
```
docker pull rvarki/moni-align
```
2. Run the moni help command. 
```
docker run rvarki/moni-align moni -h
```

> [!NOTE] 
> Make sure to include the -v option in the docker run command before attempting to run any moni commands on your data. The -v option will allow you to [mount](https://docs.docker.com/storage/bind-mounts/) your host directory into the Docker container so that moni can interact with your data. Users are encouraged to mount into `\mnt` directory in the Docker containter. 

## Singularity/Apptainer

If working on a high performance computing cluster, you might not have access to Docker due to it needing root permissions. Instead, you can use the Singularity/Apptainer version of the tool.

1. Build the SIF image of Moni-align
```
singularity pull moni.sif docker://rvarki/moni-align:latest
```
or
```
apptainer build moni.sif docker://rvarki/moni-align:latest
```
2. Run the moni help command

```
./moni.sif moni -h
```

## Building from Source

1. Currently, building from source only works on Linux machines. 

```
git clone https://github.com/maxrossi91/moni-align.git
cd moni-align
mkdir build
cd build
cmake ..
make
cmake ..
make
```

> [!WARNING]
> The first make command will fail. Continue with the rest of the steps. 

2. Run the moni help command
```
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$PWD/thirdparty/lib"
./moni -h
```

> [!IMPORTANT]
> Make to run the export LD_LIBRARY_PATH command prior to running any moni commands if building the project from source.

# Example

We provide test data located in the `data/mouse` directory to test the commands.




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