# Aligning to only reference FASTA file issue

The reason why `moni-align` currently cannot do alignments with just the FASTA file is tied to its coupling with `LevioSAM` to liftover alignments to the reference genome. The `moni-build` function will not generate the `.ldx` or `.lidx` files if a VCF file is not provided. These are `LevioSAM` files and it requires a VCF to build these files. It is a non-trivial task to decouple the use of `LevioSAM` from `moni-align`. I believe that I have found a workaround to the issue although its far from unideal.

## Workaround Idea #1 (Better idea)
Create a dummy VCF file and build the index with the reference and dummy VCF file. Then do the alignment and afterwards filter the alignments with samtools. For full alignments, you should not have to do any filtering afterwards but it would be necessary if only desiring MEMs.

1. Create a dummy 1-line VCF file (can use the script provided in the workaround directory)

```
./dummy_vcf.py -o ../data/mouse/vcf/dummy
```

2. Run `moni-build` with the FASTA file and dummy VCF file.

```
./moni build -r ../data/mouse/ref/mouse.chr19.fa.gz -v ../data/mouse/vcf/dummy.vcf.gz -o ../data/mouse/index/mouse
```

3. Run the `moni-align` command

```
./moni align -i ../data/mouse/index/mouse -1 ../data/mouse/reads/mouse.chr19.R1.fastq -2 ../data/mouse/reads/mouse.chr19.R2.fastq -o ../data/mouse/output/mouse.sam
```

4. Optional: Filter the alignments with the samtools view command if necessary.

## Workaround Idea #2
We can build the index with just the reference FASTA file and then afterwards create a dummy VCF file to create dummy `LevioSAM` files to satisfy `moni-align` and not cause it to error. In theory, these files should not actually do anything in `moni-align` since the index was only built on the reference genome and therefore the alignments should only be to the reference.

## Workaround Steps
0. Move to the build directory.

1. Run `moni-build` with only the FASTA file.

```
./moni build -f ../data/mouse/ref/mouse.chr19.fa.gz -o ../data/mouse/index/mouse 
```
This should produce the following files in the `data/mouse/index` directory: `mouse.moni.log`, `mouse.plain.slp`, `mouse.slcp`, and `mouse.thrbv.full.lcp.ms`

2. Create a dummy 1-line VCF file (can use the script provided in the workaround directory)

```
./dummy_vcf.py -o ../data/mouse/vcf/dummy
```
This will create a dummy VCF file called `dummy.vcf.gz` in the `data/mouse/vcf` directory.

3. Download and run [mvtool](https://github.com/marco-oliva/mvtool)   (Thanks Marco) (assuming downloaded in the workaround directory and using singularity version)

```
./mvtool_sif mvtool -r ../data/mouse/ref/mouse.chr19.fa.gz -v ../data/mouse/vcf/dummy.vcf.gz -w 10 -o ../data/mouse/index/mouse
```
This should produce the following files files in the `data/mouse/index` directory: `mouse.lengths` and `mouse.lifting`

4. Move to the index directory and rename the `mouse.lengths` and `mouse.lifting` files

```
mv mouse.lengths mouse.lidx
mv mouse.lifting mouse.ldx
```

This will rename the files as `mouse.lidx` and `mouse.ldx`. 

5. Move back to the build directory and run the `moni-align` command

```
./moni align -i ../data/mouse/index/mouse -1 ../data/mouse/reads/mouse.chr19.R1.fastq -2 ../data/mouse/reads/mouse.chr19.R2.fastq -o ../data/mouse/output/mouse.sam
```

All the original alignments should be to the reference genome.