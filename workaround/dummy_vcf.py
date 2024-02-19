#!/usr/bin/python3
import argparse
import os
import gzip

Description = """
Generate a dummy one-line VCF file 

   by Rahul Varki
"""

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--contig', help='Contig name to use in the VCF file (default. chr19)', type=str, dest="contig", default = "chr19")
    parser.add_argument('-p', '--position', help='Position of the variant in the VCF file (default. 10)', type=int, dest="pos", default = 10)
    parser.add_argument('-r', '--ref', help='Reference base to use in the VCF file (default. A)', type=str, dest="ref", default = "A")
    parser.add_argument('-a', '--alt', help='Alternate base to use in the VCF file (default. T)', type=str, dest="alt", default = "T")
    parser.add_argument('-s', '--sample', help='Name of the dummy sample (default. Dummy)', type=str, dest="sample", default = "Dummy")
    parser.add_argument('-o', '--output', help='Output prefix (without the .vcf.gz extension)', type=str, dest="output", required = True)  
    args = parser.parse_args()

    outdir = os.path.dirname(args.output)
    if not os.path.exists(outdir):
        print("Error: The output directory {} does not exist".format(outdir))
        return

    filename = args.output + ".vcf.gz"

    with gzip.open(filename, 'wb') as vcf_file:
        # Write the VCF header
        vcf_file.write('##fileformat=VCFv4.3\n'.encode())
        vcf_file.write('##contig=<ID={contig}>\n'.format(contig = args.contig).encode())
        vcf_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'.encode())
        vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format(args.sample).encode())

        # Write a variant record
        vcf_file.write('{contig}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t0/1\n'.format(contig=args.contig, pos=args.pos, ref=args.ref, alt=args.alt).encode())


if __name__ == '__main__':
    main()