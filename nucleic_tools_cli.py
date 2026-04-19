#!/usr/bin/env python3
'''
Command-line interface for NucleicTools.

Usage:
    python cli.py filter --input reads.fastq --output filtered.fastq
    python cli.py convert --input sequences.fasta

To make this script executable:
    chmod +x cli.py
    ./cli.py filter --help
'''

import argparse
from main import filter_fastq
from bio_files_processor import convert_multiline_fasta_to_oneline


def parse_args():
    '''Parse command-line arguments and return parsed args.'''
    parser = argparse.ArgumentParser(description="NucleicTools CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # filter
    filter_parser = subparsers.add_parser("filter", help="Filter FASTQ sequences")
    filter_parser.add_argument("--input", required=True, help="Input FASTQ file")
    filter_parser.add_argument("--output", required=True, help="Output FASTQ file name")
    filter_parser.add_argument("--gc-bounds", nargs=2, type=float, default=[0, 100],
                               metavar=("MIN", "MAX"), help="GC-content bounds in percent (default: 0 100)")
    filter_parser.add_argument("--length-bounds", nargs=2, type=int, default=[0, 2**32],
                               metavar=("MIN", "MAX"), help="Sequence length bounds (default: 0 4294967296)")
    filter_parser.add_argument("--quality-threshold", type=int, default=0,
                               help="Minimum average quality in Phred33 (default: 0)")
    filter_parser.add_argument("--overwrite", action="store_true",
                               help="Overwrite output file if it exists")

    # convert
    convert_parser = subparsers.add_parser("convert", help="Convert multiline FASTA to oneline")
    convert_parser.add_argument("--input", required=True, help="Input FASTA file")

    return parser.parse_args()


def main():
    '''Run the selected subcommand.'''
    args = parse_args()

    if args.command == "filter":
        filter_fastq(
            input_fastq=args.input,
            output_fastq=args.output,
            gc_bounds=tuple(args.gc_bounds),
            length_bounds=tuple(args.length_bounds),
            quality_threshold=args.quality_threshold,
            overwrite=args.overwrite,
        )

    elif args.command == "convert":
        convert_multiline_fasta_to_oneline(args.input)


if __name__ == "__main__":
    main()
