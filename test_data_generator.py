#!/usr/bin/env python3

# Copyright (c) 2018 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import argparse
import functools
import subprocess
import sys
from pathlib import Path


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--accession', type=str, help="The SRA accession number to download")
    parser.add_argument(
        '-n', '--number', type=int, default=10000, help="number of reads per output file", )
    parser.add_argument("-N", "--nof", default=1, type=int, help="number of output files")
    parser.add_argument("--bti", "--bowtie-index", type=str, help="the bowtie index to align to.")
    parser.add_argument("-a", "--adapter", type=str, default="AGATCGGAAGAG", help="Adapter for cutadapt")
    parser.add_argument("-A", type=str, help="Adapter for paired read cutadapt")
    return parser


def main():
    args = argument_parser().parse_args()
    run = functools.partial(subprocess.run, stdout=sys.stdout, stderr=sys.stderr)
    run(["fastq-dump", "--gzip", "--split-3", args.accession])
    read_one = Path(args.accession + "_1.fastq.gz")
    read_two = Path(args.accession + "_2.fastq.gz")
    read_unpaired = Path(args.accession + ".fastq.gz")

    if read_one.exists() and read_two.exists():

        cutadapt_out_one = read_one.with_suffix(".cutadapt.gz")
        cutadapt_out_two = read_two.with_suffix(".cutadapt.gz")

        run(["cutadapt",
             "-a", args.a,
             "-A", args.A,
             "-o", cutadapt_out_one,
             "-p", cutadapt_out_two,
             read_one, read_two])

        output_sam = Path(args.accession + ".sam")
        run(["bowtie", args.bowtie_index,
             "-1", cutadapt_out_one,
             "-2", cutadapt_out_two,
             "--sam", output_sam])

        run(["samtools"])
        aligned_sam = output_sam.with_suffix("aligned.sam")


        aligned_read_one = cutadapt_out_one.with_suffix(".aligned.gz")
        aligned_read_two = cutadapt_out_two.with_suffix(".aligned.gz")


        run(["picard", "SamToFastq",])


