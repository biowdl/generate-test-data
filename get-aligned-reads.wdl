version 1.0

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

import "tasks/bowtie.wdl" as bowtie
import "tasks/cutadapt.wdl" as cutadapt
import "tasks/samtools.wdl" as samtools
import "tasks/picard.wdl" as picard

workflow getAlignedReadsBowtie {
    input {
        File read1
        File? read2
        Array[File]+ bowtieIndex = [
            "reference/bowtie/reference.1.ebwt",
            "reference/bowtie/reference.2.ebwt",
            "reference/bowtie/reference.3.ebwt",
            "reference/bowtie/reference.4.ebwt",
            "reference/bowtie/reference.rev.1.ebwt",
            "reference/bowtie/reference.rev.2.ebwt"
        ]
        Array[String]+ adapter = ["AGATCGGAAGAG"]  # Illumina Universal Adapter
        Array[String]+? adapterRead2
    }


    call cutadapt.Cutadapt {
        input:
            adapter = adapter,
            adapterRead2 = adapterRead2,
            read1 = read1,
            read2 = read2
    }

    call bowtie.Bowtie {
        input:
            readsUpstream = [Cutadapt.cutRead1],
            readsDownstream = select_all([Cutadapt.cutRead2]),
            indexFiles = bowtieIndex
    }

    call samtools.View as samtoolsView {
        input:
            inFile = Bowtie.outputBam,
            outputBam = false,
            includeFilter = 2,  # Select all properly mapped reads
    }

    call picard.SamToFastq {
        input:
            inputBam = samtoolsView.outputFile,
            paired = defined(read2)
    }

    call SelectReads as selectR1 {
        input:
            original = read1,
            mapped = SamToFastq.read1
    }

    if (defined(read2)) {
        call SelectReads as selectR2 {
            input:
                original = select_first([read2]),
                mapped = select_first([SamToFastq.read2])
        }
    }

    output {
        File selectedRead1 = selectR1.outputFile
        File? selectedRead2 = selectR2.outputFile
    }
}

task SelectReads {
    input {
        File original
        File mapped
        String outputPath = basename(original, "\.f(ast)?q(.gz)?") + "selected.fastq.gz"
        String dockerImage = "quay.io/biocontainers/dnaio:0.3--py37h14c3975_1"
    }

    command <<<
        python3 ~{mapped} ~{original} ~{outputPath} \
        <<CODE
        import sys
        import dnaio

        with dnaio.FastqReader(sys.argv[1]) as mapped_reads:
            mapped_ids = {sequence.name for sequence in mapped_reads}
        with dnaio.FastqReader(sys.argv[2]) as original_reads:
            with dnaio.FastqWriter(sys.argv[3], "w") as output_file:
                for sequence in original_reads:
                    if sequence.name.split()[0] in mapped_ids:
                        output_file.write(sequence)
        CODE
    >>>

    output {
        File outputFile = outputPath
    }

    runtime {
        docker: dockerImage
    }
}