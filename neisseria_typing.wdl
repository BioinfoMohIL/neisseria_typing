version 1.0

# task ReadsQc {
# 	input {
# 		Array[File] reads
#         Int cpu = 8
# 		String docker = "quay.io/staphb/fastqc:0.12.1"
# 	}

#     output {
# 		Array[File] reports = glob("outputs/*.html")
# 	}

# 	command <<<
# 		mkdir results && fastqc -o results -q ~{sep=" " reads} -t ~{cpu}
#     >>>

# 	runtime {
# 		cpu: cpu
# 		docker: "${docker}"
# 	}	
	
# }

task Trimming {
    input {
        Array[File] reads
        String      docker="staphb/trimmomatic:0.39"
        Int         minlen = 50
        Int         head_crop =  20      
        Int         window_size = 4
        Int         quality_trim_score = 20
        Int         trailing = 3
        Int         crop = 265
        Int         cpus = 4
        String      memory = "8 GB"
    }

    output {
        File response = stdout()
    }


    command <<<
        ARRAY=(~{sep=" " reads})
        samples=()

        # Fetch Samplename
        for file in "${ARRAY[@]}"; do
            if [[ ${file} == *"_R1."* ]]; then
              modified=$(echo "${file}" | sed 's/_[^_]*$//')

              samples+=(${modified})
            fi
        done

        echo "${samples[@]}"


        # trimmomatic PE \
        # -phred33 \
        # -threads ~{cpus} \
        # ${read1} ${read2} \
        # -baseout ${samplename}.fastq.gz \
        # SLIDINGWINDOW:${window_size}:${quality_trim_score} \
        # MINLEN:${trimmomatic_minlen} > ${samplename}.trim.stats.txt
    >>>

}
#   output {
#     File       read1_trimmed = "${samplename}_R1.fastq.gz"
#     File       read2_trimmed = "${samplename}_R2.fastq.gz"
#     File       trimmomatic_stats = "${samplename}.trim.stats.txt"
#     String     version = read_string("VERSION")
#     String     pipeline_date = read_string("DATE")
#   }

#   runtime {
#     docker:       "${docker}"
#     memory:       "${memory}"
#     cpu:          cpus
#     preemptible:  0
#   }

# }

# task Assembly {
#   input {
#     File read1_trimmed
#     File read2_trimmed
#     String samplename
#     String docker = 'quay.io/staphb/shovill:1.1.0-2022Dec'
#     Int minlen
#   }

#   command <<<
#     shovill --version | head -l | tee VERSION
#     shovill \
#     --outdir out \
#     --R1 ~{read1_trimmed} \
#     --R2 ~{read2_trimmed} \
#     --minlen ~{minlen} \
#     mv out/contigs.fa out/~{samplename}_contigs.fasta
#   >>>

#   output {
#     File assembly_fasta = "out/~{samplename}_contigs.fasta"
#     File contig_gfa = "out/~{samplename}_contigs.gfa"
#     String shovill_version = read_string("VERSION")
#   }

#   runtime {
#     docker: "~{docker}"
#     memory: "16GB"
#     cpu: 4
#   }

# }

# task Typing {
#     input {
#         File assemblies
#     }

#     command <<<
    #     mkdir results

    #     samples=(~{sep=" " assemblies})
        
    #     for sample in ${samples[@]};do
    #         neisseria_typing --input ${sample} --output results
    #     done
    # >>>


#     runtime {
#         docker: "bioinfomoh/neisseria_typing:1"
#         memory: "8 GB"
#     }
# }

workflow NeisseriaTyping {

    input {
		Array[File] reads

	}

	# call ReadsQc {
	# 	input:
	# 		fastqs = fastqs
	# }

    call Trimming {
		input:
			reads = reads
	}

	# output {
	# 	Array[File] reports = ReadsQc.reports
	# }

	meta {
		author: "David Maimoun- MOH Israel"
	}
 
}
