#---------------project_directory-------------#
project_dir = '/home/TW/project/OA_RA/ATAC_part'

#-------------clean_data_directory------------#
clean_data_dir = '/home/TW/data/Inhouse/ATAC/OA_RA'

#----------------bowtie2_index----------------#
bowtie2_index = '/home/TW/data/index/bowtie2/hg19/hg19'

sample_list = ['OA5','RA8','1906218RA','1907607RA']

rule all:
    input:
        expand(project_dir + '/peak_results/{sample}_peaks.narrowPeak', sample = sample_list)

rule bowtie2_alignment:
    input:
        Read1= clean_data_dir + '/{sample}/clean_data/{sample}_1P.fq.gz',
        Read2= clean_data_dir + '/{sample}/clean_data/{sample}_2P.fq.gz'
    output:
        temp(project_dir + '/bowtie2_result/{sample}.sam')
    threads: 8
    shell:
        'bowtie2 -p {threads} -x {bowtie2_index} -1 {input.Read1} -2 {input.Read2} -S {output} -X 2000'

rule sam2bam:
    input:
        project_dir + '/bowtie2_result/{sample}.sam'
    output:
        project_dir + '/bowtie2_result/{sample}.bam'
    threads: 8
    shell:
        'samtools view -@ {threads} -b -S {input} -o {output}'

rule filter_bam:
    input:
        project_dir + '/bowtie2_result/{sample}.bam'
    output:
        project_dir + '/bowtie2_result/{sample}_filter.bam'
    threads:8
    shell:
        'samtools view -@ {threads} -b -f 2 -q 30 {input} -o {output}'

rule sort_bam_files:
    input:
        project_dir + '/bowtie2_result/{sample}_filter.bam'
    output:
        project_dir + '/bowtie2_result/{sample}_filter_sorted.bam'
    threads:8
    shell:
        'samtools sort -@ {threads} {input} -o {output}'

rule remove_duplication:
    input:
        project_dir + '/bowtie2_result/{sample}_filter_sorted.bam'
    output:
        project_dir + '/bowtie2_result/{sample}_filter_sorted_rmdup.bam',
        project_dir + '/bowtie2_result/{sample}_marked_dup_metrics.txt'
    shell:
        'picard MarkDuplicates I={input} O={output[0]} M={output[1]} REMOVE_DUPLICATES=True'

rule macs2_call_peak:
    input:
        project_dir + '/bowtie2_result/{sample}_filter_sorted_rmdup.bam'
    output:
        project_dir + '/peak_results/{sample}_peaks.narrowPeak'
    params:
        project_dir + '/peak_results/{sample}'
    shell:
        'macs2 callpeak -t {input} --nomodel -g hs -n {params} -q 0.05 --extsize 200 --shift -100 --keep-dup all -B --SPMR'
