process ALIGN {
	tag "ALIGN on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
      
        label "l_cpu"
        label "l_mem"
 
	
	input:
	tuple val(name), val(sample), path(reads)

	output:
    tuple val(name), val(sample), path("${name}.mdup.sorted.bam"), path("${name}.mdup.sorted.bai")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	echo ALIGN $name
	source activate bwa
	bwa mem -R ${rg} -t $task.cpus ${params.refindex}.fa $reads | samblaster |samtools view -Sb - | sambamba sort /dev/stdin -o ${name}.mdup.sorted.bam
        samtools index ${name}.mdup.sorted.bam ${name}.mdup.sorted.bai
	"""
}

process GATK {
       tag "GATK on $name"
       publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'
        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.vcf")

        script:
        """
        echo GATK $name
        source activate gatk4610
        gatk --java-options "-Xmx4g" HaplotypeCaller -R ${params.ref}.fa -I $bam -L ${params.varbed}  --dont-use-soft-clipped-bases true -A StrandBiasBySample -minimum-mapping-quality 0 --mapping-quality-threshold-for-genotyping 0 --enable-dynamic-read-disqualification-for-genotyping true --flow-filter-alleles-qual-threshold 0 -O ${name}.vcf
        """
}
