process ALIGN {
	tag "ALIGN on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
      
        label "l_cpu"
        label "l_mem"
 
	
	input:
	tuple val(name), val(sample), path(fwd), path(rev)

	output:
    tuple val(name), val(sample), path("${name}.mdup.sorted.bam"), path("${name}.mdup.sorted.bai")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	echo ALIGN $name
	source activate bwa
	bwa mem -R ${rg} -t $task.cpus ${params.refindex}.fa $fwd $rev | samblaster |samtools view -Sb - | sambamba sort /dev/stdin -o ${name}.mdup.sorted.bam
        samtools index ${name}.mdup.sorted.bam ${name}.mdup.sorted.bai
	"""
}

process CRAM {
        tag "CRAM on $name using $task.cpus CPUs and $task.memory memory"
        publishDir "${params.outDirectory}/${sample.run}/crams/", mode:'copy'
        container "staphb/samtools:1.20"
        label "l_cpu"
        label "l_mem"

   
        input:
        tuple val(name), val(sample), path(bam), path(bai)
        
        output:
        tuple val(name), val(sample), path("${name}.cram"), path("${name}.cram.crai")

        script:
        """
        echo CRAM $name
        samtools view -T ${params.refindex}.fa -C -o ${name}.cram $bam
        samtools index ${name}.cram
        """
}

process COVERAGE1 {
          tag "COVERAGE1 on $name"
       publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
         container "staphb/samtools:1.20"

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.coveragefin.txt")

        script:
        """
        echo COVERAGE1 $name
        samtools bedcov ${params.varbed1} $bam -d 20 > ${name}.COV
        awk '{print \$5/(\$3-\$2)}'  ${name}.COV >  ${name}.COV-mean
        awk '{print (\$6/(\$3-\$2))*100"%"}' ${name}.COV > ${name}-procento-nad-20
        paste ${name}.COV-mean ${name}-procento-nad-20 > vysledek
        echo "chr" "start" "stop" "name" ${name}.COV-mean ${name}-procento-nad-20 > hlavicka
        sed -i 's/ /\t/'g hlavicka
        paste ${params.varbed1} vysledek > coverage
        cat hlavicka coverage > ${name}.coveragefin.txt
        sed -i -e "s/\r//g" ${name}.coveragefin.txt
        """
}

process COMBINECOVERAGEMEAN {
    tag "COMBINECOVERAGEMEAN"

    input:
    tuple val(run_name), path(coverage_files)

    publishDir { "${params.outDirectory}/${run_name}/mapped/" }, mode: 'copy'

    output:
    path "coveragemeanALL"


    script:
    """
    i=0
    tmpfiles=""
    for file in ${coverage_files}; do
        awk '{print \$5}' "\$file" > column_\$i.txt
        tmpfiles="\$tmpfiles column_\$i.txt"
        i=\$((i+1))
    done
    echo "chr" "start" "stop" "name" > hlavicka
    sed -i 's/ /\t/'g hlavicka
    cat hlavicka ${params.varbed1} > bedshlavickou
    paste  bedshlavickou \$tmpfiles > coveragemeanALL
    """
}

process COMBINECOVERAGEPROCENTA {
    tag "COMBINECOVERAGEPROCENTA"

    input:
    tuple val(run_name), path(coverage_files)

    publishDir { "${params.outDirectory}/${run_name}/mapped/" }, mode: 'copy'

    output:
    path "coverageprocentoALL"


    script:
    """
    i=0
    tmpfiles=""
    for file in ${coverage_files}; do
        awk '{print \$6}' "\$file" > column_\$i.txt
        tmpfiles="\$tmpfiles column_\$i.txt"
        i=\$((i+1))
    done
    echo "chr" "start" "stop" "name" > hlavicka
    sed -i 's/ /\t/'g hlavicka
    cat hlavicka ${params.varbed1} > bedshlavickou
    paste  bedshlavickou \$tmpfiles > coverageprocentoALL
    """
}


workflow {
        rawfastq = Channel.fromPath("${params.homeDir}/samplesheet.csv")
    .splitCsv(header: true)
    .map { row ->
        def baseDir = new File("${params.baseDir}")
        def runDir = baseDir.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.endsWith(row.run)
            }
        })[0] //get the real folderName that has prepended date

        def fileR1 = file("${runDir}/processed_fastq/${row.name}_R1.fastq.gz", checkIfExists: true)
        def fileR2 = file("${runDir}/processed_fastq/${row.name}_R2.fastq.gz", checkIfExists: true)

                def meta = [name: row.name, run: row.run]
        [
            meta.name,
            meta,
            fileR1,
            fileR2,
                ]
    }
     . view()

aligned = ALIGN(rawfastq)
cramy = CRAM(aligned)

coverage_results = COVERAGE1(aligned)
coverage_files_collected = coverage_results
    .map { name, sample, f -> tuple(sample.run, file(f)) }
    .groupTuple() // groups by sample.run automatically!
finalcoverage = COMBINECOVERAGEMEAN(coverage_files_collected)
finalprocenta = COMBINECOVERAGEPROCENTA(coverage_files_collected)
}
