include { ALIGN; GATK } from "${params.projectDirectory}/modules"

workflow {
rawfastq = Channel.fromPath("${params.homeDir}/samplesheet.csv")
    . splitCsv( header:true )
    . map { row ->
        def meta = [name:row.name, run:row.run, bryja:row.bryja]
        def baseDir = new File("${params.baseDir}")
                def runDir = baseDir.listFiles(new FilenameFilter() {
                        public boolean accept(File dir, String name) {
                                return name.endsWith(meta.run)
                        }
                })[0] //get the real folderName that has prepended date
        [meta.name, meta, [
            file("${runDir}/processed_fastq/${meta.name}_R1.fastq.gz", checkIfExists: true),
            file("${runDir}/processed_fastq/${meta.name}_R2.fastq.gz", checkIfExists: true),
        ]]
    }
     . view()

aligned	= ALIGN(rawfastq)
varcalling = GATK(aligned)
//normalizovany = NORMALIZACE(varcalling)
//anotovany = ANOTACE(normalizovany)
//anotovany2 = VCF2TXT(anotovany)
//coverage = COVERAGE(aligned)
}
