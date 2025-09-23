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

process VARDICT {

         tag "VARDICT on $name"
         publishDir "${params.outDirectory}/${sample.run}/viktor/", mode:'copy'

         input:
         tuple val(name), val(sample), path(bam), path(bai)
         output:
         tuple val(name), val(sample), path("${name}.vcf")

         script:
         """
         echo VARDICT $name

        source activate vardict
        vardict -G ${params.refindex}.fa -f 0.05 -N ${name} -b ${bam} -c 1 -S 2 -E 3 -g 4 -U ${params.varbed3} | Rscript ${params.teststrandbias} | perl ${params.var2vcf_valid} -N ${name} -f 0.05 -A > ${name}.vcf
       """
}



process VAFaNORMALIZACE {
        tag "VAFaNORMALIZACE on $name"
        publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path(gatk)

        output:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo VAFaNORMALIZACE $name

        bcftools +fill-tags $gatk -Ob -o ${name}.pom2.bcf -- -t FORMAT/VAF
        bcftools convert -O v -o ${name}.vaf.vcf ${name}.pom2.bcf 
         
        bcftools norm -m-both -f ${params.ref}.fa -o ${name}.norm.vcf ${name}.vaf.vcf
        bgzip ${name}.norm.vcf
        tabix ${name}.norm.vcf.gz

        """
}

process ANOTACE_ACGT {
        tag "ANOTACEACGT on $name"
        publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.norm.acgt.vcf"), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")

        script:
        """
        source activate gatk4610
        echo ANOTACEACGT $name
        gatk --java-options "-Xmx4g"  VariantAnnotator   -V ${name}.norm.vcf.gz -O ${name}.norm.acgt.vcf.gz --resource:ACGT ${params.ACGT} --expression ACGT.AF   --expression ACGT.AC   --expression ACGT.AC_Hom   --expression ACGT.AC_Het   --expression ACGT.AC_Hemi

        gunzip ${name}.norm.acgt.vcf.gz
        """
}

process ANOTACE_OMIM {
        tag "ANOTACEOMIM on $name"
        publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.acgt.vcf"), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.norm.vaf.omim.vcf"), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo ANOTACEOMIM $name

bcftools annotate -a ${params.AR} -h ${params.ARheader} -c CHROM,FROM,TO,dedicnostAR -l dedicnostAR:append -m -xx ${name}.norm.acgt.vcf > ${name}.pom3
bcftools annotate -a ${params.AD} -h ${params.ADheader} -c CHROM,FROM,TO,dedicnostAD -l dedicnostAD:append -m -aa ${name}.pom3 > ${name}.pom4
bcftools annotate -a ${params.Xlinked}  -h ${params.Xlinkedheader} -c CHROM,FROM,TO,dedicnostXlinked -l dedicnostXlinked:append -m -bb ${name}.pom4 > ${name}.pom5
bcftools annotate -a ${params.Ylinked}  -h ${params.Ylinkedheader} -c CHROM,FROM,TO,dedicnostYlinked -l dedicnostYlinked:append -m -cc ${name}.pom5 > ${name}.pom6

bcftools annotate -a ${params.omimphenotyp} -h ${params.omimphenotypheader} -c CHROM,FROM,TO,fenotyp -l fenotyp:append -m -yy ${name}.pom6 | sed  's/xx/dedicnostAR=0/' | grep -v 'INFO=<ID=dedicnostAR=0' | sed  's/yy/fenotyp=0/' | grep -v 'INFO=<ID=fenotyp=0' | sed  's/aa/dedicnostAD=0/' | grep -v 'INFO=<ID=dedicnostAD=0' | sed  's/bb/dedicnostXlinked=0/' | grep -v 'INFO=<ID=dedicnostXlinked=0' | sed  's/cc/dedicnostYlinked=0/' | grep -v 'INFO=<ID=dedicnostYlinked=0' > ${name}.norm.vaf.omim.vcf

        """
}

process MetaRNN {
        tag "MetaRNN on $name"
        publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.vaf.omim.vcf"), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.MetaRNN.pom.vcf"), path("${name}.norm.vaf.omim.vcf")

        script:
        """
        echo MetaRNN $name
        source activate MetaRNN
        python ${params.MetaRNN} hg38 ${name}.norm.vcf.gz
        cat ${name}.norm.vcf.gz.indel.annotated ${name}.norm.vcf.gz.nsSNV.annotated > ${name}.MetaRNN
        grep -v "#" ${name}.MetaRNN > a
        sed -i 's/ENST/Varsome=ENST/' a
        sed -i 's/;/ /g' a
        sed -i 's/ /*/g' a
        awk -F "\t" '{print \$1, \$2, ". " \$3, \$4,  "340", "PASS", \$5, \$6}' a > aa
        sed -i 's/ /\t/' aa
        sed -i 's/ /\t/' aa
        sed -i 's/ /\t/' aa
        sed -i 's/ /\t/' aa
        sed -i 's/ /\t/' aa
        sed -i 's/ /\t/' aa
        sed -i 's/ /\t/' aa
        sed -i 's/ /*/g' aa
        sort -k1,1 -k2,2n aa > aaa
        cat ${params.hlavicka}  aaa > ${name}.MetaRNN.pom.vcf


        """
}

process ANOTACE_MetaRNN {
        tag "ANOTACEMetaRNN on $name"
        publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.MetaRNN.pom.vcf"), path("${name}.norm.vaf.omim.vcf")

        output:
        tuple val(name), val(sample), path("${name}.norm.metarnn.vcf.gz"), path("${name}.norm.metarnn.vcf.gz.tbi")

        script:
        """
        source activate gatk4610
        echo ANOTACEMetaRNN $name
        bgzip ${name}.MetaRNN.pom.vcf
        tabix ${name}.MetaRNN.pom.vcf.gz

        bgzip ${name}.norm.vaf.omim.vcf
        tabix ${name}.norm.vaf.omim.vcf.gz

        gatk --java-options "-Xmx4g"  VariantAnnotator  -V ${name}.norm.vaf.omim.vcf.gz  -O ${name}.norm.metarnn.vcf  --resource:MetaRNN ${name}.MetaRNN.pom.vcf.gz --expression MetaRNN.Varsome

        bgzip ${name}.norm.metarnn.vcf
        tabix ${name}.norm.metarnn.vcf.gz
        """
}

process ANOTACE_annovar {
       tag "ANOTACE on $name"
       //publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.metarnn.vcf.gz"), path("${name}.norm.metarnn.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.norm.metarnn.vcf.gz.hg38_multianno.vcf.gz"), path("${name}.norm.metarnn.vcf.gz.hg38_multianno.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo ANOTACE $name

        ${params.annovar} -vcfinput ${name}.norm.metarnn.vcf.gz ${params.annovardb}  -buildver hg38 -protocol refGeneWithVer,ensGene,1000g2015aug_all,1000g2015aug_eur,exac03nontcga,avsnp150,clinvar_20240917,dbnsfp41c,gnomad41_exome,gnomad41_genome,cosmic70,revel,GTEx_v8_eQTL \
        -operation gx,g,f,f,f,f,f,f,f,f,f,f,f -nastring . -otherinfo -polish -xreffile ${params.gene_fullxref.txt} -arg '-splicing 20 -exonicsplicing',,,,,,,,,,,, --remove
        bgzip ${name}.norm.metarnn.vcf.gz.hg38_multianno.vcf
        tabix ${name}.norm.metarnn.vcf.gz.hg38_multianno.vcf.gz
        """
}

process VCF2TXT {
       tag "VCF2TXT on $name"
       publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.metarnn.vcf.gz.hg38_multianno.vcf.gz"), path("${name}.norm.metarnn.vcf.gz.hg38_multianno.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.final.txt")

        script:
        """
        echo VCF2TXT $name
        source activate gatk4610
        gatk --java-options "-Xmx4g" VariantsToTable -R ${params.ref}.fa  --show-filtered  -V ${name}.norm.metarnn.vcf.gz.hg38_multianno.vcf.gz -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF SB -GF VAF -F dedicnostAR -F dedicnostAD -F dedicnostXlinked -F dedicnostYlinked -F fenotyp  -F ACGT.AF -F ACGT.AC -F ACGT.AC_Hom -F ACGT.AC_Het -F ACGT.AC_Hemi -F Func.refGeneWithVer -F Gene.refGeneWithVer -F GeneDetail.refGeneWithVer -F ExonicFunc.refGeneWithVer -F AAChange.refGeneWithVer -F 1000g2015aug_all -F 1000g2015aug_eur  -F gnomad41_exome_AF -F gnomad41_exome_AF_nfe -F gnomad41_genome_AF -F gnomad41_genome_AF_nfe -F avsnp150 -F CLNSIG -F REVEL -F MetaRNN.Varsome -F SIFT_pred -F MutationTaster_pred -F Gene_full_name.refGeneWithVer -F FATHMM_pred -F PROVEAN_pred -F Function_description.refGeneWithVer -F Disease_description.refGeneWithVer -F Tissue_specificityUniprot.refGeneWithVer -F Expression-egenetics.refGeneWithVer --output ${name}.final.txt
        """
}

process VEP {
    tag "VEP on $name"
    publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy' 
    //container "ensemblorg/ensembl-vep:release_108.0"
    container "ensemblorg/ensembl-vep:release_114.1"
  
    input:
    tuple val(name), val(sample), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")
 
    output:
    tuple val(name), val(sample), path("${name}.vep.vcf")

    script:
    """
    vep -i ${name}.norm.vcf.gz --cache --cache_version 114 --dir_cache $params.vep \
        --fasta ${params.ref}.fa --merged --offline --vcf --hgvs -o ${name}.vep.vcf \
    --plugin CADD,snv=${params.caddsnv},indels=${params.caddindel} --dir_plugins ${params.vepplugin} --force_overwrite --no_stats --plugin AlphaMissense,file=${params.alfamissense} 
    """
}

process BIOPET {
     tag "BIOPET on $name"
     publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

     input:
     tuple val(name), val(sample), path(vep)

     output:
     tuple val(name), val(sample), path("${name}.vepalpfa.anot.vcf")

     script:
     """
     source activate biopet
     biopet tool VepNormalizer -I $vep -O ${name}.vepalpfa.anot.vcf -m standard
     """
}

process VCFTOTXTVEP {

        tag "VCFTOTXTVEP on $name"
        publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path(biopet)

        output:
        tuple val(name), val(sample), path("${name}.vepalpfa.anot.txt")

        script:
        """
        source activate py36

        python ${params.vcfsimplify} SimplifyVCF -toType table -inVCF $biopet -out ${name}.vepalpfa.anot.txt
        """
}

process spojitannovarVEP {

        tag "spojitannovarVEP on $name"
        publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path(final_txt), path(vep_txt)

        output:
        path("${name}.merged.txt")

        script:
        """
        echo "Merging ${final_txt} and ${vep_txt}"
        awk '{print \$1, \$2, \$4, \$5, \$28, \$29, \$55, \$56}' ${vep_txt}  > vyber
        sed -i 's/ /\t/'g vyber
        paste ${final_txt}  vyber > ${name}.merged.txt
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
varcalling = GATK(aligned)
varcalling2 = VARDICT(aligned)
normalizovany = VAFaNORMALIZACE(varcalling)

anotovanyacgt = ANOTACE_ACGT(normalizovany)
anotovanyomim = ANOTACE_OMIM(anotovanyacgt)
metarnnskore = MetaRNN(anotovanyomim)
anotovanymetarnn = ANOTACE_MetaRNN(metarnnskore)
anotovany = ANOTACE_annovar(anotovanymetarnn)
anotovanyfin = VCF2TXT(anotovany)

vepovany = VEP(normalizovany)
biopetovany = BIOPET(vepovany)
textovany = VCFTOTXTVEP(biopetovany)

combined = anotovanyfin.join(textovany, by: [0,1])
spojitannovarVEP(combined)

coverage_results = COVERAGE1(aligned)
coverage_files_collected = coverage_results
    .map { name, sample, f -> tuple(sample.run, file(f)) }
    .groupTuple() // groups by sample.run automatically!
finalcoverage = COMBINECOVERAGEMEAN(coverage_files_collected)
finalprocenta = COMBINECOVERAGEPROCENTA(coverage_files_collected)
}
