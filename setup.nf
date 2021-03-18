#!/usr/bin/env nextflow

// Installation script

humanGTF = Channel.value('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz')
humanFasta = Channel.value('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz')
mouseGTF = Channel.value('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz')
mouseFasta = Channel.value('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz')
dbSNP = Channel.value('https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.38.gz')

// check for the cosmic files 

def cencusTable = new File("$projectDir/data/CosmicMutantExportCensus.tsv.gz")
def cosmicVCF = new File("$projectDir/data/CosmicCodingMuts.vcf.gz")

if (!cencusTable.exists() | !cosmicVCF.exists() ) {

println """
One or more of the cosmic files is missing!
Please make sure both "CosmicMutantExportCensus.tsv.gz" and "CosmicCodingMuts.vcf.gz" are present in the "RNA2CM/data" directory.
The files can be downloaded at https://cancer.sanger.ac.uk/cosmic/download
"""

} else {    

    file("$projectDir/data/GRCh").mkdir()
    file("$projectDir/data/GRCm").mkdir()

    process downloadHumanGTF {
        
        input:
        val gtf from humanGTF
        
        output:
        file "gencode.v34.primary_assembly.annotation.gtf" into gtf4Intervals, gtf4Star

        """
        wget $gtf
        gunzip gencode.v34.primary_assembly.annotation.gtf.gz
        """ }
        
    process downloadHumanFasta {
        
        publishDir "$projectDir/data", mode: 'copy' 
        
        input:
        val fasta from humanFasta
        
        output:
        file "GRCh38.primary_assembly.genome.fa" into fasta4Indexing, fasta4Dict, fasta4Star

        """
        wget $fasta
        gunzip GRCh38.primary_assembly.genome.fa.gz
        """ }

    process downloadMouseGenome {
        
        input:
        val fasta from mouseFasta
        val gtf from mouseGTF
        
        output:
        file "GRCm38.primary_assembly.genome.fa" into mouseFasta4Star
        file "gencode.vM25.primary_assembly.annotation.gtf" into mouseGTF4Star 

        """
        wget $fasta
        wget $gtf
        gunzip GRCm38.primary_assembly.genome.fa.gz gencode.vM25.primary_assembly.annotation.gtf.gz
        """ }

    process downloadDBSNP {
        
        input: 
        val vcf from dbSNP
        
        output:
        file "GCF_000001405.38.gz" into dbSNP4Rename
        file "GCF_000001405.38.gz.tbi" into dbSNPIndex4Rename 

        """
        wget --no-check-certificate $vcf
        wget --no-check-certificate ${vcf}.tbi
        """ }
        
    process generateStarIndex {
        
        cpus params.cpu
        memory '40GB'
            
        input:
        file fasta from fasta4Star
        file gtf from gtf4Star
        path genomeDir from Channel.fromPath("$projectDir/data/GRCh", type: 'dir')
        
        output:
        val 'foo' into bar
        path "*" into out
        
        """
        STAR --runThreadN $params.cpu --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $fasta --sjdbGTFfile $gtf --sjdbOverhang $params.readlength 
        """  }


    process generateStarIndexMouse {
        
        cpus params.cpu
        memory '40GB'
        
        input:
        file fasta from mouseFasta4Star
        file gtf from mouseGTF4Star
        path genomeDir from Channel.fromPath("$projectDir/data/GRCm", type: 'dir')
        val 'foo' from bar
        
        output:
        path "*" into outM
        """
        STAR --runThreadN $params.cpu --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $fasta --sjdbGTFfile $gtf --sjdbOverhang $params.readlength
        rm -rf $projectDir/data/GRCm38.primary_assembly.genome.fa $projectDir/data/gencode.vM25.primary_assembly.annotation.gtf $projectDir/data/gencode.v34.primary_assembly.annotation.gtf
        """  }

    process indexFasta {
         
        publishDir "$projectDir/data", mode: 'copy'
        
        input:
        file fasta from fasta4Indexing 
         
        output: 
        path '*' into fastaIndex

        """
        samtools faidx $fasta 
        """  }
        
    process createDictionery {
        
        publishDir "$projectDir/data", mode: 'copy'
        memory '8GB'
        
        input:
        file fasta from fasta4Dict
        file index from fastaIndex
         
        output: 
        path '*'

        """
        /gatk-4.1.3.0/gatk CreateSequenceDictionary -R $fasta 
        """  }
       
    process createIntervals {
        
        publishDir "$projectDir/data", mode: 'copy'
        
        input:
        file gtf from gtf4Intervals
         
        output: 
        file "GRCh38_exome.bed.gz" into exomeIntervals

        """
        awk '{if(\$3=="exon") {print \$1"\\t"\$4-100"\\t"\$5+100"\\t"substr(\$16,2,length(\$16)-3)}}' $gtf | sort -k 1,1 -k2,2n | bgzip > GRCh38_exome.bed.gz
        """ }

    process indexIntervals {
        
        publishDir "$projectDir/data", mode: 'copy'
        
        input:
        file intervals from exomeIntervals
        
        output: 
        path "GRCh38_exome.bed.gz.tbi"

        """
        tabix $intervals
        """ }

    process indexCosmic {
        
        publishDir "$projectDir/data", mode: 'copy'
        
        input:
        file cosmic from Channel.fromPath("$projectDir/data/CosmicCodingMuts.vcf.gz")
        
        output: 
        path "CosmicCodingMuts.vcf.gz.tbi" into CosmicIndex

        """
        tabix $cosmic
        """ }
        
    process renameDBSNP {
        
        cpus params.cpu
        publishDir "$projectDir/data", mode: 'copy' 
        
        input:
        file dbSNP from dbSNP4Rename
        file dbSNPindex from dbSNPIndex4Rename
        file remapNCBI from Channel.fromPath("$projectDir/data/remapNCBI.txt") 
        
        output: 
        path "dbSNPbuild154Renamed.vcf.gz" into dbsnpRenamed

        """
        bcftools annotate --threads $params.cpu --output-type z --rename-chrs $remapNCBI --output dbSNPbuild154Renamed.vcf.gz $dbSNP
        """ }

    process renameCosmic {
        
        cpus params.cpu
        publishDir "$projectDir/data", mode: 'copy' 
        
        input:
        file cosmic from Channel.fromPath("$projectDir/data/CosmicCodingMuts.vcf.gz")
        file index from CosmicIndex
        file remapCosmic from Channel.fromPath("$projectDir/data/remapCOSMIC.txt")
        
        output: 
        path "CosmicCodingMutsRenamed.vcf.gz" into cosmicRenamed

        """
        bcftools annotate --threads $params.cpu --output-type z --rename-chrs $remapCosmic --output CosmicCodingMutsRenamed.vcf.gz $cosmic
        """ }

    process indexRenamedVCFs {
        
        publishDir "$projectDir/data", mode: 'copy'

        input:
        file dbSNP from dbsnpRenamed
        file cosmic from cosmicRenamed
        
        output: 
        path "dbSNPbuild154Renamed.vcf.gz.tbi"
        path "CosmicCodingMutsRenamed.vcf.gz.tbi"

        """
        tabix $dbSNP
        tabix $cosmic
        rm -rf $projectDir/data/CosmicCodingMuts.vcf.gz $projectDir/data/CosmicCodingMuts.vcf.gz.tbi  
        """ }
}
    
