#!/usr/bin/env nextflow

// pipeline for identification of cancer related mutations

params.fastq = ''
params.fastq2 = false

params.cpu = 8
params.keepInter = true
params.filterMouse = true

if ( params.fastq2 == false ) {

    prefix = ("$params.fastq" =~ /(.+)\.(?:fastq.gz|fastq|fq|fq.gz)/)[0][1]
    
    if ( prefix.contains("/") ) {
	    prefix = prefix.split("/").last() }
    
    process trimmomaticSE {
        
        cpus params.cpu
	    memory '16GB'
        
        input:
        path fastq from Channel.fromPath("$params.fastq") 
        path adapters from Channel.fromPath("$projectDir/data/CommonAdapters.fa")
        
        output: 
        path "${prefix}.trimmed.fastq.gz" into trimmed, trimmed4mouse
        
        """
        java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads $params.cpu $fastq ${prefix}.trimmed.fastq.gz ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """ 
    }
      
             
    process alignToHumanGenomeSE {
          
        cpus params.cpu
		memory '48GB'
        
        input:
        path fastq from trimmed
        path genomeDir from Channel.fromPath("$projectDir/data/GRCh", type: 'dir')
        
        output:
        path "${prefix}Aligned.sortedByCoord.out.bam" into aligned2human
        val 'foo' into bar
        
        """
        STAR --runThreadN $params.cpu --genomeDir $genomeDir --readFilesIn $fastq --outFileNamePrefix ${prefix} --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outSAMattributes NM --twopassMode Basic --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.1
        """
    }


    process alignToMouseGenomeSE {
        
        cpus params.cpu
		memory '48GB'
           
        input:
        path fastq from trimmed4mouse
        path genomeDir from Channel.fromPath("$projectDir/data/GRCm", type: 'dir')
        val 'foo' from bar
        
        output:
        path "${prefix}GRCmAligned.sortedByCoord.out.bam" into aligned2mouse
        
        when:
        params.filterMouse == true
        
        """
        STAR --runThreadN $params.cpu --genomeDir $genomeDir --readFilesIn $fastq --outFileNamePrefix ${prefix}GRCm --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outSAMattributes NM --twopassMode Basic --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.1
        """ 
    }

} else {

    prefix = ("$params.fastq" =~ /(.+)_*[a-zA-Z0-9]*\.(?:fastq.gz|fastq|fq|fq.gz)/)[0][1] 

	if ( prefix.contains("/") ) {
    	prefix = prefix.split("/").last() }


    process trimmomaticPE {
        
        cpus params.cpu
        memory '16GB'
	
        input:
        path fastq1 from Channel.fromPath("$params.fastq") 
        path fastq2 from Channel.fromPath("$params.fastq2") 
        path adapters from Channel.fromPath("$projectDir/data/CommonAdapters.fa")
         
        output: 
        path "${prefix}_1.trimmed.fastq.gz" into trimmed1
        path "${prefix}_2.trimmed.fastq.gz" into trimmed2
        path "${prefix}_1.trimmed.fastq.gz" into trimmed4mouse1
        path "${prefix}_2.trimmed.fastq.gz" into trimmed4mouse2
        
        """
        java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $params.cpu $fastq1 $fastq2 ${prefix}_1.trimmed.fastq.gz unpaired1.fastq.gz ${prefix}_2.trimmed.fastq.gz unpaired2.fastq.gz ILLUMINACLIP:$adapters:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
        """ 
    }


    process alignToHumanGenomePE {
        
        cpus params.cpu
		memory '48GB'
           
        input:
        path fastq1 from trimmed1
        path fastq2 from trimmed2
        path genomeDir from Channel.fromPath("$projectDir/data/GRCh", type: 'dir')
        
        output:
        path "${prefix}Aligned.sortedByCoord.out.bam" into aligned2human
        
        """
        STAR --runThreadN $params.cpu --genomeDir $genomeDir --readFilesIn $fastq1 $fastq2 --outFileNamePrefix ${prefix} --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outSAMattributes NM --twopassMode Basic --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.1
        """ 
    }
     
        
    process alignToMouseGenomePE {
        
        cpus params.cpu
		memory '48GB'
                
        input:
        path fastq1 from trimmed4mouse1
        path fastq2 from trimmed4mouse2
        path genomeDir from Channel.fromPath("$projectDir/data/GRCm", type: 'dir')
        
        output:
        path "${prefix}GRCmAligned.sortedByCoord.out.bam" into aligned2mouse
        
        when:
        params.filterMouse == true
        
        """
        STAR --runThreadN $params.cpu --genomeDir $genomeDir --readFilesIn $fastq1 $fastq2 --outFileNamePrefix ${prefix}GRCm --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outSAMattributes NM --twopassMode Basic --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.1
        """ 
        }
}

if ( params.filterMouse == true ) {
    
    process XenofilteR {
        
        if (params.keepInter == true) {
        	publishDir "$launchDir", mode: 'copy'}
        cpus params.cpu
		memory '16GB'
        
        input:
        path bam1 from aligned2human
        path bam2 from aligned2mouse
        
        output:
        file "Filtered_bams/${prefix}_Filtered.bam" into filteredBam
		
        """
        #!/usr/bin/Rscript --save
        library("XenofilteR")
        bp.param <- SnowParam(workers = $params.cpu, type = "SOCK")
        sample.list <- matrix(c('$bam1','$bam2'),ncol=2)
        output.names <- c('${prefix}')
        XenofilteR(sample.list, destination.folder = "./", MM_threshold = 8, bp.param = bp.param, output.names)
        """ 
    }
    
    
} else {

    process skipXenofilteR {
        
		if (params.keepInter == true) {
            publishDir "$launchDir", mode: 'copy'}
        
        input:
        path bam from aligned2human
        
        output:
        path bam into filteredBam 
        
        """
        echo "--- mouse read filtering was not performed ---"
        """
    }
}

process markDuplicates {
    
    memory '16GB'
    
    input:
    path bam from filteredBam
    
    output:
    path "${prefix}marked_duplicates.bam" into markedDuplicates 

    """
    /gatk-4.1.3.0/gatk MarkDuplicates --CREATE_INDEX true --I $bam --O ${prefix}marked_duplicates.bam --VALIDATION_STRINGENCY SILENT --M ${prefix}marked_dup_metrics.txt
    """ 
}
    
    
process splitNcigar {
    
    memory '16GB'
    
    input:
    path bam from markedDuplicates
    path intervals from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz")
    path intIndex from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz.tbi")
    path reference_genome from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa")
    path index from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa.fai")
    path dict from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.dict")
    
    output:
    path "${prefix}splitN.bam" into splitN 

    """
    /gatk-4.1.3.0/gatk SplitNCigarReads -L $intervals -R $reference_genome -I $bam -O ${prefix}splitN.bam
    """ 
}
    
    
process addGroups {
    
    memory '16GB'
    
    input:
    path bam from splitN
    
    output:
    path "${prefix}.grouped.bam" into grouped4BQSR
    path "${prefix}.grouped.bai" into grouped4BQSRindex
    path "${prefix}.grouped.bam" into grouped4applyBQSR
    path "${prefix}.grouped.bai" into grouped4applyBQSRindex
    """
    /gatk-4.1.3.0/gatk AddOrReplaceReadGroups --CREATE_INDEX true --I $bam --O ${prefix}.grouped.bam --RGID rnasq --RGLB lb --RGPL illumina --RGPU pu --RGSM $prefix
    """  
}

    
process baseQualityRecalibration {
    
    memory '16GB'
    
    input:
    path bam from grouped4BQSR
    path index from grouped4BQSRindex
    path intervals from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz")
    path intIndex from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz.tbi")
    path reference_genome from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa")
    path index from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa.fai")
    path dict from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.dict")
    path dbSNP from Channel.fromPath("$projectDir/data/dbSNPbuild154Renamed.vcf.gz")
    path dbSNPindex from Channel.fromPath("$projectDir/data/dbSNPbuild154Renamed.vcf.gz.tbi")
    
    output:
    path "${prefix}.recal_data.table" into recalTable 

    """
    /gatk-4.1.3.0/gatk BaseRecalibrator -L $intervals -I $bam --use-original-qualities --disable-sequence-dictionary-validation true -R $reference_genome --known-sites $dbSNP -O ${prefix}.recal_data.table
    """ 
}
    
    
process applyBQSR {
    
    memory '16GB'
    
    input:
    path table from recalTable
    path bam from grouped4applyBQSR
    path bai from grouped4applyBQSRindex
    path intervals from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz")
    path intIndex from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz.tbi")
    path reference_genome from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa")
    path index from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa.fai")
    path dict from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.dict")
    
    output:
    path "${prefix}.recal_output.bam" into recalibrated
    path "${prefix}.recal_output.bai" into recalibratedIndex  

    """
    /gatk-4.1.3.0/gatk ApplyBQSR -L $intervals -R $reference_genome -I $bam --use-original-qualities --add-output-sam-program-record --bqsr-recal-file $table -O ${prefix}.recal_output.bam
    """ 
}



process callVariants {
    
    if (params.keepInter == true) {
        publishDir "$launchDir", mode: 'copy'}
    memory '16GB'
    
    input:
    path bam from recalibrated
    path bai from recalibratedIndex
    path intervals from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz")
    path intIndex from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz.tbi")
    path reference_genome from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa")
    path index from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa.fai")
    path dict from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.dict")
    
    output:
    file "${prefix}.output.vcf.gz" into variants
    file "${prefix}.output.vcf.gz.tbi" into variantsIndex 

    """
    /gatk-4.1.3.0/gatk HaplotypeCaller -L $intervals -R $reference_genome -I $bam -O ${prefix}.output.vcf.gz --dont-use-soft-clipped-bases --pcr-indel-model AGGRESSIVE
    """ 
}

    
process hardFilter {
    
    memory '16GB'
    
    input:
    path vcf from variants
    path vcfIndex from variantsIndex
    path reference_genome from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa")
    path index from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.fa.fai")
    path dict from Channel.fromPath("$projectDir/data/GRCh38.primary_assembly.genome.dict")
    
    output:
    path "${prefix}.hardfilter.vcf.gz" into hardfilterd


    """
    /gatk-4.1.3.0/gatk VariantFiltration --R $reference_genome --V $vcf --window 35 --cluster 3 --filter-name "FS" --filter "FS > 30.0" --filter-name "QD" --filter "QD < 2.0" -O ${prefix}.hardfilter.vcf.gz
    """ 
}    

    
process filterVariants {
    
    cpus params.cpu
    
    input:
    path hardfilterd
   
    output:
    path "${prefix}filtered.vcf.gz" into filterd

    """
    bcftools view --threads $params.cpu -i 'FILTER="PASS" && FORMAT/DP >= 10 && FORMAT/AD[:1] >= 5' --output-type z --output-file ${prefix}filtered.vcf.gz ${prefix}.hardfilter.vcf.gz
    """ } 

process IndexfilteredVariants {
        
    input:
    path filterd
   
    output:
    path "${prefix}filtered.vcf.gz" into filterd4annotation
    path "${prefix}filtered.vcf.gz.tbi" into filteredIndex  

    """
    tabix "${prefix}filtered.vcf.gz"
    """ } 
      
process geneAnnotation {
    
    cpus params.cpu
    
    input:
    path filterd4annotation
    path intervals from Channel.fromPath("$projectDir/data/GRCh38_exome.bed.gz")
    path header from Channel.fromPath("$projectDir/data/header.txt")
    path filteredIndex
   
    output:
    path "${prefix}named.vcf.gz" into gene
    path "${prefix}named.vcf.gz.tbi" into geneIndex  

    """
    bcftools annotate --threads $params.cpu -a $intervals -h $header -c CHROM,FROM,TO,Gene --output-type z --output ${prefix}named.vcf.gz ${prefix}filtered.vcf.gz
    tabix ${prefix}named.vcf.gz
    """ } 

process snpAnnotation {
    
    cpus params.cpu
    
    input:
    path gene
    path geneIndex
    path dbSNP from Channel.fromPath("$projectDir/data/dbSNPbuild154Renamed.vcf.gz")
    path dbSNPindex from Channel.fromPath("$projectDir/data/dbSNPbuild154Renamed.vcf.gz.tbi")
    
    output:
    path "${prefix}dbSNP.vcf.gz" into snp
    path "${prefix}dbSNP.vcf.gz.tbi" into snpIndex  

    """
    bcftools annotate --threads $params.cpu -a $dbSNP -c INFO/RS,INFO/COMMON --output-type z --output ${prefix}dbSNP.vcf.gz ${prefix}named.vcf.gz
    tabix ${prefix}dbSNP.vcf.gz
    """ } 

process cosmicAnnotation {
   
    cpus params.cpu
       
    input:
    path snp
    path snpIndex
    path cosmic_vcf from Channel.fromPath("$projectDir/data/CosmicCodingMutsRenamed.vcf.gz")
    path cosmicIndex from Channel.fromPath("$projectDir/data/CosmicCodingMutsRenamed.vcf.gz.tbi")
    
    output:
    path "${prefix}.annotated.vcf.gz" into cosmic
    path "${prefix}.annotated.vcf.gz.tbi" into cosmicIndex  

    """
    bcftools annotate --threads $params.cpu -a $cosmic_vcf -c ID,INFO/CNT --output-type z --output ${prefix}.annotated.vcf.gz ${prefix}dbSNP.vcf.gz
    tabix ${prefix}.annotated.vcf.gz
    """
    } 
    
process variantTable {
    
    if (params.keepInter == true) {
        publishDir "$launchDir", mode: 'copy'}
    
    input:
    path cosmic
    path cosmicIndex
   
    output:
    path "${prefix}varTable.tsv" into table  

    """
    bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/Gene\t%ID\t%INFO/CNT\t%INFO/RS\t%INFO/COMMON\t[%AD]\n' ${prefix}.annotated.vcf.gz > ${prefix}varTable.tsv
    """
    } 
    
process findCancerMutations {
    
    publishDir "$launchDir", mode: 'copy'
        
    input:
    path varTable from table
    path cosmic from Channel.fromPath("$projectDir/data/CosmicMutantExportCensus.tsv.gz")
    
    output:
    path "${prefix}_cancer_mutations.csv" into results
    
    """
    #!/usr/bin/python3

    import pandas as pd

    df = pd.read_csv('$varTable', sep='\t')

    df = df[df['[6]ID'] != '.']
    df['[7]CNT'] = df['[7]CNT'].astype("int32")
    df = df[df['[7]CNT'] >= 20]

    df = df[df['[9]COMMON'] != '1']
   
    mutations = pd.read_csv('$cosmic', sep='\t', compression='gzip', encoding='latin1')
    mutations = mutations[['Tier','GENOMIC_MUTATION_ID', 'Mutation AA', 'Mutation Description', 'FATHMM prediction','FATHMM score']]

    df = df.merge(mutations, left_on='[6]ID', right_on='GENOMIC_MUTATION_ID')

    df = df[df['Tier'] == 1]
    df = df[df['Mutation Description'] != 'Substitution - coding silent']
    df = df[df['FATHMM prediction'] == 'PATHOGENIC']

    df.drop_duplicates(subset='[6]ID', inplace=True)
    df = df[['[5]Gene', '# [1]CHROM', '[2]POS', '[3]REF', '[4]ALT', '[7]CNT', '[8]RS', 'GENOMIC_MUTATION_ID', 'Mutation AA', 'Mutation Description', 'FATHMM score','[10]$prefix:AD']]
    df.rename(columns = {'[5]Gene': 'Gene', '# [1]CHROM': 'CHROM', '[2]POS': 'POS', '[3]REF':'REF','[4]ALT': 'ALT', '[7]CNT': 'COSMIC_CNT','[8]RS': 'RS(dbSNP)', 'Mutation AA': 'Mutation_AA', 'Mutation Description': 'Mutation_Description', '[10]$prefix:AD': 'AD'}, inplace=True)

    df.to_csv('${prefix}_cancer_mutations.csv', index=False)
    """  
}
