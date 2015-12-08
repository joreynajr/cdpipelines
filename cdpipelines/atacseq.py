import os

from general import _make_dir
from general import JobScript

class ATACJobScript(JobScript):
    def star_align(
        self,
        r1_fastq, 
        r2_fastq, 
        rgpl, 
        rgpu, 
        star_index, 
        threads,
        genome_load='LoadAndRemove',
        star_path='STAR',
    ):
        """
        Align paired fastq files with STAR.
    
        Parameters
        ----------
        r1_fastq : str 
            Path to R1 fastq file.
    
        r2_fastq : str 
            Path to R2 fastq file.
    
        rgpl : str
            Read Group platform (e.g. illumina, solid). 
    
        rgpu : str
            Read Group platform unit (eg. run barcode). 

        Returns
        -------
        bam : str
            Path to output alignment bam file.
        
        log_out : str
            Path to log file.

        log_final_out : str
            Path to final log file.
        
        log_progress_out : str
        Path to progress log file.
        
        sj_out : str
            Path to output SJ.out.tab file.

        """
        lines = (' \\\n\t'.join([
            star_path, 
            '--runThreadN {}'.format(threads),
            '--alignIntronMax 1',
            '--genomeDir {}'.format(star_index), 
            '--genomeLoad {}'.format(genome_load),
            '--readFilesCommand zcat',
            '--readFilesIn {} {}'.format(r1_fastq, r2_fastq),
            '--outSAMattributes All', 
            '--outSAMunmapped Within',
            '--outSAMattrRGline ID:1 PL:{} PU:{} LB:{} SM:{}'.format(
                rgpl, rgpu, self.sample_name, self.sample_name),
            '--outFilterMultimapNmax 20', 
            '--outFilterMismatchNmax 999',
            '--outFilterMismatchNoverLmax 0.04',
            '--seedSearchStartLmax 20',
            '--outFilterScoreMinOverLread 0.1',
            '--outFilterMatchNminOverLread 0.1',
            '--outSAMtype BAM Unsorted']))
        lines += '\n\n'
        lines += 'if [ -d _STARtmp ] ; then rm -r _STARtmp ; fi\n\n'
        bam = os.path.join(
            self.tempdir, '{}.bam'.format(self.sample_name))
        log_out = os.path.join(
            self.tempdir, '{}_Log.out'.format(self.sample_name))
        log_final_out = os.path.join(
            self.tempdir, '{}_Log.final.out'.format(self.sample_name))
        log_progress_out = os.path.join(
            self.tempdir, '{}_Log.progress.out'.format(self.sample_name))
        sj_out = os.path.join(
            self.tempdir, '{}_SJ.out.tab'.format(self.sample_name))
        lines += 'mv Aligned.out.bam {}\n'.format(bam)
        lines += 'mv Log.out {}\n'.format(log_out)
        lines += 'mv Log.final.out {}\n'.format(log_final_out)
        lines += 'mv Log.progress.out {}\n'.format(log_progress_out)
        lines += 'mv SJ.out.tab {}\n'.format(sj_out)
        lines += '\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return bam, log_out, log_final_out, log_progress_out, sj_out

    def filter_tlen(
        self,
        bam,
        tlen_max=140,
        samtools_path='samtools',
    ):
        """
        Filter out reads whose tlen column in bam file is greater than tlen_max.
    
        Parameters
        ----------
        bam : str
            Path to bam file to filter. Uniquely mapped reads should have
            quality score greater than 255 (e.g. STAR alignments).

        tlen_max : str
            All reads with tlen greater than this will be filtered out.
    
        Returns
        -------
        out_bam : str
            Path to output filtered bam file.
    
        """
        root = os.path.splitext(os.path.split(bam)[1])[0]
        out_bam = os.path.join(self.tempdir, '{}_tlen_leq_{}.bam'.format(
            root, tlen_max))
        lines = ('{} view -h {} \\\n\t| awk \'function abs(x){{return '
                 '((x < 0.0) ? -x : x)}} {{if (abs($9) <= {}) {{print}} '
                 'else if (substr($1,1,1) == "@") {{print}}}}\' \\\n\t'
                 '| {} view -Sb - \\\n\t> {}\n\n'.format(
                     samtools_path, bam, tlen_max, samtools_path, out_bam))
        with open(self.filename, "a") as f:
            f.write(lines)
        return out_bam

    def filter_multi_mt_blacklist_reads(
        self,
        bam,
        encode_blacklist,
        bedtools_path='bedtools',
        samtools_path='samtools',
    ):
        """
        Filter out read pairs that are not mapped uniquely, that map to the
        mitochondrial genome, or that overlap the ENCODE blacklist regions.
    
        Parameters
        ----------
        bam : str
            Path to bam file to filter. Uniquely mapped reads should have
            quality score greater than 255 (e.g. STAR alignments).

        encode_blacklist : str
            Path to ENCODE blacklist bed file.
    
        Returns
        -------
        out_bam : str
            Path to output filtered bam file.
    
        """
        out_bam = os.path.join(self.tempdir, '{}_filtered.bam'.format(
            self.sample_name))
        lines = (
            '{} view -h -q 255 {} \\\n\t| '.format(samtools_path, bam) + 
            'awk \'{if ($3 != "chrM") {print} ' + 
            'else if (substr($1,1,1) == "@") {print}}\' \\\n\t| ' + 
            '{} view -Su - \\\n\t| '.format(samtools_path) + 
            '{} intersect -v -ubam -abam stdin -b {} \\\n\t| '.format(
                bedtools_path, encode_blacklist) + 
            '{} view -h - \\\n\t| '.format(samtools_path) + 
            'awk \'{if (substr($1,1,1) == "@") {print} ' + 
            'else if ($1 == c1) {print prev; print}  prev=$0; c1=$1}\' ' +
            '\\\n\t| {} view -Sb - \\\n\t> {}\n\n'.format(samtools_path, out_bam)
        )
        with open(self.filename, "a") as f:
            f.write(lines)
        return out_bam

    def count_unique_mt_reads(
        self,
        bam, 
        samtools_path='samtools',
        bg=False,
    ):
        """
        Count the number of uniquely mapped reads mapped to chrM.
    
        Parameters
        ----------
        bam : str
            Path to bam file to count chrM reads for.

        bg : boolean
            If true, run in the background using &.
    
        Returns
        -------
        out : str
            Path to output file with chrM counts.
    
        """
        out = os.path.join(
            self.tempdir,
            '{}_chrM_unique_read_counts.txt'.format(self.sample_name))
        lines = ('{} view -q 255 {} \\\n\t| cut -f 3 \\\n\t| grep chrM \\\n\t| '
                 'uniq -c \\\n\t> {}'.format(samtools_path, bam, out))
        if bg:
            lines += ' &'
        lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return out
    
    def bigwig_from_bedgraph(
        self,
        bedgraph,
        scale=False,
        tlen_max=None,
        web_available=True,
        write_to_outdir=False,
        bedGraphToBigWig_path='bedGraphToBigWig',
        bedtools_path='bedtools',
    ):
        """
        Make bigwig coverage file from bam file.
    
        Parameters
        ----------
        bedgraph : str
            Path to bedgraph file to create bigwig for.
        
        scale : bool
            If True, add note to trackline that data is scaled.

        tlen_max : int
            Add info about tlen filtering.

        web_available : bool
            If True, write trackline to self.links_tracklines, make softlink to
            self.linkdir, and set write_to_outdir = True.

        write_to_outdir : bool
            If True, write output files directly to self.outdir.
    
        bedGraphToBigWig_path : str
            Path to bedGraphToBigWig executable.

        bedtools_path : str
            Path to bedtools. If bedtools_path == 'bedtools', it is assumed that
            the hg19 human.hg19.genome file from bedtools is also in your path.

        Returns
        -------
        bigwig : str
            Path to output bigwig file.
    
        """
        if write_to_outdir or web_available:
            dy = self.outdir
        else:
            dy = self.tempdir
        if bedtools_path == 'bedtools':
            genome_file = '$(which human.hg19.genome)'
        else:
            genome_file = os.path.join(
                os.path.split(os.path.split(bedtools_path)[0])[0], 'genomes',
                'human.hg19.genome')
        root = os.path.splitext(os.path.split(bedgraph)[1])[0]
        bigwig = os.path.join(dy, '{}.bw'.format(root))
        lines = '{} \\\n\t{} \\\n\t{} \\\n\t{}\n\n'.format(
            bedGraphToBigWig_path, bedgraph, genome_file, bigwig)
        with open(self.filename, "a") as f:
            f.write(lines)
        if web_available:
            name = '{}_atac'.format(self.sample_name)
            desc = 'ATACseq coverage for {}.'.format(self.sample_name)
            if tlen_max:
                name += '_140'
                desc += ' Fragments <= 140.'
            if scale:
                name += '_scaled'
                desc = desc.replace('ATACseq', 'Scaled ATACseq')

            url = self.webpath + '/' + os.path.split(bigwig)[1]
            t_lines = (
                'track type=bigWig name="{}" '
                'description="{}" '
                'visibility=0 db=hg19 bigDataUrl={}\n'.format(
                    name, desc, url))
            with open(self.links_tracklines, "a") as f:
                f.write(t_lines)
            link = self.add_softlink(bigwig)

        return bigwig
     
    def _add_macs2_trackline(
        self,
        bed, 
    ):
        """
        Add trackline to macs2 bed file. 
    
        Parameters
        ----------
        bed : str
            Full path to macs2 bed file.
    
        Returns
        -------
        bed : str
            Path to output bed file with trackline added. The original bed file
            is overwritten, so this should be the same as the input bed file
            path.
    
        """
        path = os.path.split(bed)[0]
        track_type = os.path.splitext(bed)[1][1:]
        if track_type == 'bed':
            name = 'summits'
        if track_type == 'narrowPeak':
            name = 'narrowPeak'
        if track_type == 'broadPeak':
            name = 'broadPeak'
        if track_type == 'gappedPeak':
            name = 'gappedPeak'
        
        temp_bed = os.path.join(self.tempdir, 'temp.{}'.format(name))
        track_line = (
            'track type={} name=\\"{}_macs2_atac_{}\\" '
            'description=\\"macs2 ATAC-seq {} for {}\\" '
            'visibility=0 db=hg19'.format(
                track_type, self.sample_name, name, name, self.sample_name)
        )
        lines = 'cat <(echo "{}") {} > {}\n'.format(track_line, bed, temp_bed)
        lines += 'mv {} {}\n\n'.format(temp_bed, bed)
        with open(self.filename, "a") as f:
            f.write(lines)
        return bed
    
    def macs2(
        self,
        bam,
        write_to_outdir=True,
        web_available=True, 
        broad=False,
    ):
        """
        Call peaks with MACS2 for ATAC-seq data. The macs2 executable is assumed
        to be in your path which it should be if you installed it using pip
        install MACS2 into your Python environment.
    
        Parameters
        ----------
        bam : str
            Path to paired-end, coordinate sorted bam file.
    
        write_to_outdir : bool
            If True, write output files directly to self.outdir.
    
        web_available : bool
            If True, write trackline to self.links_tracklines, make softlink to
            self.linkdir, and set write_to_outdir = True.

        Returns
        -------
        excel : str
            Path to output excel file.
    
        out1 : str
            Path to output bed file. If broad == True, this is the broadPeak
            file. Otherwise this is the narrowPeak file.
    
        out2 : str
            Path to output bed file. If broad == True, this is the gappedPeak
            file. Otherwise this is the summits file.
    
        """
        if web_available:
            dy = self.outdir
        else:
            dy = self.tempdir
        # Run macs2.
        run_type = '--call-summits'
        if broad:
            run_type = '--broad'
        lines = ('macs2 callpeak --nomodel --nolambda --keep-dup all \\\n'
                 '\t{} -f BAMPE -g hs \\\n'
                 '\t-t {} \\\n'
                 '\t-n {} \\\n'
                 '\t--outdir {}\n\n'.format(
                     run_type, bam, self.sample_name, dy))
        lines += 'macs2 --version\n\n'
        
        # Output files.
        excel = os.path.join(dy, '{}_peaks.xls')
        if not broad:
            out1 = os.path.join(
                dy, '{}_peaks.narrowPeak'.format(self.sample_name))
            out2 = os.path.join(
                dy, '{}_summits.bed'.format(self.sample_name))
        elif broad:
            out1 = os.path.join(
                dy, '{}_peaks.broadPeak'.format(self.sample_name))
            out2 = os.path.join(
                dy, '{}_peaks.gappedPeak'.format(self.sample_name))

        with open(self.filename, "a") as f:
            f.write(lines)

        out1 = self._add_macs2_trackline(out1)
        out2 = self._add_macs2_trackline(out2)

        # If web_available, make softlinks to bed files and write URLs to file.
        if web_available:
            tf_lines = ''
            link = self.add_softlink(out1)
            url = self.webpath + '/' + os.path.split(out1)[1]
            tf_lines += url + '\n'
            link = self.add_softlink(out2)
            url = self.webpath + '/' + os.path.split(out2)[1]
            tf_lines += url + '\n'

            with open(self.links_tracklines, "a") as f:
                f.write(tf_lines)
        return excel, out1, out2

def pipeline(
    r1_fastqs, 
    r2_fastqs, 
    outdir, 
    sample_name, 
    star_index,
    encode_blacklist,
    find_intersecting_snps_path=None,
    filter_remapped_reads_path=None,
    gatk_fasta=None,
    linkdir=None,
    webpath_file=None,
    vcf=None,
    vcf_sample_name=None,
    is_phased=False,
    conda_env=None,
    modules=None,
    queue=None,
    star_genome_load='LoadAndRemove',
    rgpl='ILLUMINA',
    rgpu='',
    tempdir=None,
    mappability=None,
    expected_unique_pairs=20000000,
    star_path='STAR',
    picard_path='$picard',
    bedtools_path='bedtools',
    bedGraphToBigWig_path='bedGraphToBigWig',
    fastqc_path='fastqc',
    samtools_path='samtools',
    sambamba_path='sambamba',
    gatk_path='$GATK',
    bigWigAverageOverBed_path='bigWigAverageOverBed',
    bcftools_path='bcftools',
    bammarkduplicates_path='bammarkduplicates',
    featureCounts_path='featureCounts',
):
    """
    Make a SGE/shell scripts for running the entire ATAC-seq pipeline. The
    defaults are set for use on the Frazer lab's SGE scheduler on flh1/flh2.

    Parameters
    ----------
    r1_fastqs : list or str
        Either a list of paths to gzipped fastq files with R1 reads or path to a
        single gzipped fastq file with R1 reads.

    r2_fastqs : list or str
        Either a list of paths to gzipped fastq files with R2 reads or path to a
        single gzipped fastq file with R2 reads.

    outdir : str
        Directory to store shell scripts, stdout/stderr logs, and output files
        and directories.

    sample_name : str
        Sample name used for naming files etc.

    star_index : str
        Path to STAR index.

    encode_blacklist : str
        Path to ENCODE blacklist bed file.

    find_intersecting_snps_path : str
        Path to find_intersecting_snps.py from WASP.
    
    filter_remapped_reads_path : str
        Path to filter_remapped_reads.py from WASP.

    gatk_fasta : str
        Fasta file that corresponds to the fasta file used for STAR but is
        sorted karyotypically for GATK. Assumed to have associated dict and fai
        files. This is needed for testing allelic bias.

    linkdir : str
        Path to directory where softlinks should be made. Some pipeline parts
        may make softlinks output files here for display on the web.

    webpath_file : str
        File whose first line is the URL that points to linkdir. For example,
        if we make a link to the file s1_coord_sorted.bam in linkdir and
        webpath_file has http://site.com/files on its first line, then
        http://site.com/files/s1_coord_sorted.bam should be available on the
        web. If the web directory is password protected (it probably should be),
        then the URL should look like http://username:password@site.com/files.
        This is a file so you don't have to make the username/password combo
        public (although I'd recommend not using a sensitive password). You can
        just put the webpath_file in a directory that isn't tracked by git, 
        figshare, etc.

    vcf : str
        VCF file containing exonic variants used for allelic bias.
    
    vcf_sample_name : str
        Sample name of this sample in the VCF file (if different than
        sample_name). For instance, the sample name in the VCF file may be the
        sample name for WGS data which may differ from the ATAC-seq sample name.

    conda_env : str
        Conda environment to load at the beginning of the script.

    modules : str
        Comma-separated list of modules to load at the beginning of the script.

    rgpl : str
        Read Group platform (e.g. illumina, solid). 

    rgpu : str
        Read Group platform unit (eg. run barcode). 

    tempdir : str
        Directory to store temporary files.

    expected_unique_pairs : int
        Number of expected mapped READ PAIRS. This only needs to be a rough 
        estimate. You should keep this number constant for all samples across 
        an experiment/project (i.e. any samples you want to compare against 
        each other). For example, say you are doing RNA-seq for 10 samples on 
        a lane, the lane yields 100,000,000 pairs of reads, and you expect 
        90%% of the read pairs to align uniquely. Then you would want 
        expected_num = 9,000,000. In this case, you would set actual_num to 
        the actual number of uniquely mapped reads. The reason you want this 
        number close to the number of expected number of uniquely mapped 
        pairs is so the normalization does not drastically change the 
        coverage for samples near the expected number of input read pairs. 
        For instance, if expected_num=20M and a sample has exactly 20M 
        uniquely mapped read pairs, then the coverage will not be changed. If 
        expected_num=20M and a sample has only 10M input read PAIRS, it will 
        be normalized so that it looks like it had 20M input read pairs (i.e. 
        all coverages will be multiplied by 2).

    star_path : str
        Path to STAR aligner.

    picard_path : str
        Path to Picard tools.

    bedtools_path : str
        Path to bedtools.

    bedGraphToBigWig_path : str
        Path bedGraphToBigWig executable.

    Returns
    -------
    fn : str
        Path to submission shell script.

    """
    with open(webpath_file) as wpf:
        webpath = wpf.readline().strip()

    # Bash commands to submit jobs. I'll collect these as I make the jobs and
    # then write them to a file at the end.
    submit_commands = []
    
    ##### Job 1: Combine fastqs and align with STAR. #####
    job = ATACJobScript(
        sample_name, 
        job_suffix='alignment',
        outdir=os.path.join(outdir, 'alignment'), 
        threads=8, 
        memory=32,
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, 
        queue=queue, 
        conda_env=conda_env,
        modules=modules,
    )
    alignment_jobname = job.jobname
    
    # Input files.
    for fq in r1_fastqs + r2_fastqs:
        job.add_input_file(fq)

    # Combine R1 and R2 fastqs.
    if type(r1_fastqs) == str:
        r1_fastqs = [r1_fastqs]
    if type(r2_fastqs) == str:
        r2_fastqs = [r2_fastqs]
    r1_fastqs = [os.path.realpath(x) for x in r1_fastqs]
    r2_fastqs = [os.path.realpath(x) for x in r2_fastqs]
    combined_r1 = job.combine_fastqs(r1_fastqs, suffix='R1', bg=True)
    combined_r2 = job.combine_fastqs(r2_fastqs, suffix='R2', bg=True)
    with open(job.filename, "a") as f:
            f.write('\nwait\n\n')
    # We don't want to keep the fastqs indefinitely, but we need them for the
    # fastQC step later.
    combined_r1 = job.add_output_file(combined_r1)
    combined_r2 = job.add_output_file(combined_r2)

    # Align reads.
    (star_bam, log_out, log_final_out, log_progress_out, sj_out) = \
            job.star_align(combined_r1, combined_r2, rgpl, rgpu, star_index,
                            job.threads, genome_load=star_genome_load)
    outdir_star_bam = job.add_output_file(star_bam)
    log_final_out = job.add_output_file(log_final_out)
    [job.add_output_file(x) for x in [log_out, log_progress_out, sj_out]]
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    ##### Job 2: Run fastQC. ##### 
    job = ATACJobScript(
        sample_name, 
        job_suffix='fastqc', 
        outdir=os.path.join(outdir, 'qc'), 
        threads=1, 
        memory=4,
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, 
        queue=queue, 
        conda_env=conda_env,
        modules=modules, 
        wait_for=[alignment_jobname]
    )
    fastqc_jobname = job.jobname
 
    # Input files.
    job.add_input_file(combined_r1, delete_original=True)
    job.add_input_file(combined_r2, delete_original=True)

    # Run fastQC.
    fastqc_html, fastqc_zip = job.fastqc([combined_r1, combined_r2],
                                         fastqc_path)
    [job.add_output_file(x) for x in fastqc_html + fastqc_zip]
        
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    ##### Job 3: Filter, coordinate sort, mark duplicates and index bam. #####
    job = ATACJobScript(
        sample_name, 
        job_suffix = 'sort_rmdup_index',
        outdir=os.path.join(outdir, 'alignment'), 
        threads=4, 
        memory=6,
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, 
        queue=queue, 
        conda_env=conda_env,
        modules=modules,
        wait_for=[alignment_jobname]
    )
    sort_rmdup_index_jobname = job.jobname

    # Input files.
    star_bam = job.add_input_file(outdir_star_bam)

    # Count uniquely mapped mitochondrial reads.
    mito_counts = job.count_unique_mt_reads(
        star_bam, 
        samtools_path=samtools_path,
        bg=True,
    )
    job.add_output_file(mito_counts)

    # Filter out non-uniquely mapped and mitochondrial reads.
    filtered_bam = job.filter_multi_mt_blacklist_reads(
        bam=star_bam,
        encode_blacklist=encode_blacklist,
        bedtools_path=bedtools_path,
        samtools_path=samtools_path,
    )
    job.add_temp_file(filtered_bam)
    
    # Coordinate sort.
    coord_sorted_bam = job.sambamba_sort(
        filtered_bam, 
        tempdir=job.tempdir,
        sambamba_path=sambamba_path,
    )
    job.add_temp_file(coord_sorted_bam)

    # Index sorted bam file.
    temp_bam_index = job.sambamba_index(coord_sorted_bam, sambamba_path)
    job.add_temp_file(temp_bam_index)

    # Mark and remove duplicates.
    rmdup_bam, duplicate_metrics, removed_reads = job.biobambam2_mark_duplicates(
        coord_sorted_bam,
        remove_duplicates=True,
        bammarkduplicates_path=bammarkduplicates_path)
    outdir_rmdup_bam = job.add_output_file(rmdup_bam)
    job.add_output_file(duplicate_metrics)
    job.add_output_file(removed_reads)

    # Index rmdup bam file.
    bam_index = job.sambamba_index(rmdup_bam, sambamba_path)
    outdir_rmdup_bam_index = job.add_output_file(bam_index)

    # Add softlink to bam file in outdir and write URL and trackline.
    link = job.add_softlink(outdir_rmdup_bam)
    link = job.add_softlink(outdir_rmdup_bam_index)
    name = '{}_atac'.format(job.sample_name)
    desc = 'ATAC-seq alignment for {}.'.format(job.sample_name)
    url = job.webpath + '/' + os.path.split(outdir_rmdup_bam)[1]
    t_lines = (
        'track type=bam name="{}" '
        'description="{}" '
        'visibility=0 db=hg19 bigDataUrl={}\n'.format(
            name, desc, url))
    with open(job.links_tracklines, "a") as f:
        f.write(t_lines)
        f.write(url + '\n')

    # Filter out reads with fragment size greater than 140.
    tlen_bam = job.filter_tlen(
        rmdup_bam,
        tlen_max=140,
        samtools_path=samtools_path,
    )
    outdir_tlen_bam = job.add_output_file(tlen_bam)

    # Index sorted bam file.
    tlen_bam_index = job.sambamba_index(tlen_bam, sambamba_path)
    outdir_tlen_bam_index = job.add_output_file(tlen_bam_index)

    # Query name sort.
    query_sorted_bam = job.sambamba_sort(
        tlen_bam, 
        tempdir=job.tempdir,
        queryname=True,
        root=os.path.splitext(os.path.split(tlen_bam)[1])[0],
        sambamba_path=sambamba_path,
    )
    query_sorted_bam = job.add_output_file(query_sorted_bam)

    # Add softlink to bam file in outdir and write URL and trackline.
    link = job.add_softlink(outdir_tlen_bam)
    link = job.add_softlink(outdir_tlen_bam_index)
    name = '{}_atac_140'.format(job.sample_name)
    desc = 'ATAC-seq alignment for {}. Fragments <= 140.'.format(job.sample_name)
    url = job.webpath + '/' + os.path.split(outdir_tlen_bam)[1]
    t_lines = (
        'track type=bam name="{}" '
        'description="{}" '
        'visibility=0 db=hg19 bigDataUrl={}\n'.format(
            name, desc, url))
    with open(job.links_tracklines, "a") as f:
        f.write(t_lines)
        f.write(url + '\n')

    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
    
    ##### Job 4: Collect QC metrics. #####
    job = ATACJobScript(
        sample_name, 
        job_suffix='qc_metrics',
        outdir=os.path.join(outdir, 'qc'),
        threads=1, 
        memory=5, 
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, 
        queue=queue,
        conda_env=conda_env, 
        modules=modules,
        wait_for=[sort_rmdup_index_jobname],
    )
    qc_metrics_jobname = job.jobname
    
    # Input files.
    rmdup_bam = job.add_input_file(outdir_rmdup_bam)

    # Collect Picard insert size metrics.
    insert_metrics, insert_hist = job.picard_insert_size_metrics(
        rmdup_bam, 
        picard_path=picard_path,
        bg=False,
    )

    # Collect index stats.
    index_out, index_err = job.picard_bam_index_stats(
        rmdup_bam, 
        picard_path=picard_path,
        bg=False,
    )
    job.add_output_file(index_out)
    job.add_output_file(index_err)

    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    ##### Job 5: Make md5 hash for STAR bam file. #####
    job = ATACJobScript(
        sample_name, 
        job_suffix='md5',
        outdir=os.path.join(outdir, 'alignment'), 
        threads=1, 
        memory=1,
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, 
        queue=queue, 
        conda_env=conda_env,
        modules=modules,
        wait_for=[alignment_jobname],
    )
    md5_jobname = job.jobname
    
    # Input files.
    star_bam = job.add_input_file(outdir_star_bam)

    # Make md5 hashes for bam files.
    star_md5sum = job.make_md5sum(star_bam)
    job.add_output_file(star_md5sum)

    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
      
    ##### Job 6: Make bigwig files for final bam file. #####
    job = ATACJobScript(
        sample_name, 
        job_suffix='bigwig',
        outdir=os.path.join(outdir, 'alignment'), 
        threads=1, 
        memory=4,
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, 
        queue=queue, 
        conda_env=conda_env,
        modules=modules,
        wait_for=[sort_rmdup_index_jobname],
    )
    bigwig_jobname = job.jobname
        
    # Input files.
    rmdup_bam = job.add_input_file(outdir_rmdup_bam)
    tlen_bam = job.add_input_file(outdir_tlen_bam)

    # Make bigwig.
    bg = job.bedgraph_from_bam(
        rmdup_bam, 
        bedtools_path=bedtools_path,
    )
    job.add_temp_file(bg)
    bw = job.bigwig_from_bedgraph(
        bg,
        bedGraphToBigWig_path=bedGraphToBigWig_path,
        bedtools_path=bedtools_path,
    )
    job.add_output_file(bw)

    # Make scaled bigwig file. 
    scaled_bg = job.scale_bedgraph(
        bg,
        log_final_out,
        expected_unique_pairs,
    )
    job.add_temp_file(scaled_bg)
    scaled_bw = job.bigwig_from_bedgraph(
        scaled_bg,
        scale=True,
        bedGraphToBigWig_path=bedGraphToBigWig_path,
        bedtools_path=bedtools_path,
    )
    job.add_output_file(scaled_bw)

    # Make bigwig for tlen filtered bam.
    tlen_bg = job.bedgraph_from_bam(
        tlen_bam, 
        bedtools_path=bedtools_path,
    )
    # tlen_bg will actually overwrite bg.
    # job.add_temp_file(tlen_bg)
    tlen_bw = job.bigwig_from_bedgraph(
        tlen_bg,
        tlen_max=140,
        bedGraphToBigWig_path=bedGraphToBigWig_path,
        bedtools_path=bedtools_path,
    )
    job.add_output_file(tlen_bw)

    # Make scaled bigwig file for tlen filtered bam.
    tlen_scaled_bg = job.scale_bedgraph(
        tlen_bg,
        log_final_out,
        expected_unique_pairs,
    )
    # tlen_scaled_bg will actually overwrite scaled_bg
    # job.add_temp_file(tlen_scaled_bg)
    tlen_scaled_bw = job.bigwig_from_bedgraph(
        tlen_scaled_bg,
        scale=True,
        tlen_max=140,
        bedGraphToBigWig_path=bedGraphToBigWig_path,
        bedtools_path=bedtools_path,
    )
    job.add_output_file(tlen_scaled_bw)

    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
    
    ##### Job 7: Peak calling. #####
    job = ATACJobScript(
        sample_name, 
        job_suffix = 'macs2',
        outdir=os.path.join(outdir, 'macs2'), 
        threads=1, 
        memory=4,
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, 
        queue=queue, 
        conda_env=conda_env,
        modules=modules,
        wait_for=[sort_rmdup_index_jobname],
    )
    macs2_jobname = job.jobname
    
    # Input files.
    tlen_bam = job.add_input_file(outdir_tlen_bam)

    # Run macs2 peak calling.
    excel, narrow_peak, summits = job.macs2(
        tlen_bam,
        web_available=True, 
        broad=False,
    )
    job.add_output_file(excel)
    outdir_narrow_peak = job.add_output_file(narrow_peak)
    job.add_output_file(summits)

    # Motif analysis.
    homer_outdir = job.homer_motif_analysis(
        summits,
        size=200,
        mask=True,
        web_available=True,
        write_to_outdir=True,
    )

    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
    
    ##### Job 8: Count reads in peaks. #####
    job = ATACJobScript(
        sample_name, 
        job_suffix = 'counts',
        outdir=os.path.join(outdir, 'counts'),
        threads=4, 
        memory=8, 
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, queue=queue,
        conda_env=conda_env, 
        modules=modules,
        wait_for=[macs2_jobname],
    )
    counts_jobname = job.jobname
    
    # Input files.
    query_sorted_bam = job.add_input_file(query_sorted_bam)
    narrow_peak = job.add_input_file(outdir_narrow_peak)

    # Run featureCounts.
    counts, counts_summary = job.featureCounts_count(
        narrow_peak,
        query_sorted_bam,
        sort=False,
        filetype='bed',
        featureCounts_path=featureCounts_path,
    )
    job.add_output_file(counts)
    job.add_output_file(counts_summary)

    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
   
    # We'll only go through the allelic bias steps if a VCF was provided.
    if vcf:
        ##### Job 9: WASP first step. #####
        job = ATACJobScript(
            sample_name, 
            job_suffix = 'wasp_allele_swap',
            outdir=os.path.join(outdir, 'wasp'),
            threads=1, 
            memory=4, 
            linkdir=linkdir,
            webpath=webpath,
            tempdir=tempdir, 
            queue=queue,
            conda_env=conda_env, 
            modules=modules,
            wait_for=[sort_rmdup_index_jobname],
        )
        wasp_allele_swap_jobname = job.jobname
           
        # Input files.
        rmdup_bam = job.add_input_file(rmdup_bam)
        # The VCF might be large so we probably don't want to copy it ever.
        vcf = job.add_input_file(vcf, copy=False)
        # The exon bed file is small so we don't need to copy it ever.
        exon_bed = job.add_input_file(exon_bed, copy=False)

        # Run WASP allele swap.
        if not vcf_sample_name:
            vcf_sample_name = sample_name
        (snp_directory, hets_vcf, keep_bam, wasp_r1_fastq, wasp_r2_fastq,
         to_remap_bam, to_remap_num) = job.wasp_allele_swap(
             rmdup_bam, 
             find_intersecting_snps_path, 
             vcf, 
             exon_bed,
             gatk_fai=gatk_fasta + '.fai',
             vcf_sample_name=vcf_sample_name, 
             samtools_path=samtools_path,
             bcftools_path=bcftools_path,
        )
        # WASP outputs a file (keep_bam) that has reads that don't overlap
        # variants. I'm going to discard that file. I'll discard the WASP SNP
        # directory as well since I have the hets in a VCF file.
        job.add_temp_file(keep_bam)
        snp_directory = job.add_temp_file(snp_directory)
        hets_vcf = job.add_output_file(hets_vcf)
        wasp_r1_fastq = job.add_output_file(wasp_r1_fastq)
        wasp_r2_fastq = job.add_output_file(wasp_r2_fastq)
        to_remap_bam = job.add_output_file(to_remap_bam)
        to_remap_num = job.add_output_file(to_remap_num)

        job.write_end()
        if not job.delete_sh:
            submit_commands.append(job.sge_submit_command())
        
        ##### Job 10: WASP second step. #####
        job = ATACJobScript(
            sample_name, 
            job_suffix='wasp_remap',
            outdir=os.path.join(outdir, 'wasp'),
            threads=8, 
            memory=32, 
            linkdir=linkdir,
            webpath=webpath,
            tempdir=tempdir,
            queue=queue, 
            conda_env=conda_env, 
            modules=modules,
            wait_for=[wasp_allele_swap_jobname],
        )
        wasp_remap_jobname = job.jobname
        
        # Input files.
        wasp_r1_fastq = job.add_input_file(wasp_r1_fastq, delete_original=True)
        wasp_r2_fastq = job.add_input_file(wasp_r2_fastq, delete_original=True)

        # Realign allele-swapped fastqs.
        remapped_bam, log_out, log_final_out, log_progress_out, sj_out = \
                job.star_align(wasp_r1_fastq, wasp_r2_fastq, rgpl, rgpu,
                               star_index, job.threads,
                               genome_load=star_genome_load,
                               transcriptome_align=False)
        job.add_output_file(remapped_bam)
        job.add_output_file(log_out)
        job.add_output_file(log_final_out)
        job.add_output_file(log_progress_out)
        job.add_temp_file(sj_out)

        job.write_end()
        if not job.delete_sh:
            submit_commands.append(job.sge_submit_command())
        
        ##### Job 11: WASP third step. #####
        job = ATACJobScript(
            sample_name, 
            job_suffix = 'wasp_alignment_compare',
            outdir=os.path.join(outdir, 'wasp'),
            threads=1, 
            memory=5, 
            linkdir=linkdir,
            webpath=webpath,
            tempdir=tempdir, 
            queue=queue,
            conda_env=conda_env, 
            modules=modules,
            wait_for=[wasp_remap_jobname],
        )
        wasp_alignment_compare_jobname = job.jobname

        # Input files.
        to_remap_bam = job.add_input_file(to_remap_bam, delete_original=True)
        to_remap_num = job.add_input_file(to_remap_num, delete_original=True)
        remapped_bam = job.add_input_file(remapped_bam, delete_original=True)

        # Compare alignments.
        temp_filtered_bam = job.wasp_alignment_compare(
            to_remap_bam, 
            to_remap_num,
            remapped_bam, 
            filter_remapped_reads_path,
        )
        job.add_temp_file(temp_filtered_bam)
            
        # Coordinate sort and index filtered bam file.
        wasp_filtered_bam, wasp_bam_index = job.picard_coord_sort(
            temp_filtered_bam, 
            index=True,
            picard_path=picard_path,
        )
        # I'll keep this bam file and its index as a record of which reads were
        # used to calculate allelic bias. It might be useful for visualization.
        # The bam file is pretty small anyway.
        job.add_output_file(wasp_filtered_bam)
        job.add_output_file(wasp_bam_index)

        # Reorder bam file so it will work with GATK.
        reordered_bam = job.picard_reorder(
            wasp_filtered_bam, 
            fasta=gatk_fasta,
            picard_path=picard_path,
        )
        job.add_temp_file(reordered_bam)

        # Index reordered bam.
        reordered_index = job.picard_index(
            reordered_bam, 
            picard_path=picard_path,
            bg=False,
        )
        job.add_temp_file(reordered_index)

        # Get allele counts.
        allele_counts = job.count_allele_coverage(
            reordered_bam, 
            hets_vcf,
            gatk_fasta, 
            gatk_path=gatk_path,
        )
        job.add_output_file(allele_counts)

        job.write_end()
        if not job.delete_sh:
            submit_commands.append(job.sge_submit_command())
        
        ##### Job 12: Run MBASED for allelic bias. #####
        job = ATACJobScript(
            sample_name,
            job_suffix='mbased',
            outdir=os.path.join(outdir, 'mbased'),
            threads=8, 
            memory=16, 
            linkdir=linkdir,
            webpath=webpath,
            tempdir=tempdir, 
            queue=queue,
            conda_env=conda_env, 
            modules=modules,
            wait_for=[wasp_alignment_compare_jobname],
        )
        mbased_jobname = job.jobname
    
        # Input files.
        allele_counts = job.add_input_file(allele_counts)

        mbased_infile, locus_outfile, snv_outfile = job.mbased(
            allele_counts, 
            gene_bed, 
            is_phased=is_phased, 
            num_sim=1000000, 
            vcf=vcf,
            vcf_sample_name=vcf_sample_name, 
            mappability=mappability,
            bigWigAverageOverBed_path=bigWigAverageOverBed_path,
        )
        job.add_output_file(mbased_infile)
        job.add_output_file(locus_outfile)
        job.add_output_file(snv_outfile)

        job.write_end()
        if not job.delete_sh:
            submit_commands.append(job.sge_submit_command())

    ##### Submission script #####
    # Now we'll make a submission script that submits the jobs with the
    # appropriate dependencies.
    import datetime as dt
    now = str(dt.datetime.now())
    now = now.replace('-', '_').replace(' ', '_').replace(':', '_').replace('.', '_')
    submit_fn = os.path.join(outdir, 'sh', '{}_submit_{}.sh'.format(
        sample_name, now))
    with open(submit_fn, 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('\n'.join(submit_commands))

    return submit_fn


# import os
# 
# from general import _bedgraph_to_bigwig
# from general import _bigwig_files
# from general import _combine_fastqs
# from general import _cutadapt_trim
# from general import _fastqc
# from general import JobScript
# from general import _flagstat
# from general import _picard_insert_size_metrics
# from general import _make_softlink
# # from general import _pbs_header
# from general import _picard_coord_sort
# from general import _picard_query_sort
# from general import _picard_mark_duplicates
# from general import _samtools_index
# 
# def _star_align(r1_fastqs, r2_fastqs, sample, rgpl, rgpu, star_index, star_path,
#                 threads):
#     """
#     Align paired ATAC-seq fastq files with STAR.
# 
#     Parameters
#     ----------
#     r1_fastqs : str
#         Gzipped R1 fastq file(s). If multiple files, each file should be
#         separated by s space and should be ordered the same as the R2 files.
# 
#     r2_fastqs : str
#         Gzipped R2 fastq file(s). If multiple files, each file should be
#         separated by s space and should be ordered the same as the R1 files.
# 
#     sample : str
#         Sample name.
# 
#     rgpl : str
#         Read Group platform (e.g. illumina, solid). 
# 
#     rgpu : str
#         Read Group platform unit (eg. run barcode). 
# 
#     """
#     # I use threads - 2 for STAR so there are open processors for reading and
#     # writing. The parameters here are set in part based on the discussion here:
#     # https://groups.google.com/forum/#!searchin/rna-star/chip$20seq/rna-star/E_mKqm9jDm0/ZpB6yRcWi60J.
#     # --alignIntronMax 1 disables spliced alignments. --alignEndsType EndToEnd
#     # prohibits soft clipping. I may or may not want to prohibit soft clipping.
#     # Smaller-sized fragments may still have adapter sequence and won't align
#     # without soft clipping. I think that since our reads are long, we should be
#     # ok soft clipping since they should align pretty well. I'm only outputting 
#     # unique alignments right now.
#     line = (' \\\n'.join([star_path, 
#                           '\t--runThreadN {}'.format(threads - 2),
#                           '\t--alignIntronMax 1',
#                           '\t--genomeDir {}'.format(star_index), 
#                           '\t--genomeLoad LoadAndRemove', 
#                           '\t--readFilesCommand zcat',
#                           '\t--readFilesIn {} {}'.format(r1_fastqs, 
#                                                          r2_fastqs),
#                           '\t--outSAMtype BAM Unsorted', 
#                           '\t--outSAMattributes All', 
#                           '\t--outSAMunmapped Within',
#                           ('\t--outSAMattrRGline ID:1 ' + 
#                            'PL:{} '.format(rgpl) + 
#                            'PU:{} '.format(rgpu) + 
#                            'LB:{0} SM:{0}'.format(sample)), 
#                           '\t--outFilterMultimapNmax 20', 
#                           '\t--outFilterMismatchNmax 999',
#                           '\t--outFilterMismatchNoverLmax 0.04',
#                           '\t--seedSearchStartLmax 20']) + '\n\n')
#     return line
# 
# def _homer(bam, sample_name, temp_tagdir, final_tagdir, homer_path, link_dir,
#            bedtools_path, bigwig=False):
#     """
#     Make tag directory and call peaks with HOMER. Optionally make bigwig file.
# 
#     Parameters
#     ----------
#     bam : str
#         Path to paired-end, coordinate sorted bam file.
# 
#     temp_tagdir : str
#         Path to temporary tag directory that HOMER will create.
# 
#     final_tagdir : str
#         Path to final tag directory (what the temp_tagdir will be copied to).
# 
#     link_dir : str
#         Path to directory where softlinks should be made. HOMER will put the
#         bigwig file in the directory specified when setting up HOMER but a
#         softlink to a bed file with the HOMER peaks will be made here.
# 
#     bigwig : bool
#         If True, have HOMER make a bigwig file.
# 
#     Returns
#     -------
#     lines : str
#         Lines to be printed to shell/PBS script.
# 
#     """
#     temp_tagdir = os.path.realpath(temp_tagdir)
#     final_tagdir = os.path.realpath(final_tagdir)
#     name = '{}_atac_homer'.format(sample_name)
#     lines = []
#     lines.append('{}/makeTagDirectory {} {}'.format(homer_path, temp_tagdir,
#                                                     bam))
#     if bigwig:
#         lines.append('{}/makeBigWig.pl {} hg19 -name {}'
#                      ' -url www.fake.com/ -webdir {}'.format(
#                          homer_path, os.path.split(temp_tagdir)[1], name, 
#                          os.path.split(temp_tagdir)[0]))
#         lines.append('mv {}/{}_tags.ucsc.bigWig {}'.format(
#             os.path.split(temp_tagdir)[0], sample_name, temp_tagdir))
#     lines.append('{}/findPeaks {} -style histone -size 75 -minDist 75 '
#                  '-o auto'.format(homer_path, temp_tagdir))
#     # softlink_lines = []
#     posfile = os.path.join(temp_tagdir, 'regions.txt')
#     bed = os.path.join(temp_tagdir, '{}_peaks.bed'.format(name))
#     lines.append(_convert_homer_pos_to_bed(
#         posfile, bed, name, homer_path, bedtools_path))
#     bed = os.path.join(final_tagdir, '{}_peaks.bed'.format(name))
#     # softlink, name = _make_softlink(bed, name, os.path.join(link_dir, 'atac',
#     #                                                         'peak'))
#     # softlink_lines.append(softlink)
#     
#     lines = '\n'.join(lines) + '\n\n'
#     # softlink_lines = '\n'.join(softlink_lines)
#     # return lines, softlink_lines
#     return lines
# 
# def _combined_homer(input_tagdirs, combined_name, temp_tagdir, final_tagdir,
#                     homer_path, link_dir, bedtools_path, tracklines_file,
#                     web_path_file, bigwig=False):
#     """
#     Combine tag directories and call peaks with HOMER. Optionally make bigwig
#     file.
# 
#     Parameters
#     ----------
#     input_tagdirs : list 
#         A list of paths to the tagdirs to be combined for calling peaks.
# 
#     combined_name : str
#         Used for naming output files and directories.
# 
#     temp_tagdir : str
#         Path to temporary tag directory that HOMER will create.
# 
#     final_tagdir : str
#         Path to final tag directory (what the temp_tagdir will be copied to).
# 
#     link_dir : str
#         Path to directory where softlinks should be made. HOMER will put the
#         bigwig file in the directory specified when setting up HOMER but a
#         softlink to a bed file with the HOMER peaks will be made here.
# 
#     bigwig : bool
#         If True, have HOMER make a bigwig file.
# 
#     Returns
#     -------
#     lines : str
#         Lines to be printed to shell/PBS script.
# 
#     """
#     with open(web_path_file) as wpf:
#         web_path = wpf.readline().strip()
#     web_path = web_path + '/atac'
# 
#     temp_tagdir = os.path.realpath(temp_tagdir)
#     final_tagdir = os.path.realpath(final_tagdir)
#     name = '{}_combined_peak'.format(combined_name)
#     bed = os.path.join(temp_tagdir,
#                        '{}_combined_homer_peaks.bed'.format(combined_name))
#     lines = []
#     tf_lines = []
#     softlink_lines = []
#     lines.append('{}/makeTagDirectory {} -d {}'.format(homer_path, temp_tagdir, 
#                                                       ' '.join(input_tagdirs)))
#     if bigwig:
#         lines.append('{}/makeBigWig.pl {} hg19 -name {}'
#                      ' -url www.fake.com/ -webdir {}'.format(
#                          homer_path, os.path.split(temp_tagdir)[1], name, 
#                          os.path.split(temp_tagdir)[0]))
#         lines.append('mv {0}.ucsc.bigWig {0}'.format(temp_tagdir))
#         temp_link_dir = os.path.join(link_dir, 'bigwig')
#         bw = os.path.join(
#             final_tagdir,
#             '{}.ucsc.bigWig'.format(os.path.split(temp_tagdir)[1]))
#         softlink, name = _make_softlink(bw, name, os.path.join(
#             temp_link_dir, 'atac', 'bigwig'))
#         softlink_lines.append(softlink)
#         temp_web_path = web_path + '/bigwig'
#         try:
#             os.makedirs(temp_link_dir)
#         except OSError:
#             pass
#         tf_lines += ('track type=bigWig name="{0}_atac_cov" '
#                      'description="ATAC-seq '
#                      'coverage for {0}" visibility=0 db=hg19 bigDataUrl='
#                      '{1}/{0}.ucsc.bigWig\n'.format(
#                          os.path.split(temp_tagdir)[1], temp_web_path))
#     lines.append('{}/findPeaks {} -style super -size 75 '
#                  ' -o auto'.format(
#                      homer_path, temp_tagdir))
#     lines.append('{}/findPeaks {} -style histone -size 75 -minDist 75 '
#                  '-o auto'.format(homer_path, temp_tagdir))
#    
#     posfile = os.path.join(temp_tagdir, 'regions.txt')
#     name = '{}_combined'.format(combined_name)
#     bed = os.path.join(temp_tagdir, '{}_homer_atac_peaks.bed'.format(name))
#     lines.append(_convert_homer_pos_to_bed(
#         posfile, bed, name, homer_path, bedtools_path))
#     bed = os.path.join(final_tagdir, '{}_homer_atac_peaks.bed'.format(name))
#     softlink, name = _make_softlink(bed, name, os.path.join(link_dir, 'atac',
#                                                             'peak'))
#     softlink_lines.append(softlink)
#     temp_web_path = web_path + '/peak'
#     tf_lines += '{}/{}_homer_atac_peaks.bed\n'.format(temp_web_path, name)
# 
#     posfile = os.path.join(temp_tagdir, 'superEnhancers.txt')
#     name = '{}_combined_super_enhancers'.format(combined_name)
#     bed = os.path.join(temp_tagdir, '{}_homer_atac_peaks.bed'.format(name))
#     lines.append(_convert_homer_pos_to_bed(
#         posfile, bed, name, homer_path, bedtools_path))
#     bed = os.path.join(final_tagdir, '{}_homer_atac_peaks.bed'.format(name))
#     softlink, name = _make_softlink(bed, name, os.path.join(link_dir, 'atac',
#                                                             'peak'))
#     softlink_lines.append(softlink)
#     tf_lines += '{}/{}_homer_atac_peaks.bed\n'.format(temp_web_path, name)
#     lines = '\n'.join(lines) + '\n\n'
#     softlink_lines = '\n'.join(softlink_lines)
# 
#     # Write tracklines and URLs.
#     if os.path.exists(tracklines_file):
#         with open(tracklines_file) as f:
#             existing_lines = f.read()
#     else:
#         existing_lines = ''
# 
#     with open(tracklines_file, 'w') as tf:
#         tf.write(existing_lines + ''.join(tf_lines))
#     
#     return lines, softlink_lines
# 
# def _convert_homer_pos_to_bed(posfile, bed, sample_name, homer_path,
#                               bedtools_path):
#     """
#     Convert HOMER results file to bed file.
# 
#     Parameters
#     ----------
#     posfile : str
#         Full path to HOMER results file (i.e. regions.txt or
#         superEnhancers.txt).
# 
#     bed : str
#         Name of output bed file.
# 
#     sample_name : str
#         Used for naming output files.
# 
#     Returns
#     -------
#     lines : str
#         Lines to be printed to shell/PBS script.
# 
#     """
#     lines = []
#     tagdir = os.path.split(posfile)[0]
#     lines.append('{}/pos2bed.pl {} | grep -v \# > temp.bed'.format(
#         homer_path, posfile))
#     track_line = ' '.join(['track', 'type=bed',
#                            'name=\\"{}_homer_atac_peaks\\"'.format(
#                                sample_name),
#                            ('description=\\"HOMER ATAC-seq peaks for '
#                             '{}\\"'.format(sample_name)),
#                            'visibility=0',
#                            'db=hg19'])
#     lines.append('{} sort -i temp.bed > temp2.bed'.format(bedtools_path))
#     lines.append('cat <(echo {}) temp2.bed > {}'.format(track_line, bed))
#     lines.append('rm temp.bed temp2.bed')
#     return '\n'.join(lines) + '\n'
# 
# def combined_homer_peaks(
#     tagdirs, 
#     outdir, 
#     combined_name, 
#     tracklines_file,
#     link_dir,
#     web_path_file,
#     bedtools_path,
#     homer_path,
#     environment,
#     conda_env='',
#     tempdir='/scratch', 
#     threads=6, 
#     shell=False,
# ):
#     """
#     Make a PBS or shell script for combining together HOMER tag directories and
#     calling peaks on the combined tags.
# 
#     Parameters
#     ----------
#     tagdirs : list 
#         A list of paths to the tagdirs to be combined for calling peaks.
# 
#     outdir : str
#         Directory to store PBS/shell file and results.
# 
#     combined_name : str
#         Name used for naming files and directories for combined data.
# 
#     tracklines_file : str
#         Path to file for writing tracklines. The tracklines will be added to the
#         file; the contents of the file will not be overwritten. These tracklines
#         can be pasted into the genome browser upload for custom data.
# 
#     link_dir : str
#         Path to directory where softlinks for genome browser should be made.
# 
#     web_path_file : str
#         File whose first line is the URL that points to link_dir. For example,
#         if we make a link to the file s1_coord_sorted.bam in link_dir and
#         web_path_file has http://site.com/files on its first line, then
#         http://site.com/files/s1_coord_sorted.bam should be available on the
#         web. If the web directory is password protected (it probably should be),
#         then the URL should look like http://username:password@site.com/files.
#         This is a file so you don't have to make the username/password combo
#         public (although I'd recommend not using a sensitive password). You can
#         just put the web_path_file in a directory that isn't tracked by git, 
#         figshare, etc.
# 
#     bedtools_path : str
#         Path to bedtools.
# 
#     homer_path : str
#         Path to HOMER bin.
# 
#     environment : str
#         Bash file with PATH information that can be sourced. This should include
#         the paths to executables HOMER will need like bedGraphToBigWig.
# 
#     conda_env : str
#         If provided, load conda environment with this name. This will control
#         which version of MACS2 is used.
# 
#     tempdir : str
#         Directory to store files as STAR runs.
# 
#     threads : int
#         Number of threads to reserve using PBS scheduler. 
# 
#     shell : boolean
#         If true, make a shell script rather than a PBS script.
#     
#     Returns
#     -------
#     fn : str
#         Path to PBS/shell script.
# 
#     """
#     assert threads >= 3
# 
#     if shell:
#         pbs = False
#     else: 
#         pbs = True
# 
#     tempdir = os.path.join(tempdir, '{}_combined_peaks'.format(combined_name))
#     outdir = os.path.join(outdir, '{}_combined_peaks'.format(combined_name))
# 
#     # I'm going to define some file names used later.
#     local_tagdir = '{}_combined_tags'.format(combined_name)
#     temp_tagdir = os.path.join(tempdir, local_tagdir)
#     final_tagdir = os.path.join(outdir, local_tagdir)
#     
#     # Files to copy to output directory.
#     files_to_copy = [temp_tagdir]
#     # Temporary files that can be deleted at the end of the job. We may not want
#     # to delete the temp directory if the temp and output directory are the
#     # same.
#     files_to_remove = [os.path.split(os.path.realpath(x))[1] for x in tagdirs]
# 
#     try:
#         os.makedirs(outdir)
#     except OSError:
#         pass
# 
#     if shell:
#         fn = os.path.join(outdir, '{}_combined_peaks.sh'.format(combined_name))
#     else:
#         fn = os.path.join(outdir, '{}_combined_peaks.pbs'.format(combined_name))
# 
#     f = open(fn, 'w')
#     f.write('#!/bin/bash\n\n')
#     if pbs:
#         out = os.path.join(outdir,
#                            '{}_combined_peaks.out'.format(combined_name))
#         err = os.path.join(outdir,
#                            '{}_combined_peaks.err'.format(combined_name))
#         job_name = '{}_combined_peaks'.format(combined_name)
#         f.write(_pbs_header(out, err, job_name, threads))
#     
#     if conda_env != '':
#         f.write('source activate {}\n'.format(conda_env))
#     f.write('mkdir -p {}\n'.format(tempdir))
#     f.write('cd {}\n'.format(tempdir))
#     f.write('rsync -avz \\\n{} \\\n\t.\n\n'.format(
#         ' \\\n'.join(['\t{}'.format(x) for x in tagdirs])))
# 
#     # Add executables to path because HOMER expects them there.
#     f.write('source {}\n\n'.format(environment))
#     
#     # Run HOMER.
#     td = [os.path.split(os.path.realpath(x))[1] for x in tagdirs]
#     lines, softlink_lines = _combined_homer(td, combined_name, temp_tagdir,
#                                             final_tagdir, homer_path, link_dir,
#                                             bedtools_path, tracklines_file,
#                                             web_path_file, bigwig=True)
#     f.write(lines)
#     f.write('wait\n\n')
# 
#     if tempdir != outdir:
#         f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
#             ' \\\n\t'.join(files_to_copy), outdir))
# 
#     if len(files_to_remove) > 0:
#         f.write('rm -r \\\n\t{}\n\n'.format(' \\\n\t'.join(files_to_remove)))
# 
#     if tempdir != outdir:
#         f.write('rm -r {}\n'.format(tempdir))
# 
#     f.write(softlink_lines)
# 
#     f.close()
#     return fn
# 
# def _nucleoatac(bam, bed, sample_name, fasta, threads):
#     lines = ('nucleoatac run --bed {} \\\n'
#              '\t--bam {} \\\n'
#              '\t--out {} \\\n'
#              '\t--fasta {} \\\n'
#              '\t--cores {}\n\n'.format(bed, bam, sample_name, fasta, threads))
#     return lines
# 
# def nucleoatac(
#     bam, 
#     bed, 
#     sample_name, 
#     fasta, 
#     outdir,
#     tracklines_file,
#     link_dir,
#     web_path_file,
#     environment,
#     conda_env='',
#     tempdir='/scratch', 
#     threads=4, 
#     shell=False,
# ):
#     """
#     Make a PBS or shell script for estimating nucelosome occupancy using
#     ATAC-seq data.
# 
#     Parameters
#     ----------
#     bam : str 
#         Path to bam file with aligned reads to use for estimating occupancy.
# 
#     bed : str
#         Path to bed file with positions to estimate occupancy for.
# 
#     sample_name : str
#         Name used for naming files and directories.
# 
#     fasta : str
#         Path to genome fasta. Must be indexed.
# 
#     outdir : str
#         Directory to store directory containing PBS/shell file and results.
# 
#     tracklines_file : str
#         Path to file for writing tracklines. The tracklines will be added to the
#         file; the contents of the file will not be overwritten. These tracklines
#         can be pasted into the genome browser upload for custom data.
# 
#     link_dir : str
#         Path to directory where softlinks for genome browser should be made.
# 
#     web_path_file : str
#         File whose first line is the URL that points to link_dir. For example,
#         if we make a link to the file s1_coord_sorted.bam in link_dir and
#         web_path_file has http://site.com/files on its first line, then
#         http://site.com/files/s1_coord_sorted.bam should be available on the
#         web. If the web directory is password protected (it probably should be),
#         then the URL should look like http://username:password@site.com/files.
#         This is a file so you don't have to make the username/password combo
#         public (although I'd recommend not using a sensitive password). You can
#         just put the web_path_file in a directory that isn't tracked by git, 
#         figshare, etc.
# 
#     environment : str
#         Bash file with PATH information that can be sourced. This should include
#         the paths to executables HOMER will need like bedGraphToBigWig.
# 
#     conda_env : str
#         If provided, load conda environment with this name. This will control
#         which version of nucleoatac is used.
# 
#     tempdir : str
#         Directory to store files as nucleoatac runs.
# 
#     threads : int
#         Number of threads to reserve using PBS scheduler and for nucleoatac to
#         use.
# 
#     shell : boolean
#         If true, make a shell script rather than a PBS script.
#     
#     Returns
#     -------
#     fn : str
#         Path to PBS/shell script.
# 
#     """
#     if shell:
#         pbs = False
#     else: 
#         pbs = True
# 
#     tempdir = os.path.join(tempdir, '{}_nucleoatac'.format(sample_name))
#     outdir = os.path.join(outdir, '{}_nucleoatac'.format(sample_name))
# 
#     # I'm going to define some file names used later.
#     temp_bam = os.path.join(tempdir, os.path.split(bam)[1])
#     
#     # Files to copy to output directory.
#     files_to_copy = []
#     # Temporary files that can be deleted at the end of the job. We may not want
#     # to delete the temp directory if the temp and output directory are the
#     # same.
#     files_to_remove = [temp_bam]
# 
#     try:
#         os.makedirs(outdir)
#     except OSError:
#         pass
# 
#     if shell:
#         fn = os.path.join(outdir, '{}_nucleoatac.sh'.format(sample_name))
#     else:
#         fn = os.path.join(outdir, '{}_nucleoatac.pbs'.format(sample_name))
# 
#     f = open(fn, 'w')
#     f.write('#!/bin/bash\n\n')
#     if pbs:
#         out = os.path.join(outdir,
#                            '{}_nucleoatac.out'.format(sample_name))
#         err = os.path.join(outdir,
#                            '{}_nucleoatac.err'.format(sample_name))
#         job_name = '{}_nucleoatac'.format(sample_name)
#         f.write(_pbs_header(out, err, job_name, threads))
#     
#     if conda_env != '':
#         f.write('source activate {}\n'.format(conda_env))
#     f.write('mkdir -p {}\n'.format(tempdir))
#     f.write('cd {}\n'.format(tempdir))
# 
#     f.write('source {}\n\n'.format(environment))
#     f.write('rsync -avz \\\n\t{} \\\n\t.\n\n'.format(bam))
#     
#     # Run nucleoatac.
#     lines = _nucleoatac(temp_bam, bed, sample_name, fasta, threads)
#     f.write(lines)
#     f.write('wait\n\n')
# 
#     if tempdir != outdir:
#         if len(files_to_copy) > 0:
#             f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
#                 ' \\\n\t'.join(files_to_copy), outdir))
# 
#     if len(files_to_remove) > 0:
#         f.write('rm -r \\\n\t{}\n\n'.format(' \\\n\t'.join(files_to_remove)))
# 
#     if tempdir != outdir:
#             f.write('rsync -avz {} {}\n\n'.format(os.path.join(tempdir, '*'),
#                                                   outdir))
# 
#     if tempdir != outdir:
#         f.write('rm -r {}\n'.format(tempdir))
# 
#     f.close()
#     return fn
# 
# def motif_analysis(
#     bed, 
#     sample_name, 
#     outdir,
#     tracklines_file,
#     link_dir,
#     web_path_file,
#     environment,
#     mask=False,
#     conda_env='',
#     tempdir='/scratch', 
#     threads=4, 
#     shell=False,
# ):
#     """
#     Make a PBS or shell script for analyzing motifs with HOMER.
# 
#     Parameters
#     ----------
#     bed : str 
#         Bed file with positions to analyze.
# 
#     sample_name : str
#         Name used for naming files and directories.
# 
#     outdir : str
#         Directory to store directory containing PBS/shell file and results.
# 
#     tracklines_file : str
#         Path to file for writing tracklines. The tracklines will be added to the
#         file; the contents of the file will not be overwritten. These tracklines
#         can be pasted into the genome browser upload for custom data.
# 
#     link_dir : str
#         Path to directory where softlinks for genome browser should be made.
# 
#     web_path_file : str
#         File whose first line is the URL that points to link_dir. For example,
#         if we make a link to the file s1_coord_sorted.bam in link_dir and
#         web_path_file has http://site.com/files on its first line, then
#         http://site.com/files/s1_coord_sorted.bam should be available on the
#         web. If the web directory is password protected (it probably should be),
#         then the URL should look like http://username:password@site.com/files.
#         This is a file so you don't have to make the username/password combo
#         public (although I'd recommend not using a sensitive password). You can
#         just put the web_path_file in a directory that isn't tracked by git, 
#         figshare, etc.
# 
#     environment : str
#         Bash file with PATH information that can be sourced. This should include
#         the paths to executables HOMER will need like bedGraphToBigWig.
# 
#     mask : bool
#         Whether to pass the -mask parameter to HOMER.
# 
#     conda_env : str
#         If provided, load conda environment with this name. 
# 
#     tempdir : str
#         Directory to store files as motif analysis runs.
# 
#     threads : int
#         Number of threads to reserve using PBS scheduler and for HOMER to use.
# 
#     shell : boolean
#         If true, make a shell script rather than a PBS script.
#     
#     Returns
#     -------
#     fn : str
#         Path to PBS/shell script.
# 
#     """
#     if shell:
#         pbs = False
#     else: 
#         pbs = True
# 
#     tempdir = os.path.join(tempdir, '{}_motif'.format(sample_name))
#     outdir = os.path.join(outdir, '{}_motif'.format(sample_name))
# 
#     # I'm going to define some file names used later.
#     
#     # Files to copy to output directory.
#     files_to_copy = []
#     # Temporary files that can be deleted at the end of the job. We may not want
#     # to delete the temp directory if the temp and output directory are the
#     # same.
#     files_to_remove = []
# 
#     try:
#         os.makedirs(outdir)
#     except OSError:
#         pass
# 
#     if shell:
#         fn = os.path.join(outdir, '{}_motif.sh'.format(sample_name))
#     else:
#         fn = os.path.join(outdir, '{}_motif.pbs'.format(sample_name))
# 
#     f = open(fn, 'w')
#     f.write('#!/bin/bash\n\n')
#     if pbs:
#         out = os.path.join(outdir,
#                            '{}_motif.out'.format(sample_name))
#         err = os.path.join(outdir,
#                            '{}_motif.err'.format(sample_name))
#         job_name = '{}_motif'.format(sample_name)
#         f.write(_pbs_header(out, err, job_name, threads))
#     
#     if conda_env != '':
#         f.write('source activate {}\n'.format(conda_env))
#     f.write('mkdir -p {}\n'.format(tempdir))
#     f.write('cd {}\n'.format(tempdir))
# 
#     f.write('source {}\n\n'.format(environment))
# 
#     # Prepare some stuff for making softlinks and web links.
#     link_dir = os.path.join(link_dir, 'atac', 'motif')
#     with open(web_path_file) as wpf:
#         web_path = wpf.readline().strip()
#     web_path = web_path + '/atac/motif'
#     if os.path.exists(tracklines_file):
#         with open(tracklines_file) as tf:
#             tf_lines = tf.read()
#     else:
#         tf_lines = ''
#     try:
#         os.makedirs(link_dir)
#     except OSError:
#         pass
# 
#     # Run HOMER motif analysis.
#     if mask:
#         lines = ('findMotifsGenome.pl {} hg19 {} -size given -mask '
#                  '-p {}\n'.format(bed, outdir, threads))
#     else:
#         lines = 'findMotifsGenome.pl {} hg19 {} -size given -p {}\n'.format(
#             bed, outdir, threads)
#     f.write(lines)
#     f.write('wait\n\n')
# 
#     new_lines, name = _make_softlink(outdir, sample_name, link_dir)
#     f.write(new_lines)
#     tf_lines += '{}/{}\n'.format(web_path, os.path.split(outdir)[1])
# 
#     # Write tracklines and URLs.
#     with open(tracklines_file, 'w') as tf:
#         tf.write(tf_lines)
#     if tempdir != outdir:
#         if len(files_to_copy) > 0:
#             f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
#                 ' \\\n\t'.join(files_to_copy), outdir))
# 
#     if len(files_to_remove) > 0:
#         f.write('rm -r \\\n\t{}\n\n'.format(' \\\n\t'.join(files_to_remove)))
# 
#     if tempdir != outdir:
#             f.write('rsync -avz {} {}\n\n'.format(os.path.join(tempdir, '*'),
#                                                   outdir))
# 
#     if tempdir != outdir:
#         f.write('rm -r {}\n'.format(tempdir))
# 
#     f.close()
#     return fn
