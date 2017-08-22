import os  
import glob 
import itertools as it 
import datetime as dt
from cdpipelines.general import JobScript

JOBNAMES = ['bamtofastq', 'run_phlat', 'align_to_hla', 'run_vbseq']

def get_scripts(
        outdir,
        sample
    ):  
    """
    Return a glob of SGE/shell scripts for running the entire HLA pipeline. 

    Parameters
    __________
    sample : str, UUID for the sample.
    
    Returns 
    -------
    shell_fns : str,

    """
    shell_fns = glob.glob(os.path.join(outdir, sample, 'sh', '*.sh'))
    return shell_fns 



def get_logs(
        outdir,
        sample
    ):  
    """
    Return a glob of SGE/shell logs from the HLA pipeline. 

    Parameters
    __________
    sample : str, UUID for the sample.
    
    Returns 
    -------
    shell_fns : str,

    """
    log_fns = glob.glob(os.path.join(outdir, sample, 'logs', '*'))
    return log_fns 



def delete_script(
        outdir, 
        sample,
        jobname
    ):
    """
    Delete SGE/shell scripts. 

    Parameters
    __________
    jobnames : list, list of jobnames to delete (specified name are in the pipeline function)

    """
    shell_fns = glob.glob(os.path.join(outdir, sample, 'sh', '{}_{}*.sh'.format(sample, jobname)))
    try:
        for shell_fn in shell_fns:
            os.remove(shell_fn)
        return (sample, jobname, True)
    except:
        return (sample, jobname, False)



class WGS_HLAJobScript(JobScript):

    def bwa_index(
        self,
        in_fasta,
    ):
        """
        Index fasta file using bwa.

        Parameters
        ----------
        in_fasta : str
            Path to file input fasta file.

        index : str
            Path to index file to be written. If not provided, the index is
            written to the bwa default {in_bam}.bai in the current working
            directory.

        Returns
        -------
        index : str
            Path to index file for input fasta file.

        """
        indexes = [os.path.join(self.tempdir, os.path.splitext(in_fasta)[0] + ext) \
                    for ext in ['.amb', '.ann', '.bwt', '.pac', '.sa']]
        
        cmd = 'bwa index {}'.format(in_fasta)
        lines = self._add_execution_date(cmd)
        with open(self.filename, "a") as f:
            lines = '\n'.join(lines)
            f.write(lines)
        return indexes
    
    def sambamba_index(
        self,
        in_bam, 
        sambamba_path='sambamba',
        suffix=None, 
    ):
        """
        Index bam file using sambamba.
    
        Parameters
        ----------
        in_bam : str
            Path to file input bam file.
    
        bg : boolean
            Whether to run the process in the background.
    
        Returns
        -------
        index : str
            Path to output index file.
    
        """
        index = os.path.join(in_bam + '.bai')
        cmd = '{} index -t {} \\\n\t{} \\\n\t{}'.format(
            sambamba_path, self.threads, in_bam, index)
        lines = self._add_execution_date(cmd)
        with open(self.filename, "a") as f:
            lines = '\n'.join(lines)
            f.write(lines)
        return index

    def picard_downsampling(self,
        in_bam,
        random_seed,
        read_probability,
        suffix=None,
    ):
        """
        Down sample the bam file using picard.

        Parameters
        ----------
        in_bam : str
            Path to file input bam file.
            
        random_seed: int
            Random seed for random sampling. 
        
        read_probability: float
            Probability of keeping a read. 
        
        Returns
        -------
        
        out_bam : str
            Path to file output bam file.
        """
        ds_int = int(read_probability * 100)
        out_bam = os.path.basename(in_bam).replace('ds100', 'ds{}'.format(ds_int))
        out_bam = os.path.join(self.outdir, out_bam)
        cmd = 'java -Xmx2G -jar $picard DownsampleSam I={} O={} RANDOM_SEED={} P={}'.\
            format(in_bam, out_bam, random_seed, read_probability)
            
        lines = self._add_execution_date(cmd)
        with open(self.filename, "a") as f:
            lines = '\n'.join(lines)
            f.write(lines)
        return out_bam

    
    #    def samtools_extract_regions(self, in_bam, regions):
    #        
    #        Extract hits aligning to HLA regions
    #        
    #        Parameters
    #        __________
    #        in_bam : str
    #            Path to input bam file.
    #        regions : str
    #            Regions that will be extracted as specified by the samtools format. 
    #        
    #        Returns
    #        -------
    #        mhc_bam : str
    #            Path to bam with extractions file.
    #    
    #                
    #        # GRCH37 HLA regions 
    #        # Make reads directory for extraction of HLA R1 and R2 fastq files 
    #        #### Deprecated Sambamba version has some bug. ####
    #        # Version 0.5.9 and 0.6.1 produce output but switching to samtools to be sure
    #        # with open(self.filename, 'a') as f:
    #            # f.write('sambamba view -h -f bam {} {} -o {} -t {}'.format(bam, hla_regions, extraction, self.threads))
    #        mhc_bam = os.path.join(self.outdir,'{}_hla.bam'.format(self.sample_name))
    #        #### Samtools version 
    #        cmd = 'samtools view -h -b -@ {} -o {} {} {}'.format(self.threads, mhc_bam, in_bam, regions)
    #        lines = self._add_execution_date(cmd)
    #
    #        with open(self.filename, 'a') as f:
    #            f.write('\n'.join(lines))
    #        return out_bam

    def samtools_extract_bed(self, in_bam, bed_fn, loci='extract'):
        """
        Extract hits aligning to HLA regions
        
        Parameters
        __________
        in_bam : str
            Path to input bam file.
        regions : str
            Regions that will be extracted as specified by the samtools format. 
        
        Returns
        -------
        mhc_bam : str
            Path to bam with extractions file.
    
        """        
        # GRCH37 HLA regions 
        # Make reads directory for extraction of HLA R1 and R2 fastq files 
        #### Deprecated Sambamba version has some bug. ####
        # Version 0.5.9 and 0.6.1 produce output but switching to samtools to be sure
        # with open(self.filename, 'a') as f:
            # f.write('sambamba view -h -f bam {} {} -o {} -t {}'.format(bam, hla_regions, extraction, self.threads))
        out_bam = os.path.join(self.outdir,'{}_{}.bam'.format(self.sample_name, loci))
        #### Samtools version 
        cmd = 'samtools view -h -b -@ {} -L {} -o {} {}'.format(self.threads, bed_fn, out_bam, in_bam)
        lines = self._add_execution_date(cmd)

        with open(self.filename, 'a') as f:
            f.write('\n'.join(lines))
        return out_bam
    
    def bedtools_bamtofastq(self, in_bam, paired=False):
        """
        Converts bam files into R1 and R2 fastq files.
        
        Parameters
        __________
        in_bam : str
            Path to input bam file.
        paired : bool
            Whether the bam file if from a pair-end sequence. 
        
        Returns
        -------
        r1_fastq : str
            Path to r1 fastq file.
        r2_fastq : str
            Path to r2 fastq file.
        """        
        r1_fastq = os.path.join(self.outdir, '{}_R1.fastq'.format(self.sample_name))
        r2_fastq = os.path.join(self.outdir, '{}_R2.fastq'.format(self.sample_name))
        if paired:
            cmd = 'bedtools bamtofastq -i {} -fq {} -fq2 {}'.format(in_bam, r1_fastq, r2_fastq)
            lines = self._add_execution_date(cmd)
            with open(self.filename, 'a') as f:
                f.write('\n'.join(lines))
            return (r1_fastq, r2_fastq)
        else:
            cmd = 'bedtools bamtofastq -i {} -fq {}'.format(in_bam, r1_fastq)
            lines = self._add_execution_date(cmd)
            with open(self.filename, 'a') as f:
                f.write('\n'.join(lines))
            return (r1_fastq,)

    def bwa_map_to_allele(
        self,
        allele_ref,
        r1_fastq,
        r2_fastq, 
        allele_name, 
        stringent=False,
    ):  
        """
        Aligning paired end data to specific allele reference sequences using bwa.  
        
        - bwa P is used for paired-end data. 
        - samtools -F 4 is used to extract mapped reads only. 
        
        
        Parameters
        __________
        hla_ref : str
            Path to hla_ref file.
        r1_fastq : str
            Path to input r1 fastq file.
        r2_fastq : str
            Path to input r2 fastq file.
            Path to bowtie2.
        
        Returns
        -------
        out_bam : str
            Path to bam file.
        """         
        out_bam = os.path.join(self.outdir, '{}_{}.bam'.format(self.sample_name, allele_name))

        if stringent: 
            cmd = 'bwa mem -P -B 40 -O 60 -E 10 -L 10000 -t {} {} {} {} | samtools view -F 4 -b - > {}'.\
                format(self.threads,
                allele_ref,
                r1_fastq,
                r2_fastq,
                out_bam)
        else: 
            cmd = 'bwa mem -P -L 10000 -t {} {} {} {} | samtools view -F 4 -b - > {}'.\
                format(self.threads,
                allele_ref,
                r1_fastq,
                r2_fastq,
                out_bam)

        lines = self._add_execution_date(cmd)
        with open(self.filename, 'a') as f:
            f.write('\n'.join(lines))
        return out_bam     
    
    def bwa_multi_map(
        self,
        r1_fastq,
        r2_fastq, 
        hla_ref,
    ):  
        """
        Aligning paired end data multiples times using bwa. The -a option 
        which produces all the alignments, not just the primary 
        alignment.

        -L 10000 is used to avoid clipping 
        
        Parameters
        __________
        hla_ref : str
            Path to hla_ref file.
        r1_fastq : str
            Path to input r1 fastq file.
        r2_fastq : str
            Path to input r2 fastq file.
            Path to bowtie2.
        
        Returns
        -------
        sam : str
            Path to phlat sum file.
        """         
        out_bam = os.path.join(self.outdir, '{}_hla_db_multi_mapped.bam'.format(self.sample_name))
        cmd = 'bwa mem -t {} -L 10000 -P -a {} {} {} | samtools view -b - > {}'.\
            format(self.threads,
            hla_ref,
            r1_fastq,
            r2_fastq,
            out_bam)
        lines = self._add_execution_date(cmd)
        with open(self.filename, 'a') as f:
            f.write('\n'.join(lines))
        return out_bam     

    def bedtools_coverage(self, 
        bam, 
        bed, 
        bedtools_path='bedtools',
    ):
        """
        Calculate the coverage.
    
        Parameters
        ----------
        bam : str
            Bam file to calculate coverage for.

        bed : str
            Bed file to filter sites on.
    
        bedtools_path : str
            Path to bedtools. If bedtools_path == 'bedtools', it is assumed that
            the hg19 human.hg19.genome file from bedtools is also in your path.

        Returns
        -------
        coverage : str
            Path to output bedgraph file.
    
        """

        coverage = os.path.join(self.tempdir, '{}_mhc_region.coverage'.format(self.sample_name))
        cmd = '{} coverage -b {} -a {} > {}'.format(bedtools_path, bed, bam, coverage)

        lines = self._add_execution_date(cmd)
        with open(self.filename, "a") as f:
            f.write('\n'.join(lines))
        return coverage  

    def samtools_depth(self, 
        bam, 
        bed, 
        bedtools_path='samtools',
    ):
        """
        Calculate the depths at the MHC and surrounding region.
    
        Parameters
        ----------
        bam : str
            Bam file to calculate coverage for.

        bed : str
            Bed file to filter sites on the MHC region and 10 random MHC regions.
    
        bedtools_path : str
            Path to bedtools. If bedtools_path == 'bedtools', it is assumed that
            the hg19 human.hg19.genome file from bedtools is also in your path.

        Returns
        -------
        coverage : str
            Path to output bedgraph file.
    
        """

        coverage = os.path.join(self.tempdir, '{}_mhc_and_surrounding.coverage'.format(self.sample_name))
        cmd = '{} depth -b {} {} > {}'.format(bedtools_path, bed, bam, coverage)

        lines = self._add_execution_date(cmd)
        with open(self.filename, "a") as f:
            f.write('\n'.join(lines))
        return coverage  

    def samtools_mpileup(self, 
        sorted_bam, 
        allele_fasta,
        allele_bed, 
        samtools_path='samtools',
        suffix=None,
    ):
        """
        Calculate the depths at the MHC and surrounding region.
    
        Parameters
        ----------
        bam : str
            Bam file to calculate coverage for.

        bed : str
            Bed file to filter sites on the MHC region and 10 random MHC regions.
    
        samtools_path : str
            Path to bedtools. If bedtools_path == 'bedtools', it is assumed that
            the hg19 human.hg19.genome file from bedtools is also in your path.

        Returns
        -------
        coverage : str
            Path to output bedgraph file.
    
        """
        if suffix:
            mpileup = os.path.join(self.outdir, '{}_{}.mpileup'.format(self.sample_name, suffix))

        else:
            mpileup = os.path.join(self.outdir, '{}.mpileup'.format(self.sample_name))

        cmd = '{} mpileup -A -a -f {} --positions {} -o {} {}'.format(samtools_path, allele_fasta, allele_bed, mpileup, sorted_bam)

        lines = self._add_execution_date(cmd)
        with open(self.filename, "a") as f:
            f.write('\n'.join(lines))
        return mpileup  


    def phlat_typing(
        self,
        r1_fastq,
        r2_fastq, 
        bowtie2 = '/software/bowtie2-2.2.6/bowtie2',
        phlat_dir = '/frazer01/home/joreyna/software/Release/phlat-release',
    ):  
        """
        HLA type paired fastq files with PHLAT.
        
        Parameters
        __________
        r1_fastq : str
            Path to input r1 fastq file.
        r2_fastq : str
            Path to input r2 fastq file.
        bowtie2 : str
            Path to bowtie2.
        phlat_dir : str
        
        Returns
        -------
        phlat_fn : str
            Path to phlat sum file.
        """         
        phlat_fn = os.path.join(self.outdir, '{}_HLA.sum'.format(self.sample_name))
        new_phlat_fn = phlat_fn.replace('_HLA','').replace('sum', 'phlat')

        cmd = 'python -O {} -1 {} -2 {} -index {} -p {} -b2url {} -orientation "--fr" -tag {} -e {} -o {}'.\
            format(os.path.join(phlat_dir, 'dist/PHLAT.py'),
            r1_fastq,
            r2_fastq,
            os.path.join(phlat_dir, 'b2folder'),
            self.threads, 
            bowtie2, self.sample_name, phlat_dir, self.outdir)
        lines = self._add_execution_date(cmd)
        lines.append('mv {} {}'.format(phlat_fn, new_phlat_fn))
        
        with open(self.filename, 'a') as f:
            f.write('\n'.join(lines))
        return new_phlat_fn 
    
    def vbseq_typing(
        self,
        sam,
        hla_ref,
    ):
        """
        HLA type paired fastq files with HLA-VBSeq.
        
        Parameters
        __________
        sam : str
            Path to input sam file.
        hla_ref : str
            Path to the hla reference file (WGS uses the hla_gen file whereas RNAS uses the hla_nuc file.)
        
        Returns
        -------
        vbseq_fn : str
            Path to vbseq file.
        """
        
        vbseq_fn = os.path.join(self.outdir, '{}.vbseq'.format(self.sample_name))        
        cmd = 'java -Xmx2G -jar $HLA_VBSeq {} {} {} --alpha_zero 0.01 --is_paired'.format(hla_ref, sam, vbseq_fn)
        lines = self._add_execution_date(cmd)
        with open(self.filename, 'a') as f:
            f.write('\n'.join(lines))
        return vbseq_fn 
    
    def parse_vbseq_results(
        self,
        hla,
        hla_allele_list, 
        email=False,
        hla_vbseq_parse_results_path='/repos/cardips-pipelines/HLA_Typing/sources/parse_result.pl',
    ):
        """
        Parse results from HLA-VBSeq.

        Parameters
        __________
        hla : str
            Path to the vbseq file.
        hla_allele_list: str 
            Path to the HLA allele list.

        Returns
        -------
        vbseq_result_fn : str
            Path to the parsed vbseq file.
        """

        vbseq_result_fn = os.path.join(self.outdir, '{}.vbseq.avgdp'.format(self.sample_name))        
        cmd = 'perl {} {} {} > {}'.format(hla_vbseq_parse_results_path, hla_allele_list, hla, vbseq_result_fn)
        lines = self._add_execution_date(cmd)
        with open(self.filename, 'a') as f:
            f.write('\n'.join(lines))
            if email == True:
                f.write('echo Hello Mars! | mail -r joreyna@flh1.ucsd.edu -s "Jobs are complete." joreyna@live.com\n')
        return vbseq_result_fn  

    def generate_vbseq_top_alleles_fasta(self, resolution, top, hla_bed, hla_fasta, raw_vbseq):
        """
        Down sample the bam file using picard.

        Parameters
        ----------
        top: int 
            Number of top alleles to consider
        resolution: int
            HLA type resolution to consider. 
        hla_bed: str
            Bed file derived from hla_gen.fasta
        hla_fasta
            HLA reference file from IPD/IMGT HLA database, (e.g. hla_gen.fasta)
            
            
        
        Returns
        -------
        
        out_bam : str
            Path to file output bam file.
        """   
        
        output_bed = os.path.join(self.outdir, '{}_res_{}_top_{}.bed'.format(self.sample_name, resolution, top))
        output_fasta  = os.path.join(self.outdir, '{}_res_{}_top_{}.fasta'.format(self.sample_name, resolution, top))
        
        args = [self.sample_name, resolution, top, hla_bed, hla_fasta, self.outdir, raw_vbseq]
        
        from __init__ import _scripts 
        generate_top_alleles_fasta_script = os.path.join(_scripts, 'generate_vbseq_top_alleles_fasta.py')
        cmd = 'python {} \\'.format(generate_top_alleles_fasta_script)
        for i, x in enumerate(args):
            if i < len(args) - 1:
                cmd += '\n\t' + str(x) + ' \\'
            else:
                cmd += '\n\t' + str(x)
                
        lines = self._add_execution_date(cmd)
        with open(self.filename, "a") as f:
            lines = '\n'.join(lines)
            f.write(lines)
            
        return output_bed, output_fasta

    
    def bash_submit_command(self):
        """Get command to submit script."""
        if self.wait_for:
            return 'bash {} wait'.format(self.filename)
        else:
            return 'bash {}'.format(self.filename)


def pipeline(
        sample_name,
        in_bam,
        linkdir, 
        outdir,
        hla_regions_bed,
        coverage_bed,
        hla_ref,
        hla_allele_list,
        queue = None,
        webpath = None,
        email = False,
    ):  
    """
    Make SGE/shell scripts for running the entire HLA pipeline. The defaults
    are set for use on the Frazer lab's SGE schedule on flh1/flh2. 

    Parameters
    __________
    in_bam : str
    linkdir : str, 
    outdir : str,
    hla_regions_bed: str,
    coverage_bed: str,
    bowtie2 : str
    phlat_dir : str
    queue : str

    """
    
    submit_commands = [] 

    #### Job #1: Extracting alignments in the MHC loci and surrounding region into an MHC BAM file. ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='extract_mhc_data',
        threads=8,
        memory=32,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'reads'),
        queue=queue,
        conda_env='hla',
        modules='samtools/1.2')
    extract_mhc_job = job.jobname
    job.add_input_file(in_bam)
    mhc_bam = job.samtools_extract_bed(in_bam, coverage_bed, loci='mhc')
    job.add_output_file(mhc_bam)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    #### Job #2: Calculating depth for the MHC loci and surroudng region using the MHC BAM. ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='calculate_depth',
        threads=1,
        memory=4,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'reads'),
        queue=queue,
        conda_env='hla',
        modules='samtools/1.2', 
        wait_for=[extract_mhc_job])
    depth_job = job.jobname
    job.add_input_file(mhc_bam)
    depth = job.samtools_depth(mhc_bam, coverage_bed)
    job.add_output_file(depth)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    #### Job #3: Extracting alignments in the HLA loci from the MHC BAM. ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='extract_hla_data',
        threads=8,
        memory=32,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'reads'),
        queue=queue,
        conda_env='hla',
        modules='samtools/1.2',
        wait_for=[extract_mhc_job])
    extract_hla_job = job.jobname
    job.add_input_file(mhc_bam)
    hla_bam = job.samtools_extract_bed(mhc_bam, hla_regions_bed, loci='hla')
    job.add_output_file(hla_bam)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    #### Job #3: Query sorting the HLA BAM. ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='query_sort_hla_bam',
        threads=8,
        memory=32,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'reads'),
        queue=queue,
        conda_env='hla',
        modules='sambamba/0.6.1',
        wait_for=[extract_hla_job])
    qsort_job = job.jobname 
    job.add_input_file(hla_bam)
    qsort_bam = job.sambamba_sort(hla_bam, queryname=True)
    job.add_output_file(qsort_bam)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    #### Job #4: Converting the HLA BAM into R1 and R2 fastq's ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='bamtofastq',
        threads=1,
        memory=4,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'reads'),
        queue=queue,
        conda_env='hla',
        modules='bedtools/2.25.0',
        wait_for=[qsort_job])
    bamtofastq_job = job.jobname
    job.add_input_file(in_bam)
    (fastq_r1, fastq_r2) = job.bedtools_bamtofastq(qsort_bam, paired=True)
    fastq_r1 = job.add_output_file(fastq_r1)
    fastq_r2 = job.add_output_file(fastq_r2)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    # Removing the indexing. I'm going to assume that the index files are already there.
    # This is true for most coverage experiments except for 25. DJ forgot to index them 
    # or deleted the data for these files.  
    #index_bam = job.sambamba_index(in_bam)
    #job.add_temp_file(index_bam)


    #### Job #5: Generate HLA types using PHLAT #### 
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='run_phlat',
        threads=8,
        memory=32,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'hla'),
        queue=queue,
        conda_env='hla',
        modules='',
        wait_for=[bamtofastq_job])
    phlat_job = job.jobname
    job.add_input_file(fastq_r1)
    job.add_input_file(fastq_r2)
    phlat = job.phlat_typing(fastq_r1, fastq_r2)
    job.add_output_file(phlat)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command()) 
        
    ### Job #6: Preprocessing for VBSeq, aligning to HLA sequence database. ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='multi_map_mhc_reads',
        threads=8,
        memory=32,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'reads'),
        queue=queue,
        conda_env='hla',
        modules='bwa,samtools/1.2',
        wait_for=[bamtofastq_job])
    multi_map_mhc_reads_job = job.jobname
    job.add_input_file(fastq_r1)
    job.add_input_file(fastq_r2)
    multi_map_sam = job.bwa_multi_map(fastq_r1, fastq_r2, hla_ref)
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
        
    #### Job #7: Running VBSeq and parsing results. ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='run_vbseq',
        threads=1,
        memory=5,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'hla'),
        queue=queue,
        conda_env='hla',
        modules='HLA-VBSeq',
        wait_for=[multi_map_mhc_reads_job])
    vbseq_job = job.jobname
    job.add_input_file(multi_map_sam)
    raw_vbseq = job.vbseq_typing(multi_map_sam, hla_ref)
    parsed_vbseq = job.parse_vbseq_results(raw_vbseq, hla_allele_list, email)
    job.add_output_file(raw_vbseq)
    job.add_output_file(parsed_vbseq)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
        
    ##### Submission script #####
    now = str(dt.datetime.now())
    now = now.replace('-', '_').replace(' ', '_').replace(':', '_').replace('.', '_')
    submit_fn = os.path.join(outdir, 'sh/', '{}_submit_{}.sh'.format(sample_name, now))
    with open(submit_fn, 'w') as f:
            f.write('#!/bin/bash\n\n')
            f.write('\n'.join(submit_commands))   
    return submit_fn

def vbseq_variability_downsample_pipeline(
        sample_name,
        mhc_bam,
        linkdir, 
        outdir,
        hla_regions_bed,
        coverage_bed,
        hla_ref,
        hla_allele_list,
        seed,
        read_probability,
        queue = None,
        webpath = None,
        email = False,
        run_all = False
    ):  
    """
    Make SGE/shell scripts for running the entire HLA pipeline. The defaults
    are set for use on the Frazer lab's SGE schedule on flh1/flh2. 

    Parameters
    __________
    mhc_bam : str
    linkdir : str, 
    outdir : str,
    hla_regions_bed: str,
    coverage_bed: str,
    read_probability: float,
    seed: int,
    bowtie2 : str
    phlat_dir : str
    queue : str
    webpath: str
    email: bool 
    run_all: bool
        This pipeline is meant to test the variability of HLA-VBSeq
        predictions. So only some parts of the complete pipeline are
        run. If run_all is set to true then all parts will run.

    """
    
    submit_commands = [] 
   
    #### Job #1: Down samplings reads in the MHC region. ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='picard_downsample',
        threads=1,
        memory=5,
        linkdir=linkdir,
        outdir=os.path.join(outdir, 'reads'),
        queue=queue,
        conda_env='hla',
        modules='picard,sambamba')
    down_sample_job = job.jobname
    job.add_input_file(mhc_bam)
    ds_mhc_bam = job.picard_downsampling(mhc_bam, \
        seed, read_probability, suffix=None)
    ds_mhc_bam_index = job.sambamba_index(ds_mhc_bam, sambamba_path='sambamba')
    job.add_output_file(ds_mhc_bam)
    job.add_output_file(ds_mhc_bam_index)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    if run_all:
        #### Job #2: Calculating depth for the MHC loci and surroudng region using the MHC BAM. ####
        job = WGS_HLAJobScript(sample_name=sample_name,
            job_suffix='calculate_depth',
            threads=1,
            memory=4,
            linkdir=os.path.join(linkdir),
            outdir=os.path.join(outdir, 'reads'),
            queue=queue,
            conda_env='hla',
            modules='samtools/1.2', 
            wait_for=[down_sample_job])
        depth_job = job.jobname
        job.add_input_file(ds_mhc_bam)
        depth = job.samtools_depth(ds_mhc_bam, coverage_bed)
        job.add_output_file(depth)
        job.write_end()
        if not job.delete_sh:
            submit_commands.append(job.sge_submit_command())

    #### Job #3: Extracting alignments in the HLA loci from the MHC BAM. ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='extract_hla_data',
        threads=8,
        memory=32,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'reads'),
        queue=queue,
        conda_env='hla',
        modules='samtools/1.2',
        wait_for=[down_sample_job])
    extract_hla_job = job.jobname
    job.add_input_file(ds_mhc_bam)
    hla_bam = job.samtools_extract_bed(ds_mhc_bam, hla_regions_bed, loci='hla')
    job.add_output_file(hla_bam)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    #### Job #3: Query sorting the HLA BAM. ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='query_sort_hla_bam',
        threads=8,
        memory=32,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'reads'),
        queue=queue,
        conda_env='hla',
        modules='sambamba/0.6.1',
        wait_for=[extract_hla_job])
    qsort_job = job.jobname 
    job.add_input_file(hla_bam)
    qsort_bam = job.sambamba_sort(hla_bam, queryname=True)
    job.add_output_file(qsort_bam)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    #### Job #4: Converting the HLA BAM into R1 and R2 fastq's ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='bamtofastq',
        threads=1,
        memory=4,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'reads'),
        queue=queue,
        conda_env='hla',
        modules='bedtools/2.25.0',
        wait_for=[qsort_job])
    bamtofastq_job = job.jobname
    job.add_input_file(ds_mhc_bam)
    (fastq_r1, fastq_r2) = job.bedtools_bamtofastq(qsort_bam, paired=True)
    fastq_r1 = job.add_output_file(fastq_r1)
    fastq_r2 = job.add_output_file(fastq_r2)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    # Removing the indexing. I'm going to assume that the index files are already there.
    # This is true for most coverage experiments except for 25. DJ forgot to index them 
    # or deleted the data for these files.  
    #index_bam = job.sambamba_index(ds_mhc_bam)
    #job.add_temp_file(index_bam)

    if run_all:
        #### Job #5: Generate HLA types using PHLAT #### 
        job = WGS_HLAJobScript(sample_name=sample_name,
            job_suffix='run_phlat',
            threads=8,
            memory=32,
            linkdir=os.path.join(linkdir),
            outdir=os.path.join(outdir, 'hla'),
            queue=queue,
            conda_env='hla',
            modules='',
            wait_for=[bamtofastq_job])
        phlat_job = job.jobname
        job.add_input_file(fastq_r1)
        job.add_input_file(fastq_r2)
        phlat = job.phlat_typing(fastq_r1, fastq_r2)
        job.add_output_file(phlat)
        job.write_end()
        if not job.delete_sh:
            submit_commands.append(job.sge_submit_command()) 
        
    ### Job #6: Preprocessing for VBSeq, aligning to HLA sequence database. ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='multi_map_mhc_reads',
        threads=8,
        memory=32,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'reads'),
        queue=queue,
        conda_env='hla',
        modules='bwa,samtools/1.2',
        wait_for=[bamtofastq_job])
    multi_map_mhc_reads_job = job.jobname
    job.add_input_file(fastq_r1)
    job.add_input_file(fastq_r2)
    multi_map_sam = job.bwa_multi_map(fastq_r1, fastq_r2, hla_ref)
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
        
    #### Job #7: Running VBSeq and parsing results. ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='run_vbseq',
        threads=1,
        memory=5,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'hla'),
        queue=queue,
        conda_env='hla',
        modules='HLA-VBSeq',
        wait_for=[multi_map_mhc_reads_job])
    vbseq_job = job.jobname
    job.add_input_file(multi_map_sam)
    raw_vbseq = job.vbseq_typing(multi_map_sam, hla_ref)
    parsed_vbseq = job.parse_vbseq_results(raw_vbseq, hla_allele_list, email)
    job.add_output_file(raw_vbseq)
    job.add_output_file(parsed_vbseq)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
        
    ##### Submission script #####
    now = str(dt.datetime.now())
    now = now.replace('-', '_').replace(' ', '_').replace(':', '_').replace('.', '_')
    submit_fn = os.path.join(outdir, 'sh/', '{}_submit_{}.sh'.format(sample_name, now))
    with open(submit_fn, 'w') as f:
            f.write('#!/bin/bash\n\n')
            f.write('\n'.join(submit_commands))   
    return submit_fn

def vbseq_variability_pipeline(
        sample_name,
        fastq_r1,
        fastq_r2,
        hla_ref,
        hla_allele_list,
        linkdir, 
        outdir,
        queue = None,
        webpath = None,
        email = False,
    ):  
    """
    Make SGE/shell scripts for running the entire HLA pipeline. The defaults
    are set for use on the Frazer lab's SGE schedule on flh1/flh2. 

    Parameters
    __________
    multi_map_bam : str
    linkdir : str, 
    outdir : str,
    hla_regions_bed: str,
    coverage_bed: str,
    queue : str

    """
    
    submit_commands = [] 

    ### Job #1: Preprocessing for VBSeq, aligning to HLA sequence database. ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='multi_map_mhc_reads',
        threads=8,
        memory=12,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'reads'),
        queue=queue,
        conda_env='hla',
        modules='bwa,samtools/1.2')
    multi_map_mhc_reads_job = job.jobname
    job.add_input_file(fastq_r1)
    job.add_input_file(fastq_r2)
    multi_map_bam = job.bwa_multi_map(fastq_r1, fastq_r2, hla_ref)
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    #### Job #2: Running VBSeq and parsing results. ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='run_vbseq',
        threads=1,
        memory=5,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'hla'),
        queue=queue,
        conda_env='hla',
        wait_for=[multi_map_mhc_reads_job],
        modules='HLA-VBSeq'
    )
    vbseq_job = job.jobname
    job.add_input_file(multi_map_bam)
    raw_vbseq = job.vbseq_typing(multi_map_bam, hla_ref)
    parsed_vbseq = job.parse_vbseq_results(raw_vbseq, hla_allele_list, email)
    job.add_output_file(raw_vbseq)
    job.add_output_file(parsed_vbseq)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
        
    ##### Submission script #####
    now = str(dt.datetime.now())
    now = now.replace('-', '_').replace(' ', '_').replace(':', '_').replace('.', '_')
    submit_fn = os.path.join(outdir, 'sh/', '{}_submit_{}.sh'.format(sample_name, now))
    with open(submit_fn, 'w') as f:
            f.write('#!/bin/bash\n\n')
            f.write('\n'.join(submit_commands))   
    return submit_fn

def rerun_top_alleles_pipeline(
        sample_name,
        resolution,
        top, 
        hla_bed, 
        hla_ref, 
        raw_vbseq, 
        fastq_r1,
        fastq_r2,
        linkdir, 
        outdir,
        hla_allele_list,
        queue = None,
        webpath = None,
        email = False,
    ):  
    """
    Make SGE/shell scripts for running the entire HLA pipeline. The defaults
    are set for use on the Frazer lab's SGE schedule on flh1/flh2. 

    Parameters
    __________
    resolution: int
        HLA type resolution to consider. 
    top: int 
        Number of top alleles to consider
    hla_bed: str
        Bed file derived from hla_gen.fasta
    hla_ref
        HLA reference file from IPD/IMGT HLA database, (e.g. hla_gen.fasta)
    mhc_bam : str
    linkdir : str, 
    outdir : str,
    hla_regions_bed: str,
    coverage_bed: str,
    queue : str
    webpath: str
    """
    
    submit_commands = [] 
    #### Job #1: Extract the top alleles 
    job = WGS_HLAJobScript(sample_name=sample_name, 
        job_suffix='gen_top_sequences', 
        threads=1, 
        memory=5, 
        linkdir=linkdir,
        outdir = os.path.join(outdir, 'reads'), 
        queue=queue,
        conda_env='hla', 
        modules='bedtools,samtools,bwa')
    top_alleles_job = job.jobname
    job.add_input_file(raw_vbseq)
    top_alleles_bed, top_alleles_fasta = \
        job.generate_vbseq_top_alleles_fasta(resolution, top, hla_bed, hla_ref, raw_vbseq)
    job.add_output_file(top_alleles_bed)
    job.add_output_file(top_alleles_fasta)
    top_alleles_fasta_indexes = job.bwa_index(top_alleles_fasta)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    ### Job #2: Preprocessing for VBSeq, aligning to HLA sequence database. ####
    job = WGS_HLAJobScript(sample_name=sample_name,
        job_suffix='multi_map_mhc_reads',
        threads=8,
        memory=12,
        linkdir=os.path.join(linkdir),
        outdir=os.path.join(outdir, 'reads'),
        queue=queue,
        conda_env='hla',
        modules='bwa,samtools/1.2',
        wait_for=[top_alleles_job])
    multi_map_mhc_reads_job = job.jobname
    job.add_input_file(fastq_r1)
    job.add_input_file(fastq_r2)
    multi_map_bam = job.bwa_multi_map(fastq_r1, fastq_r2, top_alleles_fasta)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    #### Job #3: Running VBSeq and parsing results. ####
    job = WGS_HLAJobScript(sample_name=sample_name, 
        job_suffix='run_vbseq',
        threads=1,
        memory=5,
        linkdir=linkdir,
        outdir=os.path.join(outdir, 'hla'),
        queue=queue,
        conda_env='hla',
        modules='HLA-VBSeq',
        wait_for=[multi_map_mhc_reads_job])

    vbseq_job = job.jobname
    job.add_input_file(multi_map_bam)
    job.add_input_file(top_alleles_fasta)
    raw_vbseq = job.vbseq_typing(multi_map_bam, top_alleles_fasta)
    parsed_vbseq = job.parse_vbseq_results(raw_vbseq, hla_allele_list, email)
    job.add_output_file(raw_vbseq)
    job.add_output_file(parsed_vbseq)
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
        
    ##### Submission script #####
    now = str(dt.datetime.now())
    now = now.replace('-', '_').replace(' ', '_').replace(':', '_').replace('.', '_')
    submit_fn = os.path.join(outdir, 'sh/', '{}_submit_{}.sh'.format(sample_name, now))
    with open(submit_fn, 'w') as f:
            f.write('#!/bin/bash\n\n')
            f.write('\n'.join(submit_commands))   
    return submit_fn
