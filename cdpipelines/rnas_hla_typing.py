import os  
import glob 
import shutil
import subprocess 
import datetime as dt
from cdpipelines.general import JobScript
JOBNAMES = ['bamtofastq', 'run_phlat', 'align_to_hla', 'run_vbseq']

class RNAS_HLAJobScript(JobScript):
	
	@classmethod 
	def _add_execution_date(self, lines, cmd):
		lines.append('date 1>&2')
		lines.append('cmd="{}"'.format(cmd))
		lines.append('echo Executing: $cmd 1>&2')
		lines.append('eval $cmd')
		lines.append('date 1>&2\n\n')
		return lines 
	
	def sambamba_index(
		self,
		in_bam, 
		sambamba_path='sambamba',
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
		lines = []
		lines = self._add_execution_date(lines, cmd)
		with open(self.filename, "a") as f:
			lines = '\n'.join(lines)
			f.write(lines)
		return index
	
	def grep_extract_hla(self, in_sam):
		"""
		Extract hits aligning to HLA isoforms 
		
		Parameters
		__________
		in_bam : str
			Path to input bam file. 
		
		Returns
		-------
		extract_bam : str
			Path to bam with extractions file.
	
		"""		
		extract_bam = os.path.join(self.outdir,'{}_gencode.hla.bam'.format(self.sample_name))
		#### Samtools version 
		lines = []
		cmd = 'cat {} | grep HLA | grep -v HHLA | grep -v SCHLAP1 | '			'sambamba view -f bam -S -o {} /dev/stdin'.format(in_sam, extract_bam)
		lines = self._add_execution_date(lines, cmd)

		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))
		return extract_bam
	
	def bwa_align(
		self,
		r1_fastq,
		r2_fastq, 
		sam_suffix,
		bamify=False,
		reference='/frazer01/home/joreyna/projects/hla_typing/pipeline/supp/hla_nuc_160614.fasta',
	):  
		"""
		Align paired fastq files with bwa.
		
		Parameters
		__________
		reference : str
			Path to reference file.
		r1_fastq : str
			Path to input r1 fastq file.
		r2_fastq : str
			Path to input r2 fastq file.
		
		Returns
		-------
		sam : str
			Path to phlat sum file.
		"""		 
		sam = os.path.join(self.outdir, '{}_{}.sam'.format(self.sample_name, sam_suffix))
		lines = []
		if bamify == False:
			cmd = 'bwa mem -t {} {} {} {} > {}'				.format(self.threads,						 reference,						 r1_fastq,						 r2_fastq,
						sam)
		else:
			cmd = 'bwa mem -t {} {} {} {} | sambamba view -S -f bam -o {} /dev/stdin '				.format(self.threads,						 reference,						 r1_fastq,						 r2_fastq,
						sam.replace('.sam', '.bam'))
		lines = self._add_execution_date(lines, cmd)
		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))
		return sam   
	
	def bwa_pre_vbseq_align(
		self,
		r1_fastq,
		r2_fastq, 
		sam_suffix,
		bamify=False,
		reference='/frazer01/home/joreyna/projects/hla_typing/pipeline/supp/hla_nuc_160614.fasta',
	):  
		"""
		Align paired fastq files with bwa. Uses the -a option 
		which produces all the alignments, not just the primary 
		alignment.
		
		Parameters
		__________
		reference : str
			Path to reference file.
		r1_fastq : str
			Path to input r1 fastq file.
		r2_fastq : str
			Path to input r2 fastq file.
		
		Returns
		-------
		sam : str
			Path to phlat sum file.
		"""		 
		sam = os.path.join(self.outdir, '{}_{}.sam'.format(self.sample_name, sam_suffix))
		lines = []
		if bamify == False:
			cmd = 'bwa mem -t {} -L 10000 -P -a {} {} {} > {}'				.format(self.threads,						 reference,						 r1_fastq,						 r2_fastq,
						sam)
		else:
			sam = sam.replace('.sam', '.bam')
			cmd = 'bwa mem -t {} -L 10000 -P -a {} {} {} | sambamba view -S -f bam -o {} /dev/stdin '				.format(self.threads,						 reference,						 r1_fastq,						 r2_fastq,
						sam)
		lines = self._add_execution_date(lines, cmd)
		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))
		return sam

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
		r1_fastq = os.path.join(self.outdir, '{}_R1.hla.fastq'.format(self.sample_name))
		r2_fastq = os.path.join(self.outdir, '{}_R2.hla.fastq'.format(self.sample_name))
		lines = []
		if paired:
			cmd = 'bedtools bamtofastq -i {} -fq {} -fq2 {}'.format(in_bam, r1_fastq, r2_fastq)
			lines = self._add_execution_date(lines, cmd)
			with open(self.filename, 'a') as f:
				f.write('\n'.join(lines))
			return (r1_fastq, r2_fastq)
		else:
			cmd = 'bedtools bamtofastq -i {} -fq {}'.format(in_bam, r1_fastq)
			lines = self._add_execution_date(lines, cmd)
			with open(self.filename, 'a') as f:
				f.write('\n'.join(lines))
			return (r1_fastq,)
	
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

		lines = []
		cmd = 'python -O {} -1 {} -2 {} -index {} -b2url {} -orientation "--fr" -tag {} -e {} -o {}'			.format(os.path.join(phlat_dir, 'dist/PHLAT.py'),					 r1_fastq,					 r2_fastq,					 os.path.join(phlat_dir, 'b2folder'),					 bowtie2, self.sample_name, phlat_dir, self.outdir)		
		lines = self._add_execution_date(lines, cmd)
		lines.append('mv {} {}'.format(phlat_fn, new_phlat_fn))
		
		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))
		return new_phlat_fn 
	
	def vbseq_typing(
		self,
		sam,
		hla_ref='/frazer01/home/joreyna/projects/hla_typing/pipeline/supp/hla_nuc_160614.fasta'
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
		lines = []
		cmd = 'java -Xmx27G -jar $HLA_VBSeq {} {} {} --alpha_zero 0.01 --is_paired'.format(hla_ref, sam, vbseq_fn)
		lines = self._add_execution_date(lines, cmd)
		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))
		return vbseq_fn 
	
	def parse_vbseq_results(
		self,
		hla,
		allele_ls = '/frazer01/home/joreyna/projects/hla_typing/pipeline/supp/Allelelist_160614.txt'
	):
		"""
		Parse results from HLA-VBSeq.

		Parameters
		__________
		hla : str
			Path to the vbseq file.

		Returns
		-------
		vbseq_result_fn : str
			Path to the parsed vbseq file.
		"""

		vbseq_result_fn = os.path.join(self.outdir, '{}.vbseq.avgdp'.format(self.sample_name))		
		lines = []
		cmd = 'perl /frazer01/home/joreyna/projects/hla_typing/pipeline/supp/parse_result.pl {} {} > {}'.			format(allele_ls, hla, vbseq_result_fn)
		lines = self._add_execution_date(lines, cmd)
		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))
		return vbseq_result_fn  
	
	def bash_submit_command(self):
		"""Get command to submit script."""
		if self.wait_for:
			return 'bash {} wait'.format(self.filename)
		else:
			return 'bash {}'.format(self.filename)

	def seq2hla_typing(
		self,
		r1_fastq,
		r2_fastq, 
	):  
		"""
		HLA type paired fastq fies with seq2hla.

		Parameters
		__________
		r1_fastq : str
			Path to input r1 fastq file.
		r2_fastq : str
			Path to input r2 fastq file.

		Returns
		-------
		phlat_fn : str
			Path to seq2hla directory.
		"""		 
		seq2hla_dir = os.path.join(self.outdir, 'seq2hla')
		lines = []
		cmd = 'python $seq2hla -p {} -1 {} -2 {} -r {}'	.format(self.threads, r1_fastq, r2_fastq, os.path.join(self.outdir, 'seq2hla', self.sample_name))
		lines = self._add_execution_date(lines, cmd)		
		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))

		return (os.path.join(seq2hla_dir, 'seq2hla', self.sample_name + '-ClassI.HLAgenotype4digits'),
				os.path.join(seq2hla_dir, 'seq2hla', self.sample_name +'-ClassII.HLAgenotype4digits'))


	def filter_IMGT_HLA_sequences(self, alleles, imgt_hla_nuc_bed, imgt_hla_nuc_reference):
		"""
		filter the IMGT/HLA sequence database for the sample specific cDNA sequences.
		
		Parameters
		----------
		alleles: list
			A list of 8-digit HLA alleles.
		imgt_hla_nuc_bed: str
			Path to a IMGT/HLA bed file derived from the nuc fasta file. 
		imgt_hla_nuc_reference: str
			The nuc fasta file.
		
		Returns
		-------
		hla_alleles: str
			Path to the alleles text file.
		hla_bed: str
			Path to the sample specific hla bed file.
		hla_sequences: str
			path to the sample specific hla sequence file.
		
		"""
		
		hla_alleles = os.path.join(self.outdir, '{}.alleles'.format(self.sample_name))
		hla_bed = os.path.join(self.outdir, '{}.bed'.format(self.sample_name))
		hla_sequences = os.path.join(self.outdir, '{}.fa'.format(self.sample_name))

		lines = []
		
		# WRITE out the sequences to filter for them later in the IMGT/HLA RNA bed file
		with open(hla_alleles, 'w') as f:
			f.write('\n'.join([x.replace('*', '\*') for x in alleles]))
		
		
		# GREP for coordinates in the IMGT/HLA RNA bed file which correspond to the sample's allelic sequences 
		lines.append('grep {} -f {} > {}'.format(imgt_hla_nuc_bed, hla_alleles, hla_bed))
		
		# EXTRACT the cDNA sequences from the IMGT/HLA nuc fasta file. 
		lines.append('') 
		lines.append( \
			 'bedtools getfasta -name -fi {} -bed {} -fo {}'.\
			 format(imgt_hla_nuc_reference, hla_bed, hla_sequences))
		
		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))
		
		return hla_alleles, hla_bed, hla_sequences

	def concatenate_GENCODE_and_HLA_sequences(self, gencode_template, hla_sequences):
		"""
		
		Parameters
		----------
		gencode_template: str
			Path to the GENCODE template file, a derivative of the original GENCODE file after exluding HLA sequences.
			
		hla_sequences: str
			path to the sample specific hla sequence file.
			
			
		Returns
		-------
		gencode_reference: str
			Path to the sample specific gencode reference. 
			
		"""

		gencode_reference = os.path.join(self.outdir, '{}.gencode_w_hla_types.fa'.format(self.sample_name))
		lines = []
		lines.append('') 
		lines.append('') 
		lines.append('cat {} {} > {}'.format(hla_sequences, gencode_template, gencode_reference))
		
		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))
		
		return gencode_reference

	def rsem_prepare_references(self, bowtie2, gencode_reference):
		"""
		
		Parameters
		----------
		bowtie2: str
			Path to bowtie2.
			
		gencode_reference: str
			Path to the sample specific gencode reference. 
			
		Returns 
		-------
		rsem_reference: str
			RSEM reference location. 
			
		"""
		rsem_reference = os.path.join(self.outdir, self.sample_name)
		lines = []
		lines.append('rsem-prepare-reference --bowtie2 --bowtie2-path {} {} {}'.\
					 format(bowtie2, gencode_reference, rsem_reference))
		
		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))
		
		return rsem_reference 

	def rsem_calculate_expression(self, bowtie2, r1_fastq, r2_fastq, rsem_reference):

		"""
		
		Parameters
		----------
		bowtie2: str
			Path to bowtie2.
			
		 r1_fastq: str
			 Path to R1 fastq file. 
			 
		 r2_fastq: str
			 Path to R1 fastq file. 
			 
		rsem_reference: str
			RSEM reference location.		   
			
		Returns 
		-------
		
		rsem_gene_expression: str
			RSEM gene expression data file. 
		
		"""
		
		lines = []
		rsem_output = os.path.join(self.outdir, self.sample_name)
		rsem_gene_expression = os.path.join(rsem_output, '{}.genes.results'.format(self.sample_name))

		
		lines.append('rsem-calculate-expression -p {} --paired-end --bowtie2 --bowtie2-path {} --estimate-rspd {} {} {} {}'.\
			 format(self.threads, bowtie2, r1_fastq, r2_fastq, rsem_reference, rsem_output))
		
		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))
		
		
		return rsem_gene_expression 

# # Running for RNAS
def gene_expression_pipeline(
	sample_name = 'sample',
	alleles = None, 
	r1_fastqs = None, 
	r2_fastqs = None, 
	imgt_hla_nuc_bed = None,
	imgt_hla_nuc_reference = None, 
	gencode_template = None,
	bowtie2 = None,
	outdir = None, 
	linkdir = None,
	webpath = None, 
	queue = None, 
):
	"""
	Make SGE/shell scripts for running the entire RNA-seq pipeline. The defaults
	are set for use on the Frazer lab's SGE scheduler on flh1/flh2.
	Parameters
	----------
	sample_name: str
		The name of the sample. Used to assign file names. 
	alleles: list
		A list of 8-digit HLA alleles.
	r1_fastqs : list or str
		Either a list of paths to gzipped fastq files with R1 reads or path to a
		single gzipped fastq file with R1 reads. If you want to process SRA
		files, pass None here.
	r2_fastqs : list or str
		Either a list of paths to gzipped fastq files with R2 reads or path to a
		single gzipped fastq file with R2 reads. If you want to process SRA
		files, pass None here.
	imgt_hla_nuc_bed: str
		Path to a IMGT/HLA bed file derived from the nuc fasta file. 
	imgt_hla_nuc_reference: str
		The nuc fasta file.
	gencode_template: str
		Path to the GENCODE template file, a derivative of the original GENCODE file after exluding HLA sequences.
	bowtie2: str
		Path to bowtie2.
	outdir: str 
	linkdir: str
	webpath: str 
	queue: str 
	"""
	
	submit_commands = []

	##### Job 1: Filter the IMGT/HLA genomic bed file for the sample specific alleles. #####
	job = RNAS_HLAJobScript(
		sample_name=sample_name,
		job_suffix='rsem_prep',
		linkdir=linkdir,
		outdir= os.path.join(outdir, 'rsem_prep'),
		threads=1,
		memory=4,
		queue=queue,
		conda_env='hla',
		modules='cardips',
	)
	rsem_prep_jobname = job.jobname 


	# EXTRACT sample specific cDNA sequences 
	hla_alleles, hla_bed, hla_sequences = job.filter_IMGT_HLA_sequences(alleles, imgt_hla_nuc_bed, imgt_hla_nuc_reference)
	job.add_output_file(hla_alleles)
	job.add_output_file(hla_bed)
	job.add_output_file(hla_sequences)

	# GENERATE the sample specific GENCODE reference 
	gencode_reference = job.concatenate_GENCODE_and_HLA_sequences(gencode_template, hla_sequences)
	job.add_output_file(gencode_reference)
	job.write_end()

	if not job.delete_sh:
		submit_commands.append(job.sge_submit_command())

	##### Job 2: Run rsem-prepare-references. #####
	job = RNAS_HLAJobScript(
		sample_name=sample_name,
		job_suffix='rsem_prepare_references',
		linkdir=linkdir,
		outdir= os.path.join(outdir, 'rsem_prepare_references'),
		threads=8,
		memory=24,
		queue=queue,
		conda_env='hla',
		modules='cardips',
		wait_for=[rsem_prep_jobname]
	)
	rsem_prepare_references_jobname = job.jobname 

	# RUN rsem-prepare-references
	rsem_reference = job.rsem_prepare_references(bowtie2, gencode_reference)
	job.add_output_file(rsem_reference)

	if not job.delete_sh:
		submit_commands.append(job.sge_submit_command())

	##### Job 3: Combine fastqs and calculate gene expression using RSEM. #####
	job = RNAS_HLAJobScript(
		sample_name=sample_name,
		job_suffix='rsem_calculate_expression',
		linkdir=linkdir,
		outdir= os.path.join(outdir, 'rsem_calculate_expression'),
		threads=8,
		memory=24,
		queue=queue,
		conda_env='hla',
		modules='cardips',
		wait_for=[rsem_prepare_references_jobname]
	)

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

	# RUN rsem-calculate-expression
	rsem_gene_expression = job.rsem_calculate_expression(bowtie2, combined_r1, combined_r2, rsem_reference)
	job.add_output_file(rsem_reference)

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

# # Running for RNAS
def pipeline(
	sample_name,
	r1_fastqs, 
	r2_fastqs,
	outdir, 
	linkdir, 
	queue=None
):
	"""
	Make SGE/shell scripts for running the entire RNA-seq pipeline. The defaults
	are set for use on the Frazer lab's SGE scheduler on flh1/flh2.
	Parameters
	----------
	r1_fastqs : list or str
		Either a list of paths to gzipped fastq files with R1 reads or path to a
		single gzipped fastq file with R1 reads. If you want to process SRA
		files, pass None here.
	r2_fastqs : list or str
		Either a list of paths to gzipped fastq files with R2 reads or path to a
		single gzipped fastq file with R2 reads. If you want to process SRA
		files, pass None here.
	"""
	
	submit_commands = []
	
	##### Job 1: Combine fastqs and align with BWA. #####
	job = RNAS_HLAJobScript(
		sample_name, 
		job_suffix='transcriptome_align',
		linkdir=os.path.join(linkdir, 'reads'),
		outdir=os.path.join(outdir, 'reads'), 
		threads=8, 
		memory=32,
		queue=queue, 
		conda_env='hla',
		modules='bwa',
	)
	alignment_jobname = job.jobname

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

	# Transcriptome Alignment 
	job.add_input_file(combined_r1)
	job.add_input_file(combined_r2)
	transcriptome_sam = job.bwa_align(combined_r1, combined_r2, 
		  sam_suffix='gencode',
		  reference='/frazer01/home/joreyna/projects/hla_typing/pipeline/supp/gencode.v24lift37.transcripts.fa')
	job.add_output_file(transcriptome_sam)
	transcriptome_align = job.jobname
	if not job.delete_sh:
		submit_commands.append(job.sge_submit_command())
		
		
	##### Job 2: Isoform extraction and fastq generation. #####
	job = RNAS_HLAJobScript(
		sample_name, 
		job_suffix='bamtofastq',
		linkdir=os.path.join(linkdir, 'reads'),
		outdir=os.path.join(outdir, 'reads'), 
		threads=1, 
		memory=10,
		queue=queue, 
		conda_env='hla',
		modules='sambamba,bedtools',
		wait_for=[transcriptome_align]
	)
	# extracting for HLA isoforms.
	job.add_input_file(transcriptome_sam)
	hla_bam = extract_hla = job.grep_extract_hla(transcriptome_sam)
	job.add_output_file(hla_bam)

	# bam to fastq.
	hla_r1_fastq, hla_r2_fastq = job.bedtools_bamtofastq(hla_bam, paired=True)
	job.add_output_file(hla_r1_fastq)
	job.add_output_file(hla_r2_fastq)
	
	hla_extraction = job.jobname
	if not job.delete_sh:
		submit_commands.append(job.sge_submit_command())
		
		
	#### Job #3: Generate HLA types using PHLAT #### 
	job = RNAS_HLAJobScript(sample_name=sample_name,						 job_suffix='run_phlat',						 threads=8,						 memory=32,						 linkdir=os.path.join(linkdir, 'hla'),						 outdir=os.path.join(outdir, 'hla'),						 queue=queue,						 conda_env='hla',						 modules='',
						wait_for=[hla_extraction])
	phlat = job.phlat_typing(hla_r1_fastq, hla_r2_fastq)
	job.add_output_file(phlat)
	job.write_end()
	if not job.delete_sh:
		submit_commands.append(job.sge_submit_command()) 
		
		
	#### Job #4: Preprocessing for VBSeq, aligning to HLA fasta ####
	job = RNAS_HLAJobScript(sample_name=sample_name,
						job_suffix='pre_align_vbseq',
						threads=8,
						memory=32,
						linkdir=os.path.join(linkdir, 'reads'),
						outdir=os.path.join(outdir, 'reads'),
						queue=queue,
						conda_env='hla',
						modules='bwa,sambamba',
						wait_for=[hla_extraction])
	job.add_input_file(hla_r1_fastq)
	job.add_input_file(hla_r2_fastq)
	sam = job.bwa_pre_vbseq_align(hla_r1_fastq, hla_r2_fastq, 'imgt.hla', bamify=True)
	pre_align_vbseq = job.jobname
	if not job.delete_sh:
		submit_commands.append(job.sge_submit_command())

	#### Job #5: Running VBSeq ####
	job = RNAS_HLAJobScript(sample_name=sample_name,
						job_suffix='run_vbseq',
						threads=1,
						memory=32,
						linkdir=os.path.join(linkdir, 'hla'),
						outdir=os.path.join(outdir, 'hla'),
						queue=queue,
						conda_env='hla',
						modules='HLA-VBSeq',
						wait_for=[pre_align_vbseq])
	#job.add_temp_file(sam)
	vbseq = job.vbseq_typing(sam)
	job.add_output_file(vbseq)
	# Parsing vbseq 
	vbseq_parsed = job.parse_vbseq_results(vbseq) 
	job.add_output_file(vbseq_parsed)
	job.write_end()
	if not job.delete_sh:
		submit_commands.append(job.sge_submit_command())
		
	#### Job #5: Running seq2hla ####
	job = RNAS_HLAJobScript(sample_name=sample_name,
						job_suffix='run_seq2hla',
						threads=8,
						memory=32,
						linkdir=os.path.join(linkdir, 'hla'),
						outdir=os.path.join(outdir, 'hla'),
						queue=queue,
						conda_env='seq2hla',
						modules='seq2hla',
						wait_for=[hla_extraction])
	make_dir(os.path.join(job.outdir, 'seq2hla'))
	seq2hla_fns = job.seq2hla_typing(hla_r1_fastq, hla_r2_fastq)
	for fn in seq2hla_fns:
		job.add_output_file(fn)
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

