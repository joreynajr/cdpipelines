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
	
	def samtools_extract_regions(self, in_bam, regions):
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
		extract_bam : str
			Path to bam with extractions file.
	
		"""		
		# GRCH37 HLA regions 
		# Make reads directory for extraction of HLA R1 and R2 fastq files 
		#### Deprecated Sambamba version has some bug. ####
		# Version 0.5.9 and 0.6.1 produce output but switching to samtools to be sure
		# with open(self.filename, 'a') as f:
			# f.write('sambamba view -h -f bam {} {} -o {} -t {}'.format(bam, hla_regions, extraction, self.threads))
		extract_bam = os.path.join(self.outdir,'{}_hla.bam'.format(self.sample_name))
		#### Samtools version 
		lines = []
		cmd = 'samtools view -h -b -@ {} -o {} {} {}'.format(self.threads, extract_bam, in_bam, regions)
		lines = self._add_execution_date(lines, cmd)

		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))
		return extract_bam

	def samtools_extract_bed(self, in_bam, bed_fn):
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
		extract_bam : str
			Path to bam with extractions file.
	
		"""		
		# GRCH37 HLA regions 
		# Make reads directory for extraction of HLA R1 and R2 fastq files 
		#### Deprecated Sambamba version has some bug. ####
		# Version 0.5.9 and 0.6.1 produce output but switching to samtools to be sure
		# with open(self.filename, 'a') as f:
			# f.write('sambamba view -h -f bam {} {} -o {} -t {}'.format(bam, hla_regions, extraction, self.threads))
		extract_bam = os.path.join(self.outdir,'{}_hla.bam'.format(self.sample_name))
		#### Samtools version 
		lines = []
		cmd = 'samtools view -h -b -@ {} -L {} -o {} {}'.format(self.threads, bed_fn, extract_bam, in_bam)
		lines = self._add_execution_date(lines, cmd)

		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))
		return extract_bam
	
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
	
	def bwa_align(
		self,
		r1_fastq,
		r2_fastq, 
		hla_ref,
	):  
		"""
		Align paired fastq files with bwa. Uses the -a option 
		which produces all the alignments, not just the primary 
		alignment.
		
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
		sam = os.path.join(self.outdir, '{}_bwa.sam'.format(self.sample_name))
		lines = []
		cmd = 'bwa mem -t {} -L 10000 -P -a {} {} {} > {}'.\
			format(self.threads,
			hla_ref,
			r1_fastq,
			r2_fastq,
			sam)
		lines = self._add_execution_date(lines, cmd)
		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))
		return sam	 


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
		cmd = 'python -O {} -1 {} -2 {} -index {} -p {} -b2url {} -orientation "--fr" -tag {} -e {} -o {}'.\
			format(os.path.join(phlat_dir, 'dist/PHLAT.py'),
			r1_fastq,
			r2_fastq,
			os.path.join(phlat_dir, 'b2folder'),
			self.threads, 
			bowtie2, self.sample_name, phlat_dir, self.outdir)
		lines = self._add_execution_date(lines, cmd)
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
		lines = []
		cmd = 'java -Xmx27G -jar $HLA_VBSeq {} {} {} --alpha_zero 0.01 --is_paired'.format(hla_ref, sam, vbseq_fn)
		lines = self._add_execution_date(lines, cmd)
		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))
		return vbseq_fn 
	
	def parse_vbseq_results(
		self,
		hla,
		hla_allele_list, 
		email=False,
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
		lines = []
		cmd = 'perl /software/HLA-VBseq_d16.06.14/parse_result.pl {} {} > {}'.format(hla_allele_list, hla, vbseq_result_fn)
		lines = self._add_execution_date(lines, cmd)
		with open(self.filename, 'a') as f:
			f.write('\n'.join(lines))
			if email == True:
				f.write('echo Hello Mars! | mail -r joreyna@flh1.ucsd.edu -s "Jobs are complete." joreyna@live.com\n')
		return vbseq_result_fn  
	
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
		hla_regions,
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
	hla_regions; str,
	bowtie2 : str
	phlat_dir : str
	queue : str

	"""
	
	submit_commands = [] 

	#####PHLAT jobs#####
	#### Job #1: Extracting HLA Regions into R1 and R2 fastq's ####
	job = WGS_HLAJobScript(sample_name=sample_name,
		job_suffix=JOBNAMES[0],
		threads=8,
		memory=32,
		linkdir=os.path.join(linkdir),
		outdir=os.path.join(outdir, 'reads'),
		queue=queue,
		conda_env='hla',
		modules='sambamba/0.6.1,samtools/1.2,bedtools/2.25.0')
	job.add_input_file(in_bam)
	index_bam = job.sambamba_index(in_bam)
	job.add_temp_file(index_bam)

	# Extract HLA Regions 
	extract_bam = job.samtools_extract_bed(in_bam, hla_regions)
	job.add_output_file(extract_bam)

	# Sort bam by query (aka read name) 
	qsort_bam = job.sambamba_sort(extract_bam, queryname=True)
	job.add_output_file(qsort_bam)

	# Generate R1 and R2 Fastq files 
	(fastq_r1, fastq_r2) = job.bedtools_bamtofastq(qsort_bam, paired=True)
	fastq_r1 = job.add_output_file(fastq_r1)
	fastq_r2 = job.add_output_file(fastq_r2)
	job.write_end()
	extract_job = job.jobname

	if not job.delete_sh:
		submit_commands.append(job.sge_submit_command())

	#### Job #2: Generate HLA types using PHLAT #### 
	job = WGS_HLAJobScript(sample_name=sample_name,
		job_suffix=JOBNAMES[1],
		threads=8,
		memory=32,
		linkdir=os.path.join(linkdir),
		outdir=os.path.join(outdir, 'hla'),
		queue=queue,
		conda_env='hla',
		modules='',
		wait_for=[extract_job])
	phlat = job.phlat_typing(fastq_r1, fastq_r2)
	job.add_output_file(phlat)
	job.write_end()
	if not job.delete_sh:
		submit_commands.append(job.sge_submit_command()) 
		
	### Job #3: Preprocessing for VBSeq, aligning to HLA fasta ####
	job = WGS_HLAJobScript(sample_name=sample_name,
		job_suffix='pre_align_vbseq',
		threads=8,
		memory=32,
		linkdir=os.path.join(linkdir),
		outdir=os.path.join(outdir, 'reads'),
		queue=queue,
		conda_env='hla',
		modules='bwa',
		wait_for=[extract_job])
	job.add_input_file(fastq_r1)
	job.add_input_file(fastq_r2)
	sam = job.bwa_align(fastq_r1, fastq_r2, hla_ref)
	pre_align_vbseq = job.jobname
	if not job.delete_sh:
		submit_commands.append(job.sge_submit_command())
		
	#### Job #4: Running VBSeq and parsing results####
	job = WGS_HLAJobScript(sample_name=sample_name,
		job_suffix='run_vbseq',
		threads=1,
		memory=32,
		linkdir=os.path.join(linkdir),
		outdir=os.path.join(outdir, 'hla'),
		queue=queue,
		conda_env='hla',
		modules='HLA-VBSeq',
		wait_for=[pre_align_vbseq])
	job.add_temp_file(sam)
	# running VBSeq
	vbseq = job.vbseq_typing(sam, hla_ref)
	job.add_temp_file(vbseq)
	# parsing results
	vbseq_parsed = job.parse_vbseq_results(vbseq, hla_allele_list, email)
	job.add_output_file(vbseq_parsed)
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
