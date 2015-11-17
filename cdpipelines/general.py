import os
import subprocess

def _git_info():
    """Get current git version"""
    # Necessary to pipe to cat because git will open the result in less
    # otherwise.
    d = os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[0:-2] +
                         ['.git'])
    command = ('git --git-dir {0} log -1 --pretty=oneline --decorate | '
               'cat'.format(d))
    res = subprocess.Popen(command, shell=True,
                           stdout=subprocess.PIPE).communicate()
    res = res[0][:-1].strip() 
    return res 

def _make_dir(d):
    """Make directory d if it doesn't exist"""
    try:
        os.makedirs(d)
    except OSError:
        pass

class JobScript:
    def __init__(self, sample_name, job_suffix, outdir, threads, memory,
                 shell_fn=None, tempdir=None, queue=None, conda_env=None,
                 modules=None, wait_for=None, copy_input=False):
        """
        Create SGE/shell script object.

        Parameters
        ----------
        sample_name : str
            Sample name used for naming directories, files, etc.

        job_suffix : str
            This suffix will be used for naming directories, files, etc.

        outdir: str
            Path to directory where final output files should be stored.

        threads : int
            Number of threads to request for SGE scripts and to use for
            multi-threaded software.

        memory : int
            Amount of memory in Gb to request for SGE scripts.

        shell_fn : str
            Path to output shell script. If not provided, the script will be
            written in the output directory.

        tempdir : str
            Path to directory where temporary directory should be made. If not
            provided, the output directory will be used as the temp directory.

        queue : str
            SGE queue to use if writing SGE script. If not provided, jobs will
            go into the default week queue.

        conda_env : str
            Path to conda environment to load when job begins.

        modules : str
            Modules (separated by commas e.g. bedtools,samtools) to load at
            beginning of script.

        wait_for : list
            A list of jobnames to wait for before starting this job. This is
            accomplished using -hold_jid.

        copy_input : bool
            Whether to copy input files to temp directory. 

        """
        self.sample_name = sample_name
        self.job_suffix = job_suffix
        self.jobname = '{}_{}'.format(sample_name, job_suffix)
        self.outdir = outdir
        if tempdir:
            self.tempdir = os.path.realpath(os.path.join(tempdir, self.jobname))
        else:
            self.tempdir = self.outdir
        _make_dir(self.outdir)
        assert type(threads) is int
        self.threads = threads
        assert type(memory) is int
        self.memory = memory
        self.queue = queue
        self.conda_env = conda_env
        if modules:
            self.modules = modules.split(',')
        else:
            self.modules = None
       
        _make_dir(os.path.join(os.path.split(outdir)[0], 'logs'))
        self.out = os.path.join(os.path.split(self.outdir)[0], 'logs',
                                '{}.out'.format(self.jobname))
        self.err = os.path.join(os.path.split(self.outdir)[0], 'logs',
                                '{}.err'.format(self.jobname))
        self.wait_for = wait_for
        self.copy_input = copy_input
        self.input_files_to_copy = []
        self.output_files_to_copy = []
        self.temp_files_to_delete = []
        self.softlinks = []
        self._set_filename()
        self._write_header()

    def sge_submit_command(self):
        """Get command to submit script."""
        if self.wait_for:
            return 'qsub -hold_jid {} {}'.format(','.join(self.wait_for),
                                                 self.filename)
        else:
            return 'qsub {}'.format(self.filename)

    def _set_filename(self):
        """Make SGE/shell script filename."""
        if shell_fn:
            self.filename = shell_fn
        else:
            _make_dir(os.path.join(os.path.split(self.outdir)[0], 'sh'))
            self.filename = os.path.join(os.path.split(self.outdir)[0], 'sh',
                                         '{}.sh'.format(self.jobname))
    
    def _write_header(self):
        with open(self.filename, "a") as f:
            f.write('#!/bin/bash\n\n')
            if self.queue:
                assert self.queue in ['short', 'long']
                f.write('#$ -q {}\n'.format(self.queue))
            f.write('#$ -N {}\n'.format(self.jobname))
            f.write('#$ -l h_vmem={}\n'.format(
                self.memory / float(self.threads)))
            f.write('#$ -pe smp {}\n'.format(self.threads))
            f.write('#$ -S /bin/bash\n')
            f.write('#$ -o {}\n'.format(self.out))
            f.write('#$ -e {}\n\n'.format(self.err))
            f.write('# Git repository version:\n# {}\n\n'.format(_git_info()))
            if self.modules:
                for module in self.modules:
                    f.write('module load {}\n\n'.format(module))
            if self.conda_env:
                f.write('source activate {}\n\n'.format(self.conda_env))
            if self.tempdir:
                f.write('mkdir -p {}\n'.format(self.tempdir))
                f.write('cd {}\n\n'.format(self.tempdir))

    def add_softlink(self, target, link):
        """Add a target, link pair to the JobScript instance so that a softlink
        from target to link is made at the end of the jobscript. This happens at
        the end of the script so that the files actually exist and the softlink
        works. Useful for linking files into web directories for instance."""
        self.softlinks.append([target, link])

    def add_temp_file(self, fn, copy=False):
        """Add temporary file to self.temp_files_to_delete. If copy == True,
        add to self.output_files_to_copy. Returns temp path."""
        if copy:
            self.output_files_to_copy.append(fn)
        return self.temp_file_path(fn)

    def add_input_file(self, fn, copy=None):
        """Add input file to self.input_files_to_copy if self.copy_input or
        copy == True. Return temp path. If both self.copy_input and copy are
        False, the realpath of the input file is returned."""
        if copy:
            self.input_files_to_copy.append(fn)
        elif copy == False:
            pass 
        else:
            if self.copy_input:
                self.input_files_to_copy.append(fn)
                copy = True
        if copy:
            return self.temp_file_path(fn)
        else:
            return os.path.realpath(fn)

    def temp_file_path(self, fn):
        """Return the path to the temporary version of an input file"""
        return os.path.join(self.tempdir, os.path.split(fn)[1])

    def copy_input_files(self):
        if len(self.input_files_to_copy) > 0:
            with open(self.filename, "a") as f:
                f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format( 
                    '\\\n\t'.join(self.input_files_to_copy),
                    self.tempdir))
            # We will delete any input files we copy over.
            self.temp_files_to_delete += [
                os.path.join(self.tempdir, os.path.split(x)[1]) for x in
                self.input_files_to_copy]

    def _copy_output_files(self):
        if (len(self.output_files_to_copy) > 0 and 
            os.path.realpath(self.tempdir) != os.path.realpath(self.outdir)):
            with open(self.filename, "a") as f:
                f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format( 
                    '\\\n\t'.join(self.output_files_to_copy),
                    self.outdir))

    def _delete_temp_files(self):
        if len(self.temp_files_to_delete) > 0:
            if (os.path.realpath(self.tempdir) == os.path.realpath(self.outdir) 
                or self.tempdir is None):
                self.temp_files_to_delete = [
                    x for x in self.temp_files_to_delete if x not in
                    self.output_files_to_copy
                ]
            with open(self.filename, "a") as f:
                f.write('rm -r \\\n\t{}\n\n'.format(
                    ' \\\n\t'.join(self.temp_files_to_delete)))

    def _delete_tempdir(self):
        if self.tempdir and (os.path.realpath(self.tempdir) !=
                             os.path.realpath(self.outdir)):
            with open(self.filename, "a") as f:
                f.write('rm -r {}\n'.format(self.tempdir))

    def _make_softlinks(self):
        with open(self.filename, "a") as f:
            for p in self.softlinks:
                f.write(_softlink(p[0], p[1]))

    def write_end(self):
        self._copy_output_files()
        self._delete_temp_files()
        self._delete_tempdir()
        self._make_softlinks()

    def picard_collect_rna_seq_metrics(
        self,
        in_bam, 
        ref_flat, 
        rrna_intervals,
        picard_path='$picard',
        picard_memory=2, 
        picard_tempdir='.',
        strand_specific=True, 
        bg=False,
    ):
        """
        Collect RNA-seq metrics using Picard. The input bam file is assumed to
        be coordinate sorted.
    
        Parameters
        ----------
        in_bam : str
            Path to input coordinate sorted bam file.
    
        ref_flat : str
            Path to refFlat file with non-rRNA genes. Can be gzipped.
    
        rrna_intervals : str
            Pato to interval list file with rRNA intervals.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        picard_memory : int
            Amount of memory in Gb to give picard.
    
        picard_tempdir : str
            Path to directory to use for Picard temporary files. Default is
            current directory.
    
        bg : boolean
            Whether to run the process in the background.

        Returns
        -------
        metrics : str
            Path to output metrics file.
    
        chart : str
            Path to output PDF file.
    
        """
        metrics = os.path.join(self.tempdir,
                               '{}_rna_seq_metrics.txt'.format(self.sample_name))
        chart = os.path.join(self.tempdir,
                             '{}_5_3_coverage.pdf'.format(self.sample_name))
        if strand_specific:
            ss = 'SECOND_READ_TRANSCRIPTION_STRAND'
        else:
            ss = 'NONE'
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(picard_memory),
            '-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
            '-Djava.io.tmpdir={}'.format(picard_tempdir), 
            '-jar {} CollectRnaSeqMetrics'.format( picard_path),
            'I={}'.format(in_bam),
            'REF_FLAT={}'.format(ref_flat),
            'STRAND_SPECIFICITY={}'.format(ss),
            'RIBOSOMAL_INTERVALS={}'.format(rrna_intervals),
            'ASSUME_SORTED=TRUE',
            'CHART_OUTPUT={}'.format(chart),
            'O={}'.format(metrics)]))
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return metrics, chart

    def picard_collect_multiple_metrics(
        self,
        in_bam, 
        picard_path='$picard',
        picard_memory=2, 
        picard_tempdir='.',
        bg=False,
    ):
        """
        Collect multiple metrics using Picard. The input bam file is assumed to
        be sorted.
    
        Parameters
        ----------
        in_bam : str
            Path to input coordinate sorted bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        picard_memory : int
            Amount of memory in Gb to give picard.
    
        picard_tempdir : str
            Path to directory to use for Picard temporary files. Default is
            current directory.
    
        bg : boolean
            Whether to run the process in the background.

        Returns
        -------
        output : tuple
            Tuple of the paths to the following output files:
                alignment_summary_metrics quality_by_cycle.pdf
                base_distribution_by_cycle.pdf quality_by_cycle_metrics
                base_distribution_by_cycle_metrics quality_distribution.pdf
                insert_size_histogram.pdf quality_distribution_metrics
                insert_size_metrics.
    
        """
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(picard_memory),
            '-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
            '-Djava.io.tmpdir={}'.format(picard_tempdir), 
            '-jar {} CollectMultipleMetrics'.format(picard_path),
            'VALIDATION_STRINGENCY=SILENT',
            'ASSUME_SORTED=TRUE',
            'I={}'.format(in_bam), 
            'O={}'.format(self.sample_name)]))
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        output = [os.path.join(self.tempdir, '{}.{}'.format(self.sample_name, x))
                               for x in [
                                   'alignment_summary_metrics',
                                   'quality_by_cycle.pdf',
                                   'base_distribution_by_cycle.pdf',
                                   'quality_by_cycle_metrics',
                                   'base_distribution_by_cycle_metrics',
                                   'quality_distribution.pdf',
                                   'insert_size_histogram.pdf',
                                   'quality_distribution_metrics',
                                   'insert_size_metrics']]
        return tuple(output)
    
    def picard_index(
        self,
        in_bam, 
        picard_path='$picard',
        picard_memory=2, 
        picard_tempdir='.',
        bg=False,
    ):
        """
        Index bam file using Picard Tools.
    
        Parameters
        ----------
        in_bam : str
            Path to file input bam file.
    
        index : str
            Path to index file for input bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        picard_memory : int
            Amount of memory in Gb to give picard.
    
        picard_tempdir : str
            Path to directory to use for Picard temporary files. Default is
            current directory.
    
        bg : boolean
            Whether to run the process in the background.
    
        Returns
        -------
        index : str
            Path to output index file.
    
        """
        index = os.path.join(self.tempdir, os.path.splitext(in_bam)[0] + '.bai')
        line = (' \\\n\t'.join([
            'java -Xmx{}g -jar'.format(picard_memory),
            '-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
            '-Djava.io.tmpdir={}'.format(picard_tempdir), 
            '-jar {} BuildBamIndex'.format(picard_path),
            'I={}'.format(in_bam),
            'O={}'.format(index)]))
        if bg:
            line += ' &\n\n'
        else:
            line += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return index
    
    def picard_merge(
        self,
        bams, 
        out_bam, 
        picard_path='$picard',
        picard_memory=2, 
        picard_tempdir='.',
        bg=False,
    ):
        """
        Merge bam files using Picard. Input bam files are assumed to be
        coordinate sorted.
    
        Parameters
        ----------
        bams : str
            Bam files to merge.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        picard_memory : int
            Amount of memory in Gb to give picard.
    
        picard_tempdir : str
            Path to directory to use for Picard temporary files. Default is
            current directory.
    
        bg : boolean
            Whether to run the process in the background.
    
        Returns
        -------
        out_bam : str
            Path to output merged bam file.
    
        """
        merge_in = ''.join(['\tI={} \\\n'.format(x) for x in bams])
        lines = ['java -Xmx{}g -jar'.format(picard_memory),
                 '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                 '\t-Djava.io.tmpdir={}'.format(picard_tempdir), 
                 '\t-jar {} MergeSamFiles'.format(picard_path),
                 '\tASSUME_SORTED=TRUE',
                 '\tUSE_THREADING=TRUE']
        for bam in bams:
            lines.append('\tI={}'.format(bam))
        lines.append('\tO={}'.format(out_bam))
        lines = (' \\\n'.join(lines))
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return out_bam
    
    def samtools_index(
        self,
        in_bam, 
        bg=False,
        samtools_path='samtools',
    ):
        """
        Index bam file using samtools.
    
        Parameters
        ----------
        in_bam : str
            Path to file input bam file.
    
        index : str
            Path to index file to be written. If not provided, the index is
            written to the samtools default {in_bam}.bai in the current working
            directory.
    
        Returns
        -------
        index : str
            Path to index file for input bam file.
    
        """
        index = os.path.join(self.tempdir, os.path.splitext(in_bam)[0] + '.bai')
        line = '{} index {}'.format(samtools_path, in_bam)
        if bg:
            line += ' &\n\n'
        else:
            line += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return index
    
    def picard_mark_duplicates(
        self,
        in_bam, 
        picard_path='$picard',
        picard_memory=2, 
        picard_tempdir='.',
        remove_dups=False,
    ):
        """
        Mark and optionally remove duplicates using Picard Tools.
    
        Parameters
        ----------
        in_bam : str
            Path to input bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        picard_memory : int
            Amount of memory in Gb to give picard.
    
        picard_tempdir : str
            Path to directory to use for Picard temporary files. Default is
            current directory.

        Returns
        -------
        out_bam : str
            Path to output bam file.
    
        duplicate_metrics : str
            Path to index file for input bam file.
    
    
        """
        mdup_bam = os.path.join(
            self.tempdir, '{}_sorted_mdup.bam'.format(self.sample_name))
        duplicate_metrics = os.path.join(
            self.outdir, '{}_duplicate_metrics.txt'.format(self.sample_name))
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(picard_memory),
            '-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
            '-Djava.io.tmpdir={}'.format(picard_tempdir), 
            '-jar {} MarkDuplicates'.format(picard_path),
            'METRICS_FILE={}'.format(duplicate_metrics),
            'VALIDATION_STRINGENCY=SILENT',
            'ASSUME_SORTED=TRUE',
            'I={}'.format(in_bam), 
            'O={}'.format(out_bam)]))
        if remove_dups:
            lines += ' \\\n\tREMOVE_DUPLICATES=TRUE\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return mdup_bam, duplicate_metrics

    def picard_gc_bias_metrics(
        self,
        in_bam, 
        picard_path='$picard',
        picard_memory=2, 
        picard_tempdir='.',
        bg=False,
    ):
        """
        Collect GC bias metrics using Picard. The input bam file is assumed to
        be sorted.
    
        Parameters
        ----------
        in_bam : str
            Path to input bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        picard_memory : int
            Amount of memory in Gb to give picard.
    
        picard_tempdir : str
            Path to directory to use for Picard temporary files. Default is
            current directory.
    
        bg : boolean
            Whether to run the process in the background.

        Returns
        -------
        metrics : str
            Path to output metrics file.
    
        chart : str
            Path to output PDF chart.
    
        out : str
            Path to picard output file.
    
        """
        metrics = os.path.join(self.tempdir,
                               '{}_gc_bias_metrics.txt'.format(self.sample_name))
        chart = os.path.join(self.tempdir,
                             '{}_gc_bias.pdf'.format(self.sample_name))
        out = os.path.join(self.tempdir,
                           '{}_gc_bias_metrics_out.txt'.format(self.sample_name))
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(picard_memory),
            '-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
            '-Djava.io.tmpdir={}'.format(picard_tempdir), 
            '-jar {} CollectGcBiasMetrics'.format(picard_path),
            'VALIDATION_STRINGENCY=SILENT',
            'I={}'.format(in_bam), 
            'O={}'.format(out),
            'CHART_OUTPUT={}'.format(chart),
            'SUMMARY_OUTPUT={}'.format(metrics),
            'ASSUME_SORTED=TRUE']))
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return metrics, chart, out
    
    def picard_bam_index_stats(
        self,
        in_bam, 
        picard_path='$picard',
        picard_memory=2, 
        picard_tempdir='.',
        bg=False):
        """
        Collect bam index stats with Picard.
    
        Parameters
        ----------
        in_bam : str
            Path to input bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        picard_memory : int
            Amount of memory in Gb to give picard.
    
        picard_tempdir : str
            Path to directory to use for Picard temporary files. Default is
            current directory.
    
        bg : boolean
            Whether to run the process in the background.

        Returns
        -------
        out : str
            Path to write stdout to. This contains the index stats.
    
        err : str
            Path to write stderr to. This contains picard stuff.
    
        """
        out = os.path.join(self.outdir,
                           '{}_index_stats.txt'.format(self.sample_name))
        err = os.path.join(self.outdir,
                           '{}_index_stats.err'.format(self.sample_name))
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(picard_memory),
            '-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
            '-Djava.io.tmpdir={}'.format(picard_tempdir), 
            '-jar {} BamIndexStats'.format(picard_path),
            'VALIDATION_STRINGENCY=SILENT',
            'I={}'.format(in_bam),
            '> {}'.format(out),
            '2> {}'.format(err)]))
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return out, err
    
    def picard_insert_size_metrics(
        self,
        in_bam, 
        picard_path='$picard',
        picard_memory=2, 
        picard_tempdir='.',
        bg=False,
    ):
        """
        Collect insert size metrics using Picard. The input bam file is assumed
        to be sorted.
    
        Parameters
        ----------
        in_bam : str
            Path to input bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        picard_memory : int
            Amount of memory in Gb to give picard.
    
        picard_tempdir : str
            Path to directory to use for Picard temporary files. Default is
            current directory.
    
        bg : boolean
            Whether to run the process in the background.
   
        Returns
        -------
        metrics : str
            Path to output metrics file.
    
        hist : str
            Path to output histogram PDF.
    
        """
        metrics = os.path.join(
            self.tempdir, '{}_insert_size_metrics.txt'.format(self.sample_name))
        hist = os.path.join(self.tempdir,
                            '{}_insert_size.pdf'.format(self.sample_name))
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(picard_memory),
            '-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
            '-Djava.io.tmpdir={}'.format(picard_tempdir), 
            '-jar {} CollectInsertSizeMetrics'.format(picard_path),
            'VALIDATION_STRINGENCY=SILENT',
            'I={}'.format(in_bam), 
            'O={}'.format(metrics),
            'HISTOGRAM_FILE={}'.format(hist),
            'ASSUME_SORTED=TRUE']))
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return metrics, hist
    
    def picard_query_sort(
        in_bam, 
        picard_path='$picard',
        picard_memory=2, 
        picard_tempdir='.',
        bg=False,
    ):
        """
        Query sort using Picard Tools.
    
        Parameters
        ----------
        in_bam : str
            Path to input bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        picard_memory : int
            Amount of memory in Gb to give picard.
    
        picard_tempdir : str
            Path to directory to use for Picard temporary files. Default is
            current directory.
    
        bg : boolean
            Whether to run the process in the background.
    
        Returns
        -------
        out_bam : str
            Path to output bam file.
    
        """
        out_bam = os.path.join(self.tempdir,
                               '{}_qsorted.bam'.format(self.sample_name))
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(picard_memory),
            '-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
            '-Djava.io.tmpdir={}'.format(picard_tempdir), 
            '-jar {} SortSam'.format(picard_path),
            'VALIDATION_STRINGENCY=SILENT',
            'I={}'.format(in_bam), 
            'O={}'.format(out_bam),
            'SO=queryname']))
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return out_bam
    
    def picard_coord_sort(
        self,
        in_bam, 
        index=False,
        picard_path='$picard',
        picard_memory=2, 
        picard_tempdir='.',
    ):
        """
        Coordinate sort using Picard Tools.
    
        Parameters
        ----------
        in_bam : str
            Path to input bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        picard_memory : int
            Amount of memory in Gb to give picard.
    
        picard_tempdir : str
            Path to directory to use for Picard temporary files. Default is
            current directory.
    
        Returns
        -------
        out_bam : str
            Path to output bam file.
    
        out_index : str
            Path to output index file. Only returned if index == True.
    
        """
        out_bam = os.path.join(self.tempdir,
                               '{}_sorted.bam'.format(self.sample_name))
        if index:
            out_index = os.path.join(
                self.tempdir, 
                os.path.join(self.tempdir, os.path.splitext(in_bam)[0] +
                             '.bai'))
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(picard_memory),
            '-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
            '-Djava.io.tmpdir={}'.format(picard_tempdir), 
            '-jar {} SortSam'.format(picard_path),
            'VALIDATION_STRINGENCY=SILENT',
            'I={}'.format(in_bam), 
            'O={}'.format(out_bam),
            'SO=coordinate']))
        if index:
            lines += ' \\\n\tCREATE_INDEX=TRUE' 
            old_index = '.'.join(out_bam.split('.')[0:-1]) + '.bai'
            lines += 'mv {} {}\n\n'.format(old_index, out_index)

        with open(job.filename, "a") as f:
            f.write(lines)
        if index:
            return out_bam, out_index
        else: 
            return out_bam
    
    def cutadapt_trim(
        self,
        fastq, 
        length, 
        out, 
        bg=False,
    ):
        """
        Cut a specified number of bases from a fastq file using cutadapt.
        Cutadapt should be installed in your Python environment.
    
        Parameters
        ----------
        fastq : str
            Fastq or gzipped/bzipped fastq.
    
        length : int
            Positive or negative integer. Positive numbers remove bases at the
            front of the read and negative numbers remove bases at the end of
            the read.
    
        out : str
            Path to output (optionally gzipped/bzipped) fastq files.
    
        bg : boolean
            Whether to run the process in the background (i.e. include an
            ampersand at the end of the command).
    
        Returns
        -------
        out : str
            Path to output (optionally gzipped/bzipped) fastq files.
    
        """
        line = 'cutadapt --cut {} -o {} {}'.format(length, out, fastq)
        if bg:
            line += ' &\n\n'
        else:
            line += '\n\n'
        with open(job.filename, "a") as f:
            f.write(lines)
        return out
    
    def bedgraph_to_bigwig(
        self,
        bedgraph, 
        bedgraph_to_bigwig_path='bedGraphToBigWig',
        bedtools_path='bedtools',
    ):
        """
        Convert bedgraph file to bigwig file.
    
        Parameters
        ----------
        bedgraph : str
            Input bedgraph file.
    
        bedgraph_to_bigwig_path : str
            Path bedGraphToBigWig executable from UCSC.

        bedtools_path : str
            Path to bedtools. If bedtools_path == 'bedtools', it is assumed that
            the hg19 human.hg19.genome file from bedtools is also in your path.

        Returns
        -------
        bigwig : str
            Path to output bigwig file.
    
        """
        bigwig = os.path.join(
            job.tempdir,
            '{}.bw'.format(os.path.splitext(os.path.split(bedgraph)[1])))
        # If bedtools is in the path, I'll assume the genome file is as well.
        if bedtools_path == 'bedtools':
            bedtools_genome_path = 'human.hg19.genome'
        else:
            bedtools_genome_path = os.path.join(
                os.path.split(os.path.split(bedtools_path)[0])[0], 'genomes',
                'human.hg19.genome')
        lines =  ' \\\n\t'.join([
            '{} {}'.format(bedgraph_to_bigwig_path, bedgraph),
            '{}'.format(bedtools_genome_path),
            '{} &\n'.format(bigwig)])
        with open(job.filename, "a") as f:
            f.write(lines)
        return bigwig
    
    def flagstat(
        self,
        bam, 
        samtools_path='samtools',
        bg=False):
        """
        Run flagstat for a bam file.
    
        Parameters
        ----------
        bam : str
            Bam file to calculate coverage for.
    
        samtools_path : str
            Path to samtools executable.
    
        Returns
        stats_file : str
            File to write flagstats to.
    
        -------
        """
        stats_file = os.path.join(
            job.tempdir, '{}_flagstat.txt'.format(
                os.path.splitext(os.path.split(bam)[1])[0]))
        lines = '{} flagstat {} > {}'.format(samtools_path, bam, stats_file)
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(job.filename, "a") as f:
            f.write(lines)
        return stats_file
    
    def coverage_bedgraph(
        self,
        bam, 
        strand='.',
        bedtools_path='bedtools',
    ):
        """
        Make lines that create a coverage bedgraph file.
    
        Parameters
        ----------
        bam : str
            Bam file to calculate coverage for.
    
        bedtools_path : str
            Path to bedtools.
    
        sample_name : str
            Sample name for naming files etc.
    
        strand : str
            If '+' or '-', calculate strand-specific coverage. Otherwise,
            calculate coverage using all reads.
    
        Returns
        -------
        bedgraph : str
            Path to output bedgraph file.
    
        """
        bedgraph = os.path.join(
            self.tempdir,
            '{}.bg'.format(os.path.splitext(os.path.split(bam)[1])[0]))
        if strand == '+' or strand == '-':
            if strand == '+':
                name = '{}_plus'.format(self.sample_name)
                bedgraph = '{}_plus.bg'.format(os.path.splitext(bedgraph)[0])
            else:
                name = '{}_minus'.format(self.sample_name)
                bedgraph = '{}_minus.bg'.format(os.path.splitext(bedgraph)[0])
            lines = ' \\\n\t'.join([
                '{} genomecov -ibam'.format(bedtools_path),
                '{}'.format(bam),
                '-g hg19.genome -split -bg ',
                '-strand {} -trackline'.format(strand),
                '-trackopts \'name="{}"\''.format(name),
                '> {} &\n\n'.format(bedgraph)])
        else:
            name = job.sample_name
            lines = ' \\\n\t'.join([
                '{} genomecov -ibam'.format(bedtools_path),
                '{}'.format(bam),
                '-g hg19.genome -split -bg ',
                '-trackline'.format(strand),
                '-trackopts \'name="{}"\''.format(name),
                '> {} &\n\n'.format(bedgraph)])
        with open(job.filename, "a") as f:
            f.write(lines)
        return bedgraph
    
    def bigwig_files(
        self,
        in_bam, 
        out_bigwig, 
        sample_name, 
        out_bigwig_minus='',
        bedgraph_to_bigwig_path='bedGraphToBigWig',
        bedtools_path='bedtools',
    ):
        """
        Make bigwig coverage files.
    
        Parameters
        ----------
        in_bam : str
            Path to bam file to create bigwigs for.
    
        out_bigwig : str
            Path to output bigwig file. If out_bigwig_minus is provided,
            out_bigwig has the plus strand coverage.
    
        out_bigwig_minus : str
            Path to output bigwig file for minus strand. If out_bigwig_minus is
            not provided, the coverage is calculated using reads from both
            strands and written to out_bigwig.
    
        Returns
        -------
        lines : str
            Lines to be printed to shell script.
    
        """
        lines = ''
        if out_bigwig_minus != '':
            lines += _coverage_bedgraph(in_bam, 'plus.bg', sample_name,
                                        strand='+')
            lines += _coverage_bedgraph(in_bam, 'minus.bg', sample_name,
                                        strand='-')
            lines += ('wait\n\n')
            lines += (_bedgraph_to_bigwig('plus.bg', out_bigwig))
            lines += (_bedgraph_to_bigwig('minus.bg', out_bigwig_minus))
                                          
            lines += ('\nwait\n\n')
            lines += ('rm plus.bg minus.bg\n\n')
        
        else:
            lines = _coverage_bedgraph(in_bam, 'both.bg', sample_name)
            lines += ('wait\n\n')
            lines += (_bedgraph_to_bigwig('both.bg', out_bigwig))
            lines += ('wait\n\n')
            lines += ('rm both.bg\n\n')
        with open(job.filename, "a") as f:
            f.write(lines)
            if out_bigwig_minus:
                return  
        
    def combine_fastqs(
        self,
        fastqs,
        out_fastq,
        bg=False,
    ):
        """
        If fastqs is a string with a path to a single fastq file, make a
        softlink to out_fastq. If fastqs is a list of fastqs, cat the fastqs
        together into the file out_fastq.
    
        Parameters
        ----------
        fastqs : list or str
            Either a list of paths to gzipped fastq files or path to a single
            gzipped fastq file.
    
        out_fastq : str
            Path to single output fastq file.
    
        Returns
        -------
        out_fastq : str
            Path to single output fastq file.
    
        """
        fastqs = sorted(fastqs)
        if type(fastqs) == list:
            lines = 'cat \\\n\t{} \\\n\t> {}'.format(' \\\n\t'.join(fastqs),
                                                       out_fastq)
        elif type(fastqs) == str:
            lines = 'ln -s {} {}'.format(fastqs, out_fastq)
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(job.filename, "a") as f:
            f.write(lines)
        return out_fastq
    
    def fastqc(
        self,
        fastqs, 
        outdir, 
        threads=1,
        fastqc_path='fastqc',
    ):
        """
        Run FastQC
    
        Parameters
        ----------
        fastqs : str or list
            Path to fastq file or list of paths to fastq files.
    
        outdir : str
            Path to directory to store FastQC results to.
    
        threads : int
            Number of threads to run FastQC with.
    
        fastqc_path : str
            Path to FastQC.
    
        Returns
        -------
    
        """
        if type(fastqs) == list:
            fastqs = ' \\\n\t'.join(fastqs)
        lines = ('{} --outdir {} --nogroup --extract --threads {} \\\n'
                 '\t{}\n'.format(fastqc_path, outdir, threads, fastqs))
        with open(job.filename, "a") as f:
            f.write(lines)
        return None # I should probably figure out what fastQC outputs and
                    # provide links here
    
    def softlink(self, target, link):
        """
        Make softlink from target to link.
    
        Parameters
        ----------
        target : str
            Full path to file to make link to.
    
        link : str
            Full path to softlink.
    
        Returns
        -------
        link : str
            Full path to softlink.
    
        """
        lines = 'ln -s {} {}\n\n'.format(target, link)
        with open(job.filename, "a") as f:
            f.write(lines)
        return link
    
    def make_softlink(self, fn, sample_name, link_dir):
        """
        Make softlink for file fn in link_dir. sample_name followed by an
        underscore will be appended to the front of fn if the sample_name isn't
        in fn.
    
        Parameters
        ----------
        fn : str
            Full path to file to make link to.
    
        link_dir : str
            Path to directory where softlink should be made.
    
        Returns
        -------
        link : str
            Path to softlink.
    
        """
        if job.sample_name not in os.path.split(fn)[1]:
            name = '{}_{}'.format(sample_name, os.path.split(fn)[1])
            link = os.path.join(link_dir, name)
        else:
            name = os.path.split(fn)[1]
            link = os.path.join(link_dir, name)
        lines = 'ln -s {} {}\n'.format(fn, link)
        with open(job.filename, "a") as f:
            f.write(lines)
        return link

def _wasp_snp_directory(vcf, directory, sample_name, regions,
                        bcftools_path='bcftools'):
    """
    Convert VCF file into input files directory and files needed for WASP. Only
    bi-allelic heterozygous sites are used. Both SNPs and indels are included.

    Parameters:
    -----------
    vcf : str
        Path to VCF file.

    directory : str
        Output directory. A directory snps will be output in this directory with
        the variants for WASP.

    sample_name : str
        Use this sample name to get heterozygous SNPs from VCF file.

    regions : str
        Path to bed file to define regions of interests (e.g. exons, peaks,
        etc.). These regions should be non-overlapping.

    """
    import glob
    # First we extract all heterozygous variants for this sample.
    if not os.path.exists(directory):
        os.makedirs(directory)

    c = ('{} view -O u -m2 -M2 -R {} -s {} {} | {} view -g het | grep -v ^\\# '
         '| cut -f1,2,4,5 | '
         'awk \'{{print $2"\\t"$3"\\t"$4 >> ("{}/"$1".snps.txt")}}\''.format(
             bcftools_path, regions, sample_name, vcf, bcftools_path,
             directory))
    subprocess.check_call(c, shell=True)

    # Now we gzip the files.
    fns = glob.glob(os.path.join(directory, '*.snps.txt'))
    for fn in fns:
        subprocess.check_call('gzip {}'.format(fn), shell=True)

def wasp_allele_swap(
    bam, 
    find_intersecting_snps_path, 
    vcf, 
    sample_name,
    outdir, 
    tempdir, 
    copy_vcf=True,
    vcf_sample_name=None, 
    conda_env=None, 
    threads=6,
    samtools_path='samtools',
):
    """
    Write shell script for identifying reads in a bam file that overlap
    specified variants and switching the variant allele. This is done using
    find_intersecting_snps.py from WASP.

    Parameters
    ----------
    bam : str
        Path to input bam file.

    find_intersecting_snps_path : str
        Path to find_intersecting_snps.py script.

    vcf : str
        VCF file containing exonic SNPs.
    
    sample_name : str
        Sample name used for naming files etc.

    outdir : str
        Directory to store shell file and results.

    tempdir : str
        Directory to store temporary files.

    copy_vcf : str
        Whether to copy vcf to temp directory.

    vcf_sample_name : str
        Sample name of this sample in the VCF file (if different than
        sample_name).

    conda_env : str
        If provided, load conda environment with this name.

    threads : int
        Number of threads to request for PBS script.

    """
    job_suffix = 'wasp_allele_swap'
    job = JobScript(sample_name, job_suffix, outdir, threads, tempdir=tempdir,
                    conda_env=conda_env)
    
    if not vcf_sample_name:
        vcf_sample_name = sample_name

    # Input files.
    temp_bam = job.add_input_file(bam)
    temp_vcf = job.add_input_file(vcf, copy=copy_vcf)
    job.copy_input_files()
    
    # Files that will be created.
    temp_uniq_bam = os.path.join(
        job.tempdir, '{}_uniq.bam'.format(sample_name))
    job.temp_files_to_delete.append(temp_uniq_bam)
    
    # Files to copy to output directory.
    prefix = '{}_uniq'.format(sample_name)
    fns = [
        '{}.keep.bam'.format(prefix),
        '{}.remap.fq1.gz'.format(prefix),
        '{}.remap.fq2.gz'.format(prefix),
        '{}.to.remap.bam'.format(prefix),
        '{}.to.remap.num.gz'.format(prefix)
    ]
    job.output_files_to_copy += fns

    from __init__ import scripts
    input_script = os.path.join(scripts, 'make_wasp_input.py')

    # Run WASP to swap alleles.
    with open(job.filename, "a") as f:
        snp_directory = os.path.join(job.tempdir, 'snps')
        all_snps = os.path.join(job.outdir, 'snps.tsv')
        f.write('python {} -s \\\n\t{} \\\n\t{} \\\n\t{} \\\n\t{} & \n\n'.format(
            input_script, vcf_sample_name, temp_vcf, snp_directory, all_snps))
        f.write('{} view -b -q 255 -F 1024 \\\n\t{} \\\n\t> {}\n\n'.format(
            samtools_path, temp_bam, temp_uniq_bam))
        f.write('wait\n\n')
        f.write('python {} -s -p \\\n\t{} \\\n\t{}\n\n'.format(
            find_intersecting_snps_path, temp_uniq_bam, snp_directory))
    
    job.write_end()
    return job.filename

def wasp_alignment_compare(
    to_remap_bam, 
    remapped_bam, 
    vcf,
    fasta, 
    filter_remapped_reads_path, 
    sample_name,
    outdir, 
    tempdir, 
    picard_path='$picard',
    picard_memory=12,
    conda_env='', 
    threads=6):
    """
    Write shell script for checking original mapping position of reads against
    remapping after swapping alleles using WASP, then count allele coverage for
    each SNP.

    Parameters
    ----------
    to_remap_bam : str
        Bam file from find_intersecting_snps.py that has reads that will be
        remapped (e.g. *.to.remap.bam).

    to_remap_num : str
        Gzipped text file from find_intersecting_snps.py (e.g.
        *.to.remap.num.gz).

    remapped_bam : str
        Bam file with remapped reads.

    vcf : str
        Path to VCF file with heterozygous SNVs.

    fasta : str
        Path to fasta file used to align data.

    filter_remapped_reads_path : str
        Path to filter_remapped_reads.py script.

    sample_name : str
        Sample name used for naming files etc.

    outdir : str
        Directory to store shell file and results.

    tempdir : str
        Directory to store temporary files.

    conda_env : str
        If provided, load conda environment with this name.

    threads : int
        Number of threads to request for SGE script.

    """
    job_suffix = 'wasp_alignment_compare'
    job = JobScript(sample_name, job_suffix, outdir, threads, tempdir=tempdir,
                    conda_env=conda_env)

    # Input files.
    temp_to_remap_bam = job.add_input_file(to_remap_bam)
    temp_to_remap_num = job.add_input_file(to_remap_num)
    temp_remapped_bam = job.add_input_file(remapped_bam)
    job.copy_input_files()

    # Files that will be created.
    temp_filtered_bam = os.path.join(
        job.tempdir, '{}_filtered.bam'.format(sample_name))
    job.temp_files_to_delete.append(temp_filtered_bam)
    coord_sorted_bam = os.path.join(
        job.tempdir, '{}_filtered_coord_sorted.bam'.format(sample_name))
    job.output_files_to_copy.append(coord_sorted_bam)
    bam_index = coord_sorted_bam + '.bai'
    job.output_files_to_copy.append(bam_index)
   
    with open(job.filename, "a") as f:
        # Run WASP alignment compare.
        f.write('python {} -p \\\n\t{} \\\n\t{} \\\n\t{} \\\n\t{}\n\n'.format(
            filter_remapped_reads_path, temp_to_remap_bam, temp_remapped_bam,
            temp_filtered_bam, temp_to_remap_num))

        # Coordinate sort and index.
        lines = _picard_coord_sort(temp_filtered_bam, coord_sorted_bam,
                                   bam_index=bam_index, picard_path=picard_path,
                                   picard_memory=picard_memory,
                                   picard_tempdir=job.tempdir)
        f.write(lines)
        f.write('\nwait\n\n')
        
        # Count allele coverage.
        counts = os.path.join(job.outdir,
                              '{}_allele_counts.tsv'.format(sample_name))
        f.write('java -jar /raid3/software/GenomeAnalysisTK.jar \\\n')
        f.write('\t-R {} \\\n'.format(fasta))
        f.write('\t-T ASEReadCounter \\\n')
        f.write('\t-o {} \\\n'.format(counts))
        f.write('\t-I {} \\\n'.format(coord_sorted_bam))
        f.write('\t-sites {} \\\n'.format(vcf))
        f.write('\t-overlap COUNT_FRAGMENTS_REQUIRE_SAME_BASE \\\n')
        f.write('\t-U ALLOW_N_CIGAR_READS \n')

        f.write('\nwait\n\n')
    
    job.write_end()
    return job.filename

def wasp_remap(
    r1_fastq, 
    r2_fastq, 
    outdir, 
    sample_name, 
    star_index,
    seq_type,
    conda_env='',
    modules='',
    rgpl='ILLUMINA',
    rgpu='',
    tempdir='/scratch', 
    threads=10, 
    picard_path='$picard',
    picard_memory=15, 
    star_path='STAR',
    samtools_path='samtools',
):
    """
    Make a shell script for re-aligning reads from with variants using STAR. The
    defaults are set for use on the Frazer lab's SGE scheduler on flh1/flh2.

    Parameters
    ----------
    r1_fastq : str
        R1 reads from find_intersecting_snps.py to be remapped.

    r2_fastq : str
        R2 reads from find_intersecting_snps.py to be remapped.

    outdir : str
        Directory to store shell file and results.

    sample_name : str
        Sample name used for naming files etc.

    star_index : str
        Path to STAR index.

    star_path : str
        Path to STAR aligner.

    samtools_path : str
        Path to samtools executable.

    seq_type : str
        Type of data. Currently supports ATAC and RNA.

    conda_env : str
        If provided, load conda environment with this name. This will control
        which version of MACS2 is used.

    rgpl : str
        Read Group platform (e.g. illumina, solid). 

    rgpu : str
        Read Group platform unit (eg. run barcode). 

    tempdir : str
        Directory to store files as STAR runs.

    threads : int
        Number of threads to reserve using SGE scheduler. This number of threads
        minus 2 will be used by STAR, so this must be at least 3.

    picard_path : str
        Path to Picard tools.

    picard_memory : int
        Amount of memory (in gb) to give Picard Tools.

    Returns
    -------
    fn : str
        Path to shell script.

    """
    assert threads >= 3
    seq_types = ['ATAC', 'RNA']
    assert seq_type in seq_types, ('Only {} currently support for '
                                   'seq_type'.format(', '.join(seq_types)))
    job_suffix = 'wasp_remap'
    job = JobScript(sample_name, job_suffix, outdir, threads, tempdir=tempdir,
                    queue=queue, conda_env=conda_env, modules=modules,
                    copy_input=True)
    
    # Input files.
    temp_r1 = job.add_input_file(r1_fastq)
    temp_r2 = job.add_input_file(r2_fastq)
    job.copy_input_files()

    # Files that will be created.
    aligned_bam = os.path.join(job.tempdir, 'Aligned.out.bam')
    job.output_files_to_copy.append(aligned_bam)
    # job.temp_files_to_delete.append(aligned_bam)
    # coord_sorted_bam = os.path.join(
    #     job.tempdir, '{}_sorted.bam'.format(sample_name))
    # job.output_files_to_copy.append(coord_sorted_bam)
    job.temp_files_to_delete.append('_STARtmp')

    # Files to copy to output directory.
    job.output_files_to_copy += ['Log.out', 'Log.final.out', 'Log.progress.out',
                                 'SJ.out.tab']

    # Run WASP remapping.
    with open(job.filename, "a") as f:
        # Align with STAR and coordinate sort.
        if seq_type == 'RNA':
            from rnaseq import _star_align
            lines = _star_align(temp_r1, temp_r2, sample_name, rgpl,
                                rgpu, star_index, star_path, threads)
            f.write(lines)
            f.write('wait\n\n')
        
        elif seq_type == 'ATAC':
            from atacseq import _star_align
            lines = _star_align(temp_r1, temp_r2, sample_name, rgpl,
                                rgpu, star_index, star_path, threads)
            f.write(lines)
            f.write('wait\n\n')

        # # Coordinate sort bam file.
        # lines = _picard_coord_sort(aligned_bam, coord_sorted_bam,
        #                            picard_path, picard_memory, job.tempdir)
        # f.write(lines)
        # f.write('wait\n\n')

    job.write_end()
    return job.filename

def _mbased(
    infile, 
    bed, 
    mbased_infile, 
    locus_outfile, 
    snv_outfile, 
    sample_name,
    is_phased=False, 
    num_sim=1000000, 
    threads=1, 
    vcf=None,
    vcf_sample_name=None, 
    mappability=None,
    bigWigAverageOverBed_path='bigWigAverageOverBed',
):
    """
    Make a shell script for running MBASED to determine allelic bias from
    sequencing reads.

    Parameters
    ----------
    infile : str
        Output file from GATK's ASEReadCounter.

    bed : str
        Path to bed file for assigning heterozygous SNVs to features.
    
    mbased_infile : str
        Path to save MBASED input file to.

    locus_outfile : str
        Path to file to store locus-level results.

    snv_outfile : str
        Path to file to store SNV-level results.

    sample_name : str
        Sample name used for naming files etc.

    is_phased : bool
        Whether the input file is phased. If so, the reference alleles are
        assumed to be in phase. Note that this only matters locus by locus.

    num_sim : int
        Number of simulations for MBASED to perform.

    threads : int
        Number of threads for MBASED to use.
    
    vcf : str
        Path to gzipped, indexed VCF file with all variant calls (not just
        heterozygous calls).

    vcf_sample_name : str
        If vcf is provided, this must be provided to specify the sample name of
        this sample in the VCF file. Required if vcf is provided.

    mappability : str
        Path to bigwig file with mappability scores. A score of one should mean
        uniquely mapping.

    bigWigAverageOverBed_path : str
        Path to bigWigAverageOverBed. Required if mappability is provided.

    Returns
    -------
    lines : str
        Lines to be printed to shell script.

    """
    from __init__ import scripts
    is_phased = str(is_phased).upper()
    script = os.path.join(scripts, 'make_mbased_input.py')
    lines = 'python {} \\\n\t{} \\\n\t{} \\\n\t{}'.format(
        script, infile, mbased_infile, bed)
    if vcf:
        lines += ' \\\n\t-v {} -s {}'.format(vcf, vcf_sample_name)
    if mappability:
        lines += ' \\\n\t-m {} -p {}'.format(mappability,
                                             bigWigAverageOverBed_path)
    lines += '\n\n'
    script = os.path.join(scripts, 'mbased.R')
    lines += 'Rscript '
    lines += ' \\\n\t'.join([script, mbased_infile, locus_outfile, snv_outfile,
                             sample_name, is_phased, str(num_sim),
                             str(threads)])
    lines += '\n'
    return lines

def run_mbased(
    infile, 
    bed,
    outdir, 
    sample_name, 
    modules=None,
    conda_env=None,
    queue=None,
    is_phased=False,
    num_sim=1000000,
    threads=6, 
    vcf=None,
    vcf_sample_name=None,
    mappability=None,
    bigWigAverageOverBed_path='bigWigAverageOverBed',
):
    """
    Make a shell script for running MBASED to determine allelic bias from
    sequencing reads.

    Parameters
    ----------
    infile : str
        Tab-separated file with following columns: chrom, pos, ref_allele,
        alt_allele, locus, name, ref_count, alt_count.

    bed : str
        Path to bed file for assigning heterozygous SNVs to features.
    
    outdir : str
        Directory to store shell file and MBASED results.

    sample_name : str
        Sample name used for naming files etc.

    modules : str
        Modules (separated by commas e.g. bedtools,samtools) to load at
        beginning of script.

    is_phased : bool
        Whether the input file is phased. If so, the reference alleles are
        assumed to be in phase. Note that this only matter locus by locus.

    num_sim : int
        Number of simulations for MBASED to perform.

    threads : int
        Number of threads to reserve using SGE scheduler and for MBASED to use.

    vcf : str
        Path to gzipped, indexed VCF file with all variant calls (not just
        heterozygous calls).

    vcf_sample_name : str
        If vcf is provided, this must be provided to specify the sample name of
        this sample in the VCF file. Required if vcf is provided.

    mappability : str
        Path to bigwig file with mappability scores. A score of one should mean
        uniquely mapping.

    bigWigAverageOverBed_path : str
        Path to bigWigAverageOverBed. Required if mappability is provided.

    Returns
    -------
    fn : str
        Path to shell script.

    """
    assert threads >= 1
    job_suffix = 'mbased'
    job = JobScript(sample_name, job_suffix, outdir, threads, queue=queue,
                    modules=modules, conda_env=conda_env,
                    copy_input=True)
    
    # I'm going to define some file names used later.
    mbased_infile = os.path.join(job.outdir,
                                 '{}_mbased_input.tsv'.format(sample_name))
    locus_outfile = os.path.join(job.outdir, '{}_locus.tsv'.format(sample_name))
    snv_outfile = os.path.join(job.outdir, '{}_snv.tsv'.format(sample_name))
    
    with open(job.filename, "a") as f:
        lines = _mbased(infile, bed, mbased_infile, locus_outfile, snv_outfile,
                        sample_name, is_phased=is_phased, num_sim=num_sim,
                        threads=threads, vcf=vcf,
                        vcf_sample_name=vcf_sample_name,
                        mappability=mappability,
                        bigWigAverageOverBed_path=bigWigAverageOverBed_path)
        f.write(lines)
        f.write('wait\n\n')
    
    job.write_end()
    return job.filename

# The method below needs to be updated for SGE and JobScript. I've done a little
# already.
# def convert_sra_to_fastq(
#     sra_files, 
#     outdir, 
#     sample_name, 
#     remove_sra_files=False,
#     max_threads=32,
#     threads_per_sra=4,
#     fastq_dump_path='fastq-dump',
# ):
#     """
#     Make a shell script for converting one or more SRA files into fastq files.
#     All R1 and R2 files will be concatenated and gzipped into two single output
#     files.
# 
#     Parameters
#     ----------
#     sra_files : list
#         List of SRA files to convert.
# 
#     outdir : str
#         Directory to store shell file and gzipped fastq files.
# 
#     sample_name : str
#         Sample name used for naming files etc.
# 
#     remove_sra_files : bool
#         Whether to remove original SRA files after conversion is complete.
# 
#     max_threads : int
#         Maximum number of threads to request from SGE scheduler.
# 
#     threads_per_sra : int
#         Request this many threads per SRA input file. max_threads /
#         threads_per_sra SRA files will be converted at a time. For instance, if
#         you provide 3 SRA files, threads_per_sra=4, and max_threads=10, then the
#         first two SRA files will be converted at the same time. After they are
#         done, the last file will be converted. This is done naively using
#         background processes (&).
# 
#     fastq-dump : str
#         Path to fastq-dump from the SRA toolkit.
# 
#     Returns
#     -------
#     fn : str
#         Path to shell script.
# 
#     """
#     threads = min(threads_per_sra * len(sra_files), max_threads)
# 
#     tempdir = os.path.join(tempdir, '{}_sra_fastq'.format(sample_name))
#     outdir = os.path.join(outdir, '{}_sra_fastq'.format(sample_name))
# 
#     # I'm going to define some file names used later.
#     r1 = os.path.join(tempdir, '{}.R1.fastq.gz'.format(sample_name))
#     r2 = os.path.join(tempdir, '{}.R2.fsatq.gz'.format(sample_name))
#     temp_sra_files = [os.path.join(tempdir, os.path.split(x)[1]) for x in
#                       sra_files]
#     
#     # Files to copy to output directory.
#     files_to_copy = [r1, r2]
#     
#     # Temporary files that can be deleted at the end of the job. We may not want
#     # to delete the temp directory if the temp and output directory are the
#     # same.
#     files_to_remove = temp_sra_files
# 
#     try:
#         os.makedirs(outdir)
#     except OSError:
#         pass
# 
#     fn = os.path.join(outdir, '{}_sra_fastq.sh'.format(sample_name))
# 
#     f = open(fn, 'w')
#     f.write('#!/bin/bash\n\n')
#     if pbs:
#         out = os.path.join(outdir, '{}_sra_fastq.out'.format(sample_name))
#         err = os.path.join(outdir, '{}_sra_fastq.err'.format(sample_name))
#         job_name = '{}_sra_fastq'.format(sample_name)
#         f.write(_pbs_header(out, err, job_name, threads))
# 
#     f.write('mkdir -p {}\n'.format(tempdir))
#     f.write('cd {}\n'.format(tempdir))
#     f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
#         ' \\\n\t'.join(sra_files), tempdir))
# 
#     def chunks(l, n):
#         """Yield successive n-sized chunks from l."""
#         for i in xrange(0, len(l), n):
#             yield l[i:i+n]
#     
#     c = chunks(temp_sra_files, threads / threads_per_sra)
#     while True:
#         try:
#             n = c.next()
#             for sra in n:
#                 f.write('{} {} --split-files &\n'.format(
#                     os.path.join(fastq_dump_path), sra))
#             f.write('\nwait\n\n')
#         except StopIteration:
#             continue
#         
#     # Concatenate files, pass through awk to remove unneeded stuff, gzip.
#     f.write('cat *_1.fastq | awk \'{if (NR % 4 == 1) {print "@"$2} '
#             'if (NR % 4 == 2 || NR % 4 == 0) {print $1} '
#             'if (NR % 4 == 3) {print "+"}}\' | '
#             'gzip -c > ' + r1 + ' &\n\n')
#     f.write('cat *_2.fastq | awk \'{if (NR % 4 == 1) {print "@"$2} '
#             'if (NR % 4 == 2 || NR % 4 == 0) {print $1} '
#             'if (NR % 4 == 3) {print "+"}}\' | '
#             'gzip -c > ' + r2 + '\n\n')
#     f.write('wait\n\n')
# 
#     if remove_sra_files:
#         f.write('rm \\\n\t{}\n\n'.format(' \\\n\t'.join(sra_files)))
#             
#     f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
#         ' \\\n\t'.join(files_to_copy),
#         outdir))
#     f.write('rm \\\n\t{}\n\n'.format(' \\\n\t'.join(files_to_remove)))
# 
#     if tempdir != outdir:
#         f.write('rm -r {}\n'.format(tempdir))
#     f.close()
# 
#     return fn

def merge_bams(
    bams, 
    outdir, 
    tempdir,
    merged_name, 
    index=True,
    bigwig=False,
    copy_bams=True,
    threads=8,
    bedgraph_to_bigwig_path='bedGraphToBigWig',
    bedtools_path='bedtools',
    picard_path='$picard',
    picard_memory=2,
):
    """
    Make a shell script for combining multiple bam files using Picard.

    Parameters
    ----------
    bams : list
        List of SRA files to convert.

    outdir : str
        Directory to store shell file and merged bam file.

    merged_name : str
        Name used for output directory, files etc.

    index : bool
        Whether to index the merged bam file.

    bigwig : bool
        Whether to make bigwig file from merged bam file.

    copy_bams : bool
        Whether to copy the input bam files to the temp directory. Not
        necessary if temp directory is on the same file system as bam files.

    threads : int
        Number of threads to request from SGE scheduler.

    Returns
    -------
    fn : str
        Path to shell script.

    """
    if bigwig:
        index = True
        assert bedgraph_to_bigwig_path
        assert bedtools_path

    job_suffix = 'merged_bam'
    job = JobScript(merged_name, job_suffix, outdir, threads, tempdir=tempdir,
                    copy_input=False)
    
    # I'm going to define some file names used later.
    merged_bam = job.add_temp_file('{}_merged.bam'.format(merged_name),
                                   copy=True)
    # merged_bam = os.path.join(tempdir,
    #                             '{}_merged.bam'.format(merged_name))
    # job.output_files_to_copy.append(merged_bam)
    if index:
        merged_bam_index = job.add_temp_file(
            '{}_merged.bam.bai'.format(merged_name), copy=True)
        # merged_bam_index = os.path.join(tempdir,
        #                             '{}_merged.bam.bai'.format(merged_name))
        # job.output_files_to_copy.append(merged_bam_index)
    
    if bigwig:
        merged_bigwig = job.add_temp_file(
            '{}_merged.bw'.format(merged_name), copy=True)
    
    temp_bams = []
    for bam in bams:
        temp_bams.append(job.add_input_file(bam))
    # self.input_files_to_copy += bams

    if copy_bams:
        job.copy_input_files()

    with open(job.filename, "a") as f:
        lines = _picard_merge(temp_bams, merged_bam, picard_memory,
                              picard_path=picard_path, picard_tempdir=tempdir)
        f.write(lines)
        
        if index:
            lines = _picard_index(merged_bam, merged_bam_index,
                                  picard_path=picard_path,
                                  picard_memory=picard_memory,
                                  picard_tempdir=tempdir, bg=bigwig)
            f.write(lines)

        if bigwig:
            lines = _bigwig_files(
                merged_bam, merged_bigwig, merged_name,
                bedgraph_to_bigwig_path=bedgraph_to_bigwig_path,
                bedtools_path=bedtools_path)
            f.write(lines)

    job.write_end()
    return job.filename