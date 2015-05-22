import os
import shutil

import pytest

import pipelines as ps

class TestPicardInsertSizeMetrics:
    def test_run(self):
        """Test to make sure the function at least runs"""
        in_bam = 'test.bam'
        out_metrics = 'metrics.txt'
        out_hist = 'hist.pdf'
        picard_path = 'path/to/picard'
        picard_memory = '58G'
        tempdir = 'path/to/temp/dir'
        bg = False
        ps.general._picard_insert_size_metrics(in_bam, out_metrics, out_hist,
                                               picard_path, picard_memory,
                                               tempdir, bg=bg)

    def test_run_bg(self):
        """Test to make sure the function at least runs with bg"""
        in_bam = 'test.bam'
        out_metrics = 'metrics.txt'
        out_hist = 'hist.pdf'
        picard_path = 'path/to/picard'
        picard_memory = '58G'
        tempdir = 'path/to/temp/dir'
        bg = True
        ps.general._picard_insert_size_metrics(in_bam, out_metrics, out_hist,
                                               picard_path, picard_memory,
                                               tempdir, bg=bg)

class TestPicardQuerySort:
    def test_run(self):
        """Test to make sure the function at least runs"""
        in_bam = 'test.bam'
        out_bam = 'test_qsorted.bam'
        picard_memory = '58G'
        picard_path = 'path/to/picard'
        tempdir = 'path/to/temp/dir'
        bg = False
        ps.general._picard_query_sort(in_bam, out_bam, picard_path,
                                      picard_memory, tempdir, bg=bg)

    def test_run_bg(self):
        """Test to make sure the function at least runs with bg"""
        in_bam = 'test.bam'
        out_bam = 'test_qsorted.bam'
        picard_memory = '58G'
        picard_path = 'path/to/picard'
        tempdir = 'path/to/temp/dir'
        bg = True
        ps.general._picard_query_sort(in_bam, out_bam, picard_path,
                                      picard_memory, tempdir, bg=bg)

class TestPicardCoordSort:
    def test_run(self):
        """Test to make sure the function at least runs"""
        in_bam = 'test.bam'
        out_bam = 'test_sorted.bam'
        picard_path = 'path/to/picard'
        picard_memory = '58G'
        tempdir = 'path/to/temp/dir'
        fn = ps.general._picard_coord_sort(
            in_bam, 
            out_bam, 
            picard_path, 
            picard_memory,
            tempdir, 
            bam_index=None
        )

    def test_run_index(self):
        """Test to make sure the function at least runs when making an index"""
        in_bam = 'test.bam'
        out_bam = 'test_sorted.bam'
        picard_path = 'path/to/picard'
        picard_memory = '58G'
        tempdir = 'path/to/temp/dir'
        bam_index = 'test_sorted.bam.bai'
        fn = ps.general._picard_coord_sort(
            in_bam, 
            out_bam, 
            picard_path, 
            picard_memory,
            tempdir, 
            bam_index=bam_index,
        )

class TestCutadaptTrim:
    def test_run(self):
        """Test to make sure the function at least runs"""
        fastq = 'test.fastq.gz'
        length = 50
        out = 'test_cut.fastq.gz'
        bg = False
        fn = ps.general._cutadapt_trim(fastq, length, out, bg=bg)

    def test_run_bg(self):
        """Test to make sure the function at least runs with bg"""
        fastq = 'test.fastq.gz'
        length = 50
        out = 'test_cut.fastq.gz'
        bg = True
        fn = ps.general._cutadapt_trim(fastq, length, out, bg=bg)

    def test_run_negative(self):
        """Test to make sure the function at least runs with negative amount to
        cut"""
        fastq = 'test.fastq.gz'
        length = -50
        out = 'test_cut.fastq.gz'
        bg = False
        fn = ps.general._cutadapt_trim(fastq, length, out, bg=bg)

class TestPicardIndex:
    def test_run(self):
        """Test to make sure the function at least runs"""
        in_bam = 'test.bam'
        index = 'test.bam.bai'
        picard_memory = '58G'
        picard_path = 'path/to/picard'
        tempdir = 'path/to/temp/dir'
        lines = ps.general._picard_index(in_bam, index, picard_memory,
                                         picard_path, tempdir)

class TestWaspAlleleSwap:
    def test_run(self):
        """Test to make sure the function at least runs"""
        bam = 'test.bam'
        find_intersecting_snps_path = 'path/to/find_intersecting_snps.py'
        snp_dir = 'test/dir'
        sample_name = 'test'
        outdir = 'path/to/out'
        tempdir = 'path/to/temp/dir'
        fn = ps.general.wasp_allele_swap(
            bam, 
            find_intersecting_snps_path,
            snp_dir, 
            sample_name, 
            outdir, 
            tempdir,
            conda_env=None, 
            shell=False, 
            threads=6
        )
        shutil.rmtree('path')

class TestWaspAlignmentCompare:
    def test_run(self):
        """Test to make sure the function at least runs"""
        to_remap_bam = 'to_remap.bam'
        to_remap_num = 'to_remap.num'
        remapped_bam = 'remapped.bam'
        filter_remapped_reads_path = 'path/to/filter_remapped_reads.py'
        sample_name = 'test'
        outdir = 'path/to/out'
        tempdir = 'path/to/temp/dir'
        picard_path = 'path/to/picard'

        fn = ps.general.wasp_alignment_compare(
            to_remap_bam, 
            to_remap_num, 
            remapped_bam,
            filter_remapped_reads_path, 
            sample_name, 
            outdir, 
            tempdir,
            picard_path, 
            picard_memory=58, 
            conda_env='', 
            shell=False, 
            threads=6
        )
        shutil.rmtree('path')
