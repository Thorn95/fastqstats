# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 10:24:36 2025

@author: marti
"""
#Import
import numpy as np
import pandas as pd
import os
import statistics

#Read FASTQ file
def read_fastq(fastq_file):
    """
    Reads in a FASTQ file by line.

    Args:
        fastq_file (str): Path to a FASTQ file.

    Returns:
        lines: A list of strings where each string is a line from the file.
    """
    with open(fastq_file, "r") as file:
        lines = file.readlines()
    return lines

#Extract sequences (bases)
def array_seq(fastqdata):
    """
    Extracts all sequence lines from FASTQ data and convert them to ASCII format bases.

    FASTQ files store reads in groups of four lines where the 2nd line contains the sequence
    (index 1, 5, 9, ...). This function extracts the sequences and converts each character base
    to its ASCII interger.

    Args:
        fastqdata (list[str]): Lines from a FASTQ file where every fourth line contains a 
        sequence.
        
    Returns:
        np.ndarray: A 2D uint8 array of shape (num_reads, read_length) where every element is 
        a base in ASCII format.
    """
    seq_lines = [
        [ord(base) for base in line.strip()]#Turn bases into numbers 'A'=65, 'C'=67, 'G'=71, 'T'=84 
        for i, line in enumerate(fastqdata)
        if i % 4 == 1 #Takes every fourth lines starting from index 1 (first line with sequence)
    ]
    return np.array(seq_lines, dtype=np.uint8)       

#Extract scores
def array_scores(fastqdata):
    """
    Extracts all scores lines from FASTQ data and converts them to phred scores.

    FASTQ files store reads in groups of four lines where the 4th line contains the scoring
    in ASCII format (index 3, 7, 11, ...). This function extracts the score line and converts 
    each ASCII encoded score value to its corresponding phred score by subtractubg 33.
    
    Args:
        fastqdata (list[str]): Lines from a FASTQ file where every fourth line contains the
        ASCII-encoded base quality scores.
        
    Returns:
        np.ndarray: A 2D uint8 array of shape (num_reads, read_length) where each element is 
        the Phred quality score for a base.
    """
    scores_lines = [
        [ord(char) - 33 for char in line.strip()]#converts ASCII to phred scores
        for i, line in enumerate(fastqdata)
        if i % 4 == 3 #Takes every fourth lines starting from index 3 (first line with scores)
    ]            
    return np.array(scores_lines, dtype=np.uint8)#shape (nreads, readlen), data type 8 bit

#Base frequencies per postition
def base_freq(seq_array):
    """
    Calculates per-base frequencies  and the maximum frequency difference at each position.

    For each position in the sequence reads, this function calculates the frequency of each
    nucleotide (A, C, G, T). The difference between the maximum and minumum nucleotide 
    frequencies at each position is also calculated for quality assessment.

    Args: 
        seq_array(np.ndarray): A 2D numpy array of shape (num_reads, read_length) containing
        the nucleotide for each position in ASCII format.

    Returns:
        diffs(np.ndarray): A 1D array of length read_length, where each value is the 
        maximum difference between two nucleotides at each position in the read.
        Higher values might indicate a bias at that position.
    """
    
    bases = ['A', 'C', 'G', 'T']
    base_codes = np.array([ord(b) for b in bases], dtype=np.uint8)

    #Count bases per position
    #Create a 2D array (4 x num_positions) of counts
    base_counts = np.vstack([
        (seq_array == code).sum(axis=0)
        for code in base_codes
    ])

    #Calculate frequencies
    total_counts = base_counts.sum(axis=0)
    base_freq = base_counts / total_counts  # shape (4, num_positions)

    #Calculate per-position frequency differences
    diffs = base_freq.max(axis=0) - base_freq.min(axis=0)

    return diffs

#Adapter identifier
def adapter_find(seq_array):
    """
    Identifies and counts the pressence of common adapter sequences in sequencing reads.
    
    This function scans each read in the FASTQ file for known adapter sequences.
    A match is counted if an adapter sequence is present at the beginning (5') or en (3')
    of a read. The proportion of adapter sequences present in FASTQ file might indicate 
    adapter contaminations which might need to be removed before downstream analysis.
    
    Adapter sequences source: 
        •	Illumina Universal Adapter—AGATCGGAAGAG
        •	Illumina Small RNA 3' Adapter—TGGAATTCTCGG
        •	Illumina Small RNA 5' Adapter—GATCGTCGGACT
        •	Nextera Transposase Sequence—CTGTCTCTTATA
    
    Args: 
        seq_array(np.ndarray): A 2D numpy array of shape (num_reads, read_length) containing
        the nucleotide for each position in ASCII format.

    Returns:
        adapter_counts(dict[str, int]: A dictionary mapping each given adapter sequence
        to the number of times it appears at the start or end of a read in the FASTQ file.
    """

    adapters = ["AGATCGGAAGAG","TGGAATTCTCGG","GATCGTCGGACT","CTGTCTCTTATA"]
    adapter_counts = {adapter:0 for adapter in adapters}

    #Reads one row at a time
    for row in seq_array:
        read = "".join(map(chr, row))

        #Count number of common adapter seq that appears at start or end
        for adapter in adapters:
            if read.startswith(adapter) or read.endswith(adapter):
                adapter_counts[adapter] += 1

    return adapter_counts 

#Main function
def stats_fastq(
    fastq_path,
    base_quality_warn = 28,
    base_quality_fail = 20,
    base_content_warn = 0.10,
    base_content_fail = 0.20,
    read_quality_warn = 28,
    read_quality_fail = 20,
    adapter_warn = 0.05,
    adapter_fail = 0.10):
    """
    Computes summary statistics and quality assessment for a FASTQ file.
    
    This function reads a FASTQ file, extracts sequences and quality scores, 
    and caulculates overall and per-base metrics for quality control. Metrics 
    include per-base sequence quality, per-base nucleotide composition, 
    per-read average quality, adapter contamination, and GC content.
    Each metric is compared against user-defined warning and failure thresholds
    to producde a qualitative assessment (PASS, Warning, FAIL) per metric. The
    results are written to a summary text file named "<filename>_summary.txt"
    
    Args:
        fastq_path (str): Path to the FASTQ file to analyze.
        base_quality_warn (int, optional): Warning threshold for per-base quality scores. Defaults to 28.
        base_quality_fail (int, optional): Failure threshold for per-base quality scores. Defaults to 20.
        base_content_warn (float, optional): Warning threshold for per-base nucleotide bias. Defaults to 0.10.
        base_content_fail (float, optional): Failure threshold for per-base nucleotide bias. Defaults to 0.20.
        read_quality_warn (int, optional): Warning threshold for per-read average quality. Defaults to 28.
        read_quality_fail (int, optional): Failure threshold for per-read average quality. Defaults to 20.
        adapter_warn (float, optional): Warning threshold for adapter content. Defaults to 0.05.
        adapter_fail (float, optional): Failure threshold for adapter content. Defaults to 0.10.
    
    Returns:
        str: The filename of the summary text file generated by the function, containing
        both quantitative statistics and qualitative quality assessments.
    """
    #Extract filename
    filename = os.path.basename(fastq_path)

    #Read in and process
    fastq_data = read_fastq(fastq_path)
    seq_array = array_seq(fastq_data)
    score_array = array_scores(fastq_data) 

    #Statistics
    #For summary
    total_seq = len(seq_array)
    read_length = round((sum(len(i) for i in seq_array)) / len(seq_array))
    total_bases = (total_seq * read_length) / 1000000
    gc_content = round((np.sum((seq_array == ord('G')) | (seq_array == ord('C')))) / (seq_array.size) * 100, 2)
    #For Quality assesment
    base_quality = min(np.median(score_array, axis=0))
    base_content = base_freq(seq_array)
    read_quality = round((sum(score_array.mean(axis=1)))/len(score_array.mean(axis=1)))
    adapter_content = np.array(list(adapter_find(seq_array).values())) / total_seq
    
    
    #Results output
    results_file = filename + "_summary.txt"
    with open(results_file, "w") as r:
        #Summary fastq
        r.write(f" Summary \n")
        r.write(f"Filename: {filename}\n")
        r.write(f"Total Sequences: {total_seq}\n")
        r.write(f"Total Bases: {total_bases}Mbp\n")
        r.write(f"Average Length: {read_length}bp\n")
        r.write(f"GC content: {gc_content}%\n")


        # Per Base Sequence Quality
        if base_quality < base_quality_fail:
            base_quality_status = "FAIL"
        elif base_quality < base_quality_warn:
            base_quality_status = "Warning"
        else:
            base_quality_status = "PASS"

        r.write(f"Per Base Sequence Quality: {base_quality_status} \n")

        # Per Base Sequence Content
        max_base_content = np.max(base_content)
        if max_base_content > base_content_fail:
            base_content_status = "FAIL"
        elif max_base_content > base_content_warn:
            base_content_status = "Warning"
        else:
            base_content_status = "PASS"

        r.write(f"Per Base Sequence content: {base_content_status} \n")
        
        # Per Read Quality
        if read_quality < read_quality_fail:
            read_quality_status = "FAIL"
        elif read_quality < read_quality_warn:
            read_quality_status = "Warning"
        else:
            read_quality_status = "PASS"

        r.write(f"Per Read Quality: {read_quality_status} \n")

        # Adapter Content
        max_adapter_content = np.max(adapter_content)
        if max_adapter_content > adapter_fail:
            adapter_status = "FAIL"
        elif max_adapter_content > adapter_warn:
            adapter_status = "Warning"
        else:
            adapter_status = "PASS"

        r.write(f"Adapter content: {adapter_status} \n")

    return results_file