# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 10:54:50 2025

@author: marti
"""

#Command line script

#Imports
import argparse
from fastqstats.fastq_qcsum import stats_fastq

#Main function
def main():
    #Argument parser
    parser = argparse.ArgumentParser(description = "FASTQ quality control and summary tool",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    #Required argument
    parser.add_argument("fastq_path", help = "Path to FASTQ file")
    
    #Optional arguments for quality threshold
    parser.add_argument("--base-quality-warn", type=float, default=28,
                        help="Base quality warning threshold")
    parser.add_argument("--base-quality-fail", type=float, default=20,
                        help="Base quality fail threshold")
    parser.add_argument("--base-content-warn", type=float, default=0.10,
                        help="Base frequency warning threshold")
    parser.add_argument("--base-content-fail", type=float, default=0.20,
                        help="Base frequency fail threshold")
    parser.add_argument("--read-quality-warn", type=float, default=28,
                        help="Mean read quality warning threshold")
    parser.add_argument("--read-quality-fail", type=float, default=20,
                        help="Mean read quality fail threshold")
    parser.add_argument("--adapter-warn", type=float, default=0.05,
                        help="Adapter contamination warning threshold")
    parser.add_argument("--adapter-fail", type=float, default=0.10,
                        help="Adapter contamination fail threshold")
    
    #Parse agument from command line
    args = parser.parse_args()
    
    #Call stats_fastq (main function)
    resfile = stats_fastq(
        args.fastq_path,
        args.base_quality_warn,
        args.base_quality_fail,
        args.base_content_warn,
        args.base_content_fail,
        args.read_quality_warn,
        args.read_quality_fail,
        args.adapter_warn,
        args.adapter_fail)
    
    #Print confirmation in console
    print(f"Results written to {resfile}")
    
#Ensure correct execution
if __name__ == "__main__":
    main()