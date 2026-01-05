\# fastqstats



fastqstats is a command-line tool for basic quality control and summary statistics of FASTQ files.



---



\## Installation



Clone the repository and install it in editable mode:

```bash

git clone git@github.com:Thorn95/fastqstats.git

cd fastqstats

pip install -e .



---



\## Usage



Run the program from the comman line:



fastqstats /path/to/file.fastq



The program writes the results file to disk and prints the output path to the terminal.



---



\## Options



To see all available command-line options and default thesholds, run:



fastqstats --help



---



\## Input



* FASTQ file (Phred+33 encoded)



---



\## Output



A text file containing summary statistics and quality classifications (PASS, Warning, FAIL).



---



\## Author



Martin Törn



---

