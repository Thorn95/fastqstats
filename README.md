# fastqstats

fastqstats is a command-line tool for basic quality control and summary statistics of FASTQ files.

---

## Installation

Clone the repository and install it in editable mode:

```bash

git clone git@github.com:Thorn95/fastqstats.git

cd fastqstats

pip install -e .

```

---

## Usage

Run the program from the comman line:

```bash

fastqstats /path/to/file.fastq

```

The program writes the results file to disk and prints the output path to the terminal.

---

## Options

To see all available command-line options and default thesholds, run:

```bash

fastqstats --help

```

---

## Input

* FASTQ file (Phred+33 encoded)

---

## Output

A text file containing summary statistics and quality classifications (PASS, Warning, FAIL).

---

## Documentation

A Jupyter notebook demonstrating `fastqstats` usage is available in `documentation/`:

- [`Project_fastq.ipynb`](Fastq_qc_sum_program/documentation/Project_fastq.ipynb) - interactive environment with python
- [`Project_fastq.html`](Fastq_qc_sum_program/documentation/Project_fastq.html) - view in browser without python

---

## Author

Martin Törn

---
