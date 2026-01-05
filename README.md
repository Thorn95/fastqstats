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
Or with optional parameter:

```bash

fastqstats /path/to/file.fastq --base-quality-warn 20

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

## File manifest

- `Fastq_qc_sum_program/fastqstats/cli.py` - command-line interface
- `Fastq_qc_sum_program/fastqstats/fastq_qcsum.py` - main functions for FASTQ analysis
- `Fastq_qc_sum_program/fastqstats/Example_files` - small example FASTQ files
- `README.md` - instructions and overview
- `Fastq_qc_sum_program/documentation` - documentation notebook
- `.gitignore` - files/folders ignored by Git

---

## Documentation

A Jupyter notebook demonstrating `fastqstats` usage is available in `documentation/`:

- [`Project_fastq.ipynb`](Fastq_qc_sum_program/documentation/Project_fastq.ipynb) - interactive environment with python
- [`Project_fastq.html`](Fastq_qc_sum_program/documentation/Project_fastq.html) - view in browser without python

---

## Author

Martin Törn

---
