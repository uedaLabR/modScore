# modScore

**modScore** is a Python-based standalone tool designed to filter RNA modification calls detected by Oxford Nanopore Technologies (ONT) Dorado and modkit.  
It integrates deep learning and known modification sites to reduce false positives.

Currently supported RNA modifications:

- **m⁶A**
- **m⁵C**
- **Ψ (Pseudouridine)**
- **Inosine**

Currently supported genome :
- **hg38**
- **mm10**
---

## 🔧 Installation

### 1. Standalone (local Python environment)

Make sure you are using **Python 3.10** or later, and install the required packages with:

```bash

git clone https://github.com/uedaLabR/modScore.git

pip install --no-cache-dir \
    numpy==1.24.4 \
    tensorflow==2.15 \
    numba==0.60.0 \
    pandas==2.2.3 \
    pysam==0.22.0 \
    click==8.0.4 \
    scikit-learn==1.5.2

## 🚀 Command-line Usage
All commands are available through the MSCmd.py entry point using the Click CLI framework.

python MSCmd.py <command> [OPTIONS]

## 🔍 1. Filter RNA modification BED file
python main.py filter \
    -bed <input_bed> \
    -bed_out <filtered_output_bed> \
    -source_path <known_site_database_dir> \
    -genome <genome_version>

input_bed: BED file from Dorado/modkit
filtered_output_bed: Output BED file after filtering
source_path: Directory containing known modification site BED files
genome_version: Reference genome version (e.g., "hg38", default: "hg38")

This command also generates a corresponding *_stats.txt file summarizing the filter results.
