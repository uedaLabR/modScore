
# modScore

**modScore** is a Python-based tool for filtering RNA modification calls from Oxford Nanopore Technologies (ONT) Dorado and modkit.  
It combines **deep learning** and **known modification sites** to reduce false positives.

### 🧬 Supported RNA modifications

- **m⁶A**
- **m⁵C**
- **Ψ (Pseudouridine)**
- **Inosine**

### 🧬 Supported genome

- hg38
- mm10
---

## 🚀 Quick Start

### Option 1: Run with Python (Standalone)

Clone the repository:

```
git clone https://github.com/uedaLabR/modScore.git
cd modScore
```

Install dependencies (Python 3.10 or later recommended):

```
pip install --no-cache-dir \
    numpy==1.24.4 \
    tensorflow==2.15 \
    numba==0.60.0 \
    pandas==2.2.3 \
    pysam==0.22.0 \
    click==8.0.4 \
    scikit-learn==1.5.2
```

Run the program:

```
python main.py <command> [OPTIONS]
```

---

### Option 2: Run with Docker

Pull the prebuilt image:

```
docker pull karkinos/modscore_v01:latest
```


---

## 🛠 Available Commands

### 1. `filter` – Filter modification BED and generate statistics

```
python MSCmd.py filter \
  --bed input.bed \
  --bed_out filtered_output.bed \
  --source_path path/to/known_sites \
  --genome hg38
```

**Arguments:**

- `--bed`: Input BED file (from Dorado/modkit)
- `--bed_out`: Output filtered BED file
- `--source_path`: Directory containing known modification site BEDs
- `--genome`: Genome version (default: `hg38`)

➡️ Outputs both a filtered BED and a `*_stats.txt` summary.

---

### 2. `reflectToBam` – Update ML tags in BAM using filtered BED

```
python main.py reflectToBam \
  --bamin input.bam \
  --bamout output.bam \
  --filter_bed filtered_output.bed
```

**Arguments:**

- `--bamin`: Original BAM file
- `--bamout`: Output BAM with updated ML tags
- `--filter_bed`: BED file from the `filter` command

---

### 3. `trainSequenceClassification` – Train deep learning model

```
python main.py trainSequenceClassification \
  --source_path training_data_dir \
  --genome hg38 \
  --fp_ivtpath ivt_data_dir \
  --outhistory train_log.txt \
  --weightpath model_weights.h5
```

**Arguments:**

- `--source_path`: Directory of training data
- `--fp_ivtpath`: Directory of positive (IVT) data
- `--outhistory`: Training history log output
- `--weightpath`: Output path for model weights

---


## 🧪 Tested Environment

- Python 3.10
- Ubuntu 22.04
- TensorFlow 2.15 (GPU-enabled)
- ONT Dorado 0.9.1 + modkit 0.4.5

---

## 📦 Docker Image

- Docker Hub: [karkinos/modscore_v01](https://hub.docker.com/r/karkinos/modscore_v01)



