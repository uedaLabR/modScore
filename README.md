# modScore

**modScore** is a Python-based standalone tool designed to filter RNA modification calls detected by Oxford Nanopore Technologies (ONT) Dorado and modkit.  
It integrates deep learning and known modification sites to reduce false positives.

Currently supported RNA modifications:

- **m⁶A**
- **m⁵C**
- **Ψ (Pseudouridine)**
- **Inosine**

---

## 🔧 Installation

### 1. Standalone (local Python environment)

Make sure you are using **Python 3.10** or later, and install the required packages with:

```bash
pip install --no-cache-dir \
    numpy==1.24.4 \
    tensorflow==2.15 \
    numba==0.60.0 \
    pandas==2.2.3 \
    pysam==0.22.0 \
    click==8.0.4 \
    scikit-learn==1.5.2
