from bam_manup import BamRecalib
import cProfile
import pstats
import pysam
import matplotlib.pyplot as plt

bamfile ="/mnt/ssdnas07/nanozero/rna/rna_modscore_v1/Adipocyte_1/Adipocyte_1/Adipocyte_1_sorted.bam"

bam = pysam.AlignmentFile(bamfile, "rb")

lengths = []
cnt = 0
for read in bam.fetch(until_eof=True):
    if not read.is_unmapped:   # }bv
        lengths.append(read.query_length)
        cnt+=1
    if cnt > 100000:
        break
bam.close()

print("Total reads:", len(lengths))

# qXgO`
plt.figure(figsize=(10,6))
plt.hist(lengths, bins=100, edgecolor="black")
plt.xlabel("Read length (bp)")
plt.ylabel("Count")
plt.yscale('log')
plt.title("Read length distribution")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("/mnt/share/ueda/read_length_histogram.png", dpi=300)