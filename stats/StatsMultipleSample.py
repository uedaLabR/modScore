import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering


def stats_result(bed_files, output_dir):
    files = bed_files.split(",")

    posDict = {}  # alt -> set of "chr:pos"
    for bed_file in files:
        print("reading",bed_file)
        with open(bed_file) as f:
            for line in f:
                cols = line.strip().split("\t")
                chrom, pos = cols[0], cols[1]
                alt = cols[3]
                passfail = cols[-1]
                if passfail == "Pass":
                    key = f"{chrom}:{pos}"
                    posDict.setdefault(alt, set()).add(key)

    file_dicts = {}
    for bed_file in files:
        d = {}
        with open(bed_file) as f:
            for line in f:
                cols = line.strip().split("\t")
                chrom, pos = cols[0], cols[1]
                alt = cols[3]
                ratio = cols[10]
                key = f"{chrom}:{pos}"
                if key in posDict[alt]:

                    d.setdefault(alt, {})[key] = ratio

        file_dicts[bed_file] = d






    header = ["chrom", "pos", "alt"] + [os.path.basename(f) for f in files]

    for alt, keys in posDict.items():
        # ALTt@CSgi: m6A  m6A, A:m  A_mj
        safe_alt = re.sub(r"[^\w\-]", "_", alt)
        output_file = os.path.join(output_dir, f"stats_out.{safe_alt}.bed")
        with open(output_file, "w") as out:
            out.write("\t".join(header) + "\n")

            for key in sorted(keys):
                chrom, pos = key.split(":")
                row = [chrom, pos, alt]
                for f in files:
                    row.append(file_dicts[f].get(alt, {}).get(key, "0"))
                out.write("\t".join(row) + "\n")


def cluster_per_alt(input_dir, output_dir):
    print("start clustring")
    for filepath in glob.glob(os.path.join(input_dir, "stats_out.*.bed")):

        print("filepath",filepath)
        alt = os.path.basename(filepath).replace("stats_out.", "").replace(".bed", "")
        df = pd.read_csv(filepath, sep="\t")

        # Extract numeric data columns (starting from the 4th column)
        data_cols = df.columns[3:]
        data = df[data_cols].astype(float)

        # Remove rows that contain 2 or more zeros
        zero_counts = (data == 0).sum(axis=1)
        data_filtered = data[zero_counts < 2]
        df_filtered = df.loc[data_filtered.index]

        # Remove rows without variation (all values are identical)
        variation_mask = (data_filtered.max(axis=1) - data_filtered.min(axis=1)) > 0.3
        data_var = data_filtered[variation_mask]
        df_var = df_filtered[variation_mask]

        if data_var.shape[0] < 2 or data_var.shape[1] < 2:
            print(f"Skipping {alt}: not enough variable positions")
            continue

        # Standardize the data
        scaled = StandardScaler().fit_transform(data_var)

        print("clustering",len(scaled))
        # Perform hierarchical clustering with 3 clusters
        cluster = AgglomerativeClustering(n_clusters=3)
        labels = cluster.fit_predict(scaled)

        # Add cluster labels to the dataframe
        df_var = df_var.copy()  # avoid SettingWithCopyWarning
        df_var["cluster"] = labels

        # Save the clustered result as TSV
        output_path = os.path.join(output_dir, f"{alt}_clustered.tsv")
        df_var.to_csv(output_path, sep="\t", index=False)

        # Create and save heatmap
        plt.figure(figsize=(10, max(3, len(df_var) * 0.2)))
        sns.heatmap(
            scaled,
            yticklabels=df_var["chrom"] + ":" + df_var["pos"].astype(str) + " (Cluster " + df_var["cluster"].astype(str) + ")",
            cmap="vlag",
            center=0
        )
        plt.title(f"Clustering heatmap for {alt}")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{alt}_heatmap.png"))
        plt.close()

        print(f"Processed {alt}: {len(df_var)} variable positions")

stats_result(
    "/mnt/ssdnas/nanozero/rna/rna_modscore_v1/Adipocyte_1/Adipocyte_1/Adipocyte_1_pileup_filter.bed,"
    "/mnt/ssdnas/nanozero/rna/rna_modscore_v1/Adipocyte_2/Adipocyte_2/Adipocyte_2_pileup_filter.bed,"
    "/mnt/ssdnas/nanozero/rna/rna_modscore_v1/Adipocyte_3/Adipocyte_3/Adipocyte_3_pileup_filter.bed",
    "/share/ueda/stats_output/"
)
cluster_per_alt("/share/ueda/stats_output/", "/share/ueda/stats_output/clustering_results/")