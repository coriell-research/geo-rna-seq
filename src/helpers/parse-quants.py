from collections import defaultdict
from pathlib import Path


data_dir = Path("/home/gcalendo/data/projects/geo-rna-seq/data")
quant_files = list(data_dir.glob("**/quant.sf"))

data = defaultdict(lambda: defaultdict(float))
n_files = len(quant_files)
file_count = 1

for quant_file in quant_files:
    print(f"Processing {file_count}/{n_files}: {quant_file}")
    file_count += 1
    with open(quant_file, "r") as infile:
        next(infile)  # skip header
        for line in infile:
            l = line.strip().split("\t")
            md5 = l[0]
            count = float(l[4])
            data[md5]["N"] += 1
            data[md5]["sum"] += count
            data[md5]["sq_sum"] += count**2

print("Writing out results to file...")
outfile = Path("/home/gcalendo/data/projects/geo-rna-seq/results/normalization/data-files/transcript-statistics.tsv")
with open(outfile, "w") as f:
    f.write("md5\tN\tsum\tsum_sq\n")
    for md5 in data:
        f.write(f"{md5}\t{data[md5]['N']}\t{data[md5]['sum']}\t{data[md5]['sq_sum']}\n")
print("Finished.")
