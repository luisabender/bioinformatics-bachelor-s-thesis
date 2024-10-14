import os
import pandas as pd

kallisto_dir = "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/kallisto_output"
sample_dirs = []

# get list of sample directories
sample_dirs = [d for d in os.listdir(kallisto_dir) if os.path.isdir(os.path.join(kallisto_dir, d))]
count_matrix = pd.DataFrame()

for sample in sample_dirs:
    tsv_file = os.path.join(kallisto_dir, sample, "abundance.tsv")
    abundance_df = pd.read_csv(tsv_file, sep="\t")
    
    # Extract target ID (1st column) and tpm
    sample_counts = abundance_df[['target_id', 'est_counts']].copy()
    
    # tpm = name of sample
    sample_counts.rename(columns={'est_counts': sample}, inplace=True)
    
    if count_matrix.empty:
        count_matrix = sample_counts
    else:
        count_matrix = pd.merge(count_matrix, sample_counts, on='target_id', how='outer')

# target_id as rownames
count_matrix.set_index('target_id', inplace=True)

# save as csv file
count_matrix.to_csv("count_matrix_raw.csv")
