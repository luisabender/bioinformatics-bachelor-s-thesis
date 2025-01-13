import os
import pandas as pd

kallisto_dir_single = "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/kallisto_paired"
kallisto_dir_paired = "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/kallisto_single"


sample_dirs_single = [d for d in os.listdir(kallisto_dir_single) if os.path.isdir(os.path.join(kallisto_dir_single, d))]
count_matrix_single = pd.DataFrame()
sample_dirs_paired = [d for d in os.listdir(kallisto_dir_paired) if os.path.isdir(os.path.join(kallisto_dir_paired, d))]
count_matrix_paired = pd.DataFrame()

# create count matrix for single end reads
for sample in sample_dirs_single:
    tsv_file_single = os.path.join(kallisto_dir_single, sample, "abundance.tsv")
    abundance_df_single = pd.read_csv(tsv_file_single, sep="\t")
    
    # extract target ID (1st column) and tpm
    sample_counts_single = abundance_df_single[['target_id', 'est_counts']].copy()
    
    # tpm = name of sample
    sample_counts_single.rename(columns={'est_counts': sample}, inplace=True)
    
    if count_matrix_single.empty:
        count_matrix_single = sample_counts_single
    else:
        count_matrix_single = pd.merge(count_matrix_single, sample_counts_single, on='target_id', how='outer')

# target_id as rownames
count_matrix_single.set_index('target_id', inplace=True)

# create count matrix for paired end reads
for sample in sample_dirs_paired:
    tsv_file_paired = os.path.join(kallisto_dir_paired, sample, "abundance.tsv")
    abundance_df_paired = pd.read_csv(tsv_file_paired, sep="\t")
    
    # extract target ID (1st column) and tpm
    sample_counts_paired = abundance_df_paired[['target_id', 'est_counts']].copy()
    
    # tpm = name of sample
    sample_counts_paired.rename(columns={'est_counts': sample}, inplace=True)
    
    if count_matrix_paired.empty:
        count_matrix_paired = sample_counts_paired
    else:
        count_matrix_paired = pd.merge(count_matrix_paired, sample_counts_paired, on='target_id', how='outer')

# target_id as rownames
count_matrix_paired.set_index('target_id', inplace=True)

count_matrix = pd.merge(count_matrix_single, count_matrix_paired, left_index=True, right_index=True, how="outer")

# save as csv file
count_matrix.to_csv("/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/bachelor-thesis-luisa/count_tables/count_matrix_raw.csv")
