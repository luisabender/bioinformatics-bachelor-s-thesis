base_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa"

# generate Genome index
STAR --runMode genomeGenerate \
     --genomeDir $base_dir/genome/Arabidospis_thaliana/  \
     --genomeFastaFiles $base_dir/genome/Arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
     --sjdbGTFfile $base_dir/genome/Arabidopsis_thaliana/gtffile \
     --runThreadN #number of available cores on server node


# single end alignment
STAR --genomeDir $base_dir/genome/Arabidospis_thaliana/ \
     --readFilesIn $base_dir/sequences/ \
     --runThreadN 


