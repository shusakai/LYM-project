# activate environment
conda activate cpdb

# 1: path to seurat_obj.rds, 2: project name, 3: output directly
Rscript /path/to/cellphonedb_scripts/cellphonedb_proc.R $1 $2 $3

# shaping
python /path/to/cellphonedb_scripts/shaping_hcount_meta.py "${3}/${2}_filtered_hcount.txt" "${3}/${2}_filtered_meta.txt"

# run cellphonedb
cellphonedb method statistical_analysis "${3}/${2}_filtered_meta.txt" "${3}/${2}_filtered_hcount.txt"  --threads 30 --output-path "${3}/${2}_cellphonedb_results"
cellphonedb plot dot_plot --means-path "${3}/${2}_cellphonedb_results/means.txt" --pvalues-path "${3}/${2}_cellphonedb_results/pvalues.txt" --output-path "${3}/${2}_cellphonedb_results" --output-name dot_plot.pdf
cellphonedb plot heatmap_plot "${3}/${2}_filtered_meta.txt" --pvalues-path "${3}/${2}_cellphonedb_results/pvalues.txt" --output-path "${3}/${2}_cellphonedb_results/"
