import pandas as pd
import sys

# path to the file
path = sys.argv[1]
path_meta = sys.argv[2]
temp = pd.read_table(path).set_index("ensembl_gene_id")
temp_meta = pd.read_table(path_meta)
temp_meta = temp_meta.iloc[:-1, :]
temp.to_csv(path, sep="\t")
temp_meta.to_csv(path_meta, sep="\t", index=None)

