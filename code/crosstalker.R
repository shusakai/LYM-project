library(CrossTalkeR)

# the method always consider the first path as control: the multiple control case will be handle soon
paths <- c(
  'Pre' = 'path/to/pre_cellphonedb_results/filtered_corrected.csv',
  'During' = 'path/to/during_cellphonedb_results/filtered_corrected.csv',
  'JustFinish' = 'path/to/justfinish_cellphonedb_results/filtered_corrected.csv',
  'After' = 'path/to/after_cellphonedb_results/filtered_corrected.csv'
)


# Generating the report and the object
data <- CrossTalkeR::generate_report(
  lrpaths=paths,
  threshold=0,
  out_path='path/to/', 
  out_file='All_time.html', 
  output_fmt = "html_document", 
 # report = FALSE
)


# compare -----------------------------------------------------------------------------
# pre vs during
paths <- c(
  'Pre' = 'path/to/pre_cellphonedb_results/filtered_corrected.csv',
  'During' = 'path/to/during_cellphonedb_results/filtered_corrected.csv'
)

# Selected gene list     
genes <- c('Type II IFNR|R', "CTLA4|L", "SIRPA|R", "CD274|L", "CD28|L", "TGFB1|L", "TIGIT|L", "CXCL10|L", "CXCL9|L", "CXCL13|L")

# Generating the report and the object
data <- CrossTalkeR::generate_report(
  genes = genes,
  lrpaths=paths,
  threshold=0,
  out_path='path/to/output', 
  out_file='pre_vs_during.html', 
  output_fmt = "html_document", 
  # report = FALSE
)

# during vs justfinish
paths <- c(
  'During' = 'path/to/during_cellphonedb_results/filtered_corrected.csv', 
  'JustFinish' = 'path/to/justfinish_cellphonedb_results/filtered_corrected.csv'
)

# Selected gene list     
genes <- c('Type II IFNR|R', "CTLA4|L", "SIRPA|R", "CD274|L", "CD28|L", "TGFB1|L", "TIGIT|L", "CXCL10|L", "CXCL9|L")

# Generating the report and the object
data <- CrossTalkeR::generate_report(
  genes = genes,
  lrpaths=paths,
  threshold=0,
  out_path='path/to/output', 
  out_file='during_vs_justfinish.html', 
  output_fmt = "html_document", 
  # report = FALSE
)

# during vs justfinish
paths <- c(
  'JustFinish' = 'path/to/justfinish_cellphonedb_results/filtered_corrected.csv', 
  'After' = 'path/to/after_cellphonedb_results/filtered_corrected.csv'
)

# Selected gene list     
genes <- c('Type II IFNR|R', "CTLA4|L", "SIRPA|R", "CD274|L", "CD28|L", "TGFB1|L", "TIGIT|L", "CXCL10|L", "CXCL9|L")

# Generating the report and the object
data <- CrossTalkeR::generate_report(
  genes = genes,
  lrpaths=paths,
  threshold=0,
  out_path='path/to/output', 
  out_file='justfinish_vs_after.html', 
  #output_fmt = "html_document", 
  # report = FALSE