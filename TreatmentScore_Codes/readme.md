# Data files and code to calculate treatment score

## What is included here
__Code to calculate treatment score__ <br><br>
__Proteome reference__ : The data is from [1] but the file is provided in [2]. Some gene names have been changed and one gene was added. <br><br>
__Reference of protein-protein affinity__ : The file is provided in [2]. Some gene names have been changed and one gene pair was added. <br><br>
__Expression matrix__ : The datasets used in this study are placed. Prepare a new file when using new datasets. <br>

## How to run the code
Open "calculate_score.py". Then, set target gene(s) and datasets.

## Python packages
The code was run in Python version 3.8.5 with following packages.<br>
numpy==1.23.5 <br>
pandas==1.5.3 <br>
scipy==1.10.1 <br>

## Reference
[1] Rieckmann JC, _et al._ Social network architecture of human immune cells unveiled by quantitative proteomics. Nature Immunology 2017;18(5):583-93 [doi: 10.1038/ni.3693](https://doi.org/10.1038/ni.3693)<br>
[2] Shilts J, _et al._ A physical wiring diagram for the human immune system. Nature 2022;608(7922):397-404.  [doi: 10.1038/s41586-022-05028-x](https://doi.org/10.1038/s41586-022-05028-x)<br>


