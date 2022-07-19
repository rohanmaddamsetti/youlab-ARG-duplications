# youlab-ARG-duplications 
### by Rohan Maddamsetti and colleagues in Lingchong You Lab.
#### Contact: rohan [dot] maddamsetti [at] Duke

## Python requirements: Python 3.6+, biopython, tqdm 
## R 4.0 is required to make figures and run statistics.

Keep a copy of this repository on your local computer.
Then copy this directory to the Duke Compute Cluster to run the following pipeline.
Copy the output files (in the results directory) to your local machine (see below).
Don't worry about the files in the results/gbk-annotation directory, to save space on your
computer.

After downloading the results of the genome data pipeline, run ARG-duplication-analysis.R
to generate all figures. It's best to run this script interactively in pieces.
Some parts can take a loooong time, depending on how much memory and disk space you
have on your laptop. It's doable with 8GB of RAM and enough disk space for swap, but
you probably want at least 16GB RAM on your laptop, minimum.

Genome Data Pipeline:

If you want to re-run using the newest data available, then
download plasmids.txt and prokaryotes.txt into ../data/GENOME_REPORTS:

wget https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/plasmids.txt  
wget https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt

For replicating the results in the manuscript, skip this step, and just use the versions
of plasmids.txt and prokaryotes.txt that are provided.

Then, filter the prokaryote data for those that have either complete genomes or both
plasmids and chromosomes:

python filter-genome-reports.py > ../results/best-prokaryotes.txt

Then, fetch genome annotation for each row in best-prokaryotes.txt,
and fetch the protein-coding genes for all chromosomes and plasmids for
each row in best-prokaryotes.txt.

Both steps can be done at the same time on the Duke Compute Cluster (DCC).
And make sure these scripts are called from the src directory.
fetch-gbk-annotation and fetch-genome-and-plasmid-cds.py run overnight...

sbatch --mem=16G -t 24:00:00 --wrap="python fetch-gbk-annotation.py"  
sbatch --mem=16G -t 24:00:00 --wrap="python fetch-genome-and-plasmid-cds.py"  

Once the data has downloaded, run the following scripts. Some run
quite quickly, so no need to submit them to a partition on DCC--
just run them in an interactive session on DCC.

python make-chromosome-plasmid-table.py  
python make-gbk-annotation-table.py ## this runs for ~35 min on DCC.  

## this runs for 20 min on DCC. 
python count-cds.py > ../results/protein_db_CDS_counts.csv

## this runs for ~36h on DCC.
sbatch --mem=16G -t 48:00:00 --wrap="python tabulate-proteins.py" 

## this runs for ~36h on DCC.
sbatch --mem=16G -t 48:00:00 --wrap="python tabulate-proteins.py --ignore-singletons"

Then, copy the following files from the results/
directory onto your local machine (same directory name and file structure).

duplicate-proteins.csv  
all-proteins.csv  
protein_db_CDS_counts.csv  
gbk-annotation-table.csv  
chromosome-plasmid-table.csv  
protein_db.faa  
prokaryotes-with-plasmids.txt  

Then, run the follow scripts to annotate the genomes, and to cross-check
the computational annotation against a subset of annotations that were conducted manually.

python annotate-ecological-category.py > ../results/computationally-annotated-gbk-annotation-table.csv  

python check-ecological-annotation.py

If everything checks out, you are now ready to run ARG-duplication-analysis.R to
generate figures!

Figure 4ABC of the manuscript uses the Pluto notebook (in julia) for simulations:
duplication-linear-ODE-model-v4.jl

To run this notebook, open julia, and type:
using Pluto; Pluto.run()
Then choose this notebook and run it in the browser.

You will have to install the julia package dependencies (like Pluto) to get this to run--
see the source code for details.

Figure 4D is generated using the R script qPCR-analysis.R.

