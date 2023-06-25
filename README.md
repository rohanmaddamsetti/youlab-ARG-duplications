# youlab-ARG-duplications 
### by Rohan Maddamsetti and colleagues in Lingchong You Lab.
#### Contact: rohan [dot] maddamsetti [at] Duke

## Python requirements: Python 3.6+, biopython, tqdm 
## R 4.0 is required to make figures and run statistics.
## Julia 1.8+ is required for some stuff too.

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

## Genome Data Pipeline:

If you want to re-run using the newest data available, then
download prokaryotes.txt into ../data/GENOME_REPORTS:

wget https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt

For replicating the results in the paper, skip this step, and just use the
prokaryotes.txt provided in the data/ folder. This one is dated March 21 2021
and was used for the paper.

### See the directions in src/README.md for how to run the bioinformatics pipeline and analysis.

### You will have to download the mobileOG-db database into a folder called "../data/mobileOG-db_beatrix-1-6_v1_all" relative to the directory containing this source code.

This database can be downloaded from the following link:  
https://mobileogdb.flsi.cloud.vt.edu/entries/database_download  

### If everything checks out, you are now ready to run ARG-duplication-analysis.R to generate figures!

## Mathematical modeling results.

Figure 1ABC of the manuscript uses the Pluto notebook (in julia) for simulations:
duplication-linear-ODE-model-v5.jl

To run this notebook, open julia, and type:
using Pluto; Pluto.run()
Then choose this notebook and run it in the browser.

You will have to install the julia package dependencies (like Pluto) to get this to run--
see the source code for details.

