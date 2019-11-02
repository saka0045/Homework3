The python script phylogeny.py takes the hw3.fna fasta file as input and produces the genetic distances file (genetic_distances.txt), 
tab delimited edges file (edges.txt) and the tree in NEWICK format (tree.txt). Example usage of code:

python phylogeny.py -i hw3.fna

The output files will be produced in the same directory where the script is. This script must be executed in Linux or Mac machine.

The provided R scripts are used to visualize the tree. The older R script for edges file must be used to produce the correct tree with the edges.txt file.
This script needs the ape package installed:

Rscript hw3-plot-edges.r edges.txt hw3-tip-labels.txt

The tree with the NEWICK format can be made with the following command:

Rscript hw3-plot-newick.r tree.txt hw3-tip-labels.txt 

The output tree is named tree.pdf
