**Please visit the Lazarus webpage. . . https://dl.dropbox.com/u/23201414/software/lazarus/tutorial.html**

Lazarus is a collection of Python scripts for running the programs codeml and baseml in the PAML software suite. Lazarus allows you to run PAML for large distributions of phylogenies, and then integrate over those trees to calculate maximum a posteriori ancestral sequences. Even if you're dealing with only one tree, Lazarus is still useful as it provides an easy way to compose and launch PAML jobs. Lazarus provides an easy-to-use terminal interface that mimics the program PhyML; if you're comfortable using PhyML, then using Lazarus should be straightforward.

Given a sequence alignment (in FASTA format), a phylogenetic tree (in Newick format), and an evolutionary model (provided by PAML), Lazarus uses codeml or baseml to reconstruct the marginal distribution of ancestral states at all internal nodes of the given tree. Lazarus then parses the output from codeml or baseml, and produces a collection of text files, each containing the posterior probability ancestral state distribution for a single ancestral node.


Relevant Citations:

Maximum Likelihood Ancestral Sequence Reconstruction:
Yang et al., Genetics 1995

Marginal reconstruction:
Koshi and Goldstein, Journal of Molecular Evolution 1996

PAML, version 4:
Yang, Molecular Biology and Evolution, 2007

On Phylogenetic Uncertainty:
Hanson-Smith et al., Molecular BIology and Evolution, 2010


Caveat:

I hope you find Lazarus useful, but be aware that I do not have a large team of software developers. I am providing this software as a resource to the research community, but without the promise of support. If you have questions, email Victor Hanson-Smith (victorhansonsmith at gmail dot com), and I'll try to help you as thoroughly as possible.