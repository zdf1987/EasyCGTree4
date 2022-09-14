====== EasyCGTree ======
	Version 4.0
by Dao-Feng Zhang (Hohai University, China)

EasyCGTree is An Easy Tool for Constructing Core-Gene Tree.

EasyCGTree is a Perl script, developed to construct genome-based Maximum-Likehood (ML) phylogenetic tree, by taking microbial genomic data (proteome) in fasta format as input data. Profile hidden Markov models (HMMs) of core gene sets prepared in advance and enclosed in the package are used for homologs search by HMMER (http://hmmer.org/), and customized gene sets can also be used to build profile HMMs by EasyCGTree for homologs search. It has integrated all the steps needed between the input data and the resulted tree file into one Perl script, which would make it easier to infer a core-gene tree.

For the installation and usage of EasyCGTree, please refer to the Manual.

There is also a brief tutorial in the Manual.
-=- Change Log 4.0 -=-

1. Replace blast+ with HMMER to search homologs.
2. Add trimAl for alignment trimming.
3. Add IQ-TREE as an optional application for phylogeny inference.
4. Optimize the command line setting and ERROR report.

-=- Change Log 3.1 -=-

1. Add support of tree inference of SNPs and supertree aproaches.
2. Optimize the command line setting and ERROR report.


-=- Change Log 2.2 -=-
1.A new script (GetRepRef.pl) was developed for Linux version to screen representative query data sets. It will efficiently reduce time cost when using it before a run.

-=- Change Log 2.1 -=-
1. Add support of tree inference from amino acid sequence.
2. Optimize the command line setting.
3. Fix the bug in GetReferencFromGTDB.pl and FormatNames.pl.
