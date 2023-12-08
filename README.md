# Dimeric MultiDomain Biosensor Builder

Listed within this repository you will find the necessary code relevant to run the Dimeric MultiDomain Biosensor Builder used to produce chimera of CSF1R and VEGFR2 with TpoR in our work [TBC]. 

Our protocol is fundamentally based on the Rosetta software[1]. Specficially, it includes changes to the domain assembly protocol[2] by introducing constraints. This allows us to apply it to dimeric receptor design (the base version only works with monomers).

A lot of the individual steps are complex and involved; needing significant topological changes in input PDBs, the definition of various changable constraints, the building of arbitrarily changing fasta files and so on. We therefore provide here the means to automate most of the heavy lifting through a set of python and bash scripts. In theory, this should enable you to assemble any dimeric chimeric receptor based on any subset of domains you wish to use. In other words, pre-screen possible candidates that couple orthogonal input/output signals and engineer cellular behaviour.

The various python scripts are built to be "intelligent". They will attempt to be flexible to your needs. Therefore, the main requirement on your end are the exact PDB domains you want to assemble, and the sequences of the linkers you will be reconstructing between the domains.

In addition to the scripts needed to run the assembly protocol, we also provide a comprehensive tutorial/documentation that details exactly how you can run through an example based on SCFR with TpoR, and gives information on the various stages/flags you can use in your own applications. 

# Citation

If this repo proves useful to you, please cite us with:

Rath J.A., Rudden L.S.P., Nouraee N., Von Gunten C., Perez C., Bhugowon Y., Füeglistaler A., Chatzisouleiman A., Barth P., and Arber C., Designed allosteric biosensors for engineered T cell therapy against cancer, https://github.com/barth-lab/Dimeric_MultiDomain_Biosensor_Builder

Paper reference TBC

# Installation

You can download the Rosetta software here: https://www.rosettacommons.org/software/license-and-download

Once you have downloaded Rosetta, our changes to the assembly protocol can be included by replacing the existing source code script for the assembly protocol using the one provided in the files/ folder. This has been tested on the most recent build of Rosetta (08/12/2023), and should be compatible with versions moving backwards to 01/04/2019. 

```
cp ./files/mp_domain_assembly.cc /path/to/Rosetta/main/source/src/apps/public/membrane/mp_domain_assembly.cc
```

Where /path/to/Rosetta denotes your respective path to wherever you have it installed. We are currently in the process of merging these changes into the main branch of Rosetta, and these changes will be available on publication.

Once you have copied over the file, you can build Rosetta following the standard process: https://www.rosettacommons.org/demos/latest/tutorials/install_build/install_build

All Python scripts should work with the most up-to-date version of Python (as of December 2023), but we will assume you are using Python 3.10 to keep it compatible with the work described in our paper. We recommend first creating a safe conda environment (https://docs.conda.io/projects/conda/en/latest/index.html) to work in. For example,

```
conda create --name py310 python=3.10
```

Creates a Python 3.10 environment within GNU/Linux and other unix-based operating systems. This can then be activated:

```
conda activate py310
```

You can deactivate the environment with:

```
conda deactivate
```

And subsequently delete the environment if you desire:

```
conda remove --name py310 --all
```

The following Python packages to run the protocol should be automatically installed when you create the environment:

- numpy
- argparse
- re

The necessary package, biobox[3], performs the majority of needed macromolecular, biophysical and topological functions. It can be installed with:

```
conda install biobox -c conda-forge
```

More information on the package can be found at: https://github.com/Degiacomi-Lab/biobox. All necessary python scripts for the protocol are in either the main/ or the lib/ folder.

Since the amount of adaquate sampling required to generate receptors is significant, across numerous stages of assembly, we provide also the means to parallelize stages on a slurm-based cluster using the slurm scripts in HPC_main. These should be adaptable if you are working with a different workload manager on some other HPC. These baseline scripts will need updating, as the paths to Rosetta etc. are not currently provided.

The installation of Rosetta can take up to half a day on a standard desktop PC. Beyond that, the setting up of the conda environment etc. should take 10 minutes at most.

# Running the builder protocol

A full demo that takes you through the example of building an SCFR-TPoR chimera based on available structural models is provided through the Domain_Assembly_protocol.pdf file. This should be detailed enough to follow and apply to your own chimera, with plenty of information on the various flags needed for each stage of the assembly should you want to deviate. 

For running the example discussed in the tutorial, we provide the necessary structural examples and corresponding input scripts in the example/folder. Beyond this, all you should need to do is update the paths in the various .sh or .slurm scripts you will be using throughout. 

A simplified integration test is being built now (to be completed ~15/12/2023).

# Contact

If you have any questions or concerns about the computational side of the method, please feel free to either submit a ticket above or contact me: Lucas Rudden, lucas.rudden@epfl.ch 

# References

[1] Leaver-Fay et al., _Methods Enzymol_, 2011, (10.1016/B978-0-12-381270-4.00019-6)

[2] Koehler Leman, J. & Bonneau, R., _Biochemistry_, 2017 (10.1021/acs.biochem.7b00995)

[3] Rudden et al., _Bioinformatics_, 2022 (10.1093/bioinformatics/btab785)
