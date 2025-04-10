# Dimeric MultiDomain Biosensor Builder

Listed within this repository you will find the necessary code relevant to run the Dimeric MultiDomain Biosensor Builder used to produce chimera of CSF1R and VEGFR2 with TpoR in our work [TBC]. 

Our protocol is fundamentally based on the Rosetta software[1]. Specficially, it includes changes to the domain assembly protocol[2] by introducing constraints. This allows us to apply it to dimeric receptor design (the base version only works with monomers).

A lot of the individual steps are complex and involved; needing significant topological changes in input PDBs, the definition of various changable constraints, the building of arbitrarily changing fasta files and so on. We therefore provide here the means to automate most of the heavy lifting through a set of python and bash scripts. In theory, this should enable you to assemble any dimeric chimeric receptor based on any subset of domains you wish to use. In other words, pre-screen possible candidates that couple orthogonal input/output signals and engineer cellular behaviour.

The various python scripts are built to be "intelligent". They will attempt to be flexible to your needs. Therefore, the main requirement on your end are the exact PDB domains you want to assemble, and the sequences of the linkers you will be reconstructing between the domains.

In addition to the scripts needed to run the assembly protocol, we also provide a comprehensive tutorial/documentation that details exactly how you can run through an example based on SCFR with TpoR, and gives information on the various stages/flags you can use in your own applications. 

# CHANGELOG

14/03/2024 - Bugfixes and general improvements, particularly to documentation. You can now have a ligand be in any of the input protein domains. For the prepare_scaffold flag, you need to provide the domain corresponding to the dimeric domain the ligand is in, rather than just a flat True or False.

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

All Python scripts should work with the most up-to-date version of Python (as of December 2023), but we will assume you are using Python 3.10 to keep it compatible with the work described in our paper. You cannot use Python 2. We recommend first creating a safe conda environment (https://docs.conda.io/projects/conda/en/latest/index.html) to work in. For example,

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

You will also need the following python packages that will need to be installed separately:

- pandas

The necessary package, biobox[3], performs the majority of needed macromolecular, biophysical and topological functions. It can be installed with:

```
conda install biobox -c conda-forge
```

More information on the package can be found at: https://github.com/Degiacomi-Lab/biobox. All necessary python scripts for the protocol are in either the main/ or the lib/ folder.

Since the amount of adaquate sampling required to generate receptors is significant, across numerous stages of assembly, we provide also the means to parallelize stages on a slurm-based cluster using the slurm scripts in HPC_main. These should be adaptable if you are working with a different workload manager on some other HPC. These baseline scripts will need updating, as the paths to Rosetta etc. are not currently provided.

The installation of Rosetta can take up to half a day on a standard desktop PC. Beyond that, the setting up of the conda environment etc. should take 10 minutes at most.

# Tutorial on the builder protocol

A full demo that takes you through the example of building an SCFR-TPoR chimera based on available structural models is provided through the Domain_Assembly_protocol.pdf file. This should be detailed enough to follow and apply to your own chimera, with plenty of information on the various flags needed for each stage of the assembly should you want to deviate. 

For running the example discussed in the tutorial, we provide the necessary structural examples and corresponding input scripts in the example/folder. Beyond this, all you should need to do is update the paths in the various .sh or .slurm scripts you will be using throughout. 

# Quick Guide

A simplified version of the following is given via the integration_test.sh script in the root folder. Simply provide the script with your Rosetta build and it should run the assembly process on the example folder automatically.

```
bash integration_test.sh -R /path/to/rosetta
```

On a single desktop PC, this should take ~1 hour to complete. You may need to repeat it (or bump the numbers up for models generated) if it fails due to too many constraint violations.

The following demonstrates the general pipeline with a short explanation for each flag. For more detail, please refer to the tutorial above. We'll be running the protocol locally and therefore will keep the number of models generated slim, meaning the results will not be particularly meaningful. For proper sampling of the vast conformational space these receptors can adopt, you should generate many more samples - hence the HPC slurm scripts provided and used in the main tutorial. The various mp_assemble_stage etc. bash files we'll be using here are functionally equivalent to the slurm scripts provided in HPC_main.

You will need to choose your structural domains before running through the protocol. The most obvious source are pre-existing native domains from the PDB. However, you can also take advantage of structure prediction generative tools like AlphaFold2 (generate in browser with https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb) or RoseTTAfold (in browser with https://robetta.bakerlab.org/). You'll need to identify the linker sequence which will be rebuilt by Rosetta.

Create a test folder and copy your PDB files into there. You may also need additional linker metadata files depending on what you're doing (see the tutorial for detailed information). I'll continue as though we're using the files provided in the example folder.

```
/location/to/protocol/main/prepare_scaffold.sh -T 4 -s "D13.pdb D45.pdb D6.pdb D7.pdb" -l remove_linkers.txt -d "1 2 4" -a "1 2 3 4" -x add_linkers.txt -L True
```

| Flag        | Meaning          | 
| :-------------: |-------------|
| T | The domain in which the TM is located relative to the order given |
| s | List of domains as a string (in order from N to C termini) |
| l | The name of the remove linkers file, to remove pre-existing linkers in your input structures |
| L | True/False depending on whether a ligand is present in the input |
| d | A list of the input domains where dimerisation is occuring (in order from N to C termini)  |
| x | The name of the add linkers file, containing the sequence information on needed linkers to rebuild  |
| a | Number equivalent to domain being input. Any domain not mentioned will have their linker completely discarded |
 
This will create a new folder called run with the current UNIX time. It is in this folder that you will need to run all of the following scripts (including the slurm scripts if you're using them) besides the clustering, which need to be run from inside generated output folders. Inside this run folder you will find another folder called input_scaffold with the necessary files for assembly. We run the first stage with:

```
/location/to/protocol/main/mp_assembly_stage1.sh -T 4 -R /path/to/Rosetta -d "input_scaffold/D13_cut.pdb input_scaffold/D45_cut.pdb input_scaffold/D6_cut.pdb input_scaffold/D7_cut.pdb"
```

| Flag        | Meaning          |      
|:-------------:|-------------|
| T | The domain in which the TM is located relative to the order given (as above) |
| R | The location of your Rosetta installation |
| d | Based on your current directory location, the full name of your list of input domains |
| N | Number of output structures (6 is the default, 100 is the default in the slurm script) per assembly run |

The runtime for this first stage should be around 10 minutes. 

cd into the output_scaffold. For this quick setup with only 3 models produced it is not really necessary to cluster, but we'll include it here anyway because it is important to the actual protocol.

```
/location/to/protocol/main/remove_constraint_violations.sh
```

This will remove any models that violate the constraints calculated in the prepare_scaffold stage. Now we can cluster, though if you're running the example there is not really any purpose to this.

```
/location/to/protocol/main/run_clustering.sh -R /path/to/Rosetta -S filtered.silent
``` 

This will fail if filtered.silent is empty (i.e. all three models had violated constraints). Since we're aiming to just run through an example, we'll ignore a failure state and just move on. This will likely mean the final stage of assembly where we build linkers will give nonsense answers in terms of the linker conformation. This is another reason why a large amount of sampling is necessary.

In the event that filtered.silent was empty, just take all three models from the initial assembly stage. You'll need to extract them from the silent file for this:

```
/location/to/rosetta/source/bin/extract_pdbs.linuxgccrelease -in:file:silent_struct_type binary -in:file:silent out.silent
```

Rename the three files to c.0.0.pdb, c.1.0.pdb and c.2.0.pdb respectively - our three "clustered centres". Now we prepare for the next assembly stage. Before running this, create a new folder called 0.0 and move the "clustered" files into there. assemble_domains requires a cluster radius to take data samples from so this is just a hack to get around that. Again, this is a symptom of you needing far more sampling than we are applying here.

```
/location/to/protocol/main/assemble_domains.sh -S 2 -d "1 2 4" -s "D6.pdb c.0.0.pdb" -C 0.0 -p 3
```

| Flag        | Meaning          |
|:-------------:|-------------|
| S | Stage of assembly (2 is the default) |
| s | List of domains as a string. Includes the missing domain we're inserting, and the name of the first cluster centre |
| d | A list of the input domains where dimerisation is occuring (in order from N to C termini) |
| C | The cluster folder name based on chosen clustering radius |
| p | The position in the current order of domains where we need to insert our domains, relative to the initial inputs (e.g. 3 when D6.pdb comes after D13.pdb and D45.pdb) |

Now we can begin the actual assembly.

```
/location/to/protocol/main/mp_assemble_stage2.sh -R /path/to/Rosetta -S 2 -T 2 -D "input_scaffold_2/D6_cut.pdb"
```

| Flag        | Meaning          |
|:-------------:|-------------|
| R | The location of your Rosetta installation |
| S | Stage of assembly (2 is the default) |
| T | The domain in which the TM is located relative to where c.n.0.pdb is versus your added domain |
| D | The name of your new added domain |
| N | Number of output structures (6 is the default, 100 is the default in the slurm script) per assembly run |

This will take about 30 minutes to loop across the three "cluster centres" and generate three models each. After this assembly stage has completed, you should ideally cluster again before moving on. We won't do that for this example. Instead, we'll move on to rebuilding the linkers. Again, you should create a dummy 0.0 folder and move your "cluster centres" into there.

First we need to prepare the loop rebuilding files - the final phase of assembly.

```
/location/to/protocol/main/prepare_linkers.sh -C 0.0 -d "1 2 4" -o "3 6 7 8 4 5 1 2" -l 3 -S 2
``` 

| Flag        | Meaning          |
|:-------------:|-------------|
| C | The cluster folder name based on chosen clustering radius |
| d | A list of the input domains where dimerisation is occuring (in order from N to C termini) |
| o | The desired order of domains for the target topology based on the current order of domains - split now as a dimer |
| l | The domain to which the ligand is bound, based on the split set of domains to accomodate the dimer state |
| S | Stage of assembly (2 is the default) |

A full explanation on exactly what these new flags are doing and why they're important is provided in the tutorial.

Now we can actually build these linkers in.

```
/location/to/protocol/main/build_linkers.sh -R /path/to/Rosetta
```

| Flag        | Meaning          |
|:-------------:|-------------|
| R | The location of your Rosetta installation |
| N | Number of output structures per input (1 is the default) |

After this has completed (roughly 30 minutes), you should have your final set of models in output_loop as PDBs. 

You may get a segementation fault with some of the loop rebuilding runs. This occurs because Rosetta cannot fit the loop into the space given, this stems from avoiding culling constraint violation models in previous steps. With much greater sampling on a HPC, you won't have this problem as these violating models will not be present in the final stage.

As a final step, you should relax these models (see www.rosettacommons.org/docs/latest/application_documentation/structure_prediction/relax) to get the true energies, as the nature of this multi-step assembly process does mean you'll have many poorly packed rotamers/contacts.

If you want to use the same assessment approach as we did, you will next need to calculate the free energy of dimerisation and the coupling between residues for each of your receptors, and ideally average these across your different conformers.

# Contact

If you have any questions or concerns about the computational side of the method, please feel free to either submit a ticket above or contact me: Lucas Rudden, lucas.rudden@epfl.ch 

# References

[1] Leaver-Fay et al., _Methods Enzymol_, 2011, (10.1016/B978-0-12-381270-4.00019-6)

[2] Koehler Leman, J. & Bonneau, R., _Biochemistry_, 2017 (10.1021/acs.biochem.7b00995)

[3] Rudden et al., _Bioinformatics_, 2022 (10.1093/bioinformatics/btab785)

# Common Errors

After the running buildlinkers, you might find that your input PDB structures look weird, with residues seemingly "hopping" all over the place (ignoring the actual linker residues at the origin of your system). This means your domain ordering, the `-o` flag, was incorrect. You'll need to shuffle and play with this until it looks correct flag, was incorrect. You'll need to shuffle and play with this until it looks correct.
