# domain_assembly_constraints
Assembly procotol to construct dimeric chimeras built from diverse structural native monomeric units, coupling orthgonal input/output signals. (in construction)

Builds on previous mp_domain_assembly method within Rosetta[1], but introduces constraints to enable more complex dimeric chimera receptor constructions and linker optimisation between individual domains.

The actual mp_domain_assembly methods are with the main branch of Rosetta. Here, you will find assembly procotols to automatically construct a scaffold based on input domains, built to be "intelligent" in the sense it can perform the majority of the construction for you. The main input requirement is the exact PDB domains you want to assembly, as well as information on which linkers between domains you want to reconstruct.

Dr. Lucas S.P. Rudden, 02/2022

[1] Koehler Leman, J. & Bonneau, R. A novel domain assembly routine for creating full-length models of membrane proteins from known domain structures. Biochemistry acs.biochem.7b00995 (2017). doi:10.1021/acs.biochem.7b00995
