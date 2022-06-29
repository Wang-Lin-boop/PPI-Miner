# Motif-Dock
Motif-Dock is a protein-protein interaction (PPI) search pipeline for proteins containing **specific substructures (3D)** or **sequence motifs (2D)**. Motif-Dock can be used to search for potential PPI partners for specific receptors or protein scaffolds fusing specific motifs such as antigenic epitopes.  


Introduction
----

A number of proteins can bind to the same receptor at the same site with similar binding patterns. Most of the time, these proteins belong to the same protein family and have a high degree of sequence similarity. However, these proteins can sometimes be scattered among a diverse family of proteins, maintaining sequence conservation for a short conserved sequence or even no apparent sequence conservation. For one thing, the **short conserved sequence** between the disordered protein and the receptor is called **sequence motif**, such as the TRAF6 and its substrates (shown in Figure 1A). For another thing, some **small sub-structures**, called **structural motifs**, are no apparent sequence conservation and scattered among different protein families, such as CRBN and its substrate proteins with different 3D overall structures, as shown in Figure 1B. The development of Motif-Dock was inspired by these PPIs that are difficult to search by sequence alignment methods (e.g. Blast).  

![**Figure 1. Kinds of PPI protein motifs.**](https://user-images.githubusercontent.com/58931275/174751397-d529dfaf-f970-43f2-a0fe-0f3d99c006f7.png)  
**Figure 1. Kinds of PPI protein motifs.**  

Installation
----
**1. Install the dependent packages.**  

a. Intsall the [Rosetta](https://www.rosettacommons.org/docs/latest/build_documentation/Build-Documentation) packages and set up environment variables.  
```
echo "export rosetta_app=<path to rosetta apps>" >> ~/.bashrc
echo "export rosetta_db=<path to rosetta database>" >> ~/.bashrc
echo "export rosetta_version=<rosetta versions, mpi or static>" >> ~/.bashrc
```  

b. Intsall Julia and related packages.   

```  
wget https://julialang-s3.julialang.org/bin/linux/x64/1.5/julia-1.5.3-linux-x86_64.tar.gz
tar zxvf julia-1.5.3-linux-x86_64.tar.gz  
cd julia-1.5.3/bin
echo "export PATH=${PWD}:\$PATH" >> ~/.bashrc
source ~/.bashrc
julia
] # in Julia REPL
add BioStructures # in Julia REPL
add PDBTools # in Julia REPL
exit()
```
**2. Download Motif-Dock and configure environment variables.**   
```
git clone https://github.com/Wang-Lin-boop/Motif-Dock/scripts/
cd Motif-Dock/scripts/
echo "export PATH=${PWD}:\$PATH" >> ~/.bashrc
source ~/.bashrc
```


Database
----
1. Download the library of prepared human protein structures and disorder sequeneces.

2. Use the 3DMotif-Dock to generate your own structure libraries.

3. Use the scripts to generate your own disorder sequences libraries.

Usage of 2DMotif-Dock
----


Usage of 3DMotif-Dock
----


Citation
----


Acknowledgements
----

