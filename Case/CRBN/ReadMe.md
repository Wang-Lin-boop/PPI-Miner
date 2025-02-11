Input
===
- receptor: CRBN_CTD_Thalidomide.pdb   
- motif: BetaHairpin.pdb     

Command
===
The parameters used into our case study:   
```
3DMotif-Dock -l ${HumanProteinPDS} -f "unlimited" -3 -d 1.0 -r CRBN_CTD_Thalidomide.pdb -m BetaHairpin.pdb 
```
If you want to re-produce the CRBN substratres database locally,  we recommend:   
```
3DMotif-Dock -l ${HumanProteinPDS} -f 1.64 -3 -d 1.0 -r CRBN_CTD_Thalidomide.pdb -m BetaHairpin.pdb 
```
If you have sufficient computing resources, we recommend you relax all the human protein structures from AlphaFold DB and PDB at the frist. We believe that this will greatly improve the quality of the output structures.    


Output
===
- result table: CRBN_Substrates_Database.xlsx  
- result structure: too large to upload, you can view it at our [database](https://bailab.siais.shanghaitech.edu.cn/services/crbn-subslib).



