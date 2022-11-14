Input
===
- receptor: CRBN_CTD_Thalidomide.pdb   
- motif: BetaHairpin.pdb   
Command
===
The parameters used into our case study:   
```
3DMotif-Dock -i ${HumanProteinPDS} -f "unlimited" -3 -d 1.0 -r CRBN_CTD_Thalidomide.pdb -m BetaHairpin.pdb 
```
If you want to re-produce the CRBN substratres database locally,  we recomond:   
```
3DMotif-Dock -i ${HumanProteinPDS} -f 1.64 -3 -d 1.0 -r CRBN_CTD_Thalidomide.pdb -m BetaHairpin.pdb 
```

Output
===
- result table: CRBN_Substrates_Database.xlsx  
- result structure: too large to upload, you can view it at our [database](https://bailab.siais.shanghaitech.edu.cn/services/crbn-subslib).



