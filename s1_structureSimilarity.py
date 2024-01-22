from copy import deepcopy
from zttian.readfiles.readtable import read_table
from zttian.readfiles.writetable import write_table
import os
import pickle
import sys
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem


fileT = read_table("./input/structureTable4c.txt")
fileT2 = read_table("./input/i1_MetabolitesInfo3_adFomular.txt")
filterList = read_table("./input/reaction_pathway_network4_duplicate.txt")

file2 = []
molecularDict2 = {}
for i in fileT2[1:]:
  try:
    temp = Chem.MolFromSmiles(i[9])
    molecularDict2[i[0]] = temp
    i[1] = AllChem.GetMorganFingerprint(temp, 4, useFeatures=True)
    file2.append(i)
  except:
    pass

file = []
molecularDict = {}
for i in fileT:
  try:
    temp = Chem.MolFromSmiles(i[1])
    molecularDict[i[0]] = temp
    i[1] = AllChem.GetMorganFingerprint(temp, 4, useFeatures=True)
    file.append(i)
  except:
    pass

reactKegg = set()
for i in filterList:
  reactKegg.add(i[1])
  reactKegg.add(i[2])
libraryStructure = []
keggStructure = []
for i in file:
  if i[7] in reactKegg:
    keggStructure.append(i)
  else:
    libraryStructure.append(i)

libraryStructure2 = []
keggStructure2 = []
for i in file2:
  if i[12] in reactKegg:
    keggStructure2.append(i)
  else:
    libraryStructure2.append(i)

result = []
numT = 0
for i in libraryStructure2:
  for j in keggStructure:
    temp2 = round(DataStructs.DiceSimilarity(i[1], j[1]), 2)
    if molecularDict2[i[0]].HasSubstructMatch(molecularDict[j[0]]):
      if temp2 > 0.3:
        numT += 1
        num = "r" + str(numT).zfill(9)
        result.append([num, j[7], i[0]])
    elif temp2 > 0.7:
      numT += 1
      num = "r" + str(numT).zfill(9)
      result.append([num, j[7], i[0]])
    with open("./output/r1_jxb_Library_similarity.txt", "a", encoding="utf-8") as f:
      f.write(j[7]+"\t" + i[0]+"\t"+str(temp2)+"\n")

for ti in range(len(libraryStructure2))[:-1]:
  i = libraryStructure2[ti]
  for tj in range(len(libraryStructure2))[ti+1:]:
    j = libraryStructure2[tj]
    temp2 = round(DataStructs.DiceSimilarity(i[1], j[1]), 2)
    if molecularDict2[i[0]].HasSubstructMatch(molecularDict2[j[0]]):
      if temp2 > 0.3:
        numT += 1
        num = "r" + str(numT).zfill(9)
    elif temp2 > 0.7:
      numT += 1
      num = "r" + str(numT).zfill(9)
    with open("./output/r2_jxb_Library_similarity2.txt", "a", encoding="utf-8") as f:
      f.write(j[0]+"\t" + i[0]+"\t"+str(temp2)+"\n")
      f.write(i[0]+"\t" + j[0]+"\t"+str(temp2)+"\n")
