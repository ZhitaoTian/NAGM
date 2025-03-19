from zttian.readfiles.readtable import read_table
from zttian.readfiles.writetable import write_table
import os
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import stats
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem

#构建反应模式框架[reactionID, enzyme, reactantsList, productsList]
edges = read_table("./output/r4_all_edge.txt")
nodes = read_table("./output/r4_all_node.txt")
enzyme2ko = read_table("./input/enzyme2ko.txt")
ko2reaction = read_table("./input/ko_reaction.txt")
gene2ko = read_table("./input/i8_gene2ko.txt")
enzyme2reaction = read_table("./input/kegg_enzyme2reaction.txt")
metabolites = pd.read_table("./input/i2_Metabolites.txt")
libraryStructureT = read_table("./input/i1_MetabolitesInfo3_adFomular.txt")
keggStructureT = read_table('./input/structureTable4c.txt')
keggCompound2name = read_table("./input/kegg_compound2name.txt")


reaction2enzymeDict = {}
for i in enzyme2reaction:
  enzymeTemp = i[0].split(".")[:3]
  if enzymeTemp[0] == "2" and enzymeTemp[1] == "3":
    enzymeTemp = "2.3"
  else:
    enzymeTemp = ".".join(enzymeTemp[:3])
  try:
    reaction2enzymeDict[i[1]].add(enzymeTemp)
  except:
    reaction2enzymeDict[i[1]] = set([enzymeTemp])
print("job")
reactionsDict = {}
numberT = 0
errReaction = []
for i in edges[1:]:
  if i[3] == "kegg":
    try:
      reactionsDict[i[0]][2].add(i[1])
      reactionsDict[i[0]][3].add(i[2])
    except:
      try:
        reactionsDict[i[0]] = [
            i[0], reaction2enzymeDict[i[0]], set([i[1]]), set([i[2]]), i[3]]
      except:
        errReaction.append(i)
        pass
  else:
    numberT += 1
    modificationIdT = "M"+str(numberT).zfill(5)
    reactionsDict[modificationIdT] = [
        modificationIdT, set([i[0]]), set([i[1]]), set([i[2]]), i[3]]
print("job")
nodesDetected = set([i[0] for i in nodes if i[3] == "detected"])
reactionDetected = []
reactionKegg = []
for i, j in reactionsDict.items():
  for jj in j[3]:
    if jj in nodesDetected:
      reactionDetected.append(j)
      break
  if i.startswith("R"):
    reactionKegg.append(j)
print("job")
nodesId2SmilesDict = {i[0]: i[4] for i in nodes}


def similarityCount(a, b):
  temp1 = Chem.MolFromSmiles(a)
  temp2 = Chem.MolFromSmiles(b)
  temp1 = AllChem.GetMorganFingerprint(temp1, 4, useFeatures=True)
  temp2 = AllChem.GetMorganFingerprint(temp2, 4, useFeatures=True)
  return DataStructs.DiceSimilarity(temp1, temp2)


reactionSimilarity = []
for reactionT1 in reactionDetected:
  for reactionT2 in reactionKegg:
    if len(reactionT1[1] & reactionT2[1]) > 0:
      reactantListT = []
      productListT = []
      for reactantT1 in reactionT1[2]:
        for reactantT2 in reactionT2[2]:
          similarityT = similarityCount(
              nodesId2SmilesDict[reactantT1], nodesId2SmilesDict[reactantT2])
          reactantListT.append([reactantT1, reactantT2, similarityT])
      for productT1 in reactionT1[3]:
        for productT2 in reactionT2[3]:
          similarityT = similarityCount(
              nodesId2SmilesDict[productT1], nodesId2SmilesDict[productT2])
          productListT.append([productT1, productT2, similarityT])
      reactantListT.sort(key=lambda x: -x[2])
      productListT.sort(key=lambda x: -x[2])
      if productListT[0][2] > 0.4 and reactantListT[0][2] > 0.2:
        reactantSet = set()
        productSet = set()
        reactantSimilarityT = []
        productSimilarityT = []
        for i in reactantListT:
          if i[0] not in reactantSet and i[1] not in reactantSet:
            reactantSimilarityT.append(i[2])
            reactantSet.add(i[0])
            reactantSet.add(i[1])
        for i in productListT:
          if i[0] not in productSet and i[1] not in productSet:
            productSimilarityT.append(i[2])
            productSet.add(i[0])
            productSet.add(i[1])
        SimilarityT = ((sum(reactantSimilarityT)/len(reactantSimilarityT))
                       * (sum(productSimilarityT)/len(productSimilarityT)))**0.5
        reactionSimilarity.append([reactionT1[0], reactionT2[0], SimilarityT])


# 建立反应与ko号以及基因的字典
reaction2koDict = {}
for i in ko2reaction:
  try:
    reaction2koDict[i[1]].append(i[0])
  except:
    reaction2koDict[i[1]] = [i[0]]
ko2geneDict = {}
gene2koDict = {}
for i in gene2ko:
  try:
    ko2geneDict[i[1]].append(i[0])
  except:
    ko2geneDict[i[1]] = [i[0]]
gene2koDict = {i[0]: i[1] for i in gene2ko}
reaction2geneDict = {}
for i, j in reaction2koDict.items():
  for jj in j:
    try:
      temp = deepcopy(ko2geneDict[jj])
      try:
        reaction2geneDict[i].extend(temp)
      except:
        reaction2geneDict[i] = temp
    except:
      pass
# 根据反应相似性结果添加相应的酶号以及基因号
# 预测反应，比对的反应，反应相似性，ko号，基因号
for i in reactionSimilarity:
  try:
    i.append(",".join(reaction2koDict[i[1]]))
  except:
    i.append("")
  try:
    i.append(",".join(reaction2geneDict[i[1]]))
  except:
    i.append("")

reactionSimilarityDict = {}
for reactionSimilarityT in reactionSimilarity:
  reactionName = reactionSimilarityT[0]
  score = round(reactionSimilarityT[2], 2)
  if len(reactionSimilarityT[4]) > 0:
    genes = reactionSimilarityT[4].split(",")
    for gene in genes:
      if reactionName not in reactionSimilarityDict:
        reactionSimilarityDict[reactionName] = {}
        reactionSimilarityDict[reactionName][gene] = [
            reactionSimilarityT[1], gene2koDict[gene], score]  # 相似反应ID,KO,得分
      else:
        try:
          if score > reactionSimilarityDict[reactionName][gene][2]:
            reactionSimilarityDict[reactionName][gene] = [
                reactionSimilarityT[1], gene2koDict[gene], score]  # 相似反应ID,KO,得分
        except:
          reactionSimilarityDict[reactionName][gene] = [
              reactionSimilarityT[1], gene2koDict[gene], score]  # 相似反应ID,KO,得分


# 转录与代谢相关性计数表格整理

metabolites.index = metabolites['ID']
metabolites = metabolites.iloc[:, 1:]
for c in metabolites.columns:
  metabolites[c] = np.log2(metabolites[c]+1)


def groupMean(dataT, step):
  lenT = len(dataT.columns)
  dataT2 = pd.DataFrame()
  for i in range(lenT)[::step]:
    numT = int(i/step)
    dataT2[numT] = dataT.iloc[:, i:i+3].mean(axis=1)
  return dataT2


metabolites2 = groupMean(metabolites, 3)
metabolites2 = metabolites2.T
transcripts = pd.read_table("./input/DESeq_count.txt")
transcripts.index = transcripts.iloc[:, 0]
transcripts = transcripts.iloc[:, 1:]
transcripts = transcripts[metabolites.columns]
transcripts2 = groupMean(transcripts, 3)
transcripts2 = transcripts2.T

nodeId2metaboliteIdDict = {i[0]: i[1] for i in nodes if i[3] == "detected"}

results3 = deepcopy(reactionSimilarityDict)
for key, value in results3.items():
  reactants = reactionsDict[key][2]
  products = reactionsDict[key][3]
  for i in products:
    if i in nodeId2metaboliteIdDict:
      metaboliteT = nodeId2metaboliteIdDict[i]
      break  # 只挑选了一个产物
  for geneT, value1 in value.items():
    try:
      rT, pT = stats.spearmanr(metabolites2[metaboliteT], transcripts2[geneT])
      value1.extend([round(rT, 2), pT, round(rT*value1[2], 2)])
    except:
      value1.extend([0, 1, 0])
# results3
# key1 reactionId
# key2 geneId
# value (similarReaction ko reactionSimilarityScore gene_metabolite_r gene_metabolite_p score)
metaboliteId2nameDict = {i[0]: i[6] for i in libraryStructureT}
id2nameDict = {}
for i in keggCompound2name:
  id2nameDict[i[0]] = i[1]
for i, j in nodeId2metaboliteIdDict.items():
  id2nameDict[i] = metaboliteId2nameDict[j]


results4 = []
results5 = []
for i, j in results3.items():
  reactionT = reactionsDict[i]
  reactantName = []
  productName = []
  for ii in reactionT[2]:
    reactantName.append(id2nameDict[ii])
  for ii in reactionT[3]:
    productName.append(id2nameDict[ii])
  if reactionT[0].startswith("R"):
    reactionT[0] = '=HYPERLINK("https://www.kegg.jp/entry/' + \
        reactionT[0]+'","'+reactionT[0]+'")'
  temp1 = [reactionT[0],
           ";".join(reactionT[1]),
           ";".join(reactionT[2]),
           ";".join(reactantName),
           ";".join(reactionT[3]),
           ";".join(productName),
           reactionT[4]]
  temp1b = [reactionT[0],
            ";".join(reactionT[1]),
            ";".join(reactionT[2]),
            ";".join(reactantName),
            ";".join(reactionT[3]),
            ";".join(productName),
            reactionT[4]]
  temp2 = []
  for ii, jj in j.items():

    temp2.append([ii]+jj)
  temp2 = [iT for iT in temp2 if iT[6] >
           0.5 and iT[5] < 0.05]  # 总分大于0.5，相关性pvalue<0.05
  temp2.sort(key=lambda x: -abs(x[6]))
  numberT = 0
  for temp3 in temp2:
    numberT += 1
    temp3[1] = '=HYPERLINK("https://www.kegg.jp/entry/' + \
        temp3[1]+'","'+temp3[1]+'")'
    temp1.extend([numberT]+temp3)
    results5.append(temp1b+[numberT]+temp3)
  results4.append(temp1)
# results4.insert(0, ["反应", "酶", "底物","底物名称",  "产物","产物名称", "反应类型", "排名", "基因", "基因在kegg中的反应", "基因的ko号",
#                     "反应相似性得分", "基因与检测产物的spearman相关系数", "相关系数显著性值", "综合指数"])
results4.insert(0, ["Reaction (R0)", "Enzyme", "Reactant (r0)", "r0 name", "Product (p0)", "p0 name", "Reaction type", "Rank", "Gene (G1)", "Reaction for G1 (R1)", "KO number of the G1", "Reaction similarity between R1 and R0", "Correlation between G1 and p0", "significance of the correlation between G1 and p0",
                    "Comprehensive index for G1", "Rank", "Gene (G2)", "Reaction for G2 (R2)", "KO number of the G2", "Reaction similarity between R2 and R0", "Correlation between G2 and p0", "significance of the correlation between G2 and p0", "Comprehensive index for G2"])
print("job over")
write_table(results4, "./output/r5_candidateGene.txt")
write_table(results5, "./output/r5_candidateGeneCount.txt")
