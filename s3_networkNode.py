from zttian.readfiles.readtable import read_table
from zttian.readfiles.writetable import write_table
import os

pathway2reaction = read_table(
    "./input/reaction_pathway_network2_libraryFilter1.txt")

gene2ko = read_table("./input/i8_gene2ko.txt")
enzyme2ko = read_table("./input/enzyme2ko.txt")
ko2reaction = read_table("./input/ko_reaction.txt")
pathway2compoundT = read_table("./input/pathway2compound.txt")

compound2koT = read_table("./output/r1_jxb_Library_similarity.txt")
compound2koDict = {}
for i in compound2koT:
  try:
    if float(i[2]) > compound2koDict[i[1]][1]:
      compound2koDict[i[1]] = [i[0], i[2]]
  except:
    compound2koDict[i[1]] = [i[0], i[2]]
pathway2compound = []
for i, j in compound2koDict.items():
  pathwayT = [ii[0] for ii in pathway2compoundT if ii[1] == j[0]]
  for k in pathwayT:
    pathway2compound.append([k, i])


reaction2koDict = {}
for i in ko2reaction:
  try:
    reaction2koDict[i[1]].append(i[0])
  except:
    reaction2koDict[i[1]] = [i[0]]
ko2geneDict = {}
for i in gene2ko:
  try:
    ko2geneDict[i[1]].append(i[0])
  except:
    ko2geneDict[i[1]] = [i[0]]

reaction2gene = {}
for i, j in reaction2koDict.items():
  for jj in j:
    try:
      temp = ko2geneDict[jj]
      try:
        reaction2gene[i].extend(temp)
      except:
        reaction2gene[i] = temp
    except:
      pass


enzyme2koDict = {}
for i in enzyme2ko:
  try:
    enzyme2koDict[i[0]].append(i[1])
  except:
    enzyme2koDict[i[0]] = [i[1]]
enzyme2gene = {}
for i, j in enzyme2koDict.items():
  for jj in j:
    try:
      temp = ko2geneDict[jj]
      try:
        enzyme2gene[i].extend(temp)
      except:
        enzyme2gene[i] = temp
    except:
      pass

keggReaction = read_table(
    "./input/reaction_pathway_network2_libraryFilter2.txt")
keggLibraryReaction = read_table("./output/r3_libraryReaction2enzyme.txt")
reactionNode = set()
for i in keggReaction:
  reactionNode.add(i[1])
  reactionNode.add(i[2])
for i in keggLibraryReaction:
  reactionNode.add(i[0])
  reactionNode.add(i[1])


keggStructureT = read_table('./input/structureTable4c.txt')
keggStructureT = [i for i in keggStructureT if i[7] != "None"]
keggStructureDict = {}
for i in keggStructureT:
  if i[7] != "None":
    keggStructureDict[i[7]] = i[1]

libraryStructureT = read_table("./input/i1_MetabolitesInfo3_adFomular.txt")
libraryStructureDict = {}
nodeResult = []
includedKegg = set()
for i in libraryStructureT[1:]:
  libraryStructureDict[i[0]] = i[9]
  if i[12] != "" and i[12] in reactionNode:
    nodeResult.append([i[12], i[0], "keggId", "detected", i[9]])
    includedKegg.add(i[12])
  else:
    nodeResult.append([i[0], i[0], "libraryId", "detected", i[9]])

for i in keggStructureT:
  if i[7] not in includedKegg:
    nodeResult.append([i[7], i[7], "keggId", "notDetected", i[1]])


def selectPathway(selectPathway, name):
  reactionNum = [i[1] for i in pathway2reaction if i[0] in selectPathway]
  compoundList = [i[1] for i in pathway2compound if i[0] in selectPathway]
  selectEdge = [i for i in allEdge if i[0] in reactionNum]
  selectNode = set()
  for i in selectEdge:
    selectNode.add(i[1])
    selectNode.add(i[2])
  for i in keggLibraryReaction:
    if i[0] in selectNode:
      selectEdge.append([i[3], i[0], i[1], "putative"])
      selectNode.add(i[1])
    elif i[0] in compoundList or i[1] in compoundList:
      selectEdge.append([i[3], i[0], i[1], "putative"])
      selectNode.add(i[0])
      selectNode.add(i[1])
  selectNode = [i for i in nodeResult if i[0] in selectNode]
  selectEdge.insert(0, ["reaction", "reactant", "product", "type"])
  selectNode.insert(0, ["id", "name", "type1", "type2", "smiles"])
  write_table(selectEdge, os.path.join("./output/", "r4_"+name+"_edge.txt"))
  write_table(selectNode, os.path.join("./output/", "r4_"+name+"_node.txt"))


allEdge = [i+["kegg"] for i in keggReaction]
for i in keggLibraryReaction:
  allEdge.append([i[3], i[0], i[1], "putative"])
nodeset = set()
for i in allEdge:
  nodeset.add(i[1])
  nodeset.add(i[2])
allNode = set(nodeset)
nodeResultSet = set(i[0] for i in nodeResult)
differentNode = list(allNode - nodeResultSet)
for i in differentNode:
  nodeResult.append([i, i, "keggId", "notDetected", ""])

nodeNetwork = [i for i in nodeResult if i[0] in allNode]

allEdge.insert(0, ["reaction", "reactant", "product", "type"])
nodeResult.insert(0, ["id", "name", "type1", "type2", "smiles"])
write_table(allEdge, os.path.join("./output/", "r4_all_edge.txt"))
write_table(nodeResult, os.path.join("./output/", "r4_all_node.txt"))
# write_table(nodeNetwork, os.path.join("./output/", "r4_alledge_node.txt"))

# 选择输出特定的代谢途径的反应网络，结合Cytoscape软件(https://cytoscape.org/)进行查看
# selectPathway1 = ["map00941", "map00943", "map00944"]
# selectPathway2 = ["map00330", "map00940"]
# selectPathway3 = ["map00945"]
# selectPathway(selectPathway1, "flavonoid")
# selectPathway(selectPathway2, "Arginine")
# selectPathway(selectPathway3, "Stilbenoid_diarylheptanoid_gingerol")
