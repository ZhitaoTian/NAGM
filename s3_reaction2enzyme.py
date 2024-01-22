from zttian.readfiles.readtable import read_table
from zttian.readfiles.writetable import write_table
from rdkit import Chem
from rdkit.Chem import AllChem

keggStructureT = read_table('./input/structureTable4c.txt')
keggStructureDict = {}
for i in keggStructureT:
  if i[7] != "None":
    keggStructureDict[i[7]] = Chem.MolFromSmiles(i[1])

libraryStructureT = read_table("./input/i1_MetabolitesInfo3_adFomular.txt")
libraryStructureDict = {}
for i in libraryStructureT[1:]:
  libraryStructureDict[i[0]] = Chem.MolFromSmiles(i[9])

enzyme2reaction = read_table("./input/enzyme2reaction.txt")
for i in enzyme2reaction:
  i.append(Chem.MolFromSmarts(i[2].split(">>")[0]))
  i[2] = AllChem.ReactionFromSmarts(i[2])
enzyme2reactionDict = {"0.3": {}, "0.5": {}, "0.7": {}}
for i in enzyme2reaction:
  print(i)
  enzyme2reactionDict[i[4]][i[0]] = i

similarityT = read_table("./output/jxb_Library_similarity.txt")
similarity = []
for i in similarityT:
  i[-1] = float(i[-1])
  if i[-1] > 0.3:
    similarity.append(i)

similarityT2 = read_table("./output/jxb_Library_similarity2.txt")
similarity2 = []
for i in similarityT2:
  i[-1] = float(i[-1])
  if i[-1] > 0.3:
    similarity2.append(i)


def write_result(listT):
  with open("./output/r3_libraryReaction2enzyme.txt", "a", encoding="utf-8") as f:
    listT = [str(i) for i in listT]
    strT = "\t".join(listT)+"\n"
    f.write(strT)


result = []
for i in similarity:
  reactionOver = False
  if i[-1] > 0.3:
    reaction = keggStructureDict[i[0]]
    product = libraryStructureDict[i[1]]
    productSmile = Chem.MolToSmiles(product, isomericSmiles=False)
    for reactionT, reactionList in enzyme2reactionDict["0.3"].items():
      enzyme = reactionList[1]
      if reaction.HasSubstructMatch(reactionList[-1]):
        resT = reactionList[2].RunReactants([reaction])
        for resTT in resT:
          if product.HasSubstructMatch(resTT[0]):
            write_result(i+[enzyme])
            result.append(i+[enzyme])
            reactionOver = True
            break
      if reactionOver == True:
        break
    if i[-1] > 0.5 and reactionOver == False:
      for reactionT, reactionList in enzyme2reactionDict["0.5"].items():
        enzyme = reactionList[1]
        if reaction.HasSubstructMatch(reactionList[-1]):
          resT = reactionList[2].RunReactants([reaction])
          for resTT in resT:
            if productSmile == Chem.MolToSmiles(resTT[0], isomericSmiles=False):
              write_result(i+[enzyme])
              result.append(i+[enzyme])
              reactionOver = True
              break
        if reactionOver == True:
          break
      if i[-1] > 0.7 and reactionOver == False:
        for reactionT, reactionList in enzyme2reactionDict["0.7"].items():
          enzyme = reactionList[1]
          if reaction.HasSubstructMatch(reactionList[-1]):
            resT = reactionList[2].RunReactants([reaction])
            for resTT in resT:
              if productSmile == Chem.MolToSmiles(resTT[0], isomericSmiles=False):
                write_result(i+[enzyme])
                result.append(i+[enzyme])
                reactionOver = True
                break
          if reactionOver == True:
            break


result = []
for i in similarity2:
  reactionOver = False
  if i[-1] > 0.3:
    reaction = libraryStructureDict[i[0]]
    product = libraryStructureDict[i[1]]
    productSmile = Chem.MolToSmiles(product, isomericSmiles=False)
    for reactionT, reactionList in enzyme2reactionDict["0.3"].items():
      enzyme = reactionList[1]
      if reaction.HasSubstructMatch(reactionList[-1]):
        resT = reactionList[2].RunReactants([reaction])
        for resTT in resT:
          if product.HasSubstructMatch(resTT[0]):
            write_result(i+[enzyme])
            result.append(i+[enzyme])
            reactionOver = True
            break
      if reactionOver == True:
        break
    if i[-1] > 0.5 and reactionOver == False:
      for reactionT, reactionList in enzyme2reactionDict["0.5"].items():
        enzyme = reactionList[1]
        if reaction.HasSubstructMatch(reactionList[-1]):
          resT = reactionList[2].RunReactants([reaction])
          for resTT in resT:
            if productSmile == Chem.MolToSmiles(resTT[0], isomericSmiles=False):
              write_result(i+[enzyme])
              result.append(i+[enzyme])
              reactionOver = True
              break
        if reactionOver == True:
          break
      if i[-1] > 0.7 and reactionOver == False:
        for reactionT, reactionList in enzyme2reactionDict["0.7"].items():
          enzyme = reactionList[1]
          if reaction.HasSubstructMatch(reactionList[-1]):
            resT = reactionList[2].RunReactants([reaction])
            for resTT in resT:
              if productSmile == Chem.MolToSmiles(resTT[0], isomericSmiles=False):
                write_result(i+[enzyme])
                result.append(i+[enzyme])
                reactionOver = True
                break
          if reactionOver == True:
            break
