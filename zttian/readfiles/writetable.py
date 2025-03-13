import json

def write_table(data, filename, encoding="utf-8", sep="\t"):
  with open(filename, "w", encoding=encoding) as f:
    for i in data:
      for j in i[0:-1]:
        f.write(str(j))
        f.write(sep)
      f.write(str(i[-1]))
      f.write("\n")
