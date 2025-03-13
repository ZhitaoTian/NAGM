import json
def read_table(file_name,encoding="utf-8",sep="\t"):
    with open(file_name, 'r',encoding=encoding) as f:
        content = f.read()
    final_list = list()
    rows = content.split('\n')
    for row in rows:
        final_list.append(row.split(sep))
    del final_list[-1]
    return final_list
