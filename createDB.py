import re

# Read in rclass file
file = open("C:/Databases/KEGG/rclass/rclass", "r")
print("Reading File")
rclass = file.readlines()
rclass = rclass[0:1000]

# Find end of each record
endList = []
for line_i, line in enumerate(rclass, 1):
    if re.match('///', line):
        endList.append(line_i)

# Loop through each record and extract reaction class data
print("Parsing Data")
count2 = 0
for ref in range(0, len(endList)-1):
    count1 = count2
    count2 = endList[ref]
    r = rclass[count1:count2]
    tokenList = ['ENTR', 'DEFI', 'RPAI', 'PATH']
    tokens = {}
    for line_i, line in enumerate(r, 1):
        if line[0:4] in tokenList:
            tokens[line[0:4]] = line_i
    print(tokens)

