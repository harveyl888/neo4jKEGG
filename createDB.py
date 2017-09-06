#!/usr/bin/python

import re


def main():
    rclass = kegg_reactions("C:/Databases/KEGG/rclass/rclass")
    compounds = kegg_compounds("C:/Databases/KEGG/compound/compound")
    return 0


def kegg_reactions(filename):
    # Read in rclass file
    file = open(filename, "r")
    print("Reading Reaction File")
    rclass = file.readlines()

    # Find end of each record
    endList = []
    for line_i, line in enumerate(rclass, 1):
        if re.match("///", line):
            endList.append(line_i)

    # Loop through each record extracting tokens
    print("Parsing Reaction Data")
    count2 = 0
    rclassData = []
    for ref in range(0, len(endList)-1):
        count1 = count2
        count2 = endList[ref]
        r = rclass[count1:count2]
        tokenList = ["ENTR", "DEFI", "RPAI", "PATH"]
        tokens = dict()
        for line_i, line in enumerate(r, 1):
            if line[0:4] in tokenList:
                tokens[line[0:4]] = line_i - 1
        # Parse tokens, extracting information into a dictionary
        if "ENTR" in tokens.keys():
            record = dict()
            record["entry"] = r[tokens["ENTR"]][12:19]
            if "DEFI" in tokens.keys():
                foundEnd = False
                count = 0
                v_definition = []
                while not foundEnd:
                    currentLine = r[tokens["DEFI"] + count]
                    if currentLine[0:4] in ["DEFI", "    "]:
                        v_definition.append(currentLine[12:].rstrip())
                        count += 1
                    else:
                        foundEnd = True
                record['definition'] = v_definition
            else:
                record['definition'] = []
            if "RPAI" in tokens.keys():
                foundEnd = False
                count = 0
                v_rpair = []
                while not foundEnd:
                    currentLine = r[tokens["RPAI"] + count]
                    if currentLine[0:4] in ["RPAI", "    "]:
                        v_rpair.extend(re.findall("([C0-9_]+)", currentLine[12:]))
                        count += 1
                    else:
                        foundEnd = True
                record['rpairs'] = v_rpair
            else:
                record['rpairs'] = []
            if "PATH" in tokens.keys():
                foundEnd = False
                count = 0
                v_pathway = []
                while not foundEnd:
                    currentLine = r[tokens["PATH"] + count]
                    if currentLine[0:4] in ["PATH", "    "]:
                        v_pathway.append(currentLine[12:20])
                        count += 1
                    else:
                        foundEnd = True
                record['pathway'] = v_pathway
            else:
                record['pathway'] = []
            # Add record to list
            rclassData.append(record)
    print(len(rclassData), "reaction records created\n")
    return rclassData

def kegg_compounds(filename):
    # Read in compound file
    file = open(filename, "r")
    print("Reading Compound File")
    compound = file.readlines()

    # Find end of each record
    endList = []
    for line_i, line in enumerate(compound, 1):
        if re.match("///", line):
            endList.append(line_i)

    # Loop through each record extracting tokens
    print("Parsing Compound Data")
    count2 = 0
    compoundData = []
    for ref in range(0, len(endList)-1):
        count1 = count2
        count2 = endList[ref]
        c = compound[count1:count2]
        tokenList = ["ENTR", "PATH", "EXAC"]
        tokens = dict()
        for line_i, line in enumerate(c, 1):
            if line[0:4] in tokenList:
                tokens[line[0:4]] = line_i - 1
        # Parse tokens, extracting information into a dictionary
        if "ENTR" in tokens.keys():
            record = dict()
            record["entry"] = c[tokens["ENTR"]][12:18]
            if "EXAC" in tokens.keys():
                massText = re.search("[0-9.]+", c[tokens["EXAC"]])
                if (massText):
                    record['mass'] = float(massText.group(0))
                else:
                    record['mass'] = None
            else:
                record['mass'] = None
            if "PATH" in tokens.keys():
                foundEnd = False
                count = 0
                v_pathway = []
                while not foundEnd:
                    currentLine = c[tokens["PATH"] + count]
                    if currentLine[0:4] in ["PATH", "    "]:
                        v_pathway.append(currentLine[12:20])
                        count += 1
                    else:
                        foundEnd = True
                record['pathway'] = v_pathway
            else:
                record['pathway'] = []

            # Add record to list
            compoundData.append(record)
    print(len(compoundData), "compound records created\n")
    return compoundData


if __name__ == "__main__":
    main()
