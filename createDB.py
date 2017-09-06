#!/usr/bin/python

import re


def main():
    rclass = kegg_reactions("C:/Databases/KEGG/rclass/rclass")
    print(len(rclass), "reaction records created")
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
                        v_pathway.append(currentLine[12:20].rstrip())
                        count += 1
                    else:
                        foundEnd = True
                record['pathway'] = v_pathway
            else:
                record['pathway'] = []
        # Add record to list
        rclassData.append(record)
    return rclassData

if __name__ == "__main__":
    main()
