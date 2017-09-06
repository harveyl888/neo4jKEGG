#!/usr/bin/python

import re


def main():
    rclass = kegg_reactions("C:/Databases/KEGG/rclass/rclass")
    compounds = kegg_compounds("C:/Databases/KEGG/compound/compound")
    triples = find_triples(rclass)
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


def find_triples(rclass):
    triples = []
    for r in rclass:
        for p in r['rpairs']:
            rpair_triple = p.split("_")
            rpair_triple.append(r["entry"])
            triples.append(rpair_triple)
    return triples


# Create and populate a neo4j database
def createDB(triples, rclass, compounds):
    from py2neo import authenticate, Graph
    authenticate("localhost:7474", "neo4j", "neo4jpw")
    graph = Graph("localhost:7474/db/data/")
    graph.delete_all()

    for i in triples[0:300]:
        c = [next((item for item in compounds if item['entry'] == i[0])), next((item for item in compounds if item['entry'] == i[1]))]
        r = next((item for item in rclass if item['entry'] == i[2]))

        print(i)

        mergeText = ""
        for idx, compound in enumerate(c):
            mergeText += "MERGE (c{i1}:Compound{{id:\"{id1}\"".format(i1=idx+1, id1=compound['entry'])
            if compound['mass']:
                mergeText += ", mass:{mass}".format(mass=compound['mass'])
            if compound['pathway']:
                mergeText += ", pathways:{pathway}".format(pathway=compound['pathway'])
            mergeText += "}) "
        mergeText += "CREATE (c1)-[:REACTION{{reaction:\"{reaction}\"".format(reaction=i[2])
        if r['definition']:
            mergeText += ", definition:{definition}".format(definition=r['definition'])
        mergeText += "}]->(c2)"

        print(mergeText)
        graph.run(mergeText)

#        graph.run("MERGE (c1:Compound{id:{id1}, mass:{m1}, pathways:{p1}}) "
#                  "MERGE (c2:Compound{id:{id2}, mass:{m2}, pathways:{p2}}) "
#                  "CREATE (c1)-[:REACTION{reaction:{reaction}, definition:{definition}}]->(c2)",
#                  id1=i[0], m1=c1['mass'], p1=c1['pathway'], id2=i[1], m2=c2['mass'], p2=c2['pathway'],
#                  reaction=i[2], definition=r['definition'])



if __name__ == "__main__":
    main()


