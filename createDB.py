#!/usr/bin/python

import re
import glob
import time
import xml.etree.ElementTree as ET
from py2neo import authenticate, Graph
from pandas import DataFrame


def main():
    # Connect to neo4j
    authenticate("localhost:7474", "neo4j", "neo4jpw")
    graph = Graph("localhost:7474/db/data/")

    pathway_list = read_pathway_list("C:/Databases/KEGG/pathway.list")
    main_pathway = "Metabolism"
    use_pathways = [i[2] for i in pathway_list if i[0] in main_pathway]
    use_pathways.remove('01100')

    metabolic_reactions = read_kegg_xml("C:/Databases/KEGG/kgml/ko", use_pathways=use_pathways)

    # Read in files and create a new database
    reactions = kegg_reactions("C:/Databases/KEGG/reaction/reaction")
    enzymes = kegg_enzymes("C:/Databases/KEGG/enzyme/enzyme")
#    rclass = kegg_rclass("C:/Databases/KEGG/rclass/rclass")
    compounds = kegg_compounds("C:/Databases/KEGG/compound/compound")
#    triples = find_triples(rclass)

    create_db(metabolic_reactions, graph, reactions, enzymes, compounds)

#    create_db(triples, rclass, compounds, graph)

    # Test database
    test_database(graph)

    return 0


def read_pathway_list(filename):
    # Read in pathway list and return a hierarchy
    pathways = list()
    level_1_heading = ""
    level_2_heading = ""
    with open(filename) as f:
        for line in f:
            if line[0:1] == "#":
                if line[0:2] == "##":
                    level_2_heading = line[2:].strip()
                else:
                    level_1_heading = line[1:].strip()
            else:
                pathways.append([level_1_heading, level_2_heading, line[0:5], line[6:].strip()])
    f.close()
    return pathways


def read_kegg_xml(folder, use_pathways=None, ignore_pathways=None):
    # Read in and parse xml files
    xml_files = glob.glob(folder + '/*.xml')
    if type(use_pathways) == list:
        r = re.compile("|".join(use_pathways))
        l_keepme = map(lambda x: len(r.findall(x)) > 0, xml_files)
        keep_me = [i for i, x in enumerate(list(l_keepme)) if x]
        xml_files = [xml_files[i] for i in keep_me]
    elif type(ignore_pathways) == list:
        r = re.compile("|".join(ignore_pathways))
        l_keepme = map(lambda x: len(r.findall(x)) > 0, xml_files)
        keep_me = [i for i, x in enumerate(list(l_keepme)) if not x]
        xml_files = [xml_files[i] for i in keep_me]

    reactions = list()
    for filename in xml_files:
        tree = ET.parse(filename).getroot()
        pathway_number = tree.attrib['number']
        pathway_name = tree.attrib['title']
        reactions_xml = tree.findall("reaction")
        for r in reactions_xml:
            reaction = dict()
            reaction['id'] = r.attrib['id']
            reaction['name'] = [x[3:] for x in r.attrib['name'].split(" ")]
            reaction['type'] = r.attrib['type']
            reaction['pathway_id'] = pathway_number
            reaction['pathway_name'] = pathway_name
            for child in r:
                # Currently takes just a single substrate or product when multiple are present
                cpd = dict()
                cpd['id'] = child.attrib['id']
                if child.attrib['name'][0:2] == "g1":
                    cpd['name'] = child.attrib['name'][3:9]  # glycan
                else:
                    cpd['name'] = child.attrib['name'][4:10]  # limit attribute to 6 characters (only first cpd entry)
                reaction[child.tag] = cpd
            reactions.append(reaction)
    return reactions


def kegg_reactions(filename):
    # Read in reactions file
    file = open(filename, "r")
    print("Reading Reaction File")
    reactions = file.readlines()
    file.close()

    # Find end of each record
    end_list = []
    for line_i, line in enumerate(reactions, 1):
        if re.match("///", line):
            end_list.append(line_i)

    # Loop through each record extracting tokens
    print("Parsing Reaction Class Data")
    count2 = 0
    reaction_data = dict()
    for ref in range(0, len(end_list)-1):
        count1 = count2
        count2 = end_list[ref]
        r = reactions[count1:count2]
        token_list = ["ENTR", "NAME", "DEFI", "EQUA", "RCLA", "ENZY"]
        tokens = dict()

        for line_i, line in enumerate(r, 1):
            if line[0:4] in token_list:
                tokens[line[0:4]] = line_i - 1
        # Parse tokens, extracting information into a dictionary
        if "ENTR" in tokens.keys() and "RCLA" in tokens.keys():
            record = dict()
            # Entry token
            record["entry"] = r[tokens["ENTR"]][12:18]
            # RClass token
            found_end = False
            count = 0
            v_rclass = []
            while not found_end:
                current_line = r[tokens["RCLA"] + count]
                if current_line[0:4] in ["RCLA", "    "]:
                    v_rclass.append(re.findall("([RC0-9]+)", current_line[12:]))
                    count += 1
                else:
                    found_end = True
                record['rclass'] = v_rclass
            # Name token
            if "NAME" in tokens.keys():
                record['name'] = r[tokens["NAME"]][12:].strip()
            # Definition token
            if "DEFI" in tokens.keys():
                record['definition'] = r[tokens["DEFI"]][12:].strip()
            # Equation token
            if "EQUA" in tokens.keys():
                record['equation'] = r[tokens["EQUA"]][12:].strip()
            # Enzyme token
            if "ENZY" in tokens.keys():
                record['enzyme'] = r[tokens["ENZY"]][12:].strip()
            # Add record to list
            reaction_data[record['entry']] = record
    print(len(reaction_data), "reaction records created\n")
    return reaction_data


def kegg_enzymes(filename):
    # Read in enzyme file
    file = open(filename, "r")
    print("Reading Enzyme File")
    enzymes = file.readlines()
    file.close()

    # Find end of each record
    end_list = []
    for line_i, line in enumerate(enzymes, 1):
        if re.match("///", line):
            end_list.append(line_i)

    # Loop through each record extracting tokens
    print("Parsing Enzyme Data")
    count2 = 0
    enzyme_data = dict()
    for ref in range(0, len(end_list)-1):
        count1 = count2
        count2 = end_list[ref]
        e = enzymes[count1:count2]
        token_list = ["ENTR", "NAME", "PATH"]
        tokens = dict()

        for line_i, line in enumerate(e, 1):
            if line[0:4] in token_list:
                tokens[line[0:4]] = line_i - 1
        # Parse tokens, extracting information into a dictionary
        if "ENTR" in tokens.keys():
            record = dict()
            entry_id = re.search("EC \d+.\d+.\d+.\d+", e[tokens["ENTR"]])
            if entry_id:
                record["entry"] = entry_id.group()
            else:
                break
            # Name token
            if "NAME" in tokens.keys():
                name_text = e[tokens["NAME"]].rstrip()
                if name_text[len(name_text) - 1] == ";":
                    record['name'] = name_text[12:len(name_text) - 1]
                else:
                    record['name'] = name_text[12:]
            # Pathway token
            if "PATH" in tokens.keys():
                found_end = False
                count = 0
                v_pathway = []
                while not found_end:
                    current_line = e[tokens["PATH"] + count]
                    if current_line[0:4] in ["PATH", "    "]:
                        v_pathway.append(current_line[12:19])
                        count += 1
                    else:
                        found_end = True
                record['pathway'] = v_pathway
            # Add record to list
            enzyme_data[record['entry']] = record
    print(len(enzyme_data), "enzyme records created\n")
    return enzyme_data


def kegg_rclass(filename):
    # Read in rclass file
    file = open(filename, "r")
    print("Reading Reaction Class File")
    rclass = file.readlines()
    file.close()

    # Find end of each record
    end_list = []
    for line_i, line in enumerate(rclass, 1):
        if re.match("///", line):
            end_list.append(line_i)

    # Loop through each record extracting tokens
    print("Parsing Reaction Class Data")
    count2 = 0
    rclass_data = dict()
    for ref in range(0, len(end_list)-1):
        count1 = count2
        count2 = end_list[ref]
        r = rclass[count1:count2]
        token_list = ["ENTR", "DEFI", "RPAI", "PATH"]
        tokens = dict()
        for line_i, line in enumerate(r, 1):
            if line[0:4] in token_list:
                tokens[line[0:4]] = line_i - 1
        # Parse tokens, extracting information into a dictionary
        if "ENTR" in tokens.keys():
            record = dict()
            record["entry"] = r[tokens["ENTR"]][12:19]
            if "DEFI" in tokens.keys():
                found_end = False
                count = 0
                v_definition = []
                while not found_end:
                    current_line = r[tokens["DEFI"] + count]
                    if current_line[0:4] in ["DEFI", "    "]:
                        v_definition.append(current_line[12:].rstrip())
                        count += 1
                    else:
                        found_end = True
                record['definition'] = v_definition
            if "RPAI" in tokens.keys():
                found_end = False
                count = 0
                v_rpair = []
                while not found_end:
                    current_line = r[tokens["RPAI"] + count]
                    if current_line[0:4] in ["RPAI", "    "]:
                        v_rpair.extend(re.findall("([C0-9_]+)", current_line[12:]))
                        count += 1
                    else:
                        found_end = True
                record['rpairs'] = v_rpair
            if "PATH" in tokens.keys():
                found_end = False
                count = 0
                v_pathway = []
                while not found_end:
                    current_line = r[tokens["PATH"] + count]
                    if current_line[0:4] in ["PATH", "    "]:
                        v_pathway.append(current_line[12:20])
                        count += 1
                    else:
                        found_end = True
                record['pathway'] = v_pathway
            # Add record to list
            rclass_data[record['entry']] = record
    print(len(rclass_data), "reaction class records created\n")
    return rclass_data


def kegg_compounds(filename):
    # Read in compound file
    file = open(filename, "r")
    print("Reading Compound File")
    compound = file.readlines()
    file.close()

    # Find end of each record
    end_list = []
    for line_i, line in enumerate(compound, 1):
        if re.match("///", line):
            end_list.append(line_i)

    # Loop through each record extracting tokens
    print("Parsing Compound Data")
    count2 = 0
    compound_data = dict()
    for ref in range(0, len(end_list)-1):
        count1 = count2
        count2 = end_list[ref]
        c = compound[count1:count2]
        token_list = ["ENTR", "NAME", "FORM", "EXAC", "PATH"]
        tokens = dict()
        for line_i, line in enumerate(c, 1):
            if line[0:4] in token_list:
                tokens[line[0:4]] = line_i - 1
        # Parse tokens, extracting information into a dictionary
        if "ENTR" in tokens.keys():
            record = dict()
            record["entry"] = c[tokens["ENTR"]][12:18]
            # Name token
            if "NAME" in tokens.keys():
                name_text = c[tokens["NAME"]].rstrip()
                if name_text[len(name_text) - 1] == ";":
                    record['name'] = name_text[12:len(name_text) - 1]
                else:
                    record['name'] = name_text[12:]
            # Formula token
            if "FORM" in tokens.keys():
                formula_text = c[tokens["FORM"]].rstrip()
                record['formula'] = formula_text[12:]
            # Exact mass token
            if "EXAC" in tokens.keys():
                mass_text = re.search("[0-9.]+", c[tokens["EXAC"]])
                if mass_text:
                    record['mass'] = float(mass_text.group(0))
            # Pathways token
            if "PATH" in tokens.keys():
                found_end = False
                count = 0
                v_pathway = []
                while not found_end:
                    current_line = c[tokens["PATH"] + count]
                    if current_line[0:4] in ["PATH", "    "]:
                        v_pathway.append(current_line[12:20])
                        count += 1
                    else:
                        found_end = True
                record['pathway'] = v_pathway

            # Add record to list
            compound_data[record['entry']] = record
    print(len(compound_data), "compound records created\n")
    return compound_data


def find_triples(rclass):
    triples = []
    for r in rclass:
        for p in r['rpairs']:
            rpair_triple = p.split("_")
            rpair_triple.append(r["entry"])
            triples.append(rpair_triple)
    return triples


def flatten_dict(d):
    out = list()
    for k, v in d.items():
        if type(v) is str:
            out.append("{k}: \"{v}\"".format(k=k, v=v))
        else:
            out.append("{k}: {v}".format(k=k, v=v))
    return ", ".join(out)


# Create and populate a neo4j database
def create_db(metabolic_reactions, graph, reactions=None, enzymes=None, compounds=None):
    # clear old data
    graph.delete_all()
    # begin Cypher transaction
    tx = graph.begin()
    # iterate over each member of metabolic_reactions and add to database
    # include optional data as properties if available
    print("Processing", len(metabolic_reactions), "reactions")
    start_time = time.time()
    for index, m in enumerate(metabolic_reactions):
        # Compound 1
        if m["substrate"]["name"] in compounds.keys():
            cpd1 = "MERGE (c1:Compound {{{id}}})".format(id=flatten_dict(compounds[m["substrate"]["name"]]))
        else:
            cpd1 = "MERGE (c1:Compound {{entry: \"{id}\"}})".format(id=m["substrate"]["name"])
        # Compound 2
        if m["product"]["name"] in compounds.keys():
            cpd2 = "MERGE (c2:Compound {{{id}}})".format(id=flatten_dict(compounds[m["product"]["name"]]))
        else:
            cpd2 = "MERGE (c2:Compound {{entry: \"{id}\"}})".format(id=m["product"]["name"])
        # Mass change
        delta_mass = None
        if all([m["substrate"]["name"] in compounds.keys(), m["product"]["name"] in compounds.keys()]):
            if all(['mass' in compounds[m["substrate"]["name"]], 'mass' in compounds[m["product"]["name"]]]):
                delta_mass = round(compounds[m["product"]["name"]]['mass'] -
                                   compounds[m["substrate"]["name"]]['mass'], 4)
        # Reaction
        if m["type"] == "reversible":
            relation = "-"
        else:
            relation = "->"
        react_data = dict()
        react_data["type"] = m["type"]

        # loop through reactions
        for react_name in m["name"]:
            react_data["entry"] = react_name
            if delta_mass is not None:
                react_data["delta_mass"] = delta_mass
                react_data["abs_delta_mass"] = abs(delta_mass)
            if react_name in reactions.keys():
                if "definition" in reactions[react_name]:
                    react_data["definition"] = reactions[react_name]["definition"]
                if "equation" in reactions[react_name]:
                    react_data["equation"] = reactions[react_name]["equation"]
                if "name" in reactions[react_name]:
                    react_data["reaction_name"] = reactions[react_name]["name"]
                if "rclass" in reactions[react_name]:
                    react_data["reaction_class"] = [x[0] for x in reactions[react_name]["rclass"]]
                if "enzyme" in reactions[react_name]:
                    enz = reactions[react_name]["enzyme"]
                    react_data["enzyme"] = enz
                    if enz in enzymes.keys():
                        if "name" in enzymes[enz]:
                            react_data["enzyme_name"] = enzymes[enz]["name"]
                        if "pathway" in enzymes[enz]:
                            react_data["pathway"] = enzymes[enz]["pathway"]
            # generate MERGE for cypher transaction
            react = "[:REACTION {{{r}}}]".format(r=flatten_dict(react_data))
            merge_text = "{cpd1} {cpd2} MERGE (c1)-{react}{relation}(c2)"\
                .format(cpd1=cpd1, react=react, relation=relation, cpd2=cpd2)

    #        print("{i}    {m}".format(i=index, m=react))

            tx.append(merge_text)
        # Create simple connection between pairs of compounds
        if delta_mass is not None:
            merge_text = "{cp1} {cp2} MERGE (c1)-[:CONNECTION {{abs_delta_mass: \"{m1}\"}}]-(c2)"\
                .format(cp1=cpd1, cp2=cpd2, m1=abs(delta_mass))
        else:
            merge_text = "{cp1} {cp2} MERGE (c1)-[:CONNECTION]-(c2)" .format(cp1=cpd1, cp2=cpd2)
        tx.append(merge_text)

    end_time = time.time()
    print("Time to create transaction =", int(end_time - start_time), "seconds")
    start_time = time.time()
    tx.commit()
    end_time = time.time()
    print("Time to upload transaction =", int(end_time - start_time), "seconds")
    return


# Create and populate a neo4j database
def create_db_from_triples(triples, rclass, compounds, graph):
    # clear old data
    graph.delete_all()
    # begin Cypher transaction
    tx = graph.begin()
    # iterate over each triple and add to database
    print("Processing", len(triples), "relationships")
    start_time = time.time()
    for index, t in enumerate(triples):
        # Lookup each value in triple.  If does not exist then return a record with just the Entry id
        # Compound 1
        try:
            c1 = next(item for item in compounds if item['entry'] == t[0])
        except StopIteration:
            c1 = {'entry': t[0], 'name': None, 'formula': None, 'mass': None, 'pathway': []}
        # Compound 2
        try:
            c2 = next(item for item in compounds if item['entry'] == t[1])
        except StopIteration:
            c2 = {'entry': t[1], 'name': None, 'formula': None, 'mass': None, 'pathway': []}
        c = [c1, c2]
        # Reaction
        try:
            r = next((item for item in rclass if item['entry'] == t[2]))
        except StopIteration:
            r = {'entry': t[2], 'definition': [], 'rpairs': [], 'pathway': []}
        # Delta mass
        try:
            delta_mass = round(c2['mass'] - c1['mass'], 4)
        except TypeError:
            delta_mass = None

        merge_text = ""
        for idx, compound in enumerate(c):
            merge_text += "MERGE (c{i1}:Compound{{id:\"{id1}\"".format(i1=idx+1, id1=compound['entry'])
            if compound['name']:
                merge_text += ", name:\"{name}\"".format(name=compound['name'])
            if compound['formula']:
                merge_text += ", formula:\"{formula}\"".format(formula=compound['formula'])
            if compound['mass']:
                merge_text += ", mass:{mass}".format(mass=compound['mass'])
            if compound['pathway']:
                merge_text += ", pathways:{pathway}".format(pathway=compound['pathway'])
            merge_text += "}) "
        merge_text += "CREATE (c1)-[:REACTION{{reaction:\"{reaction}\"".format(reaction=t[2])
        if r['definition']:
            merge_text += ", definition:{definition}".format(definition=r['definition'])
        if delta_mass:
            merge_text += ", delta_mass:{delta_mass}, abs_delta_mass:{abs_delta_mass}".\
                format(delta_mass=delta_mass, abs_delta_mass=abs(delta_mass))
        merge_text += "}]->(c2)"
        tx.append(merge_text)
    end_time = time.time()
    print("Time to create transaction =", int(end_time - start_time), "seconds")
    start_time = time.time()
    tx.commit()
    end_time = time.time()
    print("Time to upload transaction =", int(end_time - start_time), "seconds")
    return


# A few tests to make sure we can execute cypher queries
def test_database(graph):

    tests = list()
    tests.append({"title": "Return first 10 compound IDs",
                  "query": "MATCH (c:Compound) RETURN c.id LIMIT 10", "dataframe": False})
    tests.append({"title": "Return first 10 relationships",
                  "query": "MATCH p=()-[r:REACTION]->() RETURN p LIMIT 10", "dataframe": False})
    tests.append({"title": "Using dataframes",
                  "query": "MATCH p=(c1:Compound)-[r:REACTION]->(c2:Compound) RETURN c1.id, r.reaction, c2.id LIMIT 25",
                  "dataframe": True})
    ppm_limit = 100
    mass_search = 18.0105
    tests.append({"title": "Search for mass match - dehydration reactions",
                  "query": "MATCH p=(c1:Compound)-[r:REACTION]->(c2:Compound) WHERE "
                           "1E6*abs({mass_search} - r.abs_delta_mass)/{mass_search} < {ppm_limit} RETURN "
                           "c1.id, r.reaction, c2.id LIMIT 25".format(mass_search=mass_search, ppm_limit=ppm_limit),
                  "dataframe": True})

    for t in tests:
        print("\n" + t['title'])
        print(t['query'])
        if t['dataframe']:
            out = DataFrame(graph.data(t['query']))
        else:
            out = graph.data(t['query'])
        print(out)
    return


if __name__ == "__main__":
    main()
