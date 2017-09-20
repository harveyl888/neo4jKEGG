# neo4jKEGG
Python functions to parse KEGG reaction class and compound files and push the derived reactions into a neo4j database.  

The databases have been downloaded from the KEGG ftp site:
ftp://ftp.pathway.jp/kegg/ligand/compound.tar.gz
ftp://ftp.pathway.jp/kegg/ligand/enzyme.tar.gz
ftp://ftp.pathway.jp/kegg/ligand/rclass.tar.gz
ftp://ftp.pathway.jp/kegg/ligand/reaction.tar.gz

In addition, kgml files from the following are required
ftp://ftp.pathway.jp/kegg/xml/kgml/metabolic/ko.tar.gz
ftp://ftp.pathway.jp/kegg/pathway/pathway.list

For more information see [https://harveyl888.github.io/blog/KEGG-neo4j/](https://harveyl888.github.io/blog/KEGG-neo4j/)
