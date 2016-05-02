import sys
import urllib2
import xml.etree.ElementTree as ET


url = 'http://www.rcsb.org/pdb/rest/search'

queryText = """

<?xml version="1.0" encoding="UTF-8"?>

<orgPdbQuery>

<version>B0907</version>

<queryType>org.pdb.query.simple.UniprotGeneNameQuery</queryType>

<description>UniProt Gene Name: </description>

<newsortfield> Residue%20Count%20Descending </newsortfield>

<query>"""+sys.argv[1]+"""</query>

</orgPdbQuery>


"""


print "query:\n", queryText

print "querying PDB...\n"

req = urllib2.Request(url, data=queryText)

f = urllib2.urlopen(req)

queryresult = f.read()


if queryresult:
    print queryresult
    print "Found number of PDB entries:", queryresult.count('\n')
else:
    print "Failed to retrieve results"

results = [i.split(":")[0] for i in queryresult.strip().split("\n")]
i=0
ids=[]
for result in results:
    pdbdesc = urllib2.urlopen("http://www.rcsb.org/pdb/rest/describePDB?structureId="+result)
    pdbtree = ET.parse(pdbdesc)
    pdbroot = pdbtree.getroot()
    descript = pdbroot.find("PDB")
    ids.append(descript.attrib['structureId'])
    print str(i) + ": " + descript.attrib['structureId'] + " - " + descript.attrib['title']
    i+=1
ind=int(input("Which PDB file?  (use index number on left side)\n"))
pid=ids[ind]
pdb = urllib2.urlopen("http://files.rcsb.org/download/"+pid+".pdb")
with open(pid+'.pdb','wb') as output:
  output.write(pdb.read())
protstruct = urllib2.urlopen("http://www.rcsb.org/pdb/rest/describeMol?structureId="+pid)
prottree = ET.parse(protstruct)
protroot = prottree.getroot()
"""
entity type="protein" use length as actual protein length if not a fragment.  look for fragment desc and residues if it is a fragment, regex e.g. "residues 3-35"
"""
i=0
chains=[]
length=[]
for i,(entry,desc) in enumerate(zip(protroot.iter('polymer'), protroot.iter('polymerDescription'))):
# if there are multiple chain ids for the same structure, it will just take the first one, which is what we want
    chain = entry.find('chain').attrib['id']
    descript = desc.attrib['description']
    if entry.attrib['type'] == 'protein' and entry.find('fragment') != None and 'residues' in entry.find('fragment').attrib['desc'].lower():
        s=entry.find('fragment').attrib['desc'].lower()
        if len(s.lower().strip(")").split(",")) == 1:
            length.append(s.strip(")").replace(" ","").split("residues")[-1])
        else:
            t = s.strip(")").split(",")
            for j in range(0,len(t)):
                length.append(t[j].replace(" ","").split("residues")[-1])
    elif entry.attrib['type'] == 'protein' and entry.find('fragment') == None:
        length.append("1-"+entry.attrib['length'])
    print str(i) + ": " + chain + " - " + descript + ", " + ", ".join(length)
    chains.append([chain,length])
    length = []
ind=int(input("Which chain to use?  (use index number on left side)\n")) # the printout will not print the chains for DNA or non-proteins
chain=chains[ind]
print chain
ans=input("Would you like to see the Pfam structure? (y or n)\n").lower()
if ans == 'yes' or 'y':
    pfamstructure = urllib2.urlopen("http://www.rcsb.org/pdb/rest/hmmer?structureId="+pid)

