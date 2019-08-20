import xml.etree.ElementTree as ET

file1 = '../../projects/ralstonia/pre-pool-files/ReutrophaH16Nottingham_transposon_coord_to_feature_lookup_table_gene.xml'
file2 = '../../projects/ralstonia/pre-pool-files/ReutrophaH16Nottingham_transposon_coord_to_feature_lookup_table_gene2.xml'

xml1 = open(file1, 'r')
xml2 = open(file2, 'r')

tree1 = ET.parse(xml1)
tree2 = ET.parse(xml2)

root1 = tree1.getroot()
root2 = tree2.getroot()

coords1 = root1.findall('coord')
coords2 = root2.findall('coord')

tot = len(coords1)
off1 = []
off2 = []
coord_lsts = {}
for coord in coords1:
    co = int(coord.attrib['coord'])
    loci = coord.findall('locus')
    coord_lsts[co] = set()
    for loc in loci:
        tag = loc.attrib['locus']
        if tag in coord_lsts[co]:
            print('wtf')
        coord_lsts[co].add(tag)

for coord in coords2:
    co = int(coord.attrib['coord'])
    loci = coord.findall('locus')
    for loc in loci:
        tag = loc.attrib['locus']
        if tag in coord_lsts[co]:
            coord_lsts[co].remove(tag)
        else:
            off2.append((tag, co))
    for tag in coord_lsts[co]:
        off1.append((tag,co))

print(off1)
print(off2)
