import requests
import pdb

ORGANISM = "eco"

pathways = requests.get('http://rest.kegg.jp/list/pathway/' + ORGANISM)

for line in pathways.content.decode('utf-8').split('\n'):
  pathwayid = line.split('\t')[0].replace('path:','')
  kgml = requests.get('http://rest.kegg.jp/get/'+pathwayid+"/kgml")
  my_content = kgml.content.decode('utf-8')
  f = open(pathwayid + '.xml','w')
  f.write(my_content)
  f.close

