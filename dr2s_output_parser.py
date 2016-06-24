"""
This script extracts sequences and identified alleles from d2rs output folders
Expected Folder structure : <rootFolder>/<barcode>/checked/
Expected File Names : [A|B]*.json (from which the allele name/s is extracted, not json parsing as we only need 1 line from the json
		      [A|B]*.fa (from which the sequence is extracted 

The <rootFolder>s are supplied as cmd line arguments to the script
If there is no "checked" folder found, the script will output "No Data Found" for that barcode

Output is a csv file written to the same location as where this script is 
with the structure : "RootFolder\tBC\tHapA Allele\tHapA Seq\tHapB Allele\tHapB Seq" where BC is the barcode
"""


from sys import argv
from glob import glob
from os import path
from Bio import SeqIO

mainPath = "/home/vineeth/bioinf/DPB1-output"
hapNames = ["A","B"]

def getSeqFromFasta(fastaFile):
	parser = SeqIO.parse(open(fastaFile), "fasta")
	for record in parser: return str(record.seq)

def getTypeLoaderAllele(jsonFile):
	jsonText = open(jsonFile).read()
	return jsonText.split("\"closest_allele\":")[1].split(",")[0].replace("\"","") 

sequences = {}
for rootFolder in argv[1:]:
	sequences[rootFolder] = {}
	lbcPath = path.join(mainPath,rootFolder,"*")
	lbcFolders = glob(lbcPath)

	for lbcFolder in lbcFolders:
		lbcId = lbcFolder.split("/")[-1]
		sequences[rootFolder][lbcId] = {}
		
		# dr2s idiosyncracy : for barcodes starting with ID, dr2s seems to output a single subfolder with name form - HLA-DPB1.ref.alt
		if lbcId.find("ID") != -1: lbcFolder = glob(path.join(lbcFolder,"*"))[0] 
	
		for hapName in hapNames:
		
			try:	
				hapJson = glob(path.join(lbcFolder,"checked","%s*.json" % hapName))[0]
				hapFasta = glob(path.join(lbcFolder, "checked", "%s*.fa" % hapName))[0]
				sequences[rootFolder][lbcId][hapName] = {"allele" : getTypeLoaderAllele(hapJson), \
									 "seq": getSeqFromFasta(hapFasta)}
			except IndexError:
				continue 

		
seqFile = open("dr2s_output.txt","w")
seqFile.write("RootFolder\tBC\tHapA Allele\tHapA Seq\tHapB Allele\tHapB Seq\n")

for rootFolder in sequences.keys():
	lbcIds = sequences[rootFolder].keys()
	for lbcId in lbcIds:
		seqFile.write("%s\t%s" % (rootFolder,lbcId))
		for hapName in hapNames:
			try:
				seqFile.write("\t%s" % sequences[rootFolder][lbcId][hapName]["allele"])
				seqFile.write("\t%s" % sequences[rootFolder][lbcId][hapName]["seq"])
			except KeyError:
				seqFile.write("\tNo Data Found")
				seqFile.write("\tNo Data Found")
		seqFile.write("\n")
seqFile.close()
		
		
			
	
