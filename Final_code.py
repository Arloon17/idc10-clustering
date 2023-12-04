import sys
# need to first import gzip when using .gz compressed files
import gzip
import matplotlib.pyplot as plt
from wordcloud import WordCloud, STOPWORDS
import networkx as nx
from displaygenelocations import *
#-------------------------------------------------------
# Preliminary analysis: ordering UK biobank medical information
#-------------------------------------------------------
def sort_medinfo(medinfofile):
	# List containing UK biobank medical information
	medinfo = []
	# Open the UK Biobank file
	with open(medinfofile,"r") as f:
		# GO over each file of the UK Biobank file
		for i, line in enumerate(f):
			new_line = line.rstrip()
			split_line = new_line.split("\t")
			medinfo.append(split_line)
	# Sort the ICD10 diagnoses with respect ot the counts
	medinfo = sorted(medinfo, key=lambda x: -int(x[1]))
	# Open a file in which the sorted list is saved
	with open("Sorted_" + medinfofile, "w") as save:
		for i in range(len(medinfo)):
			save.write("\t".join(medinfo[i][:3]) + "\n")
			save.flush()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ("sortmedinfo" in sys.argv):
	sort_medinfo( "Alzheimer_dict_nz_for_icd10.corrected.modular.txt")
#-------------------------------------------------------
# First step of the project: record the list of disease names contained
# in the Percha and Altman gene-disease network. The list of disease names are saved
# in a file named "PA_diseases.txt".
#-------------------------------------------------------
# INPUT: a string [the_string]
# OUTPUT: a lower-case version of the input string in which dashes
# underscores and apostrophes where  replaced with spaces
def normalization_Text(the_string):
	# Change upper-case letters to lower-case ones
	lower_string = the_string.lower()
	# Replace dashes with spaces
	lower_string = lower_string.replace("-", " ")
	# Replace underscores with spaces
	lower_string = lower_string.replace("_", " ")
	# Replace apostrophes with spaces
	lower_string = lower_string.replace("'", " ")
	# Return the resulting string
	return lower_string
#-------------------------------------------------------
# INPUTS: a list [the_list] of items and an item [element] to be added to the list
# OUTPUT: a boolean value (true or false) indicating whether the item [element]
# needed to be added to the list or not
def nonrepeating_Append(the_list, element):
	# Variable indicating if [element] is already contained in [the_list]
	found = False
	# GO over all the items of [the_list]
	for ele in the_list:
		# Check whether [the_list] contains [element]
		if ele == element:
			# Indicate that [element] is already contained in [the_list]
			found = True
			#If found, the loop is stopped
			break
	# Check whether [element] was found or not during the loop
	if found == False:
		# If not found, the item [element] is appended to the list
		the_list.append(element)
	# Return whether the item [element] needed to be added or not
	return found
#-------------------------------------------------------
# INPUT: a string [filename] referring to the name of a gene-disease network from the
# Percha and Altman database
# OUTPUT: the list of diseases contained in the network.
def record_Diseases(filename):
	# The list of diseases
	diseases = []
	# Open the file containing the gene-disease network.
	# Gzip is only used when you open the file
	with gzip.open(filename, 'r') as f:
		# Go through each line of the file
		for i, line in enumerate(f):
			#.decode() is used to decode binary to ASCII
			new_line = line.decode().rstrip()
			split_line = new_line.split("\t")
			# Add disease names without repetitions
			nonrepeating_Append(diseases, normalization_Text(split_line[4]))
			# Message on the prompt to check the state of the computation
			print("computing...", split_line[2][0], split_line[4])
	# Return the list of diseases contained in the network
	return diseases
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ("recdis" in sys.argv):
	# The  first step of the project consisted in collecting the list of disease names contained
	# in the gene-disease network of Percha and Altman. To do so, the function record_Disease was applied
	# of the file "part-ii-dependency-paths-gene-disease-sorted-with-themes.txt.gz" to obtain the list
	# of disease name contained in the file. Then each disease name is saved on a line of a file
	# named "PA_diseases.txt".
	output = record_Diseases("part-ii-dependency-paths-gene-disease-sorted-with-themes.txt.gz")

	with open("PA_diseases.txt",'w') as f:
		for i in range(len(output)):
			f.write(str(i)+"\t"+output[i]+"\n")
			f.flush()
#-------------------------------------------------------
#-------------------------------------------------------
#Second step of the project: Count how many times a word occurs
# in the 3rd column of a table containing ICD10 descriptions
#-------------------------------------------------------
# INPUTS: a string [filename] referring to the name of a file containing a single word at each
# of its lines.
# OUTPUTS: The list of words contained in the input file
# DESCRIPTION: looping over the elements of the list and checking whether there
#are elements of this list that are substrings of the string
def get_Exclusions(filename):
	excluded_words = []
	with open(filename,"r") as f:
		for i, line in enumerate(f):
			line = line.rstrip()
			s = line.split("\t")
			nonrepeating_Append(excluded_words,s[0])
	return excluded_words
#-------------------------------------------------------
# INPUTS: a list [the_list] of items and an item [element] to be added to the list
# OUTPUT: a boolean value (true or false) indicating whether the item [element]
# was found in the list before being added
def counting_Append(the_list, element):
	# Variable indicating if [element] is already contained in [the_list]
	found = False
	# GO over all the items of [the_list]
	for i in range(len(the_list)):
		# Check whether [the_list] contains [element]
		if the_list[i][0] == element:
			# Increment the count associated with the item [element]
			the_list[i][1] += 1
			# Indicate that [element] is already contained in [the_list]
			found = True
			# If found, the loop is stopped
			break
	# Check whether [element] was found or not during the loop
	if found == False:
		# If not found, the item [element] is appended to the list
		the_list.append([element,1])
	# Boolean value indicating whether the item [element] needed to be added or not
	return found
# -------------------------------------------------------
#INPUTS: the function takes two strings:
# - a string [medinfofile] referring to the name of a file containing a table
# of 3 columns such that its 1st column contains ICD10 codes, its 2nd column contains numerical values
# and its 3rd column contains descriptions of ICD10 codes.
# - a string [exclusionfile] referring to a file whose lines each contain a single word
# OUTPUT: a list containing lists of the form [w,c] where
# - w is a word extracted from one of the ICD10 descriptions contained in the file [medinfofile]
# and not present in the text of [exclusionfile] and
# - c is an integer counting the number of occurrences of that word in [medinfofile]
def record_Words(medinfofile,exclusionfile):
	# The list of words output by the function
	wordlist = []
	# The list of words to not be included in [wordlist]
	exclusions = get_Exclusions(exclusionfile)
	# Open the file [medinfofile] which contains the ICD10 code descriptions
	with open(medinfofile, 'r') as g:
		# Go over each line of the file
		for i, line in enumerate(g):
			new_line = line.rstrip()
			split_line = new_line.split("\t")
			# separate each word of the english description of an ICD10
			# code in the list [UKB_description]
			UKB_description = split_line[2].split(" ")
			# Go over every word of the english description
			for word in UKB_description[1:]:
				word = word.lower()
				# If the lower-case version of the word is not in the excluded word
				# then add that word to the list [wordlist]. Next to each word is the count for
				# how many times this word was found in the text of the 3rd column of [medinfofile]
				if not (word in exclusions):
					counting_Append(wordlist, word)
			#print(UKB_description)
	# Order the list with respect to their associated counts and their lengths as words
	wordlist = sorted(wordlist, key=lambda x: (-int(x[1]), len(x[0])))
	# Return the list of words
	return wordlist
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ("recword" in sys.argv):

	wordlist = record_Words("Alzheimer_dict_nz_for_icd10.corrected.modular.txt","excluded_words.txt")

	with open("wordlist.txt", "w") as savefile:
		for i in range(len(wordlist)):
			savefile.write("\t".join(map(str, wordlist[i])) + "\n")
			savefile.flush()
#-------------------------------------------------------
#-------------------------------------------------------
# Third step of the project: Translate the UK Biobank ICD10 diagnoses into a file that associates
# ICD10 code to lists of related diseases.
#-------------------------------------------------------
# INPUT: a string [diseasefile] referring to the name of a file whose lines are numbered and contain
# the name of a disease. Specifically, in this project, this file is the file named
# "PA_diseases.txt" created in the first step of the project
# OUTPUT: a list containing all the disease names of the file
def list_Diseases(diseasefile):
	# List containing the disease names
	disease_list = []
	# Open the file [diseasefile], which contains disease names
	with open(diseasefile, 'r') as f:
		# Go over each line of the file
		for i, line in enumerate(f):
			new_line = line.rstrip()
			split_line = new_line.split("\t")
			# each line is numbered, so the second column of the file contains the disease name,
			# which is added to the list disease_list
			disease_list.append(split_line[1])
			# Message printed on the prompt to verify the state of the computation
			print(i,split_line[1])
	# Return the list of disease names
	return disease_list
#--------------------------------------------------------
# INPUT: a string [word] containing a word from an ICD10 description, a list of strings [list_of_diseases]
# containing disease names and a list of string [exclusions] that contains words to not consider
# OUTPUT: a list of strings, which all belong to the list [list_of_diseases]
def match_String(word, list_of_diseases,exclusions):
	# List containing disease names from [list_of_diseases] that are longer than 3 letters,
	# that are not in the list [exclusions] of excluded words and that
	# - either contain the string [word]
	# - or is contained in the string [word]
	matches = []
	# Only consider matches with words that are not in exclusions
	if not(word in exclusions):
		# Turn the item [word] into a lower-case string where "_", "-" and "'" are
		# replaced with spaces
		normalized_word = normalization_Text(word)
		# Go over every disease in [list_of_diseases]
		for disease in list_of_diseases:
			# The string [disease] is considered to match if its length is greater than 3 and if it
			# - either contains the string [word]
			# - or is contained in the string [word]
			if len(disease) > 3 and len(normalized_word) > 3 and (disease in normalized_word or normalized_word in disease):
				matches.append(disease)
	# Return the list of disease names associated with the string [word]
	return matches
#-------------------------------------------------------
#INPUTS: the function takes two strings:
# - a string [medinfofile] referring to the name of a file containing a table
# of 3 columns such that its 1st column contains ICD10 codes, its 2nd column contains numerical values
# and its 3rd column contains descriptions of ICD10 codes.
# - a string [diseasefile] referring to a file whose lines each contain a disease name
# OUTPUT: no output -- the function directly writes its outputs in a file "DICT_icd10_PA_diseases.txt"
# containing a table of 2 columns such that its first column contains ICD10 codes and its second column
# contains disease names related to the ICD10 code diagnosis
def disease_Dictionary(medinfofile,diseasefile):
	# List of disease names stored in the file named as [diseasefile]
	disease_list = list_Diseases(diseasefile)
	# Open the file in which the function will write
	savefile = open("DICT_icd10_PA_diseases.txt","w")
	# List of words to not consider when parsing the ICD10 code descriptions of [medinfofile]
	exclusions = get_Exclusions("excluded_words.txt")
	# Open the file named [medinfofile]
	with open(medinfofile, 'r') as g:
		# Go over every line of [medinfofile]
		for i, line in enumerate(g):
			new_line = line.rstrip()
			split_line = new_line.split("\t")
			# [UKB_description] is the diagnosis description associated with the ICD10 code
			UKB_description = split_line[2].split(" ")
			# [saveline] is a list starting with the ICD10 code of the current line. This list will be completed
			# with a list of disease names associated with that specific ICD10 code
			saveline = [split_line[0]]
			# For every word in the ICD10 description, collection in [saveline] the disease names
			# that are not in [exclusions] and are associated with the word specifically
			# (see match_String)
			for word in UKB_description[1:]:
				saveline = saveline + match_String(word.lower(),disease_list,exclusions)
			# Message on the prompt to verify the state of the computation
			print(UKB_description)
			# Write the content of [saveline] on a new line of the file "DICT_icd10_PA_diseases.txt"
			# The line consists of an ICD10 code and a list of disease names
			savefile.write("\t".join(saveline)+"\n")
			savefile.flush()
	# Close the file [savefile]
	savefile.close()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ("disdict" in sys.argv):
	disease_Dictionary("Alzheimer_dict_nz_for_icd10.corrected.modular.txt","PA_diseases.txt")
#-------------------------------------------------------
#-------------------------------------------------------
# Fourth step of the project: Translate the file produced in step 3 into a file that associates
# ICD10 code to lists of related genes.
#-------------------------------------------------------
# INPUTS: a string [disgenfile] referring to the name of a gene-disease network from the
# Percha and Altman database and a string [disease] referring to the name of a disease
# OUTPUT: a list of gene names associated with the input disease in the gene-disease network passed
# in the first argument.
def get_Genes(disgenfile,disease):
	# List containing genes
	genes = []
	# Open the disease-gene network. Gzip is used with gz-files (compressed files)
	with gzip.open(disgenfile, 'r') as f:
		# Go over each line of the network
		for i, line in enumerate(f):
			#.decode() is used in combination with gzip to decode binary to ASCII
			new_line = line.decode().rstrip()
			split_line = new_line.split("\t")
			# If the normalization of the disease name is equal to the normalization
			# of one of the disease in the network, then the associated gene in the network is added to
			# the list [genes]
			if normalization_Text(disease) == normalization_Text(split_line[4]):
				nonrepeating_Append(genes, normalization_Text(split_line[2]))
	# Return the list of genes
	return genes
#-------------------------------------------------------
# INPUTS: the function takes two inputs:
# - a string [dictionary] referring to the name of a file containing a table of 2 columns such that its
# first column contains ICD10 codes and its second column contains lists of diseases associated
# with the ICD10 code
# - a string [disgenfile] referring to the name of a gene-disease network from the
# Percha and Altman database
# OUTPUT: no output -- the function directly writes its outputs in a file
# "DICT_icd10_PA_genes.txt" containing a table of 2 columns such that its first column
# contains ICD10 codes and its second column contains gene names related to the disease names
# associated with the corresponding ICD10 code in [dictionary]
def gene_Dictionary(dictionary,disgenfile):
	# Open the file in which the translation of [dictionary] will be saved
	savefile = open("DICT_icd10_PA_genes.txt","w")
	# Open the file associating ICD10 codes with lists of disease names
	with open(dictionary, "r") as f:
		# Go over each line of the table contained in [dictionary]
		for i, line in enumerate(f):
			# Message on the prompt to verify the state of the computation
			print("line",i)
			new_line = line.rstrip()
			split_line = new_line.split("\t")
			# List containing all the disease names associated with the current line
			diseases = split_line[1:]
			# [saveline] is a list starting with the ICD10 code of the current line. This list will be
			# completed with the translation of the elements of [diseases] into gene names
			saveline = [split_line[0]]
			# For every disease name in [diseases], the loop uses get_Genes and the Percha & Altman
			# network to translate the disease name [diseases[j]] into a list of associated genes
			for j in range(len(diseases)):
				# Message on the prompt to verify the state of the computation
				print(str(j) + ") computing...[" + split_line[0] + "] " + diseases[j])
				# Store in a temporary list [the_genes] all the genes associated with [diseases[j]]
				# in the Percha and Altman network
				the_genes = get_Genes(disgenfile,diseases[j])
				# Message on the prompt to see how many genes have been found for
				# the disease [diseases[j]]
				print("--->",len(the_genes))
				# For every gene names in the list [the_genes], use [nonrepeating_Append] to append the
				# genes of [the_genes] to the list [saveline] without any repetition
				for k in range(len(the_genes)):
					nonrepeating_Append(saveline,the_genes[k])
			# Message on the prompt to see the beginning of the information that will be saved in
			# the file "DICT_icd10_PA_genes.txt"
			print(saveline[:4])
			# Write the content of [saveline] on a new line of the file "DICT_icd10_PA_genes.txt"
			# The line consists of an ICD10 code and a list of related genes
			savefile.write("\t".join(saveline)+"\n")
			savefile.flush()
	# Close the file [savefile]
	savefile.close()
#-------------------------------------------------------
# INPUTS: the function takes two inputs:
# - a string [dictionary] referring to the name of a file containing a table of 2 columns such that its
# first column contains ICD10 codes and its second column contains lists of diseases associated
# with the ICD10 code
# - a string [disgenfile] referring to the name of a gene-disease network from the
# Percha and Altman database
# - a list [codes] containing strings referring to ICD10 codes
# - a string [name] appended to the name of the output file created by the function
# OUTPUT: no output -- the function directly writes its outputs in a file
# "DICT_icd10_PA_genes_"+name+".txt" containing a table of 2 columns such that its first column
# contains ICD10 codes in the list [codes] and its second column contains gene names related to the disease names
# associated with the corresponding ICD10 code in [dictionary]
def gene_Dictionary2(dictionary,diseasegenefile,codes,name):
	# Go over each ICD10 contained in the list [codes]
	for u in range(len(codes)):
		# Open the file in which the translation of [dictionary] will be saved
		savefile = open("DICT_icd10_PA_genes_"+name+".txt","a")
		# Open the file associating ICD10 codes with lists of disease names
		with open(dictionary, "r") as f:
			# Go over each line of the table contained in [dictionary]
			for i, line in enumerate(f):
				# Message on the prompt to verify the state of the computation
				print("line", i)
				new_line = line.rstrip()
				split_line = new_line.split("\t")
				# continue makes the loop go to the next iteration without executing
				# the rest of the code
				if not(split_line[0] == codes[u]):
					continue
				# List containing all the disease names associated with the current line
				diseases = split_line[1:]
				# Write the ICD10 at the beginning of a new line of the file "DICT_icd10_PA_genes_"+name+".txt"
				# The line will consist of an ICD10 code and a list of related genes
				savefile.write(split_line[0])
				savefile.flush()
				# The content of the list [gene_list] will be added progressively to the current
				# line in the file [savefile]
				gene_list = []
				for j in range(len(diseases)):
					# Message on the prompt to verify the state of the computation
					print(str(j) + ") computing...["+codes[u]+"] "+ diseases[j])
					# Store in a temporary list [the_genes] all the genes associated with [diseases[j]]
					# in the Percha and Altman network
					the_genes = get_Genes(diseasegenefile,diseases[j])
					# Message on the prompt to see how many genes have been found for
					# the disease [diseases[j]]
					print("--->",len(the_genes))
					# For every gene names in the list [the_genes], use [nonrepeating_Append] to append the
					# genes of [the_genes] to the list [gene_list] without any repetition and succesively
					# add the gene name to the current line in the file
					for k in range(len(the_genes)):
						found = nonrepeating_Append(gene_list,the_genes[k])
						if found == False:
							savefile.write("\t" +the_genes[k])
							savefile.flush()
			# Add an end of line character to end the line in [savefile]
			savefile.write("\n")
			savefile.flush()
			# Message on the prompt to visualize the beginning of the line added to [savefile]
			print(split_line[0], gene_list[:3])
	# Close the file [savefile]
	savefile.close()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ("gendict" in sys.argv):
	gene_Dictionary("DICT_icd10_PA_diseases.txt","part-ii-dependency-paths-gene-disease-sorted-with-themes.txt.gz")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ("gendict2" in sys.argv):
	codes = []
	codes_filename = sys.argv[2]
	with open(codes_filename+".txt", "r") as f:
		for i, line in enumerate(f):
			s = line.rstrip().split("\t")
			codes.append(s[0])
	gene_Dictionary2("DICT_icd10_PA_diseases.txt",
					 "part-ii-dependency-paths-gene-disease-sorted-with-themes.txt.gz",
					 codes,
					 codes_filename)
#-------------------------------------------------------
#-------------------------------------------------------
# Fifth step of the project: Reorganize a dataset extracted from BioMart into a new file containing
# a table of 2 columns such that its first column contains gene names and its second column contains
# information regarding the location of this gene in the genome. The gene location are formatted to be
# compatible with the graphic library "displaygenelocations.py"
#-------------------------------------------------------
# INPUT:  a string [filename] referring to the name of a file extracted from BioMart and containing
# a table whose columns give the following information:
# Gene stable ID,
# Gene start (bp),
# Gene end (bp),
# Chromosome/scaffold name,
# Gene description,
# Gene name
# OUTPUT: a list containing lists of the form [g,l] where
# - g is a string referring to a gene name from the input file [filename]
# - l is a concatenation of strings of the form c + "_" + s + "_" + e where
#		- c is the chromosome number of the chromosome on which g is located
#		- s is the gene start location of the gene g (in basepair)
#		- e is the gene end location of the gene g (in basepair)
def organize_GeneLocations(filename):
	# List [outpout] will contain lists of the form [g,l] where g is a gene name and l encodes
	# the location of that gen in the genome
	output = []
	# Open the file (extracted from BioMart) that contains gene and location information
	with open(filename,"r") as f:
		# Go over each file of the file
		for i, line in enumerate(f):
			# Do not parse the first line since it contains the titles
			if i!=0:
				s = line.rstrip().split(",")
				#Some gene descriptions (in s[5]) start with a space
				#these descriptions are too vague and needs to be ignored
				if s[5][0] != " ":
					#s[5] is the gene name while s[3] is the chromosomal location and
					# s[1] and s[2] are the start and end locations of the gene s[5], respectively
					output.append([normalization_Text(s[5]),s[3]+"_"+s[1]+"_"+s[2]])
	# The list output is sorted in alphabetic order with respect to the gene names (elements indexed by 0)
	output = sorted(output, key = lambda x : x[0])
	# Return the list [output]
	return output
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ("orgloc" in sys.argv):
	output = organize_GeneLocations("BIOMART_gene_name_location.txt")
	with open("DICT_gene_locations.txt",'w') as f:
		for i in range(len(output)):
			f.write("\t".join(output[i])+"\n")
			f.flush()
#-------------------------------------------------------
#-------------------------------------------------------
# Sixth step of the project: Translate the file produced in step 4 into a file that associates
# ICD10 code to lists of chromosomal locations using the file "DICT_gene_locations.txt" created in step 5.
#-------------------------------------------------------
# INPUT: the function takes two inputs:
# - a string [genlocfile] referring to a file containing a table of two columns such that its first
# column contains gene names and its second column contains a chromosomal location
# in the format "chromosome_start_end"
# - a string [gene] referring to the name of a gene
# OUTPUT: a string containing the chromosomal location of the input gene
def get_Location(genlocfile,gene):
	# The location associated with a gene name is returned as a string
	# using the format "chromosome_start_end"
	location = ""
	# Open the file containing gene names and their genomic locations
	with open(genlocfile, 'r') as f:
		# Go over each line of the file
		for i, line in enumerate(f):
			# .decode() is used to decode binary to ASCII
			new_line = line.rstrip()
			split_line = new_line.split("\t")
			# If the gene name [split_line[0]] of the current line is equal to the beginning of the
			# string [gene] (from index 0 to index len(split_line[0])-1), then store the associated
			# location in the variable [location] and exit the loop.
			if split_line[0] == gene[0:len(split_line[0])]:
				location = split_line[1]
				break
	# Return the location found
	return location
#-------------------------------------------------------
# INPUT: the function takes two inputs:
# - a string [dictionary] referring to a file containing table of two columns such that the first column
# contains ICD10 codes and the second column contains a list of associated genes
# - a string [genlocfile] referring to a file containing a table of two columns such that its first
# column contains gene names and its second column contains a chromosomal location
# in the format "chromosome_start_end"
# OUTPUT: no output -- the function directly writes its outputs in a file
# "DICT_icd10_GRCH38_locations.txt" containing a table of 2 columns such that its first column
# contains ICD10 codes and its second column contains chromosomal locations related to the gene names
# associated with the corresponding ICD10 code in [dictionary]
def location_Dictionary(dictionary,genlocfile):
	# Open file in which ICD10 and their asdsociated chromosomal locations will be saved
	savefile = open("DICT_icd10_GRCH38_locations.txt","w")
	# Open file containing the icd10 and their associated genes to be translated into locations
	with open(dictionary,"r") as f:
		# Go over each line of the file
		for i, line in enumerate(f):
			# Message on the prompt to verify the state of the computation
			print("line",i)
			new_line = line.rstrip()
			split_line = new_line.split("\t")
			# List containing all the gene names associated with the current line
			genes = split_line[1:]
			# [saveline] is a list starting with the ICD10 code of the current line. This list will be
			# completed with the translation of the elements of [genes] into locations
			saveline = [split_line[0]]
			# For every disease name in [genes], the loop uses get_Location and the BioMart database
			# to translate the gene name [genes[j]] into its associated chromosomal location
			for j in range(len(genes)):
				# Message on the prompt to verify the state of the computation
				print("computing... "+split_line[0]+".. "+ genes[j])
				the_location = get_Location(genlocfile,genes[j])
				if the_location != "":
					saveline.append(the_location)
					#print("--->",the_location)
			# Sort the location in alphabetic order
			saveline[1:] = sorted(saveline[1:])
			# Message on the prompt to see the beginning of the information that will be saved in
			# the file "DICT_icd10_GRCH38_locations.txt"
			print(saveline[:4])
			# Write the content of [saveline] on a new line of the file "DICT_icd10_GRCH38_locations.txt"
			# The line consists of an ICD10 code and a list of related genes
			savefile.write("\t".join(saveline) + "\n")
			savefile.flush()
	# Close the file [savefile]
	savefile.close()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ("locdict" in sys.argv):
	location_Dictionary("DICT_icd10_PA_genes.txt","DICT_gene_locations.txt")
#-------------------------------------------------------
#-------------------------------------------------------
# Seventh step of the project: Define a weighted network of the ICD10 codes with respect to
# their proximal chromosomal locations with a margin of 100,000 bp.
#-------------------------------------------------------
# INPUT: a string [ukbgenfile] referring to the name of file containing a table of 2 columns such
# that its first column contains ICD10 codes and its second column contains chromosomal locations
# and a string [icd10_code] containing an ICD10 code
# OUTPUT: a list of strings that contain the chromosomal locations associated with the input ICD10 code
def get_Locations(ukbgenfile,icd10_code):
	# List containing the chromosomal locations associated with the ICD10 code [icd10_code]
	locations = []
	# Open the file containing gene names and their genomic locations
	with open(ukbgenfile,"r") as f:
		# Go over each line of the file
		for i, line in enumerate(f):
			new_line = line.rstrip()
			split_line = new_line.split("\t")
			# If the ICD10 [split_line[0]] of the current line is equal to the input [icd10_code],
			# then store the associated locations in the variable [locations] and exit the loop.
			if split_line[0] == icd10_code:
				locations = split_line[1:]
				break
	# Return the locations found
	return locations
#-------------------------------------------------------
# INPUTS: two strings containing gene location in the format "chromosome_start_end"
# OUTPUT: Boolean value (True or False) indicating whether location1 is near location2 within a
# 100,000 basepair interval
def nearby_Locations(location1,location2):
	# Splits the location information into lists containing
	# - the chromosome number
	# - start of the gene location
	# - end of the gene location
	loc1 = location1.split("_")
	loc2 = location2.split("_")
	# By default, the two locations are considered far away from each other
	nearby = False
	# The locations are compared only if they are on the same chromosome
	# Below we compare the two location with respect to their interval [a1,b1] and [a2,b2] where
	# - a1 is the start of loc1
	# - b1 is the end of loc1
	# - a2 is the start of loc2
	# - b2 is the end of loc2
	try:
		if loc1[0] == loc2[0]:
			# Checks intersection of the form a1...a2...b1
			if int(loc1[1]) <= int(loc2[1]) <= int(loc1[2]):
				nearby = True
			# Checks intersection of the form a1...b2...b1
			if int(loc1[1]) <= int(loc2[2]) <= int(loc1[2]):
				nearby = True
			# Checks intersection of the form b1...(100000 bp)...a2
			if 0 <= int(loc2[1]) - int(loc1[2]) <= 100000:
				nearby = True
			# Checks intersection of the form b2...(100000 bp)...a1
			if 0 <= int(loc1[1]) - int(loc2[2]) <= 100000:
				nearby = True
	except:
		pass
	return nearby
#-------------------------------------------------------
# INPUTS: two lists whose elements are strings containing gene locations in the format
# "chromosome_start_end"
# OUTPUT: list of locations from the first input [locations1] that are nearby locations from the second
# input [locations2]
def intersect_Locations(locations1, locations2):
	# List of locations from [locations1] that are nearby locations present in [locations2]
	intersection = []
	# Go over each location in [locations1]
	for i in range(len(locations1)):
		# Go over each location in [locations2]
		for j in range(len(locations2)):
			# Check whether the location [locations1[i]] is nearby the location [locations2[j]]. If
			# this the case, then [locations1[i]] is stored in the list [intersection]
			if nearby_Locations(locations1[i],locations2[j]):
				intersection.append(locations1[i])
				break
	# Return the locations of [locations1] are were successfully found nearby locations in [locations2]
	return intersection
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ("intersect" in sys.argv):
	save = open("NETW_icd10_GRCH38_locations.txt","w")
	f1 = open("Alzheimer_dict_nz_for_icd10.corrected.modular.txt","r")
	for i, line1 in enumerate(f1):
		f2 = open("Alzheimer_dict_nz_for_icd10.corrected.modular.txt", "r")
		for j, line2 in enumerate(f2):
			s1 = line1.rstrip().split("\t")
			s2 = line2.rstrip().split("\t")
			# Query the lists of chormosomal locations associated with
			# every pair of ICD10 codes (contained in the variables [s1[0]] and [s2[0]])
			a = get_Locations("DICT_icd10_GRCH38_locations.txt",s1[0])
			b = get_Locations("DICT_icd10_GRCH38_locations.txt",s2[0])
			intersection = intersect_Locations(a, b)
			# Compute the cardinal of the intersection
			r = float(len(intersection))
			# Compute the cardinal of the union
			q = float(len(a)) + float(len(b)) - r
			# The Jaccard index of [a] and [b] is contained in the variable [J]
			J = 0
			if q != 0:
				J = r/q
			save.write("\t".join([s1[0], s2[0],str(J)]+intersection) + "\n")
			save.flush()
			print("\t".join([s1[0], s2[0],str(J)]))
	f2.close()
	f1.close()
	save.close()
#-------------------------------------------------------
#-------------------------------------------------------
# Height-th step of the project: Threshold and visualize clusters of ICD10 in a network of ICD10 codes
#-------------------------------------------------------
# INPUTS: a weighted network file [netwfile] between icd10 codes, a list [icd10_codes]
# of all ICD10 codes in the network, a list [lengths] containing the numbers of locations
# associated with each ICD10 code, a float number [theshold]
# OUTPUT: a list consisting of lists of the form [a,b,w] where
# # - a and b are two strings referring to ICD10 codes
# # - w is a  weight of the network between each pair of ICD10 codes whose
# numbers of locations are not 0. The weights are set to 0 if they are less than the value [threshold].
def threshold_Network(netwfile,icd10_codes,lengths,threshold):
	# Thresholded network
	network = []
	# Open the network file
	with open(netwfile,"r") as f:
		# The value of [current_icd10_code] defines a line in the table
		# GO over each line of the network file
		for i, line in enumerate(f):
			s = line.rstrip().split("\t")
			# Message displayed on the prompt to verify the state of the computation
			print("threshold_network", s[0])
			# The number of locations associated with the ICD10 code [s[0]]
			t = lengths[icd10_codes.index(s[0])]
			# Do not display the ICD10 scores that are associated with no locations
			if t == 0:
				continue
			# Weight between the ICD10 codes [s[1]] and [s[2]]
			x = float(s[2])
			# The number of locations associated with the ICD10 code [s[0]]
			u = lengths[icd10_codes.index(s[1])]
			# Do not display the ICD10 scores that are associated with no locations
			if u != 0:
				if x >= threshold:
					network.append([s[0],s[1],x])
				else:
					network.append([s[0],s[1],0])
	# Return the network
	return network

# -------------------------------------------------------
# INPUTS: a list consisting of lists of the form [a,b,w] where
# - a and b are two strings referring to ICD10 codes
# - w is a float value
# OUTPUT: the function has two outputs:
# - a table representing the weights of the network between each pair of ICD10 codes whose
# numbers of locations are not 0.
# - the list of averages weights for every row of the output table
def present_Network(network):
	# The output table
	table = []
	# Variable containing ICD10 codes and average scores of the row associated with the ICD10 code
	averages = []
	# Variable used to store the rows of the table before being added to the table [table]
	row = []
	# The value of [current_icd10_code] defines a line in the table
	current_icd10_code = None
	# GO over each edge of the network
	for i in range(len(network)):
		s = network[i]
		# Message displayed on the prompt to verify the state of the computation
		print("table", s[0])
		# If the first ICD10 code of the line changes, a new line is created in the table
		if current_icd10_code != s[0]:
			current_icd10_code = s[0]
			# Do not include empty rows in the tabler (this can happen at the first iteration)
			if row != []:
				table.append(row)
				averages.append([sum(row)/float(len(row)),sum(row),len(row),current_icd10_code])
			# After adding [row] to [table], start a new row in the table
			row = []
		# Add weight of the edge to the table
		row.append(s[2])
	# Sort ICD10 codes from greatest to lowest averages
	averages = sorted(averages, key = lambda x : -float(x[0]))
	# Return the table
	return table, averages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ("showtable" in sys.argv):
	icd10_codes = []
	lengths = []
	with open("Alzheimer_dict_nz_for_icd10.corrected.modular.txt","r") as f:
		for i, line in enumerate(f):
			s = line.rstrip().split("\t")
			icd10_codes.append(s[0])
			locations = get_Locations("DICT_icd10_GRCH38_locations.txt", s[0])
			lengths.append(len(locations))
			print("lengths", s[0],len(locations))

	network = threshold_Network("NETW_icd10_GRCH38_locations.txt",icd10_codes,lengths,0.7)
	table, averages = present_Network(network)

	network = sorted(network,key = lambda x : -float(x[2]))

	for i in range(0,40):
		print(i,averages[i])

	fig, ax = plt.subplots()
	ax.imshow(table)
	plt.show()
#-------------------------------------------------------
#-------------------------------------------------------
# Ninth step of the project: visualize chromosomal locations using the graphic library
# "displaygenelocations.py" sent to me by a postdoc
#-------------------------------------------------------
if ("displaygen" in sys.argv):
	display_locations(sys.argv[2], "DICT_icd10_GRCH38_locations.txt", chr_lengths, chr_names, genome, chr_colors,1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ("displaynet" in sys.argv):
	display_locations(sys.argv[2] + "\t" + sys.argv[3], "NETW_icd10_GRCH38_locations.txt", chr_lengths, chr_names, genome, chr_colors,2)
#-------------------------------------------------------
# ------------------------------------------------------------------------------

if "shownet" in sys.argv:


	def load_file(filename):
		output = []
		with open(filename,"r") as f:
			for i, s in enumerate(f):
				s = s.rstrip().split("\t")
				output.append(s)
		return output

	edges = load_file("NETW_icd10_GRCH38_locations.txt")
	nodes = load_file("Alzheimer_dict_nz_for_icd10.corrected.modular.txt")

	G = nx.Graph()

	e = []
	n = []

	for i in range(len(edges)):
		if float(edges[i][2]) > 0.8:
			nonrepeating_Append(e,(edges[i][0], edges[i][1]))
			nonrepeating_Append(n,edges[i][0])
			nonrepeating_Append(n,edges[i][1])

	G.add_nodes_from(n)
	G.add_edges_from(e)

	print("Nodes of graph: ")
	print(G.nodes())
	print("Edges of graph: ")
	print(G.edges())

	#nx.draw(G, node_color="orange", edgecolors="black")
	nx.draw(G, node_size = 900, node_color="orange", with_labels = True, edgecolors = "black")
	plt.savefig("simple_path.png")  # save as png
	plt.show()  # display