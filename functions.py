#!/bin/python3
from Bio import pairwise2
from Bio import SeqIO
import pathlib
import subprocess
import json
import re
import csv
import io
import sys
import threading
import queue
import time
import glob
from subprocess import check_output

#Checked required imports
import copy
import errno
import itertools
import multiprocessing
import os
import pickle
import re
import requests
import shutil
import tempfile
import warnings

from Bio import pairwise2
from Bio import PDB
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.Seq import Seq
from pprint import pprint
from xml.dom import minidom

from user_settings import *

class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKCYAN = '\033[96m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'

def xml_parser(filename):
	pdb_structures_byname = {}
	xml_doc = minidom.parse(filename)
	#Iterate through all the structures in the file
	for xml_structure in xml_doc.getElementsByTagName('structure'):
		#As a PDB file can contain more than a single protein the dictionary names
		#are organized formed by PROTEIN_NAME and then pdb_ID to make them unique.
		#We can cope with multiple structure blocks for the different proteins in the xml file
		#as well as for the same protein and different chains.
		#To avoid problems with case-sensitivity and files, PROTEIN_NAME and pdb_IDs will be
		#always uppercase
		
		#Initialize a dictionary for the PROTEIN_NAME and pdb_ID.
		entry_name = pdb_structures_byname.setdefault(xml_structure.getAttribute("name").upper(), {}).setdefault(xml_structure.getAttribute("pdb_id").upper(),{})
		#Load the information of the structure into the dictionary
		entry_name["subfamily"] = xml_structure.getAttribute("subfamily")
		entry_name["resolution"] = float(xml_structure.getAttribute("resolution"))
		entry_name["src_org"] = xml_structure.getAttribute("src_org")
		entry_name["uniprotkb"] = xml_structure.getAttribute("uniprotkb")
		#Initialize a dictionary for the chains
		entry_name.setdefault("chains", {})
		#Iterate through the chains inside the structure
		for xml_chain in xml_structure.getElementsByTagName('chain'):
			chain_ID = xml_chain.getAttribute("id")
			#Initialize a list to contain the TM ranges. Most of the times there will be a single range
			entry_name["chains"][chain_ID] = []
			#Iterate through the TM ranges inside the chain
			for xml_tmrange in xml_chain.getElementsByTagName('tmrange'):
				entry_name["chains"][chain_ID].append((int(xml_tmrange.getAttribute("start")), int(xml_tmrange.getAttribute("end"))))
	
	#pprint.pprint(pdb_structures_byname)
	return pdb_structures_byname


#It expects a 4 digit string with the PDB ID and UPPERCASE and returns True if succesful
def download_structure(pdb_ID, local_structures, savedir):
	if os.path.exists(savedir + pdb_ID + '.pdb'):
		return True
	else:
		#First tries to get structure from OPM
		url = 'https://opm-assets.storage.googleapis.com/pdb/' + pdb_ID.lower() + '.pdb'
		res = requests.get(url, allow_redirects = True)
		if not res.status_code == 200:
			#If structure not found in OPM, look for user-provided structure file
			print("Warning: Found no record in OPM for %s. Checking %s." % (pdb_ID, local_structures))
			try:
				shutil.copyfile(local_structures + pdb_ID + '.pdb', savedir + pdb_ID + '.pdb')
			except:
				# If no user-provided structure file, try to get structure from PDB
				print("Warning: Found no provided structure file for %s in %s. Checking PDB." %(pdb_ID, local_structures))
				pdb_url = 'https://files.rcsb.org/download/' + pdb_ID + '.pdb'
				pdb_res = requests.get(pdb_url, allow_redirects = True)
				if not pdb_res.status_code == 200:
					# If structure not found in PDB, print warning and skip structure
					print("Warning: found no record in OPM, %s, or PDB for %s, so it will be ignored." % (local_structures, pdb_ID))
					return False
				else:
					open(savedir + pdb_ID + '.pdb', 'wb').write(pdb_res.content)
					return True
			return True
		else:
			open(savedir + pdb_ID + '.pdb', 'wb').write(res.content)
			return True

def filter_pdb_residues(protein_name, pdb_ID, pdb_structures_byname, source_struct_dir, dest_struct_dir):
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', PDBConstructionWarning)
		#Filter out everything but the amino acids in the TM ranges in the selected chains
		class ResSelect(PDB.Select):
			#If you want to write out a part of the structure, make use of the Select class (also in PDBIO). Select has four methods:
			#  accept_model(model)
			#  accept_chain(chain)
			#  accept_residue(residue)
			#  accept_atom(atom)
			#
			#By default, every method returns 1 (which means the model/chain/residue/atom is included in the output).
			#By subclassing Select and returning 0 when appropriate you can exclude models, chains, etc. from the output.
			#Cumbersome maybe, but very powerful.
			def accept_chain(self, chain):
				#Only the chain we are looping through will be accepted
				for chain_ID in pdb_structures_byname[protein_name][pdb_ID]["chains"].keys():
					if chain.get_id() == chain_ID:
						return True
				else:
					return False
			def accept_residue(self, residue):
				#Only the residues that are not HETATM, in the ranges of the TM and
				#in the right chain will be accepted.
				for residue_range in pdb_structures_byname[protein_name][pdb_ID]["chains"][residue.parent.id]:
					if(
						residue.id[0] == " " and
						residue.id[1] >= int(residue_range[0]) and
						residue.id[1] <= int(residue_range[1])
					):
						return True
				else:
					return False
		
		#Saving function here to get each full structure filtered
		pdb_parser = PDB.PDBParser()
		pdb_writer = PDB.PDBIO()
		pdb_structure = pdb_parser.get_structure(pdb_ID, source_struct_dir + pdb_ID + ".pdb")
		pdb_writer.set_structure(pdb_structure)
		pdb_writer.save(dest_struct_dir + pdb_ID + "." + protein_name + ".pdb", ResSelect())
		
		return True

#pdb_IDN refers to a string combining pdb_ID and protein name: PRD_ID.PROTEIN_NAME
def run_frtmalign(stat_pdb_chain_IDN, stat_path, mobi_pdb_chain_IDs_byname, mobi_path, out_path = None, out_file = False, sup_file = False, rot_mat = False):
	
	#Check if a working path has been selected
	if out_path:
		#Create the working folder if it doesn't exist
		try:
			os.makedirs(out_path, exist_ok = True)
		except:
			raise SystemExit(out_path + ' could not be created.')
		
		#Readjust the paths relative to frtmalign working directory. Absolute paths will crash.
		mobi_path = os.path.relpath(os.path.abspath(mobi_path), out_path) + "/"
		stat_path = os.path.relpath(os.path.abspath(stat_path), out_path) + "/"
	
	
	#Dictionary to contain the results after paralellizing
	results = {}
	for mobi_protein_name, mobi_pdb_chain_IDs in mobi_pdb_chain_IDs_byname.items():
		for mobi_pdb_ID, mobi_chain_IDs in mobi_pdb_chain_IDs.items():
			for mobi_chain_ID in mobi_chain_IDs:
				#Generate stat_pdb_IDN to simplify the rest of the code
				mobi_pdb_chain_IDN = mobi_pdb_ID + "_" + mobi_chain_ID + "." + mobi_protein_name
				
				print("\t Aligning " + mobi_pdb_chain_IDN + " to " + stat_pdb_chain_IDN)
				
				#Check if a working path has been selected
				if out_path:
					#Check wether to generate a sup file
					if sup_file:
						sup_arg = mobi_pdb_chain_IDN + "_" + stat_pdb_chain_IDN + ".sup"
					else:
						sup_arg = "/dev/null"
					
					#Check if the alignment has been previously generated
					if os.path.exists(os.path.join(out_path, mobi_pdb_chain_IDN + "_" + stat_pdb_chain_IDN + ".frtxt")): 
						#Load the previously generated alignment
						with open(os.path.join(out_path, mobi_pdb_chain_IDN + "_" + stat_pdb_chain_IDN + ".frtxt"), "rt", encoding="utf-8") as out_filehandle:
							out_list = out_filehandle.read().splitlines()
							
							#If the output format of Fr-TM-Align ever changes, it needs to be readjusted
							if not len(out_list) == 20:
								raise SystemExit('The file ' + os.path.join(out_path, mobi_pdb_chain_IDN + "_" + stat_pdb_chain_IDN + ".frtxt") + " is invalid. Check it and restart again")
					else:
						#Run Fr-TM-Align saving results to files
						out = check_output([input_paths["frtmalign"], mobi_path + mobi_pdb_chain_IDN + ".pdb", stat_path + stat_pdb_chain_IDN + ".pdb", "-o", sup_arg, "-m ", str(int(rot_mat))], cwd = out_path, encoding="UTF-8")
						
						#If enabled, for each mobile file, a transformation matrix ("trf.mat") will be created.
						#Rename the file "trf.mat" in output_path and save for future use
						if rot_mat:
							if os.path.exists(out_path + 'trf.mat'): 
								shutil.copy2(out_path + 'trf.mat', out_path + mobi_pdb_chain_IDN + '_' + stat_pdb_chain_IDN + '.mat')
							else:
								raise SystemExit('trf.mat does not exist.')
						
						with open(out_path + mobi_pdb_chain_IDN + "_" + stat_pdb_chain_IDN + ".frtxt", 'wt', encoding = "UTF-8") as out_filehandle:
							out_filehandle.write(out)
						
						out_list = out.splitlines()
						if len(out_list) == 20:
							if out_file:
								with open(out_path + mobi_pdb_chain_IDN + "_" + stat_pdb_chain_IDN + ".fasta", 'w+') as out_filehandle:
									out_filehandle.write(">" + mobi_pdb_chain_IDN + "\n")
									out_filehandle.write(out_list[16] + "\n")
									out_filehandle.write("\n")
									out_filehandle.write(">" + stat_pdb_chain_IDN + "\n")
									out_filehandle.write(out_list[18] + "\n")
						else:
							raise SystemExit(out)
				else:
					#Run Fr-TM-Align without saving results to files
					out = check_output([input_paths["frtmalign"], mobi_path + mobi_pdb_chain_IDN + ".pdb", stat_path + stat_pdb_chain_IDN + ".pdb", "-m ", "0"], encoding="UTF-8")
					out_list = out.splitlines()
				if len(out_list) == 20:
					results.setdefault(stat_pdb_chain_IDN, {})["length"] = int(out_list[11].split("=")[1].split("(")[0])
					results[stat_pdb_chain_IDN].setdefault(mobi_pdb_chain_IDN, {})["length"] = int(out_list[10].split("=")[1])
					results[stat_pdb_chain_IDN][mobi_pdb_chain_IDN]["alignment_length"] = int(out_list[13].split(",")[0].split("=")[1])
					results[stat_pdb_chain_IDN][mobi_pdb_chain_IDN]["rmsd"] = float(out_list[13].split(",")[1].split("=")[1])
					results[stat_pdb_chain_IDN][mobi_pdb_chain_IDN]["tm_score"] = float(out_list[13].split(",")[2].split("=")[1])
					results[stat_pdb_chain_IDN][mobi_pdb_chain_IDN]["n_ident"] = float(out_list[13].split(",")[3].split("=")[1])
					results[stat_pdb_chain_IDN][mobi_pdb_chain_IDN]["alignment"] = {"target":{"name":stat_pdb_chain_IDN, "seq":out_list[18].strip().upper()}, "query":{"name":mobi_pdb_chain_IDN, "seq":out_list[16].strip().upper()}}
				else:
					raise SystemExit(out)
				
				
	#Fr-TM-Align returns a dictionary of IDNs because we need to be totally specific and keep a simple structure for the data
	return results


def remove_duplicates(protein_name, pdb_chain_IDs, pdb_structures_byname, structures_folder):
	#Note: This function assumes that all structures belong to the same protein
	
	try:
		os.makedirs(structures_folder + "duplicates/", exist_ok = True)
	except OSError as exc:
		if exc.errno != errno.EEXIST:
			raise SystemExit(structures_folder + ' could not be created.')
		pass
	
	#List which will contain all the IDNs to check against each other
	#In order to loop properly and do the combinations, we can't have
	#any depth layers.
	IDNs_list = []
	
	for pdb_ID, chain_IDs in pdb_chain_IDs.items():
		#In order to make the funtion work for both PDB structures and individual chains
		for chain_ID in chain_IDs:
			#Build the pdbIDN to simplify the code
			IDNs_list.append([protein_name, pdb_ID, chain_ID])
		
	index_target = 0
	index_query = 0
	while index_target < len(IDNs_list):
		
		#Build the pdbIDN to simplify the code
		IDN_target = IDNs_list[index_target][1] + "_" + IDNs_list[index_target][2] + "." + IDNs_list[index_target][0]
		
		#Check wether the item was previously deduped
		if not os.path.exists(structures_folder + IDN_target + ".pdb"):
			if os.path.exists(structures_folder + "duplicates/" + IDN_target + ".pdb"):
				#File previously deduped and probably process aborted before finishing, reusme from here
				print("\t " + protein_name + "(" + str(len(IDNs_list)) + "): " + str(index_target) + "," + str(index_query), "Aligned previously and deduped")
				del IDNs_list[index_target]
				#Go go the next loop iteration of target without increasing the index counter, as there is
				#a new item in that index.
				continue
			else:
				raise SystemExit("Previous run corrupted. Please start over")
		
		#Reset index_query for the next iteration of index_target. As 1,2 is equal to 2,1 we won't recheck past
		#values of index_target or the alignment against itself.
		index_query = index_target + 1
		
		while index_query < len(IDNs_list):
			#Build the pdbIDN to simplify the code
			IDN_query = IDNs_list[index_query][1] + "_" + IDNs_list[index_query][2] + "." + IDNs_list[index_query][0]
			
			if not os.path.exists(structures_folder + IDN_query + ".pdb"):
				if os.path.exists(structures_folder + "duplicates/" + IDN_query + ".pdb"):
					#file previously deduped and probably process aborted before finishing, reusme from here
					print("\t " + protein_name + "(" + str(len(IDNs_list)) + "): " + str(index_target) + "," + str(index_query), "Aligned previously and deduped")
					del IDNs_list[index_query]
					break
				else:
					raise SystemExit("Previous run corrupted. Please start over")
			
			frtmalign_results = run_frtmalign(IDN_target, structures_folder, {IDNs_list[index_query][0]:{IDNs_list[index_query][1]:[IDNs_list[index_query][2],]}}, structures_folder)
			
			if (
				frtmalign_results[IDN_target][IDN_query]["alignment_length"] == frtmalign_results[IDN_target]["length"] == frtmalign_results[IDN_target][IDN_query]["length"] and
				frtmalign_results[IDN_target][IDN_query]["rmsd"] < 1.00 and
				frtmalign_results[IDN_target][IDN_query]["tm_score"] > 0.99 and
				frtmalign_results[IDN_target][IDN_query]["n_ident"] == 1.00
			):
				#The two structures are considered identical
				print("\t " + protein_name + "(" + str(len(IDNs_list)) + "): " + str(index_target) + "," + str(index_query), "Aligned length =", frtmalign_results[IDN_target][IDN_query]["alignment_length"], "RMSD =", frtmalign_results[IDN_target][IDN_query]["rmsd"], "TM-score =", frtmalign_results[IDN_target][IDN_query]["tm_score"], "ID =", frtmalign_results[IDN_target][IDN_query]["n_ident"], "<- Considering them identical, omitting")
				#print("Considering them identical, popping it out")
				
				#Find the structure with the best resolution
				#print("Finding the best resolution structure")
				
				if pdb_structures_byname[IDNs_list[index_target][0]][IDNs_list[index_target][1]]["resolution"] <= pdb_structures_byname[IDNs_list[index_query][0]][IDNs_list[index_query][1]]["resolution"]:
					#If the best structure is in index_target, del index_query and just continue
					#Note: The next iteration is in the same index index_query as item in this position and there is index_target new item in place.
					#Note2: It is quicker and more efficient when index_query is the item removed
					#print("Best is index_target")
					del IDNs_list[index_query]
					
					#Move the file to folder/duplicates
					try:
						shutil.move(os.path.join(structures_folder, IDN_query + ".pdb"), os.path.join(structures_folder, "duplicates", IDN_query + ".pdb"))
					except:
						raise SystemExit("File couldn't be moved")
					
					#print(len(pdb_IDs))
				elif pdb_structures_byname[IDNs_list[index_target][0]][IDNs_list[index_target][1]]["resolution"] > pdb_structures_byname[IDNs_list[index_query][0]][IDNs_list[index_query][1]]["resolution"]:
					#If the best structure is in index_query, del index_target and break the loop.
					#print("Best is index_query")
					del IDNs_list[index_target]
					
					#Move the file to folder/duplicates
					try:
						shutil.move(os.path.join(structures_folder, IDN_target + ".pdb"), os.path.join(structures_folder, "duplicates", IDN_target + ".pdb"))
					except:
						raise SystemExit("File couldn't be moved")
					
					#print(len(pdb_IDs))
					
					#We break the loop so index_target doesn't increase as there is index_target new item in place and we need to restart the loop in index_query
					break
				else:
					raise SystemExit('Problem reading the resolution in the structures')
			else:
				print("\t " + protein_name + "(" + str(len(IDNs_list)) + "): " + str(index_target) + "," + str(index_query), "Aligned length =", frtmalign_results[IDN_target][IDN_query]["alignment_length"], "RMSD =", frtmalign_results[IDN_target][IDN_query]["rmsd"], "TM-score =", frtmalign_results[IDN_target][IDN_query]["tm_score"], "ID =", frtmalign_results[IDN_target][IDN_query]["n_ident"])
				#If no identity, increase index_query to go to the next iteration
				index_query += 1
		else:
			#If while index_query ends by exhaustion increase index_target and go to the next iteration
			index_target += 1
		#print(len(pdb_IDs))
		#print("-------------------------------------------------------")
	#print(len(pdb_IDs))
	
	#Return a pdb_chain_IDs dictionary
	#We need to return all the info including protein_name as the parent function might be parallelized
	deduped_pdb_chain_IDs_byname = {}
	
	for IDN_list in IDNs_list:
		deduped_pdb_chain_IDs_byname.setdefault(IDN_list[0], {}).setdefault(IDN_list[1], []).append(IDN_list[2])
	
	return deduped_pdb_chain_IDs_byname


def split_structure(protein_name, pdb_ID, pdb_structures_byname, source_struct_dir, dest_chain_dir):
		#Save the chains individually
		chain_IDs = []
		
		for chain_ID in pdb_structures_byname[protein_name][pdb_ID]["chains"].keys():
			class ChainSelect(PDB.Select):
				def accept_chain(self, chain):
					if chain.get_id() == chain_ID:
						return True
					else:
						return False 
			pdb_parser = PDB.PDBParser()
			pdb_writer = PDB.PDBIO()
			pdb_structure = pdb_parser.get_structure(pdb_ID + "." + protein_name, source_struct_dir + pdb_ID + "." + protein_name + ".pdb")
			pdb_writer.set_structure(pdb_structure)
			pdb_writer.save(dest_chain_dir + pdb_ID + "_" + chain_ID + "." + protein_name + ".pdb", ChainSelect())
			chain_IDs.append (chain_ID)
		
		#Remove the chains that are identical
		deduped_pdb_chain_IDs = remove_duplicates(protein_name, {pdb_ID:chain_IDs}, pdb_structures_byname, dest_chain_dir)
		
		#Returns a dict deduped_pdb_chain_IDs
		return deduped_pdb_chain_IDs



def frtmalign_2_tables(frtmalign_results, out_folder, file_basename = ""):
	
	
	#summary_table = [[None, None, None, None, None, None, None, None] for row_nr in range(1 + ((len(table_header_column)-1)*(len(table_header_column)-1)))]
	summary_table = [["Stationary IDN", "Stationary length", "Mobile IDN", "Mobile length", "Alignment length", "RMSD", "TM-score", "Seq ident"]]

	#Lists to dimension the table for the TSV files to disk with the results
	#The item in the lists is for cell 0,0 in the table
	table_col_names = set()
	table_row_names = set()
	
	#It is important to use sorted() to guarantee the order of the items in the dict
	for stat_pdb_IDN, stat_data in sorted(frtmalign_results.items()):
		#Add stat_pdb_IDN to the header row for the table
		table_col_names.add(stat_pdb_IDN)
		
		for mobi_pdb_IDN, mobi_data in sorted(stat_data.items()):
			if not mobi_pdb_IDN == "length":
				#Add mobi_pdb_IDN to the header column for the table
				table_row_names.add(mobi_pdb_IDN)
				#Add to the summary table
				summary_table.append([mobi_pdb_IDN, mobi_data["length"], stat_pdb_IDN, stat_data["length"], mobi_data["alignment_length"], mobi_data["rmsd"], mobi_data["tm_score"], mobi_data["n_ident"]])
	
	
	table_col_names = [''] + sorted(table_col_names)
	table_row_names = [''] + sorted(table_row_names)
	
	#print(table_col_names)
	#print(table_row_names)
	
	#Create template table
	results_table = [[None for col_nr in range(len(table_col_names))] for row_nr in range(len(table_row_names))]
	
	#Write headers
	results_table[0][0] = "Mobile\Stationary"
	for col_nr in range(0, len(table_col_names)):
		results_table[0][col_nr] = table_col_names[col_nr]
	for row_nr in range(0, len(table_row_names)):
		results_table[row_nr][0] = table_row_names[row_nr]

	
	#Create tables
	RMSD_table = copy.deepcopy(results_table)
	TM_score_table = copy.deepcopy(results_table)
	n_ident_table = copy.deepcopy(results_table)
	alignment_length_table = copy.deepcopy(results_table)
	
	#Fill the values in the table. We start the range at 1 so we don't overwrite the headers
	#As the headers are PDB_IDNs we split them in order to access the data in the dictionary
	for row_nr in range(1, len(table_row_names)):
		for col_nr in range(1, len(table_col_names)):
			RMSD_table[row_nr][col_nr] = frtmalign_results[table_col_names[col_nr]][table_row_names[row_nr]]["rmsd"]
			TM_score_table[row_nr][col_nr] = frtmalign_results[table_col_names[col_nr]][table_row_names[row_nr]]["tm_score"]
			n_ident_table[row_nr][col_nr] = frtmalign_results[table_col_names[col_nr]][table_row_names[row_nr]]["n_ident"]
			alignment_length_table[row_nr][col_nr] = frtmalign_results[table_col_names[col_nr]][table_row_names[row_nr]]["alignment_length"]
			
	#Save the files
	with open(out_folder + file_basename + 'alignment_rmsd.tsv', 'w', newline='') as csvfile:
		table_writer = csv.writer(csvfile, dialect = "excel-tab")
		table_writer.writerows(RMSD_table)
		
	with open(out_folder + file_basename + 'alignment_TM_score.tsv', 'w', newline='') as csvfile:
		table_writer = csv.writer(csvfile, dialect = "excel-tab")
		table_writer.writerows(TM_score_table)
		
	with open(out_folder + file_basename + 'alignment_identity.tsv', 'w', newline='') as csvfile:
		table_writer = csv.writer(csvfile, dialect = "excel-tab")
		table_writer.writerows(n_ident_table)
		
	with open(out_folder + file_basename + 'alignment_length.tsv', 'w', newline='') as csvfile:
		table_writer = csv.writer(csvfile, dialect = "excel-tab")
		table_writer.writerows(alignment_length_table)
		
	with open(out_folder + file_basename + 'alignment_output.tsv', 'w', newline='') as csvfile:
		table_writer = csv.writer(csvfile, dialect = "excel-tab")
		table_writer.writerows(summary_table)


def best_target(frtmalign_results, out_folder, file_basename = ""):
	#print(frtmalign_results_byname)
	
	#Sets to dimension the table for the TSV files to disk with the results
	#The item in the lists is for cell 0,0 in the table
	table_col_names = set()
	
	#Table row names will also be used as an internal counter for how many proteins of each kind were
	
	
	#Calculate the average value for each protein name.
	#As there are different amount of crystals for each protein, we need to
	#find the average in order not to bias the final results.
	#It is important to use sorted() to guarantee the order of the items in the dict
	stat_pdb_IDNs_vs_mobi_protein_names = {}
	for stat_pdb_IDN, stat_data in sorted(frtmalign_results.items()):
		#Add stat_pdb_IDN to the header row for the table
		table_col_names.add(stat_pdb_IDN)
		#Delete the only attribute that is not a pdb entry so we can iterate
		table_row_names = {}
		for mobi_pdb_IDN, mobi_data in sorted(stat_data.items()):
			if mobi_pdb_IDN != "length":
				mobi_protein_name = mobi_pdb_IDN.split(".")[1]
				table_row_names.setdefault(mobi_protein_name,0)
				table_row_names[mobi_protein_name] += 1
				stat_pdb_IDNs_vs_mobi_protein_names.setdefault(stat_pdb_IDN,{}).setdefault(mobi_protein_name, {}).setdefault("rmsd", 0.00)
				stat_pdb_IDNs_vs_mobi_protein_names[stat_pdb_IDN][mobi_protein_name].setdefault("tm_score", 0.00)
				stat_pdb_IDNs_vs_mobi_protein_names[stat_pdb_IDN][mobi_protein_name].setdefault("n_ident", 0.00)
				stat_pdb_IDNs_vs_mobi_protein_names[stat_pdb_IDN][mobi_protein_name].setdefault("alignment_length", 0.00)
				
				#Sum all the values from all the mobile entries
				stat_pdb_IDNs_vs_mobi_protein_names[stat_pdb_IDN][mobi_protein_name]["rmsd"] += mobi_data["rmsd"]
				stat_pdb_IDNs_vs_mobi_protein_names[stat_pdb_IDN][mobi_protein_name]["tm_score"] += mobi_data["tm_score"]
				stat_pdb_IDNs_vs_mobi_protein_names[stat_pdb_IDN][mobi_protein_name]["n_ident"] += mobi_data["n_ident"]
				stat_pdb_IDNs_vs_mobi_protein_names[stat_pdb_IDN][mobi_protein_name]["alignment_length"] += mobi_data["alignment_length"]
		for protein_name, number in table_row_names.items():
			#Average all the values from all the mobile entries
			stat_pdb_IDNs_vs_mobi_protein_names[stat_pdb_IDN][protein_name]["rmsd"] /= number
			stat_pdb_IDNs_vs_mobi_protein_names[stat_pdb_IDN][protein_name]["tm_score"] /= number
			stat_pdb_IDNs_vs_mobi_protein_names[stat_pdb_IDN][protein_name]["n_ident"] /= number
			stat_pdb_IDNs_vs_mobi_protein_names[stat_pdb_IDN][protein_name]["alignment_length"] /= number
	#print("To the table")
	table_col_names = [''] + sorted(table_col_names)
	table_row_names = [''] + sorted(table_row_names)
	#pprint.pprint(table_col_names)
	#pprint.pprint(table_row_names)
	
	#Create template table
	results_table = [[None for col_nr in range(len(table_col_names))] for row_nr in range(len(table_row_names))]
	
	#Write headers
	results_table[0][0] = "Mobile\Stationary"
	#print(table_col_names)
	#print(table_row_names)
	for col_nr in range(0, len(table_col_names)):
		results_table[0][col_nr] = table_col_names[col_nr]
	for row_nr in range(0, len(table_row_names)):
		results_table[row_nr][0] = table_row_names[row_nr]
	
	#Create tables
	RMSD_table = copy.deepcopy(results_table)
	TM_score_table = copy.deepcopy(results_table)
	n_ident_table = copy.deepcopy(results_table)
	alignment_length_table = copy.deepcopy(results_table)
	
	#Fill the values in the table. We start the range at 1 so we don't overwrite the headers
	for row_nr in range(1, len(table_row_names)):
		for col_nr in range(1, len(table_col_names)):
			RMSD_table[row_nr][col_nr] = stat_pdb_IDNs_vs_mobi_protein_names[table_col_names[col_nr]][table_row_names[row_nr]]["rmsd"]
			TM_score_table[row_nr][col_nr] = stat_pdb_IDNs_vs_mobi_protein_names[table_col_names[col_nr]][table_row_names[row_nr]]["tm_score"]
			n_ident_table[row_nr][col_nr] = stat_pdb_IDNs_vs_mobi_protein_names[table_col_names[col_nr]][table_row_names[row_nr]]["n_ident"]
			alignment_length_table[row_nr][col_nr] = stat_pdb_IDNs_vs_mobi_protein_names[table_col_names[col_nr]][table_row_names[row_nr]]["alignment_length"]
			
	
	#For the average, we can append a row to the existing tables and add the data
	RMSD_table.append(['Average',])
	TM_score_table.append(['Average',])
	n_ident_table.append(['Average',])
	alignment_length_table.append(['Average',])
	
	#Calculate the average of each value
	stat_pdb_IDNs_averages = {}
	for stat_pdb_IDN, mobi_protein_names in stat_pdb_IDNs_vs_mobi_protein_names.items():
		
		#Initialize the variables for the mathematical operations
		#If setdefault is used directly it crashes
		stat_pdb_IDNs_averages.setdefault(stat_pdb_IDN,{})["rmsd"] = 0
		stat_pdb_IDNs_averages[stat_pdb_IDN]["tm_score"] = 0
		stat_pdb_IDNs_averages[stat_pdb_IDN]["n_ident"] = 0
		stat_pdb_IDNs_averages[stat_pdb_IDN]["alignment_length"] = 0
		#Average all the values from the proteins
		for mobi_protein_name, mobi_avg_data in mobi_protein_names.items():
			stat_pdb_IDNs_averages[stat_pdb_IDN]["rmsd"] += mobi_avg_data["rmsd"]/len(mobi_protein_names)
			stat_pdb_IDNs_averages[stat_pdb_IDN]["tm_score"] += mobi_avg_data["tm_score"]/len(mobi_protein_names)
			stat_pdb_IDNs_averages[stat_pdb_IDN]["n_ident"] += mobi_avg_data["n_ident"]/len(mobi_protein_names)
			stat_pdb_IDNs_averages[stat_pdb_IDN]["alignment_length"] += mobi_avg_data["alignment_length"]/len(mobi_protein_names)
		
	#Add the values using the column names again to guarantee the same order
	for col_nr in range(1,len(table_col_names)):
		RMSD_table[len(RMSD_table)-1].append(stat_pdb_IDNs_averages[table_col_names[col_nr]]["rmsd"])
		TM_score_table[len(TM_score_table)-1].append(stat_pdb_IDNs_averages[table_col_names[col_nr]]["tm_score"])
		n_ident_table[len(n_ident_table)-1].append(stat_pdb_IDNs_averages[table_col_names[col_nr]]["n_ident"])
		alignment_length_table[len(alignment_length_table)-1].append(stat_pdb_IDNs_averages[table_col_names[col_nr]]["alignment_length"])
	
	#To create the table for the scores, we mostly copy the rows from the previous tables
	#and then add the score row to the end of it.
	stat_pdb_IDNs_score_table = []
	
	stat_pdb_IDNs_score_table.append(RMSD_table[0])
	
	stat_pdb_IDNs_score_table.append(RMSD_table[len(RMSD_table)-1])
	stat_pdb_IDNs_score_table[len(stat_pdb_IDNs_score_table)-1][0] = "RSMD average"
	
	stat_pdb_IDNs_score_table.append(TM_score_table[len(TM_score_table)-1])
	stat_pdb_IDNs_score_table[len(stat_pdb_IDNs_score_table)-1][0] = "TM score average"
	
	stat_pdb_IDNs_score_table.append(n_ident_table[len(n_ident_table)-1])
	stat_pdb_IDNs_score_table[len(stat_pdb_IDNs_score_table)-1][0] = "Identity average"
	
	stat_pdb_IDNs_score_table.append(alignment_length_table[len(alignment_length_table)-1])
	stat_pdb_IDNs_score_table[len(stat_pdb_IDNs_score_table)-1][0] = "Alignment length average"
	
	
	#Sum up all the values for the diffrerent stationary IDNs to
	#normalize the values and choose the best target
	rmsd_max = 0
	tm_score_max = 0
	n_ident_max = 0
	alignment_length_max = 0
	
	for stat_averages in stat_pdb_IDNs_averages.values():
		rmsd_max = max(rmsd_max, stat_averages["rmsd"])
		tm_score_max = max(tm_score_max, stat_averages["tm_score"])
		n_ident_max = max(n_ident_max, stat_averages["n_ident"])
		alignment_length_max = max(alignment_length_max, stat_averages["alignment_length"])
	
	#Normalize the 4 parameters and calculate the score
	#The scoring works as follows:
	# 1. Values are normalized 0 to 1 by dividing with the max value
	# 2. The value for the RMSD is inverted (1-value) so higher is better like the rest
	# 3. The values are normalized to give the final score. All values have the same weight.
	stat_pdb_IDNs_score = {}
	for stat_pdb_IDN, stat_averages in sorted(stat_pdb_IDNs_averages.items()):
		stat_pdb_IDNs_score[stat_pdb_IDN] = ((1-(stat_averages["rmsd"]/rmsd_max)) + stat_averages["tm_score"]/tm_score_max + stat_averages["n_ident"]/n_ident_max + stat_averages["alignment_length"]/alignment_length_max)/4
	
	#print(stat_pdb_IDNs_score)
	
	#Add the values using the column names again to guarantee the same order
	stat_pdb_IDNs_score_table.append(["Score"])
	for col_nr in range(1,len(table_col_names)):
		stat_pdb_IDNs_score_table[len(stat_pdb_IDNs_score_table)-1].append(stat_pdb_IDNs_score[table_col_names[col_nr]])
	
	#Save the files
	with open(out_folder + file_basename + 'protein_avg_rmsd.tsv', 'w', newline='') as csvfile:
		table_writer = csv.writer(csvfile, dialect = "excel-tab")
		table_writer.writerows(RMSD_table)
		#print(RMSD_table)
		
	with open(out_folder + file_basename + 'protein_avg_TM_score.tsv', 'w', newline='') as csvfile:
		table_writer = csv.writer(csvfile, dialect = "excel-tab")
		table_writer.writerows(TM_score_table)
		
	with open(out_folder + file_basename + 'protein_avg_identity.tsv', 'w', newline='') as csvfile:
		table_writer = csv.writer(csvfile, dialect = "excel-tab")
		table_writer.writerows(n_ident_table)
		
	with open(out_folder + file_basename + 'protein_avg_length.tsv', 'w', newline='') as csvfile:
		table_writer = csv.writer(csvfile, dialect = "excel-tab")
		table_writer.writerows(alignment_length_table)
		
	with open(out_folder + file_basename + 'scores_as_target.tsv', 'w', newline='') as csvfile:
		table_writer = csv.writer(csvfile, dialect = "excel-tab")
		table_writer.writerows(stat_pdb_IDNs_score_table)
	
	#choose the candidate with the highest score as target for the rest of the analysis
	max_score = max(stat_pdb_IDNs_score, key=stat_pdb_IDNs_score.get)
	
	#return protein_name, pdb_ID, score
	return (max_score.split(".")[1], max_score.split(".")[0], stat_pdb_IDNs_score[max_score])
	


def fasta2list(source):
	
	if isinstance(source, io.IOBase):
		print("Source is an open file")
		fasta_lines = source.readlines()
	elif isinstance(source, pathlib.Path):
		if source.is_file():
			with open(source, 'r') as fasta_file:
				fasta_lines = fasta_file.read().splitlines()
	elif isinstance(source, str):
		try:
			pathlib.Path(source).is_file()
		except:
			#print("Source is a string")
			fasta_lines = source.splitlines()
		else:
			if pathlib.Path(source).is_file():
				#print("Source is a file pointer string")
				with open(source, 'r') as fasta_file:
					fasta_lines = fasta_file.read().splitlines()
			else:
				#print("Source is a string")
				fasta_lines = source.splitlines()
	seq_list = []
	sequence_id = ""
	sequence = ""
	for line in fasta_lines:
		if line.startswith(">"):
			if sequence != "":
				#Check if sequence_id contains start and end delimitersin in
				#JalView style and append previous sequence to list
				if len(sequence_id.rsplit("/")) == 22: #FIXME with regex
					if len(sequence_id.rsplit("/")[1].split("-")) == 2:
						seq_list.append({
							"name":sequence_id.rsplit("/")[0], 
							"seq":sequence, 
							"start": int(sequence_id.rsplit("/")[1].split("-")[0]), 
							"end": int(sequence_id.rsplit("/")[1].split("-")[1])
							})
				else:
					seq_list.append({"name":sequence_id, "seq":sequence, "start": "0", "end": "0"})
				#Load the header of the next sequence and empty
				#the sequence buffer
				sequence = ""
				sequence_id = line.strip()[1:]
			else:
				sequence_id = line.strip()[1:]
		else:
			sequence += line.strip()
	else: #When the loop is over, in order to add the last sequence stored
		if sequence != "" and sequence_id != "":
			#Here code with last of the sequences
			if len(sequence_id.rsplit("/")) == 22: #FIXME with regex
				if len(sequence_id.rsplit("/")[1].split("-")) == 2:
					seq_list.append({
						"name":sequence_id.rsplit("/")[0], 
						"seq":sequence, 
						"start": int(sequence_id.rsplit("/")[1].split("-")[0]), 
						"end": int(sequence_id.rsplit("/")[1].split("-")[1])
						})
			else:
				seq_list.append({"name":sequence_id, "seq":sequence, "start": 0, "end": 0})

	return seq_list

def list2fasta(seq_list, boundaries = None):
	seq_fasta = ""
	if len(seq_list) == 2:
		if type(seq_list["name"]) is str and type(seq_list["seq"]) is str:
			seq_list = [seq_list,]
	for fasta_entry in seq_list:
		if boundaries:
			seq_fasta += ">" + fasta_entry["name"] + "\n" + fasta_entry["seq"][boundaries[0]:boundaries[1]] + "\n\n"
		else:
			seq_fasta += ">" + fasta_entry["name"] + "\n" + fasta_entry["seq"] + "\n\n"
	return seq_fasta

def sequence_download(database, accession_nr, new_header = None):
	sequence = ""
	
	if database.lower() == "uniprot":
		#Download the sequence for the UniProtKB ID
		res = requests.get("https://rest.uniprot.org/uniprotkb/" + accession_nr + ".fasta", allow_redirects = True)
		if not res.status_code == 200:
			print(bcolors.FAIL + "   Error: Couldn't communicate with the UniProtKB server" + bcolors.ENDC)
		if res.text == "":
			print(bcolors.FAIL + "   Error: Couldn't download the sequence. Check the ID" + bcolors.ENDC)
		else:
			sequence = fasta2list(res.text)[0]["seq"]
			if not new_header:
				new_header = accession_nr
	
	elif database.lower() == "jgi":
		#Fetch the local copy of the JGI proteome using the genome identifier tag
		loc_prot = fasta2list(input_paths["local_proteomes"] + accession_nr.split("|")[0] + ".fasta")
		
		#Look in all the FASTA entries for the right protein.
		for prot in loc_prot:
			#As one number can be the substring to a larger number, we need to add the delimiters
			#on both sides of the number to be sure it is the full number. We check for all possible
			#delimiters used in JGI genomes as well as for the genome identifier tag.
			if ("=" + accession_nr.split("|")[1] + " " in prot["name"] or
			"|" + accession_nr.split("|")[1] + "|" in prot["name"] or
			"|" + accession_nr.split("|")[1] + "." in prot["name"] or
			"|" + accession_nr.split("|")[1] + "\n" in prot["name"]
			):
				sequence = prot["seq"]
				if not new_header:
					new_header = accession_nr 
				break
	
	elif database.lower() == "ccds":
		
		#Download the mapping table from NCBI
		res = requests.get('https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS2Sequence.current.txt', allow_redirects = True)
		if not res.status_code == 200:
			print(bcolors.FAIL + "   Error: Couldn't communicate with the NCBI server" + bcolors.ENDC)
		else:
			print("   Mapping CCDS to NCBI proteins")
			#The file has the following structure:
			#ccds	original_member	current_member	source	nucleotide_ID	protein_ID	status_in_CCDS	sequence_status
			#In order to traverse quickly we will do a dictionary of dictionaries
			csv_file = io.StringIO(res.text)
			reader = csv.reader(csv_file, delimiter="\t")
			for row in reader:
				if accession_nr in row[0]:
					if row[2] == "1":
						if row[3] =="NCBI":
							print("   Fetching the NCBI sequence")
							sequence = sequence_download("NCBI", row[5])["seq"]
							if not new_header:
								new_header = accession_nr + " (" + row[5] + ")"
							break
	
	elif database.lower() == "ectsi":
		res = requests.get('https://bioinformatics.psb.ugent.be/orcae/annotation/Ectsi/current/' + accession_nr, allow_redirects = True)
		#If the server response is not OK (200)
		if not res.status_code == 200:
			print(bcolors.FAIL + "   Error: Couldn't communicate with the Ectsi server" + bcolors.ENDC)
		else:
			prot = re.search(r'<textarea name="data\[Protein\]\[sequence\]".*?>(.*?)</textarea>', res.text, re.M|re.I)
			if prot:
				sequence = prot.group(1).strip()
				if not new_header:
					new_header = accession_nr 
	
	elif database.lower() == "flybase":
		#Documentation at https://flybase.github.io/api/swagger-ui/
		url = 'https://api.flybase.org/api/v1.0/sequence/id/' + accession_nr
		res = requests.get(url)
		#If the server response is not OK (200)
		if not res.status_code == 200:
			print(bcolors.FAIL + "   Error: Couldn't communicate with the Flybase REST API" + bcolors.ENDC)
		else:
			sequence = json.loads(res.text)["resultset"]["result"][0]["sequence"]
			if not new_header:
				new_header = accession_nr 
	
	elif database.lower() == "ncbi":
		res = requests.get('https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=protein&report=fasta&id=' + accession_nr, allow_redirects = True)
		#If the server response is not OK (200)
		if not res.status_code == 200:
			print(bcolors.FAIL + "   Error: Couldn't communicate with the NCBI REST API" + bcolors.ENDC)
		else:
			sequence = fasta2list(res.text)[0]["seq"]
			if not new_header:
				new_header = accession_nr
	
	elif database.lower() == "wormbase":
		#Documentation at http://rest.wormbase.org/index.html
		url = 'http://rest.wormbase.org/rest/field/protein/' + accession_nr + '/sequence'
		res = requests.get(url)
		#If the server response is not OK (200)
		if not res.status_code == 200:
			print(bcolors.FAIL + "   Error: Couldn't communicate with the WormBase REST API" + bcolors.ENDC)
		else:
			sequence = json.loads(res.text)["sequence"]["data"]["sequence"]
			if not new_header:
				new_header = accession_nr
	
	else:
		loc_prot = fasta2list(input_paths["local_proteomes"] + database + ".fasta")
		for prot in loc_prot:
			if accession_nr in prot["name"]:
				sequence = prot["seq"]
				if not new_header:
					new_header = accession_nr 
	
	if sequence != "":
		#print("   " + accession_nr + " ---> " + new_header)
		return {"name":new_header, "seq":sequence.strip("*").upper()}
	
	else:
		print(bcolors.WARNING + "   Warning: No sequence for", accession_nr + bcolors.ENDC)
		return None

def pairwise2msa (ref_name, ref_seq, alignments):
	#It does not align gappy regions, if there is a gap in one of the sequences
	#from the pairwise alignment, the residue from the other won't align with
	#anything else. but a gap. A following refinement step will take care of it.
	ref_seq = ref_seq.upper()
	results_seqs = [list(ref_seq),]
	results_names = [ref_name,]
	for alignment in alignments:
		
		alignment["target"]["seq"] = alignment["target"]["seq"].upper()
		alignment["query"]["seq"] = alignment["query"]["seq"].upper()
		
		#Reference sequence does not have to be the full length sequence but one has to be contained in the other to work
		#If one sequence is contained in the other multiple times, only the first repetition will be accounted and aligned for
		#This pre-processing is necessary because as letters repeat, we need to ensure
		#we start comparing at the right place, otherwise the alignment would start at the wrong place
		if ref_seq == alignment["target"]["seq"].replace("-",""):
			target_list = list(alignment["target"]["seq"])
			query_list = list(alignment["query"]["seq"])
		elif ref_seq in alignment["target"]["seq"].replace("-",""):
			
			target_list = list(alignment["target"]["seq"])
			query_list = list(alignment["query"]["seq"])
			
			#Ungap the target and find where the match starts in the ungapped ref_seq.
			#This way we can align the beginnings which is all that matters
			ref_start = alignment["target"]["seq"].index(ref_seq.replace("-",""))
			ref_end = target_start + len(ref_seq.replace("-",""))
			
			#Add the missing sequence to the reference
			results_seqs[0].insert(0, list(alignment["target"]["seq"][:ref_start]))
			results_seqs[0].append(list(alignment["target"]["seq"][ref_end:]))
			
			#Update all the sequences already aligned with a gap to compensate for the
			#inserted sequence in both ends of the reference to keep parity in the alignment
			for i in range(1, len(results_seqs)):
				for _ in range(ref_start):
					results_seqs[i].insert(0, "-")
				for _ in range(len(ref_seq) - target_end):
					results_seqs[i].append("-")
		elif alignment["target"]["seq"].replace("-","") in ref_seq:
			#Ungap the target and find where the match starts in the ungapped ref_seq.
			#This way we can align the beginnings which is all that matters
			target_start = ref_seq.index(alignment["target"]["seq"].replace("-",""))
			target_end = target_start + len(alignment["target"]["seq"].replace("-",""))
			
			#Add the missing sequence of the reference to the target in both ends
			target_list = list(ref_seq[:target_start]) + list(alignment["target"]["seq"]) + list(ref_seq[target_end:])
			#Add a gap to compensate for the inserted sequence in both ends of the
			#query to keep parity in the pairwise alignment
			query_list = ["-"]*target_start + list(alignment["query"]["seq"]) + ["-"]*(len(ref_seq) - target_end)
		else:
			print(alignment["target"]["seq"].replace("-",""))
			print(ref_seq)
			raise Exception('References are not the same! Check sequences')
			
			
		i = 0
		while (i < len(results_seqs[0])) or (i < len(target_list)):
			#When ref_seq is longer than target, gaps after target is finished
			#We should not be able to reach the end of results before the target ends
			#as we add gaps for it every time we loop through the target, that's also
			#the reason to loop with while and not with for ..range()
			#print(i, len(results_seqs[0]), len(target_list))
			if i >= len(results_seqs[0]):
				for j in range(len(results_seqs)):
					results_seqs[j].append("-")
			elif i >= len(target_list):
				target_list.append("-")
				query_list.append("-")
			#If ref_seq has a gap, add a gap in the sequences. That will move
			#all the items one position for the next iteration.
			elif results_seqs[0][i] == "-":
				#print("Gap ref_seq")
				target_list.insert(i, "-")
				query_list.insert(i, "-")
			#If target has a gap, add a gap in all the sequences already aligned,
			#including in the reference sequence.
			elif target_list[i] == "-":
				#print("Gap_target")
				for j in range(len(results_seqs)):
					results_seqs[j].insert(i, "-")
			#If anything bizarre happens, stop
			elif target_list[i] != results_seqs[0][i]:
				#print("".join(target_list))
				#print("".join(results_seqs[0]))
				raise Exception('Unexpected behavior! Check sequences at pos:' + str(i))
			#If the residues match , let them be
			#else:
				#print("Match", target_list[i])
				#pass
			i += 1
		#When we append the values, pass a slice of the size of the ref_seq. 
		#In the case that there were gap characters after the last residue of the
		#ref_seq, that would cause a malformed FASTA file.
		results_seqs.append(query_list)
		results_names.append(alignment["query"]["name"])
	
	results_dict = {"target":{"name": results_names.pop(0),"seq": "".join(results_seqs.pop(0))},"queries":[]}
	for res_name, res_seq in zip(results_names, results_seqs):
		results_dict["queries"].append({"name": res_name, "seq": "".join(res_seq)})
	
	return results_dict


def pairwise2reference(pairwise_alignments, pdb_structures_byname):
	#This function takes a pairwise alignment, aligns the reference sequence with the 
	#target sequence of the pairwise alignment, then aligns the query sequence with the
	#the re-aligned target and then replaces the target for the reference used
	#It is thought to be used with sequences coming from structural files, which can be
	#missing residues from the loops and other parts missing in the model.
	#The pairwise alignments are expected in the following format:
	#list[
	#	{"target":{
	#		"name": str,
	#		"seq": str
	#	},
	#	"query":{
	#		"name": str,
	#		"seq": str
	#	},...
	#]
	
	adj_pairwise_alignments = []
	
	#Dictionary to act as a cache for downloaded sequences
	ref_sequences = {}
	for pairwise_alignment in pairwise_alignments:
	
		#Check if the UniProtKB sequence is in ref_sequences and if not, download it
		if pdb_structures_byname[pairwise_alignment["target"]["name"].split(".")[1]][pairwise_alignment["target"]["name"].split("_")[0]]["uniprotkb"] in ref_sequences:
			target_ref = ref_sequences[pdb_structures_byname[pairwise_alignment["target"]["name"].split(".")[1]][pairwise_alignment["target"]["name"].split("_")[0]]["uniprotkb"]]
		else:
			target_ref = sequence_download("uniprot", pdb_structures_byname[pairwise_alignment["target"]["name"].split(".")[1]][pairwise_alignment["target"]["name"].split("_")[0]]["uniprotkb"])["seq"]
			ref_sequences[pdb_structures_byname[pairwise_alignment["target"]["name"].split(".")[1]][pairwise_alignment["target"]["name"].split("_")[0]]["uniprotkb"]] = target_ref
		
		if pdb_structures_byname[pairwise_alignment["query"]["name"].split(".")[1]][pairwise_alignment["query"]["name"].split("_")[0]]["uniprotkb"] in ref_sequences:
			query_ref = ref_sequences[pdb_structures_byname[pairwise_alignment["query"]["name"].split(".")[1]][pairwise_alignment["query"]["name"].split("_")[0]]["uniprotkb"]]
		else:
			query_ref = sequence_download("uniprot", pdb_structures_byname[pairwise_alignment["query"]["name"].split(".")[1]][pairwise_alignment["query"]["name"].split("_")[0]]["uniprotkb"])["seq"]
			ref_sequences[pdb_structures_byname[pairwise_alignment["query"]["name"].split(".")[1]][pairwise_alignment["query"]["name"].split("_")[0]]["uniprotkb"]] = query_ref
		
		print("\t " + pairwise_alignment["query"]["name"], end="\r")
		#print("Reference:			" + target_ref)
		#print("Target:				" + pairwise_alignment["target"]["seq"])
		#print("Query:				" + pairwise_alignment["query"]["seq"])
		#print()
		
		#####################################################################################
		#### STEP 1: Replace gapped target sequence for ungapped sequence from UniProtKB ####
		#####################################################################################
		
		#Aligning reference to the target. The target is a gapped fragment of the reference from the crystal structure
		#In order to align properly the original gaps at the beginning and end of the sequence, we will use
		#the sign "+" for the gaps we insert to differentiate them.
		
		alignments = pairwise2.align.globalms(Seq(target_ref), Seq(pairwise_alignment["target"]["seq"]), 2, -1, -0.5, -0.1, one_alignment_only = True, gap_char = "+")
		
		#Take the results of the alignnment as adjusted reference and target
		adj_target_ref = alignments[0][0]
		adj_target_seq = alignments[0][1]
		
		#print("Adj. target ref:		" + adj_target_ref)
		#print("Adj. target:			" + adj_target_seq)
		#print()
		
		#Use the adjusted target to align the query. Residues in the reference
		#but missing in the target will become gaps in both sequences. Gaps
		#pre-existing in the original alignment will be kept as they were, 
		#maintaining full parity between original target-query and the adjusted one.
		
		target_list = list(pairwise_alignment["target"]["seq"])
		query_list = list(pairwise_alignment["query"]["seq"])
		
		#As adj_target_seq is immutable, we can access the index of the string characters directly
		#Use the target as guide to incorporate the gaps in the query
		for i in range(len(adj_target_seq)):
			#If the target is shorter than the reference, fill the rest with gaps
			if i >= len(target_list):
				target_list.append("+")
				query_list.append("+")
			#If it's not a match, proceed
			elif adj_target_seq[i] != target_list[i]:
				#If the mismatch is due to a gap, insert the gap
				if adj_target_seq[i] == "+":
					#print("Gap reference")
					target_list.insert(i, "+")
					query_list.insert(i, "+")
				#If there is a mismatch in the sequence, abort as the sequences should be the same.
				else:
					raise Exception('Unexpected behavior! Check sequences at pos:' + str(i))
			#print(target_list[i] + "<>" + query_list[i])
		
		#print("Adj. target:			" + "".join(target_list))
		#print("Adj. query:			" + "".join(query_list))
		#print()
		
		#Pair the adjusted reference with the adjusted query
		#print("New target:			" + adj_target_ref)
		#print("Adj. query:			" + "".join(query_list))
		#print()
		
		#Trim the gaps in both ends (original and not) and leave only the core aligned part.
		#For that we just have to look at the position of the first and last characters in
		#the query and trim the reference accordingly
		
		pattern = re.compile(r'[^+]')
		query_start = int(pattern.search("".join(query_list)).start())
		query_end = len(query_list) - pattern.search("".join(query_list[::-1])).start()
		
		#Pair the adjusted reference with the adjusted query
		#print("New target:			" + adj_target_ref[query_start:query_end])
		#print("New query:			" + "".join(query_list[query_start:query_end]))
		#print()
		
		#####################################################################################
		#### STEP 2: Reassign the target and the query to the new adjusted sequences     ####
		#####################################################################################
		new_target = adj_target_ref[query_start:query_end]
		new_query = "".join(query_list[query_start:query_end])
		
		
		#####################################################################################
		#### STEP 3: Replace gapped query sequence for ungapped sequence from UniProtKB  ####
		#####################################################################################
		
		alignments = pairwise2.align.globalms(Seq(query_ref), Seq(new_query), 2, -1, -0.5, -0.1, one_alignment_only = True, gap_char = "*")
		
		#Take the results of the alignnment as adjusted reference and target
		adj_query_ref = alignments[0][0]
		adj_query_seq = alignments[0][1]
		
		new_target_list = list(new_target)
		new_query_list = list(new_query)
		
		#print("Adj. new query ref:		" + adj_query_ref)
		#print("Adj. new query:			" + adj_query_seq)
		#print()
		
		
		#As adj_query_seq is immutable, we can access the index of the string characters directly
		#Use the target as guide to incorporate the gaps in the query
		for i in range(len(adj_query_seq)):
			#If the target is shorter than the reference, fill the rest with gaps
			if i >= len(new_query_list):
				new_target_list.append("*")
				new_query_list.append("*")
			#If it's not a match, proceed
			elif adj_query_seq[i] != new_query_list[i]:
				#If the mismatch is due to a gap, insert the gap
				if adj_query_seq[i] == "*":
					#print("Gap reference")
					new_target_list.insert(i, "*")
					new_query_list.insert(i, "*")
				#If there is a mismatch in the sequence, abort as the sequences should be the same.
				else:
					print(adj_query_seq[i] + "<>" + new_query_list[i])
					raise Exception('Unexpected behavior! Check sequences at pos:' + str(i))
			
		#print("Adj. new target:		" + "".join(new_target_list))
		#print("Adj. new query:			" + "".join(new_query_list))
		#print()
		
		
		pattern = re.compile(r'[^*]')
		new_target_start = int(pattern.search("".join(new_target_list)).start())
		new_target_end = len(new_target_list) - pattern.search("".join(new_target_list[::-1])).start()
		
		#Pair the adjusted reference with the adjusted query
		#print("Final ref target:		" + "".join(new_target_list[new_target_start:new_target_end]))
		#print("Final ref query:		" + adj_query_ref[new_target_start:new_target_end])
		#print()
		
		#The last step is to replace + for - in the new gaps and create the alignment object
		adj_pairwise_alignments.append({
			"target":{
				"name": pairwise_alignment["target"]["name"],
				"seq": "".join(new_target_list[new_target_start:new_target_end]).replace("+", "-").replace("*", "-")
			},
			"query":{
				"name": pairwise_alignment["query"]["name"],
				"seq": adj_query_ref[new_target_start:new_target_end].replace("+", "-").replace("*", "-")
			}
		})
		
	return adj_pairwise_alignments


def tm_seq_boundaries(query_seq, reference_seqs):
	cobalt_params = [
		input_paths["cobalt_bin"],
		"-gapopen", "15",
		"-gapextend", "1",
		"-end_gapopen", "5",
		"-end_gapextend", "1",
		"-iter", "T",
		"-rps_evalue", "1",
		"-blast_evalue", "0.1",
		"-clusters", "F",
		#"-k", "4",
		#"-max_dist", "0.8",
		#"-alph", "regular",
		"-rpsdb", input_paths["cobalt_db"],
		"-outfmt", "mfasta",
		"-i", "-"
	]
	cobalt_input = query_seq + "\n" + reference_seqs
	print("   Aligning with reference sequences")
	
	cobalt_result = check_output(cobalt_params, input = cobalt_input, encoding="UTF-8")
	
	#with open(results_dir + query_seq.split("\n")[0].strip(">").replace("|","-_-").replace("/","-") + ".fasta", "wt", encoding = "UTF-8") as alignment_file:
	#	alignment_file.write(cobalt_result)
	
	aligned = fasta2list(cobalt_result)
	
	#The start of the core protein will be given by the first residue in the reference and
	#the end of the core protein will be given by the last residue in the reference
	start_with_gaps = None
	end_with_gaps = None
	print("   Determining TM domain boundaries")
	#We will loop through the whole length of the query, which is the first sequence
	for aa_idx in range(len(aligned[0]["seq"])):
		
		#If the beginning and the end are already found, escape the loop
		if start_with_gaps != None and end_with_gaps != None:
			break
		else:
			#We start looping the sequences from 1 as 0 is the query sequence
			for seq_idx in range(1,len(aligned)):
				if aligned[seq_idx]["seq"][aa_idx] != "-" and start_with_gaps == None:
						start_with_gaps = aa_idx
				
				if aligned[seq_idx]["seq"][len(aligned[seq_idx]["seq"])-aa_idx-1] != "-" and end_with_gaps  == None:
						end_with_gaps = len(aligned[seq_idx]["seq"])-aa_idx-1
				
				if start_with_gaps != None and end_with_gaps != None:
					break
	#start = start_with_gaps - aligned[0]["seq"][0:start_with_gaps].count("-")
	#end = end_with_gaps - aligned[0]["seq"][start_with_gaps:end_with_gaps].count("-")
	print(aligned[0]["name"] + " TM start: " + str(start_with_gaps) + ", end: " + str(end_with_gaps))
	#print(">" + aligned[0][0] + "\n")
	#print(aligned[0][1][start_with_gaps:end_with_gaps] + "\n\n")
	#print(aligned[0][1].replace("-","")[start:end])
	
	return [start_with_gaps,end_with_gaps]

