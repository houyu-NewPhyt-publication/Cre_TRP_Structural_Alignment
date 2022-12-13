#!/bin/python3

import gzip
from subprocess import check_output, check_call

#import the user-defined variables
from user_settings import *
from functions import *

def prepare_paths():
	for input_name, input_path in input_paths.items():
		if input_name != "cobalt_db":
			if not os.path.exists(input_path):
				raise SystemExit(input_path + ' does not exist. Check user_settings.py')
	
	for  output_path in output_paths.values():
		try:
			os.makedirs(output_path, exist_ok = True)
		except OSError as exc:
			if exc.errno != errno.EEXIST:
				raise SystemExit(output_path + ' could not be created.')
			pass


###############################################################################
#                                                                             #
#                                   MAIN CODE                                 #
#                                                                             #
###############################################################################

#Execution of the main code starts here
if __name__ != '__main__':
	exit()

print('Info: Reading paths and setting folders')
prepare_paths()

print('Info: Reading xml file with structure info')
pdb_structures_byname = xml_parser(input_paths["structs_info"])

if os.path.exists(output_paths["savepoints"] + "selected_chains.pkz"):
	with gzip.open(output_paths["savepoints"] + "selected_chains.pkz", 'rb') as f:
		deduped_pdb_chain_IDs_byname = pickle.load(f)
else:
	
	print('Info: Downloading the structures')
	with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
		results = []
		for trp_name, trp_structures in pdb_structures_byname.items():
			for pdb_id in trp_structures.keys():
				#download_structure(pdb_id, input_paths["provided_struct"], output_paths["original_struct"])
				results.append(pool.apply_async(download_structure, args=(pdb_id, input_paths["provided_struct"], output_paths["original_struct"])))
		while (True):
			n_ready = 0
			for result in results:
				if result.ready() == True:
					n_ready += 1
			if n_ready == len(results):
				break
			else:
				print ("\t Waiting to be ready (" + str(n_ready) + "/" + str(len(results)) + ")", end="\r")
				time.sleep(1)

	print('Info: Processing the structures')
	with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
		results = []
		for trp_name, trp_structures in pdb_structures_byname.items():
			for pdb_id in trp_structures.keys():
				results.append(pool.apply_async(filter_pdb_residues, args=(trp_name, pdb_id, pdb_structures_byname, output_paths["original_struct"], output_paths["clean_struct_dir"])))
		while (True):
			n_ready = 0
			for result in results:
				if result.ready() == True:
					n_ready += 1
			if n_ready == len(results):
				break
			else:
				print ("\t Waiting to be ready (" + str(n_ready) + "/" + str(len(results)) + ")", end="\r")
				time.sleep(1)

	print('Info: Splitting structures into chains')
	pdb_chain_IDs_byname = {}
	with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
		results = []
		for trp_name, trp_structures in pdb_structures_byname.items():
			for pdb_id in trp_structures.keys():
				#split_structure(trp_name, pdb_id, pdb_structures_byname, output_paths["clean_struct_dir"], output_paths["clean_chain_dir"])
				results.append(pool.apply_async(split_structure, args=(trp_name, pdb_id, pdb_structures_byname, output_paths["clean_struct_dir"], output_paths["clean_chain_dir"])))
		while (True):
			n_ready = 0
			for result in results:
				if result.ready() == True:
					n_ready += 1
			if n_ready == len(results):
				for result in results:
					for trp_name, trp_structures in result.get().items():
						pdb_chain_IDs_byname.setdefault(trp_name, {}).update(trp_structures)
				break
			else:
				print ("\t Waiting to be ready (" + str(n_ready) + "/" + str(len(results)) + ")", end="\r")
				time.sleep(1)

	print('Info: Removing dupliates')
	deduped_pdb_chain_IDs_byname = {}
	with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
		results = []
		for trp_name, trp_structures in pdb_chain_IDs_byname.items():
			#remove_duplicates(trp_name, trp_structures, pdb_structures_byname, output_paths["clean_chain_dir"])
			results.append(pool.apply_async(remove_duplicates, args=(trp_name, trp_structures, pdb_structures_byname, output_paths["clean_chain_dir"])))
		while (True):
			n_ready = 0
			for result in results:
				if result.ready() == True:
					n_ready += 1
			if n_ready == len(results):
				for result in results:
					for trp_name, trp_structures in result.get().items():
						deduped_pdb_chain_IDs_byname.setdefault(trp_name, {}).update(trp_structures)
				break
			else:
				print ("\t Waiting to be ready (" + str(n_ready) + "/" + str(len(results)) + ")", end="\r")
				time.sleep(1)
	
	with gzip.open(output_paths["savepoints"] + "selected_chains.pkz", 'wb') as f:
		pickle.dump(deduped_pdb_chain_IDs_byname, f, pickle.HIGHEST_PROTOCOL)

if os.path.exists(output_paths["savepoints"] + "frtmalign_results.pkz"):
	with gzip.open(output_paths["savepoints"] + "frtmalign_results.pkz", 'rb') as f:
		frtmalign_results = pickle.load(f)
else:
	print('Info: Running Fr-TM-Align to align all structures pairwise')
	frtmalign_results = {}
	with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
		results = []
		for trp_name, trp_structures in deduped_pdb_chain_IDs_byname.items():
			for pdb_ID, chain_IDs in trp_structures.items():
				#In order to make the funtion work for both PDB structures and individual chains
				for chain_ID in chain_IDs:
					#Build the pdbIDN to simplify the code
					IDN_stat = pdb_ID + "_" + chain_ID + "." + trp_name
					stat_dir = output_paths["aligned_chain_dir"] + IDN_stat + "/"
					try:
						os.makedirs(stat_dir, exist_ok = True)
					except OSError as exc:
						if exc.errno != errno.EEXIST:
							raise SystemExit(output_path + ' could not be created.')
						pass
					
					#run_frtmalign(IDN_stat, output_paths["clean_chain_dir"], deduped_pdb_chain_IDs_byname, output_paths["clean_chain_dir"], out_path = stat_dir, out_file = True, sup_file = True, rot_mat = True)
					results.append(pool.apply_async(run_frtmalign, args=(IDN_stat, output_paths["clean_chain_dir"], deduped_pdb_chain_IDs_byname, output_paths["clean_chain_dir"], stat_dir, True, True, True)))
		while (True):
			n_ready = 0
			for result in results:
				if result.ready() == True:
					n_ready += 1
			if n_ready == len(results):
				for result in results:
					frtmalign_results.update(result.get())
				break
			else:
				print ("\t Waiting to be ready (" + str(n_ready) + "/" + str(len(results)) + ")", end="\r")
				time.sleep(1)
	
	with gzip.open(output_paths["savepoints"] + "frtmalign_results.pkz", 'wb') as f:
		pickle.dump(frtmalign_results, f, pickle.HIGHEST_PROTOCOL)
#pprint(frtmalign_results)

print('Info: Processing Fr-TM-Align results')
frtmalign_2_tables(frtmalign_results, output_paths["aligned_chain_dir"], file_basename = "")

bt = best_target(frtmalign_results, output_paths["aligned_chain_dir"], file_basename = "")

print('Info: Downloading the reference sequences and replacing them')

#pprint(frtmalign_results[bt[1] + "." + bt[0]])
pairwise_alignments = []
for name, alignment_values in frtmalign_results[bt[1] + "." + bt[0]].items():
	if name != "length":
		pairwise_alignments.append(alignment_values["alignment"])

adj_pairwise_alignments = pairwise2reference(pairwise_alignments, pdb_structures_byname)

print('Info: Generating structural MSA')
ref_name = adj_pairwise_alignments[0]["target"]["name"]
ref_seq = adj_pairwise_alignments[0]["target"]["seq"].replace("-","")

struct_msa = pairwise2msa(ref_name, ref_seq, adj_pairwise_alignments)
with open(output_paths["msa"] + "msa_struct_gappy.fasta", "wt", encoding="UTF-8") as file:
	file.write(list2fasta([struct_msa["target"],]))
	file.write(list2fasta(struct_msa["queries"]))

print('Info: Aligning gappy regions')
mafft_params = [
		input_paths["mafft_bin"],
		"--thread", "60",
		"--threadtb", "30",
		"--threadit", "0",
		"--inputorder",
		"--maxiterate", "1000",
		"--retree", "100",
		"--genafpair",
		"--nomemsave",
		"--seed", output_paths["msa"] + "msa_struct_gappy.fasta",
		"-"
	]
	
mafft_result = check_output(mafft_params, input = "", encoding="UTF-8")
with open(output_paths["msa"] + "msa_struct_aligned.fasta", "wt", encoding="UTF-8") as file:
	file.write(mafft_result.replace("_seed_",""))


print("Info: Downloading the TRP sequences and extracting TM domains")
#In this case the moment that the results are generated, they are written in the
#fasta file and not stored in memory because we don't need to access them from python

mafft_result = open(output_paths["msa"] + "msa_struct_aligned.fasta", "rt", encoding = "UTF-8").read()

def down_and_get(a, b, c, fasta_refs):
	
	full_seq = sequence_download(a,b,c)
	
	if full_seq:
		tm_boundaries = tm_seq_boundaries(list2fasta(full_seq), fasta_refs)
	return full_seq, tm_boundaries


with open(input_paths["trp_seqs_table"], 'r', newline='') as seqs_table_file:
	#The expected format of the table is a csvreader object with 3 columns DB	AccNR	Desired_header
	reader = csv.reader(seqs_table_file, delimiter="\t")
	next(reader, None)  # skip the headers
	with open(output_paths["sequences"] + "sequences_full.fasta", "w") as full_seqs_file:
		with open(output_paths['sequences'] + "sequences_tmonly.fasta", "w") as tm_seqs_file:
			with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
				results = []
				for row in reader:
					print(" Processing " + row[1] + "...")
					#resu = down_and_get(row[0],row[1],row[2], fasta_ref_seqs)
					#print(resu)
					#full_seqs_file.write(list2fasta(resu[0]))
					#tm_seqs_file.write(list2fasta(resu[0], resu[1]))
					#down_and_get(row[0],row[1],row[2], mafft_result.replace("-",""))
					results.append(pool.apply_async(down_and_get, args=(row[0],row[1],row[2], mafft_result.replace("-",""))))
				while (True):
					n_ready = 0
					for result in results:
						if result.ready() == True:
							n_ready += 1
					if n_ready == len(results):
						for result in results:
							if result.get():
								full_seqs_file.write(list2fasta(result.get()[0]))
								tm_seqs_file.write(list2fasta(result.get()[0], result.get()[1]))
						break
					else:
						print ("\t Waiting to be ready (" + str(n_ready) + "/" + str(len(results)) + ")", end="\r")
						time.sleep(1)


print('Info: Building multiple sequence alignment')
mafft_params = [
		input_paths["mafft_bin"],
		"--thread", "60",
		"--threadtb", "30",
		"--threadit", "0",
		"--inputorder",
		"--maxiterate", "1000",
		"--retree", "100",
		"--genafpair",
		"--nomemsave",
		"--dash",
		output_paths['sequences'] + "sequences_tmonly.fasta"
	]

mafft_result = check_output(mafft_params, encoding="UTF-8")
with open(output_paths["mafft_alignment"] + "alignment.fasta", "wt", encoding="UTF-8") as file:
	file.write(mafft_result.replace("_seed_",""))


print("Info: Running IQ-Tree")
iqtree_params = [
	input_paths["iqtree_bin"],
	#General parameters
	"-T", "60", #Number of threads to use
	"-v", #Verbose mode
	"-s", output_paths["mafft_alignment"] + "alignment.fasta", #MSA file
	"-pre", output_paths["iqtree_results"] + "iqtree_results", #Location and prefix for the IQTree files
	"-st", "AA", #Specify sequence type
	
	#Parameters for the Model Finder
	"-msub", "nuclear", #Restrict to those AA models designed for specified source
	"-m", "MFP", #Automatic model selection plus tree search
	"-madd", "C10,C20,C30,C40,C50,C60,EX2,EX3,EHO,UL2,UL3,EX_EHO,LG4M,LG4X,CF4", #Include this protein mixture models
	"-mwopt", #Turn on optimizing weights of mixture models
	
	
	#Parameters", "for the Tree search"
	"--runs", "1", #Number of independent rounds 
	"-nm", "1000", #Specify maximum number of iterations for UFBoot to stop
	"-bb", "1000", #Specify number of UF bootstrap replicates
	"-bnni", #Perform an additional step to further optimize UFBoot trees by nearest neighbor interchange (NNI) based directly on bootstrap alignments
	"-alrt", "1000", #Specify number of replicates (>=1000) to perform SH-like approximate likelihood ratio test (SH-aLRT)
	"-nstop", "100", #Specify number of unsuccessful iterations to stop.
	"-ninit", "100", #Specify number of initial parsimony trees
	"-ntop", "100",
	"-pers", "0.2",
	"-wbtl",
]

iqtree_results = check_call(iqtree_params)
