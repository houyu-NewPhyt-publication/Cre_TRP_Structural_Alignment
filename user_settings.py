#!/bin/python3

#Prevent standalone execution
if __name__ == '__main__':
	exit()

working_dir = "./"

#Define the working folders for the script either absolute or relative to working_dir
input_paths = {
	#Locations and files that MUST be adjusted
	"provided_struct"	: working_dir + "00_provided_structures/", #Folder containing the proteomes (in FASTA) for JGI and local sequences
	"structs_info"		: working_dir + "00_provided_files/struct_new.xml", #XML file containing which chains to use from which PDB
	"trp_seqs_table"		: working_dir + "00_provided_files/seqs_table.tsv", #Tab-delimited table containing the sequences to download
	"local_proteomes"	: working_dir + "00_local_proteomes/", #Folder containing the proteomes (in FASTA) for JGI and local sequences

	#Locations of the binaries that MUST be adjusted in absolute paths
	"cobalt_bin"		: "/media/data/compiled_software/ncbi-cobalt-2.1.0/bin/cobalt", #NCBI COBALT binary
	"cobalt_db"			: "/media/data/cdd_clique_0.75/cdd_clique_0.75", #NCBI COBALT DB (Downloadable at: ftp.ncbi.nlm.nih.gov/pub/cobalt/db)
	"mafft_bin"			: "/usr/bin/mafft", #MAFFT binary
	"iqtree_bin"		: "/media/data/compiled_software/iqtree-2.1.2-Linux/bin/iqtree2", #IQ-TREE2 binary
	"frtmalign"			: "/media/data/compiled_software/frtmalign/frtmalign",
}

output_paths = {
	#Default directories for results that MAY be adjusted
	"temp"				: working_dir + "temp/",
	"original_struct"	: working_dir + "01_original_structures/",
	"clean_struct_dir"	: working_dir + "02_clean_structures/",
	"clean_chain_dir"	: working_dir + "03_clean_chains/",
	"aligned_chain_dir"	: working_dir + "05_aligned_chains/",
	"msa"				: working_dir + "06_MSA/",
	"sequences"			: working_dir + "07_TRP_sequences/",
	"mafft_alignment"	: working_dir + "09_mafft_alignment/",
	"iqtree_results"	: working_dir + "10_iqtree_results/",
	"savepoints"		: working_dir + "98_savepoints/",
}
