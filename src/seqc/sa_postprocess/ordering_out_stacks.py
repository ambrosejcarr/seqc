from seqc.sa_postprocess import load_fasta as lf
import time
import pickle as pickle
from seqc.sa_postprocess import set_classes as sc

from scipy.stats import poisson

def create_lookup(lamb, max_val):
	lookup = {} 
	for count in range(max_val):
		lookup[count] = poisson.pmf(count, lamb)
	return lookup

lk = create_lookup(450, 1001)

def create_spots(jelly_file):
	in_file = open(jelly_file)
	in_read = in_file.readlines()
	in_file.close()

	gene_spots = {}

	print("Reading through file, creating gene spots.")
	line_count = 0
	for line in in_read:
		lsp = line.strip().split("\t")
		genes = [lsp[x] for x in range(0, len(lsp), 2)]
		positions = [int(lsp[x]) for x in range(1, len(lsp), 2)]
		for p, g in zip(positions, genes):
			if g not in gene_spots:
				gene_spots[g] = []
			gene_spots[g].append((p, line_count))
		line_count += 1

	return gene_spots, in_read

def determine_probabilities(gene_spots, in_read):
	print("Beginning to examine gene-level information.")
	gene_start = time.time()
	big_dictionary = {}
	giant_dictionary = {}
	gene_count = 0
	for gene in list(gene_spots.keys())[:]:
		sorted_by_position = sorted(gene_spots[gene])
		multimappings = {}
		big_dictionary[gene] = {}
		giant_dictionary[gene] = {}
		giant_multimappings = {}

		for collision in sorted_by_position:
			position, line_number = collision
			line = in_read[line_number]
			lsp = line.strip().split("\t")

			giant_genes = [lsp[x] for x in range(0, len(lsp), 2)]
			giant_positions = [int(lsp[x]) for x in range(1, len(lsp), 2)]
			for (p, g) in zip(giant_positions, giant_genes):
				if g not in giant_multimappings:
					giant_multimappings[g] = []
				giant_multimappings[g].append(position)

			genes = tuple(sorted([a for a in giant_genes if a != gene]))
			if genes not in multimappings:
				multimappings[genes] = []
			multimappings[genes].append(position)

		for colliding_gene in giant_multimappings:
			prob_to = 0
			for col_pos in giant_multimappings[colliding_gene]:
				prob_to += lk[col_pos]
			giant_dictionary[gene][colliding_gene] = prob_to

		for colliding_gene in multimappings:
			prob_to = 0
			for col_pos in multimappings[colliding_gene]:
				prob_to += lk[col_pos]
			big_dictionary[gene][colliding_gene] = prob_to #(prob_to, len((multimappings[colliding_gene])))
		gene_count += 1
	return big_dictionary, giant_dictionary

def pickle_dictionary(big_dictionary, out_file_name):
	with open(out_file_name, "wb") as p_out:
		pickle.dump(big_dictionary, p_out, protocol=4)
	p_out.close()

def determine_sc_groups(big_dictionary, giant_dictionary, original_seq_file="../preprocess/csf_labeled_fasta.fa", out_path="/data/nums_to_transcripts.txt"):
	all_genes = []
	all_seq_file = open(original_seq_file)
	seq_line = all_seq_file.readline()

	while seq_line:
		if seq_line[0] == ">":
			all_genes.append(seq_line[1:].strip())
		seq_line = all_seq_file.readline()

	all_seq_file.close()

	all_genes = set(all_genes)

	print("Done reading file")

	dj = sc.DisjointSet()
	set_dictionary = {}
	for gene in all_genes:
		new_aset = sc.ASet(gene)
		set_dictionary[gene] = new_aset
		dj.makeSet(new_aset)

	print("joining up transcripts")
	for gene in sorted(giant_dictionary.keys()):
		for pair in giant_dictionary[gene]:
			#if we only want to merge features if they belong to same gene, need the next if. If we want to just reduce features, just do if True:
			if True:
			#if pair.split("_")[1] == gene.split("_")[1]: #then they are from same gene, consider merging.
				their_sum = giant_dictionary[gene][pair] + giant_dictionary[pair][gene]
				if their_sum > 1.99:
					try:
						dj.Union(set_dictionary[gene], set_dictionary[pair])
					except KeyError:
						pass
						#print "Apparently missing one of", gene, pair

	print("Collecting info")
	raf, leaders = dj.collect() #representatives_and_families, leaders
	overlap_probabilities = {}

	out_file = open(out_path, "w")
	rep_num = 1 

	sc_translations = {}
	repres = []

	for rep in raf:
		repres.append(rep)
		for member in raf[rep]:
			out_file.write("SC%i, %s\n" %(rep_num, member))
			sc_translations[member] = rep_num
		rep_num += 1
	out_file.close()
	return sc_translations, repres

def create_sc_spots(jelly_file, sc_translations, repres):
	in_file = open(jelly_file)
	in_read = in_file.readlines()
	in_file.close()
	repres = set(repres)

	gene_spots = {}

	print("Reading through file, creating gene spots.")
	line_count = 0
	for line in in_read:
		lsp = line.strip().split("\t")
		genes = [lsp[x] for x in range(0, len(lsp), 2)]
		positions = [int(lsp[x]) for x in range(1, len(lsp), 2)]
		for p, g in zip(positions, genes):
			if g not in repres:
				continue #we don't care if they are not the representative of their sc_cluster
			else:
				g = sc_translations[g]
			if g not in gene_spots:
				gene_spots[g] = []
			gene_spots[g].append((p, line_count))
		line_count += 1

	return gene_spots, in_read

def determine_probabilities_after_sc(sc_gene_spots, in_read, sc_translations):
	print("Beginning to examine gene-level information.")
	big_dictionary = {}
	gene_count = 0
	for gene in list(sc_gene_spots.keys())[:]:
		sorted_by_position = sorted(sc_gene_spots[gene])

		multimappings = {}
		big_dictionary[gene] = {}

		positions_that_overlap = []

		for collision in sorted_by_position:
			position, line_number = collision
			line = in_read[line_number]
			lsp = line.strip().split("\t")
			genes = tuple(set(sorted([a for a in [sc_translations[lsp[x]] for x in range(0, len(lsp), 2)] if a != gene])))
			if not genes:
				continue
			if genes not in multimappings:
				multimappings[genes] = []
			multimappings[genes].append(position)
			positions_that_overlap.append(position)

		for colliding_gene in multimappings:
			prob_to = 0
			for col_pos in set(multimappings[colliding_gene]):
				prob_to += lk[col_pos]
			big_dictionary[gene][colliding_gene + (gene,)] = min(1.0, prob_to) #(prob_to, len((multimappings[colliding_gene])))

		unique_to_gene = 0
		positions_that_overlap = set(positions_that_overlap)
		for upos in range(0, 951):
			if upos not in positions_that_overlap:
				unique_to_gene += lk[upos]

		big_dictionary[gene][(gene,)] = min(1.0, unique_to_gene) #don't want probabilites over 1 but they can occur due to floating-point stuff.

		gene_count += 1

	for gene in sc_translations:
		tran_name = sc_translations[gene]
		if tran_name not in sc_gene_spots:
			big_dictionary[tran_name] = {}
			big_dictionary[tran_name][(tran_name,)] = 1.0

	return big_dictionary

def standard_run(input_file, nums_to_transcripts_path, labeled_fasta, pickle_destination):
	start_time = time.time()

	gene_spots, in_read = create_spots(input_file)
	print("Done with create_spots()")
	big_dictionary, giant_dictionary = determine_probabilities(gene_spots, in_read)
	print("Done with determine_probabilities()")
	gene_spots = {}
	sc_translations, repres = determine_sc_groups(big_dictionary, giant_dictionary, labeled_fasta, out_path=nums_to_transcripts_path)
	sc_gene_spots, in_read = create_sc_spots(input_file, sc_translations, repres)
	overlap_probabilities = determine_probabilities_after_sc(sc_gene_spots, in_read, sc_translations)

	print("Now just have to pickle...")
	pickle_dictionary(overlap_probabilities, pickle_destination)
	print(time.time() - start_time, "seconds to complete")

if __name__ == "__main__":
	standard_run("jelly_output_stack_50.txt", "nums_to_transcripts_99.txt", "my_labeled_fasta.fa", "p_coalignment.pckl")
