def load_fasta(file_lines, first_word=False):
	"""Reads in a fasta file such that the >-lines of each fasta record
	is stored as a key in the dictionary with value of sequence."""
	cur_seqs = {}
	cur_head = ""
	cur_seq = ""
	for line in file_lines:
		if line[0] == ">":
			if cur_head != "":
				cur_seqs[cur_head] = cur_seq
			if not first_word:
				cur_head = line[1:].strip()
			else:
				cur_head = line[1:].strip().split()[0]
			cur_seq = ""
		else:
			cur_seq += line.strip()

	cur_seqs[cur_head] = cur_seq
	return cur_seqs
