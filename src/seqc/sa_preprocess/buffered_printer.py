def buffered_printer(seq, count):
	"""Just outputs a fasta-compatible sequence from seq with width count.
	Just here so fasta files get line breaks at regular intervals for human readability."""
	seq = seq.replace("\n", "").replace(" ", "")
	out_seq = ""
	cur_count = 0
	for a in seq:
		if cur_count != 0 and cur_count % count == 0:
			out_seq += "\n" + a
		else:
			out_seq += a
		cur_count += 1
	return out_seq
