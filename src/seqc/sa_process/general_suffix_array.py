#This is the main program, combining work from basic_suffix_array.py and lcp.py into a single object. For sc_seq, shouldn't have to run anything else.

from scseq.sa_process.lcp import get_lcp
import scseq.sa_process.basic_suffix_array as msaw
#import child_tab as ct
import scseq.sa_process.load_fasta as lf
#import new_finish as nf
import random
import time
import subprocess

def binary_search(l, r, x, target):
	"""Goal is to find the number in the list that is closest to the query but less than or equal to it.
	Will return -1 if no such number exists in the list (query < min(list))."""
	if r < l:
		return -1
	else:
		mid = (l + r) // 2
		if x[mid] < target:
			#in production mode, pass len(x) as parameter, dumb to check each time.
			if (mid < (len(x) - 1) and x[mid + 1] > target) or mid == len(x) - 1:
				return mid
			else:
				return binary_search(mid + 1, r, x, target)
		elif x[mid] > target:
			return binary_search(l, mid - 1, x, target)
		else:
			return mid

class ConvenienceCounter:
	"""A silly class to let me keep track of sums easily."""
	def __init__(self):
		self.x = 0
	def __getitem__(self, num):
		self.x += num
		return self.x

class GeneralSuffixArray:
	"""Class that will be used to build suffix array. Can further have lcp array and child-tab array added to it for performance enhancement."""
	def __init__(self, file_name="/home/tb/dpeer_work/cterm_mappability/transcripts_with_unique_ends.fa", test_limited=False):
		"""__init__(file_name) ==> transforms the inputted fasta file into numerical array that makes sense to the suffix array constructor. Object construction does not initiate suffix array construction. That requires calling construct_SA()."""
		in_file = open(file_name)
		in_read = in_file.readlines()
		in_file.close()

		self.sa = None
		self.transput = None
		
		genes = lf.load_fasta(in_read)
		gene_keys = sorted(genes.keys()) #try to make it give same values each time...
		random.shuffle(gene_keys)
		self.words = []
		starter_code = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 5}
		if test_limited:
			#If we're doing a limited number of genes for test purposes...
			self.gene_names = gene_keys[:2000]
		else:
			self.gene_names = gene_keys[:]
		for g in self.gene_names:
			gene_g = genes[g]
			#if len(gene_g) > 1000:
				#gene_g = gene_g[(len(gene_g) - 1000):]
			self.words.append(gene_g)
		self.transform_input(starter_code)
		self.total_len = sum(len(word) for word in self.words) + len(self.words)
		print(self.total_len)
		print("Done constructing object. Now you can explicitly build the SA and the LCP array and such.")
		
	def transform_input(self, code=None):
		"""Takes a list of words, translates the nucleotides to integers,
       adds negative numbers between words, returns original words as new
       flat list."""
		num_words = len(self.words)
		starter = 1
		if code is None:
			letters = sorted(set("".join(self.words)))
			code = {letters[k]:(k + len(self.words) + 1) for k in range(len(letters))}
		else:
			starter = max(code.values()) + 1
		cc = ConvenienceCounter()
		gene_ends = [cc[len(word) + 1] for word in self.words] #by virtue of custom ConvenienceCounter class, this keeps track of ends... they're actually 1-based ends, though from what I can tell.
		gene_ends = [x - 1 for x in gene_ends] #to make it zero-based which is probably preferable.
		translated = [item for sublist in [[code[i] for i in self.words[x]] + [(starter + x)] for x in range(num_words)] for item in sublist]
		#print translated
		semi_translated = "".join([item for sublist in [[i for i in x] + ["/"] for x in self.words] for item in sublist])

		#for ge in gene_ends:
			#print semi_translated[ge]

		self.gene_ends = gene_ends
		self.semi_translated = semi_translated
		self.transput = translated
		self.final_code = code

	def construct_SA(self):
		"Builds suffix array. Takes no arguments, only needs data already constructed based on inputted fasta."""
		print("Beginning construction of suffix array.")
		start = time.time()
		self.sa = msaw.suffix_array(self.transput + [0, 0, 0], [0] * self.total_len, self.total_len, len(self.final_code) + len(self.words) + 1)
		print(time.time() - start)
		print("Done constructing suffix array.")
	
	def construct_lcp(self):
		"""Builds lcp array. Also takes no arguments, using only suffix array and constructed input."""
		print("Beginning construction of LCP.")
		start = time.time()
		if self.sa is not None and self.transput is not None:
			self.lcp_tab = get_lcp(self.sa, self.transput)
			print(time.time() - start)
			print("Done construction of LCP.")
		else:
			print("Can't construct LCP array without initializing and constructing SA.")

	def enhance_SA(self):
		"""Adds the child tab and interval tree stuff that transforms the SA into an enhanced SA. Good for searching."""
		#Going to be up, down, next_index
		print("Now adding child tab to suffix array.")
		self.cld_tab = [[None, None, None] for x in range(len(self.lcp_tab))]
		#self.whole is really self.transput
		start = time.time()
		self.build_esa()
		print(time.time() - start)
		print("Suffix array is now an enhanced suffix array.")

	def build_esa(self):
		"""Constructs esa by making up/down table, and then next_index table."""
		self.esa_up_down()
		self.esa_next_index()
		
	def esa_up_down(self):
		"""Constructs up/down table of enhanced suffix array."""
		last_index = -1
		stack = list(range(len(self.lcp_tab)))
		stack.append(0)
		n = len(self.lcp_tab)
		for i in range(1, n):
			while self.lcp_tab[i] < self.lcp_tab[stack[-1]]:
				last_index = stack.pop()
				if self.lcp_tab[i] <= self.lcp_tab[stack[-1]] and self.lcp_tab[stack[-1]] != self.lcp_tab[last_index]:
					self.cld_tab[stack[-1]][1] = last_index
			if self.lcp_tab[i] >= self.lcp_tab[stack[-1]]:
				if last_index != -1:
					self.cld_tab[i][0] = last_index
					last_index = -1
				stack.append(i)

	def esa_next_index(self):
		"""Constructs next_index table of enhanced suffix array."""
		stack = list(range(len(self.lcp_tab)))
		stack.append(0)
		n = len(self.lcp_tab)
		for i in range(1, n):
			while self.lcp_tab[i] < self.lcp_tab[stack[-1]]:
				stack.pop()
			if self.lcp_tab[i] == self.lcp_tab[stack[-1]]:
				last_index = stack.pop()
				self.cld_tab[last_index][2] = i
			stack.append(i)

	def esa_first_pass(self, i, j, char):
		"""Is used for the search, the first pass of the search is slightly different than subsequent searches. Needs to find sub-intervals beginning with letter char."""
		interval_list = []
		num_low = 0
		nex = self.cld_tab[num_low][2] #.next_index
		if self.transput[self.sa[num_low]] == char:
			return (num_low, nex - 1) 
		#interval_list.append((0, first_interval - 1))
		while nex is not None:
			interval_list.append((num_low, nex - 1))
			num_low = nex
			nex = self.cld_tab[num_low][2] #.next_index
			if self.transput[self.sa[num_low]] == char:
				if nex is not None:
					return (num_low, nex - 1)
				else:
					#Not 100% sure why appropriate to add 1...
					#It may not be appropriate.
					return (num_low + 1, j-1)
		return None

	def esa_get_interval(self, i, j, char):
		"""Used for search to get intervals matching char."""
		min_i = -1
		max_j = -1
		the_l = self.esa_get_lcp(i, j)
		if i < self.cld_tab[j+1][0] <= j:
			i_lower = self.cld_tab[j+1][0]
		else:
			i_lower = self.cld_tab[i][1]

		if self.transput[self.sa[i] + the_l] == char:
			return i, i_lower - 1

		while self.cld_tab[i_lower][2] is not None:
			i_next = self.cld_tab[i_lower][2]
			if self.transput[self.sa[i_lower] + the_l] == char:
				return i_lower, i_next - 1
			i_lower = i_next
		
		if self.transput[self.sa[i_lower] + the_l] == char:
			return i_lower, j

		if min_i == -1 and max_j == -1:
			return None, None #I don't know what to return if not found...
		else:
			return min_i, max_j

	def esa_get_lcp(self, i, j):
		"""Gets lcp for interval-tree defined by i, j."""
		if i < self.cld_tab[j+1][0] <= j:
			return self.lcp_tab[self.cld_tab[j+1][0]]
		else:
			return self.lcp_tab[self.cld_tab[i][1]]

	def esa_find(self, pattern):
		"""Wraps functionality of search in enhanced suffix array. Input pattern is the substring to search for."""
		c = 0
		query_found = True
		len_pat = len(pattern)
		n = len(self.lcp_tab) - 1
		result = self.esa_first_pass(0, n, pattern[c])
		#print "first timer:", result
		if result is not None:
			i, j = result
		else:
			i = None; j = None
		while i is not None and j is not None and c < len_pat and query_found:
			if i != j:
				low = self.esa_get_lcp(i, j) #CAT, i,j = 6, 7, low is 2.
				mini = min(low, len_pat)
				query_found = self.transput[self.sa[i] + c: self.sa[i] + mini - 1] == pattern[c:mini - 1]
				#print "query_found", query_found, mini
				
				c = mini
				if c < len_pat:
					i, j = self.esa_get_interval(i, j, pattern[c])
				#print i, j, "new get interval"
				#for index in range(i, j+1):
					#match = self.whole[self.sa[index]:self.sa[index] + len_pat]
					#print match, match == pattern				
			else:
				query_found = self.transput[self.sa[i] + c: self.sa[i] + len_pat - 1] == pattern[c:len_pat - 1] #it would be helpful to understand how to make this else condition work properly
				#print query_found
				#print self.whole[self.sa[i] + c: self.sa[i] + len_pat - 1]
				#print pattern[c:mini]
				#print pattern[c:len_pat - 1]
				#print "QUERY FOUND", query_found
				break

		if i is not None:
			#print i, j
			matching_indices = []
			#It's fine if I have to check the first one...
			#It would be nice if I could just check the first one and if that's correct,
			#then I output all. If I have to check each one then this code is O(m + w)
			#where w > z is some new number of checks to make independent of how many times
			#z is in the text.
			attempted = 0
			failed = 0
			for index in range(i, j+1):
				match = self.transput[self.sa[index]:self.sa[index] + len_pat]
				attempted += 1
				if match == pattern:
					matching_indices.append(self.sa[index])
				else:
					failed += 1
					#print index, match
					#print "In the industry, we call this a truly unfortunate occurrence."
			if failed != 0:
				print(failed, attempted, "failed, attempted")
			
			return sorted(matching_indices)
		else:
			#print "pattern", pattern, "not found"
			return []

	def deal_with_patterns(self, min_k=50, output_dir="/home/tb/misc_output/", jelly_bin="/home/tb/dpeer_work/jellyfish-2.2.0/bin/jellyfish", jelly_processes=16, jelly_s_str="400M"):
		end_list = self.gene_ends
		gene_list = self.gene_names
		end_len = len(end_list)
		all_patterns = []
		k = min_k
		sorted_indices = sorted(list(range(len(self.lcp_tab))), key=lambda indx: -self.lcp_tab[indx])
		patterns = []

		print("Adding patterns")

		tmp_out = "" #open(output_dir + "/tmp_seq.fa", "w")
		p_count = 0
		for biggest_overlap in sorted_indices:
			if self.lcp_tab[biggest_overlap] < k:
				break
			pattern = tuple(self.transput[self.sa[biggest_overlap]:self.sa[biggest_overlap] + self.lcp_tab[biggest_overlap]])
			str_pat = "".join(['?ACGTN'[x] for x in pattern]) #str_pat = 1/2 memory of the native list
			tmp_out += ">seq" + str(p_count) + "\n" + str_pat.strip() + "\n"
			p_count += 1

		print(p_count, "patterns to investigate")

		sorted_indices = []

		all_patterns = []

		out_text = ""

		command = jelly_bin + " count -m " + str(k) + " -o /dev/stdout /dev/stdin -t " + str(jelly_processes) + " -s " + jelly_s_str
		jelly_process = subprocess.Popen(command.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE)
		print("Sending jelly process!")
		out, err = jelly_process.communicate(bytes(tmp_out, 'utf-8'))
	
		if err is not None:
			print(err)
			raise RuntimeError("Something went wrong while calling jellyfish count.")

		tmp_out = ""

		print("Dumping...")
		command = jelly_bin + " dump -c -t /dev/stdin -o /dev/stdout"
		dump_process = subprocess.Popen(command.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE)

		print("About to communicate...")

		output, error = dump_process.communicate(out)
		output = output.decode('utf8')

		if error is not None:
			print(error)
			raise RuntimeError("Something went wrong while calling jellyfish dump.")

		out = ""

		print("Done with subprocesses...")	

		for kmer_read in output.split("\n"):
			if kmer_read.strip(): #not just a blank line...
				all_patterns.append(kmer_read.split("\t")[0])
		     
		output = "" 
		#wondering how much storage space each of these variables is taking...
		#there's a tradeoff between writing these guys to memory or to a file and I choose the memory route for now. Would be nice if I could stream through it a little more vs. keeping whole thing in memory, though.

		print(len(all_patterns), "patterns to search for...")
		#print "Need to sort patterns."

		all_patterns.sort() #python will sort by the values in the list, without rearranging each list

		cur_prefix = ""
		current_collection = []
		match_lines = []
		new_out = open(output_dir + "/jelly_output_stack_" + str(min_k) + ".txt", 'w')
		pattern_count = 0
		start_time = time.time()
		transl = {'?': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 5}
		first_end = end_list[0]
		for pattern in all_patterns:
			pattern = [transl[px] for px in pattern]
			first_letters = pattern[:k]
			if pattern_count % 1000000 == 0:
				print(pattern_count, time.time() - start_time)
			if first_letters != cur_prefix:

				#process current collection of patterns
				resulting_indices = self.esa_find(pattern)
				match_line = ""
				for matching_index in resulting_indices:
					if matching_index < first_end: #a bit of a special case since we don't have a gene index higher than this...
						gene_start = matching_index + 1
						gene_index = -1 #which is fine, we add one.
					else:
						gene_index = binary_search(0, end_len, end_list, matching_index)
						gene_start = matching_index - end_list[gene_index]
					matching_gene = gene_list[gene_index + 1]
					match_line += matching_gene + "\t" + str(gene_start) + "\t"
				new_out.write(match_line.strip() + "\n")
				#match_lines.append(match_line.strip())
		                     
				cur_prefix = first_letters
				current_collection = []
			current_collection.append(pattern)
			pattern_count += 1

		new_out.close()
		return True


def standard_run(fasta_input="sample_fasta.fa", output_dir="/home/tb/misc_output/", jellyfish_exe="/home/tb/dpeer_work/jellyfish-2.2.0/bin/jellyfish", k_vals=[60], jelly_num_process=16, jelly_s_str="400M", is_test_limited=False):
	start = time.time()
	gsa = GeneralSuffixArray(fasta_input, test_limited=is_test_limited)
	gsa.construct_SA()
	gsa.construct_lcp()
	gsa.enhance_SA()
	print("Finishing it off.")
	for k in k_vals:
		print("Querying with", k)
		gsa.deal_with_patterns(k, output_dir, jellyfish_exe, jelly_num_process, jelly_s_str)
	print("Finished after", time.time() - start, "seconds")

if __name__ == "__main__":
	standard_run(fasta_input="../sa_preprocess/csf_labeled.fa", output_dir="/mnt/jelly/", jellyfish_exe="/root/home/jellyfish-2.2.0/bin/jellyfish", k_vals=[60])
	print("Done all work.")
