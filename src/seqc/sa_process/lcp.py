#Algorithm is taken from Kasai...

def get_lcp(suffix_array, text):
	rank = [0,] * len(suffix_array)
	#not sure if first lcp thing should be 0 or None...
	height = [0] + [0,] * (len(suffix_array) - 1)

	for i in range(len(text)):
		rank[suffix_array[i]] = i

	h = 0
	for i in range(len(text)):
		if rank[i] > 0:
			k = suffix_array[rank[i] - 1]
			while text[i+h] == text[k+h]:
				h += 1
			height[rank[i]] = h
			if h > 0:
				h = h - 1
	return height[:]

if __name__ == "__main__":
	pos = [6, 5, 3, 1, 0, 4, 2]
	a = ['B', 'A', 'N', 'A', 'N', 'A', -1]
	#pos = [3, 7, 11, 9, 1, 5, 8, 0, 10, 2, 6, 4]
	#a = ['C', 'A', 'N', -3, 'P', 'A', 'N', -2, 'C', 'A', 'D', -1]
	#a = ['S', 'H', 'A', 'N', 'D', 'Y', -3, 'R', 'A', 'N', 'D', 'Y', -2, 'M', 'A', 'N', -1]
	#pos = [6, 12, 16, 14, 2, 8, 4, 10, 1, 13, 15, 3, 9, 7, 0, 5, 11]
	
	lcp = get_lcp(pos, a)
	print(lcp)
