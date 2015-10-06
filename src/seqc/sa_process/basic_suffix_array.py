#Python implementation of suffix array based 100% off http://www.cs.helsinki.fi/u/tpkarkka/publications/jacm05-revised.pdf

def get_I(SA12, t, n0):
	if SA12[t] < n0:
		return SA12[t] * 3 + 1
	else:
		return (SA12[t] - n0) * 3 + 2

def leq_2(a1, a2, b1, b2):
	return (a1 < b1) or (a1 == b1 and a2 <= b2)

def leq_3(a1, a2, a3, b1, b2, b3):
	return (a1 < b1) or (a1 == b1 and leq_2(a2, a3, b2, b3))

def radixPass(a, b, r, n, K):
	counter = [0,] * (K+1)
	i = 0
	for i in range(n):
		counter[r[a[i]]] += 1
	
	the_sum = 0
	for i in range(K + 1):
		t = counter[i]
		counter[i] = the_sum
		the_sum += t

	for i in range(n):
		b[counter[r[a[i]]]] = a[i]
		counter[r[a[i]]] += 1
		#b[counter[r[a[i]]]] += 1
		#b[c[r[a[i]]]

	return b

def suffix_array(T, SA, n, K):
	n0 = (n+2) // 3; n1 = (n+1) // 3; n2 = (n//3); n02 = n0 + n2
	#print n0, "n0 and n2:", n2
	#s12 = [0,] * (n02 + 3)
	s12 = []
	SA12 = [0,] * (n02 + 3)
	#s0 = [0,] * n0
	SA0 = [0,] * n0

	#Step 0
	s12 = [i for i in range(n + n0 - n1) if i % 3 != 0]
	s12.extend([0, 0, 0])

	#Step 1
	#print "1"
	#	s12, SA12, T = 
	radixPass(s12, SA12, T[2:], n02, K)
	#print "2"
	#SA12, s12, T = 
	radixPass(SA12, s12, T[1:], n02, K)
	#print "3"
	#s12, SA12, T = 
	radixPass(s12, SA12, T, n02, K) 
	#I think we actually need to sort this, though... it's done in place in original

	name = 0
	c0, c1, c2 = -1, -1, -1
	for i in range(n02):
		if T[SA12[i]] != c0 or T[SA12[i] + 1] != c1 or T[SA12[i] + 2] != c2:
			name += 1
			#print T, SA12, i
			c0 = T[SA12[i]]
			c1 = T[SA12[i] + 1]
			c2 = T[SA12[i] + 2]
		if SA12[i] % 3 == 1:
			s12[SA12[i] // 3] = name
		else:
			s12[SA12[i] // 3 + n0] = name

	if name < n02:
		suffix_array(s12, SA12, n02, name)
		i = 0
		while i < n02:
		#for i in range(n02):
			s12[SA12[i]] = i + 1
			i += 1
	else:
		i = 0
		while i < n02:
		#for i in range(n02):
			SA12[s12[i] - 1] = i
			i += 1

	#Step 2: sort nonsample suffixes

	s0 = [SA12[i] * 3 for i in range(n02) if SA12[i] < n0]

	radixPass(s0, SA0, T, n0, K)

	p = j = k = 0

	#j = 0
	#i = 0
	#for i in range(n02):
		#if SA12[i] < n0:
			#s0[j] = 3 * SA12[i]
			#j += 1
	#print "4"
	#print SA0, "before"
	#s0, SA0, T = radixPass(s0, SA0, T, n0, K)
	#print SA0, "after"

	#Step 3: Merge
	t = n0 - n1
	p = 0
	#for k in range(n):
	while k < n:
		i = get_I(SA12, t, n0)
		#print SA0, p
		j = SA0[p] if p < n0 else 0

		if SA12[t] < n0:
			litmus_test = leq_2(T[i], s12[SA12[t] + n0], T[j], s12[j//3])
		else:
			litmus_test = leq_3(T[i], T[i+1], s12[SA12[t] - n0 + 1], T[j], T[j+1], s12[(j//3) + n0])
		if litmus_test:
			SA[k] = i
			t += 1
			if t == n02:
				#for p in range(n0):
				k += 1
				while p < n0:
					#for k in [3, 1, 4, 1, 5, 9, 2, 6, 5, 3]: #oh lord, what's that thing...?:
					SA[k] = SA0[p]
					p += 1
					k += 1
		else:
			SA[k] = j
			p += 1
			if p == n0:
				k += 1
				while t < n02:
				#for k in [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7]: #oh lord...)
					SA[k] = get_I(SA12, t, n0)
					t += 1
					k += 1
		k += 1

	return SA

#BANANA$
#look into why we need the [6]... not sure how this will behave will molto more words.

#print [6] + suffix_array([2, 1, 3, 1, 3, 1, 0, 0, 0], [0] * 6, 6, 3)

if __name__ == "__main__":
	import random
	import time
	start = time.time()
	big_list = []
	num_seqs = 1500
	len_seqs = 1200
	n = num_seqs * len_seqs
	for count in range(num_seqs):
		for bp in range(len_seqs):
			big_list.append(int(random.random() * 4) + 1)
	big_list.append(0)
	big_list.append(0)
	big_list.append(0)

	suffix_array(big_list, [0] * n, n, 4)

	print(time.time() - start)
