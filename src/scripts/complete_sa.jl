##All code written by Tom Blazejewski, 2015.

### CODE originally from julia_sa.jl, copied and pasted for convenience.
function get_I(SA12, t, n0)
	if SA12[t + 1] < n0
		return SA12[t + 1] * 3 + 1
	else
		return (SA12[t + 1] - n0) * 3 + 2
	end
end

function leq_2(a1::Int, a2::Int, b1::Int, b2::Int)
	return (a1 < b1) || (a1 == b1 && a2 <= b2)
end

function leq_3(a1::Int, a2::Int, a3::Int, b1::Int, b2::Int, b3::Int)
	return (a1 < b1) || (a1 == b1 && leq_2(a2, a3, b2, b3))
end

function zeroPass(r::Array{Int, 1}, a::Array{Int, 1}, K::Int, n::Int)
	counter = zeros(Int, K+1)
	for i=1:n
		counter[r[a[i] + 1] + 1] += 1
	end
	return counter
end

function firstPass(K::Int, counter::Array{Int,1})
	the_sum::Int = 0
	t = 0
	for i in 1:(K+1)
		t = counter[i]
		counter[i] = the_sum
		the_sum += t
	end
	return counter
end

function secondPass(b::Array{Int, 1}, r::Array{Int, 1}, counter::Array{Int, 1}, a::Array{Int, 1}, n::Int)
	for i in 1:n
		b[counter[r[a[i] + 1] + 1] + 1] = a[i]
		counter[r[a[i] + 1] + 1] += 1
	end
	return b
end

function radixPass(a::Array{Int,1}, b::Array{Int,1}, r::Array{Int,1}, n::Int, K::Int)
	counter = zeroPass(r, a, K, n) #zeros(Int, K+1)

	counter = firstPass(K, counter)
	
	b = secondPass(b, r, counter, a, n)
	return b
end

function make_S12(n::Int, n0::Int, n1::Int)
	s12::Array{Int, 1} = filter(x -> x % 3 != 0, [i for i in 0:(n + n0 - n1 - 1)])
	append!(s12, zeros(Int, 3))
	return s12
end

function suffix_array(T::Array{Int,1}, SA::Array{Int,1}, n::Int, K::Int)
	n0::Int = div((n + 2), 3)
	n1::Int = div((n + 1) , 3)
	n2::Int = div(n, 3)
	n02::Int = n0 + n2

	s12 = make_S12(n, n0, n1)

	SA12 = zeros(Int, (n02 + 3))
	SA0 = zeros(Int, n0)

	radixPass(s12, SA12, T[3:end], n02, K)
	radixPass(SA12, s12, T[2:end], n02, K)
	radixPass(s12, SA12, T, n02, K)
	
	name::Int = 0
	c0::Int = -1
	c1::Int = -1
	c2::Int = -1
	i::Int = 0

	for i in 1:n02
		if T[SA12[i] + 1] != c0 || T[SA12[i] + 2] != c1 || T[SA12[i] + 3] != c2
			name += 1
			c0 = T[SA12[i] + 1]
			c1 = T[SA12[i] + 2]
			c2 = T[SA12[i] + 3]
		end
		if SA12[i] % 3 == 1
			s12[div(SA12[i], 3) + 1] = name
		else
			s12[div(SA12[i], 3) + n0 + 1] = name
		end
	end

	if name < n02
		suffix_array(s12, SA12, n02, name)
		i = 0
		while i < n02
			s12[SA12[i + 1] + 1] = i + 1
			i += 1
		end
	else
		i = 0
		while i < n02
			SA12[s12[i + 1]] = i
			i += 1
		end
	end

	#s0 = filter(ex -> ex < n0, [ind for ind in SA12]) * 3
	s0 = Int[]
	for ind in 1:n02
		if SA12[ind] < n0
			push!(s0, SA12[ind] * 3)
		end
	end

	radixPass(s0, SA0, T, n0, K)
	p::Int = 0
	j::Int = 0
	k::Int = 0

	t::Int = n0 - n1
	p = 0

	while k < n
		i = get_I(SA12, t, n0)
		if p < n0
			j = SA0[p + 1]
		else
			j = 0
		end

		#j += 1

		litmus_test::Bool

		if SA12[t + 1] < n0
			litmus_test = leq_2(T[i + 1], s12[SA12[t + 1] + n0 + 1], T[j + 1], s12[div(j + 1, 3) + 1])
		else
			#println("T: ", T, "\ni: ", i, "\nSA12: ", SA12, "\ns12: ", s12)
			litmus_test = leq_3(T[i + 1], T[i + 2], s12[SA12[t + 1] - n0 + 2], T[j + 1], T[j + 2], s12[div(j + 1,3) + 1 + n0])
		end

		if litmus_test
			SA[k + 1] = i
			t += 1
			if t == n02
				k += 1
				while p < n0
					SA[k + 1] = SA0[p + 1]
					p += 1
					k += 1
				end
			end

		else
			SA[k + 1] = j
			p += 1
			if p == n0
				k += 1
				while t < n02
					SA[k + 1] = get_I(SA12, t, n0)
					t += 1
					k += 1
				end
			end

		end
		k += 1

	end
	return SA
end

###END of julia_sa.jl

function read_fasta(file_name::String)
        in_file = open(file_name)
        in_read = readlines(in_file)
        close(in_file)

        seqs = Dict{String, String}()
        seq_head = ""
        seq = ""
        for line in in_read
                if '>' in line
                        if seq_head != ""
                                seqs[seq_head] = seq
                        end
                        seq = ""
                        seq_head = strip(line[2:end])
                else
                        seq *= strip(line)
                end
        end
        seqs[seq_head] = seq
        return seqs
end

#require("julia_sa.jl")

type GSA
	SA::Array{Int, 1}
	LCP::Array{Int, 1}
	transput::Array{Int,1}
	gene_names::Array{String}
end

function get_lcp(suffix_array::Array{Int,1}, text::Array{Int, 1})
	rank = zeros(Int, length(suffix_array))
	height = zeros(Int, length(suffix_array))
	for i in 1:length(text)
		rank[suffix_array[i]] = i
	end
	
	h = 0
	for i in 1:length(text)
		if rank[i] > 1
			k = suffix_array[rank[i] - 1]
			while text[i+h] == text[k+h]
				h += 1
			end
			height[rank[i]] = h
			if h > 0
				h = h - 1
			end
		end
	end
	return height	
end

function binary_search(l::Array{Int,1}, value::Int)
	#not sure if should be -1 on result... this is giving the index <= value.
	low = 1
	mid = 0
	high = length(l)
	while low <= high
		mid = int((low+high)/2)
		if l[mid] < value
			if (mid < (length(l) - 1) && l[mid+1] > value) || mid == length(l) - 1
				return mid
			else
				low = mid+1
			end
		elseif l[mid] > value
			high = mid-1
		else
			return mid
		end
	end
	return 0
end

function up_down(gsa::GSA, cld_tab::Array{Int,2})
  last_index = -1
  stack = collect(1:length(gsa.LCP))
  push!(stack, 1)
  n = length(gsa.LCP)
  for i in 1:n
    while gsa.LCP[i] < gsa.LCP[stack[end]]
      last_index = pop!(stack)
      if gsa.LCP[i] <= gsa.LCP[stack[end]] && gsa.LCP[stack[end]] != gsa.LCP[last_index]
        cld_tab[stack[end], 2] = last_index
      end
    end
    if gsa.LCP[i] >= gsa.LCP[stack[end]]
      if last_index != -1
        cld_tab[i, 1] = last_index
        last_index = -1
      end
      push!(stack, i)
    end
  end
	return cld_tab
end

function next_index(cld_tab::Array{Int,2}, gsa::GSA)
  stack = collect(1:length(gsa.LCP))
  push!(stack, 1)
  n = length(gsa.LCP)
  for i in 1:n
    while gsa.LCP[i] < gsa.LCP[stack[end]]
      pop!(stack)
    end
    if gsa.LCP[i] == gsa.LCP[stack[end]]
      last_index = pop!(stack)
      cld_tab[last_index, 3] = i
    end
    push!(stack, i)
  end
	return cld_tab
end

function esa_first_pass(gsa::GSA, cld_tab::Array{Int,2}, i::Int, j::Int, num_char::Int)
	interval_list = Tuple[]
	num_low = 1
	nex = cld_tab[num_low, 3]
	if gsa.transput[gsa.SA[num_low] + 1] == num_char
		return (num_low, nex - 1)
	end

	while nex != 0
		push!(interval_list, (num_low, nex-1))
		num_low = nex
		nex = cld_tab[num_low, 3]
		if gsa.transput[gsa.SA[num_low] + 1] == num_char
			if nex != 0
				return (num_low, nex - 1)
			else
				return (num_low + 1, j - 1)
			end
		end
	end
	return (0, 0)
end

function esa_get_lcp(gsa::GSA, cld_tab::Array{Int,2}, i::Int, j::Int)
	if i < cld_tab[j+1, 1] <= j
		return gsa.LCP[cld_tab[j+1, 1]]
	else
		return gsa.LCP[cld_tab[i, 2]]
	end
end

function esa_get_interval(gsa::GSA, cld_tab::Array{Int,2}, i::Int, j::Int, num_char::Int)
	min_i = -1
	max_j = -1
	the_l = esa_get_lcp(gsa, cld_tab, i, j) #this is optimistic/knowably non-functional
	if i < cld_tab[j+1, 1] <= j
		i_lower = cld_tab[j+1, 1]
	else
		i_lower = cld_tab[i, 2]
	end

	if gsa.transput[gsa.SA[i] + the_l + 1] == num_char
		return i, i_lower - 1
	end

	while cld_tab[i_lower, 3] != 0
		i_next = cld_tab[i_lower, 3]
		if gsa.transput[gsa.SA[i_lower] + the_l + 1] == num_char
			return i_lower, i_next - 1
		end
		i_lower = i_next
	end

	if gsa.transput[gsa.SA[i_lower] + the_l + 1] == num_char
		return i_lower, j
	end

	if min_i == -1 && max_j == -1
		return (0, 0)
	else
		return min_i, max_j
	end
end

#don't really expect the above to work but am just prototyping this all out for a bit.

function make_lines(gsa::GSA, gene_ends::Array{Int, 1}, parent::Int, k::Int, curly_l::Int, i::Int, j::Int, gene_names::Array{String, 1})
	first_end = gene_ends[1]
	if parent < k <= curly_l
		matching_line = String[]
		m_i = 0
		for index in i:j
			m_i = gsa.SA[index] + 1 #I think +1...
			#println(m_i, " ", first_end, " ", gene_ends[end])
			if m_i < first_end	
				gene_start = m_i + 1
				gene_index = 0
			else
				gene_index = binary_search(gene_ends, m_i)
				gene_start = m_i - gene_ends[gene_index]
			end
			push!(matching_line, gene_names[gene_index + 1])
			push!(matching_line, string(gene_start))
		end
		return join(matching_line, "\t") * "\n" #hoping this is a little more efficient...
	else
		return ""
	end
end
	#out_strs[k] *= strip(matching_line) * "\n"
#else
	#break

function inner_deal(gsa::GSA, cld_tab::Array{Int,2}, k_value::Int, i::Int, j::Int, gene_ends::Array{Int, 1}, gene_names::Array{String, 1})
	to_visit = Tuple[]
	parent_lcp = Int[]
	matching_line = ""
	lines = String[]
	if i > 0 && j > 0
		push!(to_visit, (i, j))
		push!(parent_lcp, 0)
		while length(to_visit) > 0
			next, parent = pop!(to_visit), pop!(parent_lcp)
			i, j = next
			if (i < 1 || j < 1) || (i == j)
				continue
			end
			curly_l = esa_get_lcp(gsa, cld_tab, i, j)
			matching_line = make_lines(gsa, gene_ends, parent, k_value, curly_l, i, j, gene_names)
			#print(matching_line)
			push!(lines, matching_line)
			for next_char in [1, 2, 3, 4]
				another_one = esa_get_interval(gsa, cld_tab, i, j, next_char)
				#println("Want to add ", another_one)
				#if another_one != (0, 0)
				push!(to_visit, another_one)
				push!(parent_lcp, curly_l)
				#end
			end
		end
	end
	return lines
end

function deal_with_patterns(gsa::GSA, cld_tab::Array{Int,2}, gene_ends::Array{Int, 1}, k_value::Int, gene_names::Array{String, 1}, out_file)
	sub_intervals = {0 => (0, 0)} #don't know how to define without dummy values?
	n = length(gsa.LCP)
	first_end = gene_ends[1]
	curly_l = 0
	matching_indices = Int[]
	lines = String[]
	for num_char in [1, 2, 3, 4]
		result = esa_first_pass(gsa, cld_tab, 1, n, num_char)
		sub_intervals[num_char] = result
	end
	for num_char in keys(sub_intervals)
		if num_char == 0
			continue
		end
		i, j = sub_intervals[num_char]
		write(out_file, join(inner_deal(gsa, cld_tab, k_value, i, j, gene_ends, gene_names)))
	end
	return lines
end

function give_lens(x::GSA)
	println(length(x.SA))
	println(length(x.LCP))
end

function main(fasta_in::String, text_out::String, k_value::Int)
	out_file = open(text_out * "/jelly_output_stack_" * string(k_value) * ".txt", "w")
	genes = read_fasta(fasta_in)

	starter_code = {'A' => 1, 'C' => 2, 'G' => 3, 'T' => 4, 'N' => 5}

	gene_ends = Int[]
	total_sum = 0
	for g in genes
		total_sum += length(g[2]) + 1
		push!(gene_ends, total_sum)
	end

	g_count = 0
	starter = 6 #starter is not 0 or 1 because we want to encode "special characters" to separate words, and first five letters taken for A, C, T, G, N.
	translated = Int[]
	semi_translated = String[]
	total_len = 0
	gene_names = String[]

	#gene_num_file = open(text_out * "/nums_to_transcripts.txt", "w")

	g_count = 1
	for g in genes
		push!(gene_names, string(g[1]))
		#println(gene_names)
		#write(gene_num_file, string(g_count) * ", " * string(g[1]) * "\n")
		for c in genes[g[1]]
			push!(translated, starter_code[c])
		end
		push!(translated, starter)
		#println(typeof(genes[g[1]]))
		#println(genes[g[1]])
		#semi_translated *= genes[g[1]]
		push!(semi_translated, genes[g[1]] * "/")
		total_len += length(genes[g[1]]) + 1
		starter += 1
		g_count += 1
	end

	#close(gene_num_file)

	string_semi = ""
	string_semi = join(semi_translated)

	final = append!(translated, [0, 0, 0])
	#println("Starting to make sa.")
	#println(total_len)
	sa = suffix_array(final, zeros(Int, total_len), total_len, 5 + length(genes) + 1)
	lcp = get_lcp(sa + 1, translated[1:end-3])

	a_gsa = GSA(sa, lcp, translated, gene_names)

	cld_tab = up_down(a_gsa, zeros(Int, (total_len, 3)))
	cld_tab = next_index(cld_tab, a_gsa)

	deal_with_patterns(a_gsa, cld_tab, gene_ends, k_value, gene_names, out_file)

	#println(cld_tab .- 1)
	close(out_file)
	#println("Done making sa.")
end

#println(ARGS)
if length(ARGS) == 3
	#println(ARGS)
	#Args will be fasta_file, out_directory, and k-value (int passed as str).
	main(ARGS[1], ARGS[2], int(ARGS[3]))
end
