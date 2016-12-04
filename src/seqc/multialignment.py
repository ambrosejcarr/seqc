import numpy as np
import itertools
import time
import seqc


class UnionFind:
    """Union-find data structure.

    Each unionFind instance X maintains a family of disjoint sets of
    hashable objects, supporting the following two methods:

    - X[item] returns a name for the set containing the given item.
      Each set is named by an arbitrarily-chosen one of its members; as
      long as the set remains unchanged it will keep the same name. If
      the item is not yet part of a set in X, a new singleton set is
      created for it.

    - X.union(item1, item2, ...) merges the sets containing each item
      into a single larger set.  If any item is not yet part of a set
      in X, it is added to X as one of the members of the merged set.
    """

    def __init__(self):
        """Create a new empty union-find structure."""
        self.weights = {}
        self.parents = {}

    def __getitem__(self, obj):
        """Find and return the name of the set containing the object."""

        # check for previously unknown object
        if obj not in self.parents:
            self.parents[obj] = obj
            self.weights[obj] = 1
            return obj

        # find path of objects leading to the root
        path = [obj]
        root = self.parents[obj]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]

        # compress the path and return
        for ancestor in path:
            self.parents[ancestor] = root
        return root

    def __iter__(self):
        """Iterate through all items ever found or unioned by this structure."""
        return iter(self.parents)

    def union(self, *objects):
        """Find the sets containing the objects and merge them all."""
        roots = [self[x] for x in objects]
        heaviest = max([(self.weights[r], r) for r in roots])[1]
        for r in roots:
            if r != heaviest:
                self.weights[heaviest] += self.weights[r]
                self.parents[r] = heaviest

    def union_all(self, iterable):
        for i in iterable:
            self.union(*i)

    def find_all(self, vals):
        vals = [self.find_component(v) for v in vals]
        unique = set(vals)
        reindex = dict(zip(unique, range(len(unique))))
        set_membership = np.array([reindex[v] for v in vals])
        sets = np.array(list(reindex.values()))
        return set_membership, sets

    def find_component(self, iterable):
        """Return the set that obj belongs to

        If the iterable contains items that have been unioned, then any entry in the
         iterable will be sufficient to identify the set that obj belongs to. Use the
         first entry, and return the set associated with iterable.

        If the iterable has not been entered into the structure, this method can yield
         incorrect results
        """
        return self[next(iter(iterable))]

def intersection(set_l):
    res = set_l[0]
    for s in set_l:
        res = set(set(res) & set(s))
    return res


# # Some constants
# NO_DISAMBIGUATION = 0
# RESOLVED_GENE = 1
# NO_GENE_RESOLVED = 2
# MULTIPLE_MODELS = 3


# def reduce_coalignment_array(arr, threshold = 0.0001):
#     res = {}
#     for g in arr:
#         temp = {}
#         for k in arr[g]:
#             if arr[g][k] < threshold:
#                 continue
#             temp[tuple(sorted(k))] = arr[g][k]
#         if len(temp)>0:
#             res[g] = temp
#     return res
    
# #def strip(genes):
# #    return tuple(sorted([int(g[2:]) for g in genes]))
# def strip(genes):
#     return tuple(sorted(genes))
# def strip_model(mod):
#     res = {}
#     for k in mod:
#         res[tuple(sorted(k))]=mod[k]
#     return res

# def split_to_disjoint(obs):
#     res = []
#     uf = UnionFind()
#     uf.union_all(obs.keys())
#     set_membership, sets = uf.find_all(obs.keys())
    
#     for s in sets:
#         d = {}
#         for k in np.array(list(obs.keys()))[set_membership == s]:
#             d[tuple(k)] = obs[tuple(k)]
#         res.append(d)
#     return res

# def get_indices(inds, obs_subset):
#     res = []
#     for genes in obs_subset:
#         res += inds[genes]
#     return res

# def model_to_gene(model):
#     for g in model:
#         if model[g]==1:
#             return g
            
    
# def get_combinations(l):
#     res = []
#     for i in range(len(l)):
#         res += itertools.combinations(l,i+1)
#     return res
            
# # rank the different possible models by their scores
# def best_fit_model(obs_s, coalignment_mat):
#     #obs_s = strip_model(obs)
#     gene_l = single_gene_list(obs_s)  # From the list of observation create a list of unique single genes from which different models can be inferred
      
    
#     if len(obs_s) == 1:
#         if len(list(obs_s.keys())[0]) == 1:
#             return [{gene_l[0]:1}], NO_DISAMBIGUATION
    
#     possible_genes = intersection(list(obs_s.keys()))
    
#     #There is one gene that resolve the disambiguation
#     if len(possible_genes) == 1:
#         model = {}
#         for g in gene_l:
#             model[g] = 0
#         model[list(possible_genes)[0]] = 1
#         return [model], RESOLVED_GENE
    
#     #There is more than one gene that can explain it, no model can be decided
#     if len(possible_genes) > 1:
#         return [], NO_GENE_RESOLVED

#     #There are multiple competing models. For now we don't decide bewteen them
#     return [], MULTIPLE_MODELS
# #    mod_score_list = []     
# #    for mod in get_combinations(gene_l):
# #        model = {}
# #        for k in gene_l: 
# #            if k in mod:
# #                model[k] = 1
# #            else:
# #                model[k] = 0
# #        score = model_score(model, obs_s, coalignment_mat)
# #        mod_score_list.append((model,score))
        
#     #Here to decide if there is one model that's obviously better
# #    return mod_score_list, MULTIPLE_MODELS

# # get a model and returns its likelihood score comparing the expected number of reads and the observed
# # model is basically just a bool dic of all the unique genes with flags of wether or not they're in model
# # observed is a dictionary of all gene combinations and their expected proportion
# # coalignment_mat is the coalignment matrix used to calculate the expected number of reads
# # eg:
# #   model - {A:1, B:0}
# #   observed - {A: 100 B:50, AB: 30 }
# #
# def model_score(model, observed, coalignment_mat):
#     exp = {}
#     tot = {}
#     for gene in model:
#         # patch for SC000
#         if gene==0:
#             tot[gene] = model[gene]*observed[gene,]
#         # Theres a common edge case where a gene A will only be aligned with other genes as well, in this case we update our observation vector to include A:0
#         elif (gene, ) not in observed:
#             tot[gene] = 0
#         elif gene not in coalignment_mat:
#             raise KeyError('{} not found in coalignment matrix'.format(gene))
#         elif (gene, ) not in coalignment_mat[gene]:
#             tot[gene] = 0
#         else:
#             tot[gene] = model[gene]*(observed[gene,]/coalignment_mat[gene][gene,])

#     keys = get_combinations(model.keys())   #get a list of all possible molecule combinations
    
#     # key is a set of genes and the expected number of reads for it is the sum of expected reads from all genes shared by the key,
#     # these in turn are the total reads for a gene (extrapoletaed from the uniqely mapped) multiplied by the coalignment factor (present in the coalignment matrix)
#     # e.g. if A has 20% coalignment with B and there are 80 reads mapped uniquely to A, we expect 80/0.8 * 0.2 = 20 reads to be mapped to AB from A (and more from B)
#     for k in keys:  
#         k = tuple(sorted(k))
#         sum = 0
#         for gene in k:
#             #Patch for SC000
#             if gene==0: 
#                 if k==(0,):
#                     sum=1
#                 else:
#                     sum = 0
#             #####
#             elif k in coalignment_mat[gene]:
#                 sum += tot[gene]*coalignment_mat[gene][k]
#         exp[k] = sum
    
#     score = calc_score(observed, exp)
#     return score

# def calc_score(obs, exp):
#     sum = 0
#     for k in obs:
#         if k not in exp:
#             print(k)
#             k = tuple(sorted(k))
#             print ('bad key')
#         diff = (obs[k]-exp[k])**2
#         if exp[k]!=0:
#             diff /= exp[k]
#         sum += diff
#     return sum

# #Get a dictionary of observations per gene/s and return a list of single unique genes
# def single_gene_list(obs):
#     l = []
#     for genes in obs:
#         for g in genes:
#             l.append(g)
#     return list(set(l))


