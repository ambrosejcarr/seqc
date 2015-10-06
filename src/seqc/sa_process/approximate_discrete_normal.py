#You know what, http://www.biomedcentral.com/1471-2164/13/42#B1 says that fragment length can be modeled by a Poisson, and to be honest, Poisson sorta looks like a Gaussian, but is discrete. So we're using Poisson.

from scipy.stats import poisson

def create_lookup(lamb, max_val):
	lookup = {}
	for count in range(max_val):
		lookup[count] = poisson.pmf(count, lamb)
	return lookup

def create_kde_lookup(max_val=950):
	import pickle
	pk = open("pickled_kde.pckl")
	kde = pickle.load(pk)
	lookup = {}
	for count in range(max_val + 1):
		lookup[max_val - count] = kde.integrate_box_1d(count, count+1)
	return lookup

if __name__ == "__main__":
	lk = create_kde_lookup(951)

	for k in list(lk.keys())[:5]:
		print(k, lk[k])

	nums = list(lk.values())
	print(1 - sum(nums))
