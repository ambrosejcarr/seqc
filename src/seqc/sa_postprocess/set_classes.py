class ASet:
	def __init__(self, name):
		self.parent, self.rank = None, 0
		self.name = name
		self.offspring_obj = []

class DisjointSet:
	def __init__(self):
		self.rentals = []
		self.items = []
		self.representatives = set([])

	def makeSet(self, new_addition): #where new_addition is of type ASet
		self.items.append(new_addition)
		self.rentals.append(new_addition)
		new_addition.parent = new_addition
		
	def Find(self, item):
		if item.parent == item:
			return item
		else:
			return self.Find(item.parent)
	
	def Union(self, item1, item2):
		item1_root = self.Find(item1)
		item2_root = self.Find(item2)
		if item1_root == item2_root:
			return
		if item1_root.rank < item2_root.rank:
			item1_root.parent = item2_root #loses item1_root
			
		elif item1_root.rank > item2_root.rank:
			item2_root.parent = item1_root
		else:
			item2_root.parent = item1_root
			item1_root.rank += 1

	def collect(self):
		"""Returns dictionary of representatives and the members of these representative's sets."""
		memberships = {}
		leaders = {}
		for dude in self.items:
			parent = self.Find(dude).name
			if parent not in memberships:
				memberships[parent] = []
			memberships[parent].append(dude.name)
			leaders[dude.name] = parent
		return memberships, leaders
		
	def __str__(self):
		outtie = ""	
		for dude in self.items:
			#for one in dude.offspring_obj:
			outtie += str(dude.name) + ", " + str(self.Find(dude).name) + "\n"
		return outtie
