

class Tree:

    def __init__(self, id, left=None, right=None, dist=None):
        self.id = id
        self.left = left
        self.right = right
        self.dist = dist

    def __repr__(self):
        return '<Node id=%s, left=%a, right=%a, dist=%a>' % (
            self.id,
            self.left.id if self.left is not None else None,
            self.right.id if self.left is not None else None,
            self.dist if self.dist is not None else None)

    @classmethod
    def from_linkage(cls, Z):
        current_id = Z.shape[0] * 2
        tree = {}
        for (left, right, dist, n_children) in Z[::-1]:
            tree[left] = Tree(id=left)
            tree[right] = Tree(id=right)
            if current_id not in tree:
                tree[current_id] = Tree(id=current_id, left=tree[left], right=tree[right], dist=dist)
            else:
                tree[current_id].left = tree[left]
                tree[current_id].right = tree[right]
                tree[current_id].dist = dist
            current_id -= 1
        return tree[max(tree.keys())]

    def is_leaf(self):
        return True if self.left is None and self.right is None else False

    @staticmethod
    def nodes2labels(nodes):
        return [n.id for n in nodes]

    def get_daughter(self, id_):
        for daughter in self.dfs():
            if daughter.id == id_:
                return daughter
        return None

    def has_daughter(self, id_):
        for daughter in self.dfs():
            if daughter.id == id_:
                return True
        return False

    def dfs(self):
        visited, stack = [], [self]
        while stack:
            vertex = stack.pop()
            yield vertex
            if vertex not in visited:
                visited.append(vertex)
                if vertex.left is not None:
                    stack.append(vertex.left)
                if vertex.right is not None:
                    stack.append(vertex.right)

    def bfs(self):
        visited, queue = [], [self]
        while queue:
            vertex = queue.pop(0)
            yield vertex
            if vertex not in visited:
                visited.append(vertex)
                if vertex.left is not None:
                    queue.append(vertex.left)
                if vertex.right is not None:
                    queue.append(vertex.right)
