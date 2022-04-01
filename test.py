import bte
print("Loading tree...")
t = bte.MATree("public-latest.all.masked.pb.gz")
print("Tree successfully loaded.")
parsimony = t.get_parsimony_score()
print("Parsimony score:", parsimony)
print("Traversing tree...")
nodes = t.depth_first_expansion()
print("{} nodes traversed.".format(len(nodes)))
print("Checking the tenth node...")
testnode = nodes[9]
print("Name:",testnode.id)
children = testnode.children
print("{} children.".format(len(children)))
parent = testnode.parent
print("Parent:",parent.id)
mutations = testnode.mutations
print("Mutations: ",mutations)
print("Done!")