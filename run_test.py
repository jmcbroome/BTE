import bte
import unittest

#load our test case as a universal object, outside of the test class.
t = bte.MATree("test.pb")

class TestMat(unittest.TestCase):
    def test_ps(self):
        ps = t.get_parsimony_score()
        self.assertTrue(ps > 0)

    def test_dfe(self):
        nodes = t.depth_first_expansion()
        self.assertTrue(len(nodes) > 0)

    def test_node(self):
        root = t.get_root()
        children = root.children
        self.assertTrue(len(children) > 0)
        for child in children:
            parent = child.parent
            self.assertTrue(parent.id == root.id)
            mutations = child.mutations
            self.assertTrue(len(mutations) > 0)

    