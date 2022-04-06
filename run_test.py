import bte
import unittest

#load our test case as a universal object, outside of the test class.
t = bte.MATree()
basic_newick = ""

class TestMat(unittest.TestCase):
    @classmethod
    def setUpClass(TestMat):
        t.from_pb("test.pb")
        basic_newick = t.write_newick()

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

    def test_mutation_update(self):
        minimuts = {'node_1':['A1234G']}
        t.apply_mutations(minimuts)
        n = t.get_node('node_1')
        self.assertIn('A1234G', n.mutations)

    def test_regex(self):
        #the test protobuf includes a few USA samples.
        testr = "USA/.*"
        usasub = t.get_regex(testr)
        leaves = usasub.get_leaves_ids()
        for l in leaves:
            self.assertRegex(l,testr)

    def test_mutation_accumulation(self):
        minimuts = {'node_1':['A1234G'],'node_2':['G1234A']}
        t.apply_mutations(minimuts)
        haplotypes = t.count_haplotypes()
        for h,c in haplotypes.items():
            #we fail if we have both- one or the other is fine.
            self.assertTrue(not (('A1234G' in h) and ('G1234A' in h)))

    def test_newick_load_write(self):
        basic_newick = t.write_newick()
        nomut_t = bte.MATree(nwk_string = basic_newick)
        ntl = nomut_t.get_leaves_ids()
        self.assertTrue(len(ntl) > 0)
        with self.assertRaises(Exception) as context:
            nomut_t.count_haplotypes()
        self.assertTrue('Tree does not contain explicit mutation information and this function cannot be used. You can add mutations with apply_mutations or load from a vcf or pb.' in str(context.exception))
