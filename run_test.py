import bte
import unittest
import os
#load our test case as a universal object, outside of the test class.
t = bte.MATree()

def check_tree_struct(t1,t2):
    return (t1.write_newick() == t2.write_newick())

class TestMat(unittest.TestCase):
    @classmethod
    def setUpClass(TestMat):
        t.from_pb("test.pb")

    def test_ps(self):
        ps = t.get_parsimony_score()
        self.assertTrue(ps > 0)

    def test_dfe(self):
        nodes = t.depth_first_expansion()
        self.assertTrue(len(nodes) > 0)

    def test_node(self):
        children = t.root.children
        self.assertTrue(len(children) > 0)
        for child in children:
            parent = child.parent
            self.assertTrue(parent.id == t.root.id)
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

    def test_json_load_write(self):
        t.write_json("test.json")
        t2 = bte.MATree(json_file = "test.json")
        os.remove("test.json")

    def test_newick_vcf_load_write(self):
        nwk = t.write_newick()
        with open("test.nwk",'w+') as f:
            f.write(nwk)
        t.write_vcf("test.vcf")
        t2 = bte.MATree(nwk_file = "test.nwk", vcf_file = "test.vcf")
        self.assertTrue(check_tree_struct(t,t2))
        os.remove("test.nwk")
        os.remove("test.vcf")

    def test_node_manipulation(self):
        t.create_node("node_X",t.root.id)
        self.assertTrue(t.get_node("node_X").parent.id == t.root.id)
        t.create_node("node_Y",t.root.id)
        t.move_node("node_X","node_Y")
        self.assertTrue(t.get_node("node_X").parent.id == "node_Y")
        t.remove_node("node_Y")
