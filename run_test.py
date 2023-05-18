import bte
import unittest
import os
import sys
#load our test case as a universal object, outside of the test class.
t = bte.MATree()

def check_tree_struct(t1,t2):
    return (t1.get_newick() == t2.get_newick())

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
        basic_newick = t.get_newick()
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
        nwk = t.get_newick()
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
        self.assertTrue("node_X" in [c.id for c in t.root.children])
        t.create_node("node_Y",t.root.id)
        t.move_node("node_Y","node_X")
        self.assertTrue(t.get_node("node_Y").parent.id == "node_X")
        t.remove_node("node_X")
        self.assertTrue("node_X" not in [c.id for c in t.root.children])

    def test_branch_length(self):
        t.create_node("node_Z",t.root.id, branch_length = 0.5)
        self.assertTrue(t.get_node("node_Z").branch_length == 0.5)
        t.apply_mutations({'node_Z':['A1234G']},True)
        self.assertTrue(t.get_node("node_Z").branch_length == 1)
        minfo = t.get_node("node_Z").get_mutation_information()[0]
        self.assertTrue(minfo['position'] == 1234)
        self.assertTrue(minfo['par_nuc'] == 'A')
        self.assertTrue(minfo['mut_nuc'] == 'G')
        t.remove_node("node_Z")

    def test_annotation(self):
        to_apply = {t.root.children[0].id:['ann1','ann2']}
        t.apply_node_annotations(to_apply)
        dump = t.dump_node_annotations()
        self.assertTrue(to_apply[t.root.children[0].id] == dump[t.root.children[0].id])

    def test_lca(self):
        self.assertTrue(t.LCA(['node_4','node_2']) == 'node_1')