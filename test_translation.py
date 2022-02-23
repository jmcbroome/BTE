import mat
print("Loading tree...")
t = mat.MATree("public-latest.all.masked.pb.gz")
print("Tree successfully loaded.")
print("Attempting to translate...")
trand = t.translate_tree("ncbiGenes.gtf",'NC_045512v2.fa')
print("Translation successful. {} nodes translated.".format(len(trand)))