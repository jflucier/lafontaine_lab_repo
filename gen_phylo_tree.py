import pandas as pd
from ete3 import NCBITaxa, AttrFace, faces, TreeStyle

def generate_tree(f):
    res = pd.read_csv(f, sep='\t', index_col="accession.1")

    x = res.groupby('accession.1')
    print(f"{x}")
    for acc, row in x:
        print(f"{acc}")
        print(f"{row['tax id']}")
        tax_tmp = row['tax id'].dropna()
        # tax_tmp = row[row['tax id'].notnull()]

        tx_lst = []
        for z in tax_tmp:
            tx_lst.append(int(z))

        ncbi = NCBITaxa()
        tree = ncbi.get_topology(tx_lst)

        # custom layout: adds "rank" on top of branches, and sci_name as tip names
        def my_layout(node):
            if getattr(node, "rank", None):
                rank_face = AttrFace("sci_name", fsize=15, fgcolor="indianred")
                node.add_face(rank_face, column=0, position="branch-top")
            if node.is_leaf():
                sciname_face = AttrFace("sci_name", fsize=15, fgcolor="steelblue")
                node.add_face(sciname_face, column=0, position="branch-right")

        ts = TreeStyle()
        ts.layout_fn = my_layout
        ts.show_leaf_name = False
        ts.mode = "c"
        ts.title.add_face(faces.TextFace(f"{acc}", fsize=60), 0)
        # ts.arc_start = -180 # 0 degrees = 3 o'clock
        # ts.arc_span = 180

        #print(tree.get_ascii(attributes=["sci_name", "rank"]))

        tree.render(f"/storage/Documents/service/biologie/lafontaine/20230920_riboswitch_eukaryotes/tax_tree/{acc}.png", tree_style=ts)
        # tree.render("/home/jflucier/tmp/test.svg")


if __name__ == '__main__':

    f = "/storage/Documents/service/biologie/lafontaine/20230920_riboswitch_eukaryotes/all.found.taxid.tsv"

    generate_tree(f)
