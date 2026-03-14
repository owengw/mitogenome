#!/usr/bin/env python3

import sys
from ete3 import Tree
from ete3.treeview import TreeStyle, NodeStyle

treefile=sys.argv[1]
outfile=sys.argv[2]

t=Tree(treefile)

ts=TreeStyle()
ts.mode="c"
ts.show_leaf_name=True
ts.scale=100

t.render(outfile,tree_style=ts)