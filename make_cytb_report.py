#!/usr/bin/env python3

import sys
import pandas as pd

treefile=sys.argv[1]
html=sys.argv[2]

summary=None
if len(sys.argv)>3:
    summary=sys.argv[3]

html_content="""

<html>
<head>
<title>CYTB Phylogeny Report</title>
</head>

<body>

<h1>CYTB Phylogenetic Analysis</h1>

<h2>Rectangular Tree</h2>
<img src="tree_rectangular.svg" width="900">

<h2>Circular Tree</h2>
<img src="cytb_tree_circular.svg" width="900">

"""

if summary:
    df=pd.read_csv(summary,sep="\t")
    html_content+="<h2>Dataset Summary</h2>"
    html_content+=df.to_html(index=False)

html_content+="</body></html>"

open(html,"w").write(html_content)