import pysam
import pandas as pd
import os

onlyfiles = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]
onlyfiles = [f for f in onlyfiles if not ".bai" in f ]
onlyfiles = [f for f in onlyfiles if  ".bam" in f ]

a = {"name":[], "mapped":[], "unmapped":[]}
for filename in onlyfiles:
    a["name"].append(filename.split(".")[0])
    a["mapped"].append(reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[2]) for l in pysam.idxstats(filename) ]))
    a["unmapped"].append(reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[3]) for l in pysam.idxstats(filename) ]))
b = pd.DataFrame(a)
b.to_csv("counts.csv", sep = ",", index = False)