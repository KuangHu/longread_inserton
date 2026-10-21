import argparse, json, random, collections, collections.abc
from pathlib import Path
import numpy as np
# Compatibility only: upstream uses the NumPy 1.x alias removed in NumPy 2.
if not hasattr(np, 'Inf'):
    np.Inf = np.inf
if not hasattr(collections, 'Iterable'):
    collections.Iterable = collections.abc.Iterable
from biodive import bio
p=argparse.ArgumentParser()
p.add_argument('reads')
p.add_argument('out',type=Path)
a=p.parse_args()
a.out.mkdir(parents=True,exist_ok=True)
random.seed(42);np.random.seed(42)
c=bio.Config(outdir=str(a.out)+'/',annot_fasta=[])
(a.out/'config.json').write_text(json.dumps(vars(c),indent=2)+'\n')
bio.biodive_single_sample_analysis(a.reads,c)
