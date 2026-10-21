#!/usr/bin/env python3
"""Independent exact certificate check; Python standard library and a C++17 compiler.

No NumPy, SCIP, Boost, stored mask, neural network, or search output is trusted.
The conflict mask is derived from exact lattice metadata, then all labels are
checked. Output binds the verified inputs and checker sources.
"""
import argparse,contextlib,hashlib,json,subprocess,time
from fractions import Fraction
from pathlib import Path
from independent_audit import audit
ROOT=Path(__file__).resolve().parent

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--metadata',type=Path,default=ROOT/'artifacts/candidate2.json')
 p.add_argument('--labels',type=Path,default=ROOT/'artifacts/candidate2.labels.bin')
 p.add_argument('--out',type=Path,default=ROOT/'artifacts/final_verification')
 a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True);started=time.monotonic()
 meta=json.loads(a.metadata.read_text());shape=meta['shape']
 assert len(shape)==3 and all(type(n)==int and n>0 for n in shape)
 assert meta['real_colors']==14 and meta['removed_label']==14
 assert 'half-open' in meta['boundary']
 before={str(p):digest(p) for p in (a.metadata,a.labels)}
 if 'labels_sha256' in meta and meta['labels_sha256'] != before[str(a.labels)]:
  raise ValueError('Labels do not match metadata labels_sha256')
 with (a.out/'geometry.log').open('w') as log,contextlib.redirect_stdout(log):
  mask=a.out/'exact_mask.csv';geometry=audit(a.metadata,mask)
 binary=a.out/'independent_labels'
 command=['c++','-O3','-std=c++17',str(ROOT/'independent_labels.cpp'),'-o',str(binary)]
 subprocess.run(command,check=True)
 with (a.out/'labels.log').open('w') as log:
  proc=subprocess.run([str(binary),*map(str,shape),str(a.labels),str(mask)],stdout=subprocess.PIPE,stderr=log,text=True,check=True)
 labels=json.loads(proc.stdout);assert labels['verified']
 fraction=Fraction(labels['removed_count'],labels['total_cells'])
 assert labels['removed_count']==meta['removed_count'] and labels['total_cells']==meta['total_cells']
 assert fraction==Fraction(meta['removed_fraction'])
 assert before=={str(p):digest(Path(p)) for p in before},'Inputs changed during verification'
 sources=[ROOT/n for n in ('verify.py','independent_audit.py','geometry_exact.py','independent_labels.cpp')]
 result={'verified':True,'real_colors':14,'boundary':'half-open cells, checked using their closures','removed_fraction':str(fraction),'removed_percent':float(fraction*100),'geometry':geometry,'labels':labels,'input_sha256':before,'checker_sha256':{str(p):digest(p) for p in sources},'compile_command':command,'seconds':time.monotonic()-started}
 (a.out/'certificate.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ('verified','removed_fraction','removed_percent','seconds')},indent=2))
if __name__=='__main__':main()
