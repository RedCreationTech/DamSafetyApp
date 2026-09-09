#!/usr/bin/env python3
"""Run real GoogleTest core tests with isolated CSV/Real adapters, not MOOSE tests."""
import argparse
import pathlib
import subprocess
import json
import hashlib
import sys

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--gtest-root',required=True,type=pathlib.Path)
p.add_argument('--output',required=True,type=pathlib.Path)
p.add_argument('--cxx',default='clang++')
a=p.parse_args()
root=pathlib.Path(__file__).resolve().parents[2]
gtest=a.gtest_root.resolve()/'googletest'
out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
sources=sorted((root/'src/utils').glob('*.C'))
tests=sorted((root/'unit/src').glob('AbaqusCDP*Test.C'))
command=[a.cxx,'-std=c++17','-O2','-pthread','-I',str(root/'include/utils'),'-I',str(root/'unit/standalone/include'),'-I',str(gtest/'include'),'-I',str(gtest),*map(str,sources),*map(str,tests),str(gtest/'src/gtest-all.cc'),str(gtest/'src/gtest_main.cc'),'-o',str(out/'core-unit-tests')]
(out/'build-command.json').write_text(json.dumps(command,indent=2)+'\n')
with (out/'build.log').open('w') as log:
 subprocess.run(command,cwd=root,stdout=log,stderr=subprocess.STDOUT,check=True)
with (out/'tests.log').open('w') as log:
 result=subprocess.run([str(out/'core-unit-tests'),'--gtest_output=json:'+str(out/'tests.json')],cwd=root,stdout=log,stderr=subprocess.STDOUT)
identity={str(x.relative_to(root)):hashlib.sha256(x.read_bytes()).hexdigest() for x in [*sources,*tests,*sorted((root/'unit/standalone/include').glob('*.h'))]}
(out/'source-hashes.json').write_text(json.dumps(identity,indent=2)+'\n')
print((out/'tests.log').read_text()[-1600:])
print('Scope: core software regression using GoogleTest; test-only Real/CSV adapters; no MOOSE framework, TestHarness or finite-element Job.')
sys.exit(result.returncode)
