"""Explicit one-shot LIMS submission and read-only artifact collection. No cancellation/restart."""
import argparse, contextlib, hashlib, json, subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote
import requests

ROOT = Path(__file__).resolve().parent
BASE = 'http://127.0.0.1:8200/api/sim'
SHA = 'fecadd3e13cdb7003d8c7fc7a0365bd588254730'
TERMINAL = {'succeeded','failed','canceled'}

def save(name, value):
    (ROOT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')

def get(path):
    r=requests.get(BASE+path,timeout=(5,30));r.raise_for_status();return r.json()

def job_id():
    body=json.loads((ROOT/'submission-response.json').read_text())['body']
    return body.get('job_id') or body.get('job',{}).get('job_id')

def download(item):
    name=item['path'];p=Path(name)
    assert not p.is_absolute() and '..' not in p.parts
    target=ROOT/p
    if target.exists() and item.get('sha256') and hashlib.sha256(target.read_bytes()).hexdigest()==item['sha256']:
        return {'path':name,'sha256':item['sha256'],'bytes':target.stat().st_size,'cached':True}
    target.parent.mkdir(parents=True,exist_ok=True)
    url=BASE+'/jobs/'+job_id()+'/files/'+quote(name,safe='/')
    with requests.get(url,stream=True,timeout=(10,120)) as r:
        r.raise_for_status();h=hashlib.sha256();count=0
        with target.open('wb') as f:
            for chunk in r.iter_content(1024*1024): f.write(chunk);h.update(chunk);count+=len(chunk)
    if item.get('sha256'): assert h.hexdigest()==item['sha256'],name
    if name!='task.md': assert count==item['size'],(name,count,item['size'])
    return {'path':name,'sha256':h.hexdigest(),'bytes':count,'matched_server_sha':bool(item.get('sha256'))}

def main():
    global ROOT, SHA
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','submit','status','files','download-all'])
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--solver-sha', default=SHA)
    a=p.parse_args(); ROOT=a.root.resolve(); SHA=a.solver_sha
    assert len(SHA)==40 and all(c in '0123456789abcdef' for c in SHA)
    if a.action=='prepare':
        assert not (ROOT/'submission-attempt.json').exists()
        raw=subprocess.check_output(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8','demo-process-121',
             'cat /home/kevin/DamSafetyApp-cdp-releases/'+SHA+'/release-audit.json'])
        audit=json.loads(raw);assert audit['solver_sha']==SHA
        (ROOT/'input/release-audit.json').write_bytes(raw)
        s=json.loads((ROOT/'submission.json').read_text())
        s['extra_files']=[f for f in s['extra_files'] if f['name']!='release-audit.json']
        s['extra_files'].append({'name':'release-audit.json','sha256':hashlib.sha256(raw).hexdigest()})
        save('submission.json',s);print('Prepared locked release audit; no submission');return
    if a.action=='submit':
        health = get('/health')
        assert health.get('status') == 'ok' and health.get('dam_safety_app_sha') == SHA
        mark=ROOT/'submission-attempt.json';assert not mark.exists(),'Already attempted; inspect response, never auto-repeat POST'
        agent = get('/jobs?limit=1')['agent']
        assert agent['reachable']
        r=requests.get(agent['base_url'].rstrip('/')+'/healthz',timeout=10);r.raise_for_status()
        assert r.json()['dam_safety_app_sha']==SHA
        for state in ('running','queued'):
            d=get('/jobs?state='+state+'&limit=10');assert d['agent']['reachable'] and not d['jobs']
        s=json.loads((ROOT/'submission.json').read_text());declared={s['input_file']:s['input_sha256']}
        declared.update({f['name']:f['sha256'] for f in s['mesh_files']+s['extra_files']})
        assert 'uniaxial_compression.csv' not in declared
        with contextlib.ExitStack() as stack:
            files=[]
            for name,digest in declared.items():
                path=ROOT/'input'/name;assert hashlib.sha256(path.read_bytes()).hexdigest()==digest
                files.append(('files',(name,stack.enter_context(path.open('rb')),'application/octet-stream')))
            save('submission-attempt.json',{'operation':'POST /api/sim/jobs','do_not_auto_retry':True,'health':r.json()})
            r=requests.post(BASE+'/jobs',data={'manifest':json.dumps(s)},files=files,timeout=(10,120))
            save('submission-response.json',{'status_code':r.status_code,'body':r.json()})
            r.raise_for_status();print(json.dumps({'job_id':job_id(),'state':r.json().get('state'),'case':s['case_name']},ensure_ascii=False))
        return
    jid=job_id();assert jid
    if a.action=='status':
        m=get('/jobs/'+jid);e=get('/jobs/'+jid+'/execution-status')
        save('manifest.json',m);save('execution-status.json',e)
        print(json.dumps({'job_id':jid,'state':m.get('state'),'progress':m.get('progress'),'failure':m.get('failure')},ensure_ascii=False));return
    listing=get('/jobs/'+jid+'/files');save('files.json',listing)
    if a.action=='files':
        print(json.dumps(listing,ensure_ascii=False));return
    assert listing['state'] in TERMINAL,'Do not snapshot all files while solver writes'
    with ThreadPoolExecutor(max_workers=3) as pool:
        results=list(pool.map(download,listing['files']))
    save('download-checksums.json',results)
    print(json.dumps({'files':len(results),'bytes':sum(x['bytes'] for x in results),'server_hash_matches':sum(x.get('matched_server_sha',False) or x.get('cached',False) for x in results)}))

if __name__=='__main__': main()
