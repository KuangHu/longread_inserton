"""Run an upstream command with bounded wall time and record exit/resource status."""
import argparse, datetime, json, os, resource, signal, subprocess, time
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument('--out',type=Path,required=True)
p.add_argument('--timeout',type=int,default=1200)
p.add_argument('--memory-gb',type=float,default=6)
p.add_argument('command',nargs=argparse.REMAINDER)
a=p.parse_args()
cmd=a.command[1:] if a.command[0]=='--' else a.command
a.out.mkdir(parents=True,exist_ok=True)
env=os.environ.copy()
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    env[k]='1'
def limits():
    resource.setrlimit(resource.RLIMIT_AS,(int(a.memory_gb*2**30),int(a.memory_gb*2**30)))
t=time.monotonic()
rec={'command':cmd,'cwd':os.getcwd(),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'timeout_seconds':a.timeout,'address_space_limit_gib':a.memory_gb}
(a.out/'started.json').write_text(json.dumps(rec,indent=2)+'\n')
with (a.out/'run.log').open('w') as log:
    proc=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True,preexec_fn=limits)
    try:
        rec['exit_code']=proc.wait(timeout=a.timeout)
        rec['status']='completed' if rec['exit_code']==0 else 'failed'
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid,signal.SIGTERM)
        try: proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid,signal.SIGKILL);proc.wait()
        rec.update(status='wall_time_limit',exit_code=proc.returncode)
rec['wall_seconds']=round(time.monotonic()-t,3)
rec['max_rss_kib']=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
(a.out/'run.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec,indent=2))
