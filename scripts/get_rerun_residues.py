import numpy as np
from glob import glob                                    
import os                                                

bashCommand = "sacct -n -X --format jobname,jobid -s RUNNING > .running_jobs.txt"
os.system(bashCommand)

running = []
with open('.running_jobs.txt', 'r') as f:
    for line in f:
        running.append(line.strip().split()[0])

os.chdir('basicrta-7.0/')                               
dirs = np.array(glob('?[0-9]*'))                                        

if len(running)>0:
    run_inds = np.array([np.where(dirs==res)[0][0] for res in running if res in dirs])
    if len(run_inds)>0:
        dirs = np.delete(dirs, run_inds)

rerundirs = []                                           
for adir in dirs:                                        
    if os.path.exists(f'{adir}/gibbs_110000.pkl') or os.path.exists(f'{adir}/.dataset_too_small'): 
        pass
    else:
        rerundirs.append(adir)                           
                                                         
#rerunids = [adir[:-1] for adir in rerundirs]             
rerunids = [adir for adir in rerundirs]             
with open('../rerun_residues.csv', 'w') as w:            
    for res in rerunids:                                 
        w.write(f'{res},')                               

