#!/usr/bin/bash                                                                 
if [ -z "$1" ]; then 
	echo "Provide cutoff"
	exit 1
fi
cutoff=$1
                                                                                
if [[ $2 = "--rerun" ]];                                                        
then                                                                            
    echo "submitting rerun jobs"                                                
    #./get_rerun_residues.py                                                    
    IFS=,                                                                       
    read -r -a residues < rerun_residues_$cutoff.csv                                    
else                                                                            
    echo "submitting all jobs"                                                  
    IFS=,                                                                       
    read -r -a residues < residue_list_$cutoff.csv                                      
fi                                                                              
                                                                                
                                                                                
for res in ${residues[@]}; do                                                   
    sed -r "s/RESIDUE/${res}/g;s/RESID/${res:1}/g;s/CUTOFF/$cutoff/g" submit_tmp.slu > submit.slu;
    sbatch submit.slu                                                           
    done                                                                        

