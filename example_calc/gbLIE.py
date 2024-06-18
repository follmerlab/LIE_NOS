import subprocess
import os
import pandas as pd
from datetime import datetime

def dir_exists(directory_name,verbose=False):
    """
    Checks if a directory exists in the current directory. If it does not exist, creates it.

    Parameters:
    directory_name (str): The name of the directory to check/create.
    """
    current_directory = os.getcwd()
    directory_path = os.path.join(current_directory, directory_name)
    
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        if verbose: print(f"Directory '{directory_name}' created.")
    else:
        if verbose: print(f"Directory '{directory_name}' already exists.")
        
def run_parmed(parmin, selection, ligand_name, file_type, charge=False):
    charge_command = f"change charge {selection} 0" if charge else f"strip {selection}"
    commands = f"""
    parmed <<EOF
    parm {parmin}
    {charge_command}
    parmout {file_type}{ligand_name}.parm7
    go
    EOF
    """
    subprocess.run(commands, shell=True, check=True)

def run_cpptraj(parmin, trajin, selection, charge=False):
    strip_selection = f"strip {selection}\n" if not charge else ""
    commands = f"""
    cpptraj <<EOF
    parm {parmin}
    trajin {trajin}    
    trajout temp_com.crd restart
    run
    clear all
    parm {parmin}
    trajin {trajin}      
    {strip_selection}trajout temp_rec.crd restart
    run
    clear all
    parm {parmin}
    trajin {trajin}       
    strip !{selection}
    trajout temp_lig.crd restart
    run
    exit
    EOF
    """
    subprocess.run(commands, shell=True, check=True)
    
def run_mm_gbsa(ligand_name, frame_suffix,proc=1):
    temp_gb = f"temp_gb{ligand_name}.in"
    with open('mmpbsa.in', 'r') as file:
        filedata = file.read()
    filedata = filedata.replace("XXX", ligand_name)
    filedata = filedata.replace("QQQ", f"{proc}")
    with open(temp_gb, 'w') as file:
        file.write(filedata)

    subprocess.run(f"mm_pbsa.pl {temp_gb} > bind.out", shell=True, check=True)
    subprocess.run(f"mv temp_statistics.out stats{ligand_name}_{frame_suffix}", shell=True, check=True)
    subprocess.run("rm temp*", shell=True, check=True)
    
def eel_vdw(file1, file2,alpha,beta,gamma):
    # Define column names
    column_names = ["Type", "Complex_MEAN", "Complex_STD", "Receptor_MEAN", "Receptor_STD", "Ligand_MEAN", "Ligand_STD"]

    # Reading the files with predefined column names
    data1 = pd.read_csv(file1, delim_whitespace=True, comment='#', names=column_names, skiprows=3)
    data2 = pd.read_csv(file2, delim_whitespace=True, comment='#', names=column_names, skiprows=3)

    # Extracting relevant values from the files
    ele_com = data1[data1['Type'] == 'ELE']['Complex_MEAN'].values[0]
    ele_rec = data1[data1['Type'] == 'ELE']['Receptor_MEAN'].values[0]
    gb_lig_bound = data2[data2['Type'] == 'GB']['Receptor_MEAN'].values[0]
    gb_lig_free = data1[data1['Type'] == 'GB']['Ligand_MEAN'].values[0]
    vdw_com = data1[data1['Type'] == 'VDW']['Complex_MEAN'].values[0]
    vdw_rec = data1[data1['Type'] == 'VDW']['Receptor_MEAN'].values[0]
    vdw_lig = data1[data1['Type'] == 'VDW']['Ligand_MEAN'].values[0]
    gbsur_lig_free = data1[data1['Type'] == 'GBSUR']['Ligand_MEAN'].values[0]

    # Calculating Delta EEL
    delta_eel_com_rec = ele_com - ele_rec
    delta_eel = delta_eel_com_rec + 2 * gb_lig_bound - 2 * gb_lig_free

    # Calculating Delta VDW
    delta_vdw = vdw_com - vdw_rec - gbsur_lig_free 

    # Final calculation
    #dG_lie = 0.50 * delta_eel + 0.18 * delta_vdw
    #α= 0.1005, β= 0.0336, and γ= -0.9083 

    dG_lie = alpha * delta_eel + beta * delta_vdw + gamma 
    # Output results
    print(f"dGlie = {dG_lie} for alpha={alpha}, beta={beta} and gamma={gamma}")
    print(f"dEEL and dVDW = {delta_eel} {delta_vdw}")
    return delta_eel, delta_vdw, dG_lie
    
def gblie(path,folder_list,parm7,traj,ligID, alpha , beta, gamma,proc=1,out="out"):
    now = datetime.now()
    date_ext = now.strftime("%Y_%m_%d_%H_%M")
    table = []
    for folder in folder_list:
        os.chdir(folder)
        isoform,inhibitor = folder.split('_')
        print(f"Running GB-LIE on for {inhibitor} bound to {isoform}.")
        #Trajectory and parm files
        print(os.listdir('.'))
        rundir = f"run_{date_ext}"
        dir_exists(rundir)
        os.popen(f'cp {path}/mmpbsa.in ./{rundir}')

        os.chdir(rundir)
        #run analysis

        # Initial setup
        subprocess.run(f"cp ../{parm7} complex.parm7", shell=True, check=True)
        subprocess.run(f"cp ../{traj} complex.trj.gz", shell=True, check=True)

        # Generate Receptor Parm File
        run_parmed("complex.parm7", f":{ligID}", ligID, "receptor")
        # Generate Ligand Parm File
        run_parmed("complex.parm7", f"!:{ligID}", ligID, "ligand")

        run_cpptraj("complex.parm7", "complex.trj.gz", f":{ligID}")
        run_mm_gbsa(ligID, out,proc)
        dir_exists("strip")
        subprocess.run("mv *.parm7 strip/", shell=True, check=True)

        # Generate Charge
        subprocess.run("cp ../strip.parm7 complex.parm7", shell=True, check=True)
        # Generate Receptor Parm File with charge change
        run_parmed("complex.parm7", f"!:{ligID}", ligID, "receptor", charge=True)
        # Generate Receptor Ligand File
        run_parmed("complex.parm7", f"!:{ligID}", ligID, "ligand")

        run_cpptraj("complex.parm7", "complex.trj.gz", f":{ligID}", charge=True)
        run_mm_gbsa(ligID, f"protein0_{out}",proc)
        dir_exists("charge")
        subprocess.run("mv *.parm7 charge/", shell=True, check=True)
        
        file1 = f"stats{ligID}_{out}"
        file2 = f"stats{ligID}_protein0_{out}"
        
        dEEL, dVDW, dGlie = eel_vdw(file1, file2,alpha,beta,gamma)
        table.append(
        {
            'Isoform': isoform,
            'Inhibitor': inhibitor,
            'dEEL': dEEL,
            'dVDW': dVDW,
            'dGlie': dGlie
        }
    )
        #Go to home directory
        os.chdir(path)
        data = pd.DataFrame(table)
    return data