import math
import os
import shutil
import subprocess
import pickle
import pandas as pd
import numpy as np
import plotly.io as pio
import plotly.graph_objects as go
import base64
from Bio import PDB
from config import RAMA_PREFERENCES

RAMA_PREF_VALUES = None
ACCEPTABLE_CUTOFF_MPNN = 0.03
OPTIMAL_CUTOFF_MPNN = 0.20
PRO_COLOR = "#6F8FAF"
ACCEPTABLE_COLOR = "#FAFA33"
OPTIMAL_COLOR = "#5ed35e"
BAD_COLOR = "#D3D3D3"
PRO_GRAPH_COLOR = '#517394'

def _cache_RAMA_PREF_VALUES():

    cache_file_path = "/www/ProScan/rama_pref_values_cache.pkl"

    # Check if the cache file exists
    if os.path.exists(cache_file_path):
        # If the cache file exists, load the values from the file
        with open(cache_file_path, 'rb') as cache_file:
            RAMA_PREF_VALUES = pickle.load(cache_file)

    #If cache fails, then have to recache item from the data files
    else:      
        f_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        RAMA_PREF_VALUES = {}
        for key, val in RAMA_PREFERENCES.items():
            RAMA_PREF_VALUES[key] = np.full((360, 360), 0, dtype=np.float64)
            with open(os.path.join(f_path, val["file"])) as fn:
                for line in fn:
                    if line.startswith("#"):
                        continue
                    else:
                        x = int(float(line.split()[1]))
                        y = int(float(line.split()[0]))
                        RAMA_PREF_VALUES[key][x + 180][y + 180] \
                            = RAMA_PREF_VALUES[key][x + 179][y + 179] \
                            = RAMA_PREF_VALUES[key][x + 179][y + 180] \
                            = RAMA_PREF_VALUES[key][x + 180][y + 179] \
                            = float(line.split()[2]) 
        
        # Save the computed values to the cache file
        with open(cache_file_path, 'wb') as cache_file:
            pickle.dump(RAMA_PREF_VALUES, cache_file)

    return RAMA_PREF_VALUES

#Takes PDB file name and ramachandran type to produce as strings. Writes file
#with angles of each residue and stores residue angles in dictionaries to plot
#on Ramachandran of type ramaType.
def calc_ramachandran(pdb_file_name, file_path):
    
    global RAMA_PREF_VALUES         
    RAMA_PREF_VALUES = _cache_RAMA_PREF_VALUES()
   
    extra_options_check = check_extra_options(file_path) #Check if Rosetta will be run and if scan is a glycine scan instead of proline
    run_ros = extra_options_check[0]
    run_gly = extra_options_check[1]
    get_hbonds(file_path, pdb_file_name, run_gly)     
    output_path = os.path.join(file_path, pdb_file_name + "_ProScan_results_all.txt")
    rama_list = open(output_path,"w")  
    
    #Get MPNN scores, produce reordered pdb for visualization on results page
    pdb_file = os.path.join(file_path,pdb_file_name)
    MPNN_scores = run_proteinMPNN(pdb_file_name,file_path, run_gly) 
    reorder_pdb(pdb_file_name, file_path, MPNN_scores, run_gly) 
    
    #Get basic angle classification info. Get secondary structure info
    res_info_dict = get_res_info(pdb_file_name,file_path, run_gly)
    res_info_dict = get_secondary_struct(pdb_file,res_info_dict)
    sorted_keys = sorted(res_info_dict.keys(), key=sort_res_dict)

    if run_ros == 'True' and run_gly == 'False':
        rama_list.write("Residue\tID\tChain\tPhi\tPsi\tSec_Struct\tProteinMPNN_Prob\tPro_Angle\tPre-Pro_Angle\tDDG_Pred\tNotes\n")
    elif run_ros == 'False' and run_gly == 'False':
        rama_list.write("Residue\tID\tChain\tPhi\tPsi\tSec_Struct\tProteinMPNN_Prob\tPro_Angle\tPre-Pro_Angle\tNotes\n")     
    elif run_ros == 'True' and run_gly == 'True':
        rama_list.write("Residue\tID\tChain\tPhi\tPsi\tSec_Struct\tProteinMPNN_Prob\tGly_Angle\tDDG_Pred\tNotes\n")
    else:
        rama_list.write("Residue\tID\tChain\tPhi\tPsi\tSec_Struct\tProteinMPNN_Prob\tGly_Angle\tNotes\n")

    if not MPNN_scores: #If MPNN failed, return ?
        num_residues = len(sorted_keys)
        MPNN_scores = ['0'] * num_residues
                       
    for key, MPNN_score in zip(sorted_keys, MPNN_scores):
        res_info_dict[key].append(MPNN_score)
    
    #If user checked the option, add Rosetta scores
    if run_ros == 'True':  
        execute_ros(pdb_file_name,file_path, run_gly)      
        res_info_dict = parse_ros(file_path,res_info_dict)
         
    res_info_dict = append_disrupted_hbonds(file_path, res_info_dict)

    #Write out contents to text output file
    if run_ros == 'True' and run_gly == 'False':
        for key in sorted_keys: 
            curr_res = res_info_dict[key]
            rama_list.write("\t".join([str(curr_res[0]), str(curr_res[1]), str(curr_res[2]), str(round(float(curr_res[3]),3)), str(round(float(curr_res[4]),3)), 
                           str(curr_res[7]), str(round(float(curr_res[8]),4)), str(curr_res[5]), str(curr_res[6]), str(round(float(curr_res[9]),3)), str(curr_res[10]) + "\n"]))     
    elif run_ros == 'False' and run_gly == 'False':
        for key in sorted_keys: 
            curr_res = res_info_dict[key]
            rama_list.write("\t".join([str(curr_res[0]), str(curr_res[1]), str(curr_res[2]), str(round(float(curr_res[3]),3)), str(round(float(curr_res[4]),3)), 
                        str(curr_res[7]), str(round(float(curr_res[8]),4)), str(curr_res[5]), str(curr_res[6]), str(curr_res[9]) + "\n"]))           
    elif run_ros == 'False' and run_gly == 'True':
        for key in sorted_keys: 
            curr_res = res_info_dict[key]  
            rama_list.write("\t".join([str(curr_res[0]), str(curr_res[1]), str(curr_res[2]), str(round(float(curr_res[3]),3)), str(round(float(curr_res[4]),3)), 
                           str(curr_res[6]), str(round(float(curr_res[7]),4)), str(curr_res[5]), str(curr_res[8]) + "\n"]))
    else:
        for key in sorted_keys: 
            curr_res = res_info_dict[key]  
            rama_list.write("\t".join([str(curr_res[0]), str(curr_res[1]), str(curr_res[2]), str(round(float(curr_res[3]),3)), str(round(float(curr_res[4]),3)), 
                           str(curr_res[6]), str(round(float(curr_res[7]),4)), str(curr_res[5]), str(round(float(curr_res[8]),3)), str(curr_res[9]) + "\n"]))

    rama_list.close() 
    
    return None         

def sort_res_dict(res_key):

    res_chain = res_key[0]
    res_id = int(res_key[1:])
    return (res_chain,res_id)

def get_res_info(pdb_file_name, file_path, run_gly):

    res_info_dict = {}
    pdb_file = os.path.join(file_path,pdb_file_name)
    
    if pdb_file[-4:] == ".pdb":
        structure = PDB.PDBParser().get_structure("input_structure","%s" % pdb_file)    
    elif pdb_file[-4:] == ".cif":
        structure = PDB.PDBParser().get_structure("input_structure","%s" % pdb_file[:-4] + '.pdb')
    
    model = structure.child_list[0]

    if run_gly == "False": #Scan for pro/pre-pro angle preferences
        for chain in model.get_chains():        
        
            polypeptides = PDB.PPBuilder().build_peptides(chain)
            prev_prepro_angle = "No_Angle"
        
            for poly_index, poly in enumerate(polypeptides):            
                phi_psi = poly.get_phi_psi_list()
            
                for res_index, residue in enumerate(poly):                        
                    phi, psi = phi_psi[res_index] #store phi and psi angle of current residue
                    
                    #If no phi or psi angle calculated (first/last res of chain)               
                    if not phi or not psi:
                          
                        curr_res_info = [residue.get_resname(), residue.id[1], chain.get_id(), float(0), float(0), 'No_Angle', 'No_Angle']
                        curr_res_key = str(curr_res_info[2]) + str(curr_res_info[1])
                        res_info_dict[curr_res_key] = curr_res_info              
                    elif phi and psi: #If residue has phi and psi angle

                        pro_angle = assessAngle(phi, psi, "PRO")
                        pre_angle = assessAngle(phi, psi, "PRE-PRO")
                        curr_res_info = [residue.get_resname(), residue.id[1], chain.get_id(), math.degrees(phi), math.degrees(psi), pro_angle, prev_prepro_angle]
                        curr_res_key = str(curr_res_info[2]) + str(curr_res_info[1])
                        res_info_dict[curr_res_key] = curr_res_info 
                        prev_prepro_angle = pre_angle   
    else:
        for chain in model.get_chains():  #Scan for glycine angle preferences      
        
            polypeptides = PDB.PPBuilder().build_peptides(chain)
        
            for poly_index, poly in enumerate(polypeptides):            
                phi_psi = poly.get_phi_psi_list()
            
                for res_index, residue in enumerate(poly):                        
                    phi, psi = phi_psi[res_index] #store phi and psi angle of current residue
                    
                    #If no phi or psi angle calculated (first/last res of chain)               
                    if not phi or not psi:
                          
                        curr_res_info = [residue.get_resname(), residue.id[1], chain.get_id(), float(0), float(0), 'No_Angle',]
                        curr_res_key = str(curr_res_info[2]) + str(curr_res_info[1])
                        res_info_dict[curr_res_key] = curr_res_info 
                                       
                    elif phi and psi: #If residue has phi and psi angle

                        pro_angle = assessAngle(phi, psi, "GLY")
                        curr_res_info = [residue.get_resname(), residue.id[1], chain.get_id(), math.degrees(phi), math.degrees(psi), pro_angle]
                        curr_res_key = str(curr_res_info[2]) + str(curr_res_info[1])
                        res_info_dict[curr_res_key] = curr_res_info 
                    
    return res_info_dict

def check_extra_options(file_path):

    chain_info_file = os.path.join(file_path,'chain_info.txt')
    extra_options = []

    with open(chain_info_file,'r') as chain_info:
        chain_info.readline()
        chain_info.readline()
        ros_check = chain_info.readline().strip()
        run_gly = chain_info.readline().strip()

    extra_options.append(ros_check)
    extra_options.append(run_gly)

    return extra_options

def clean_run_files(file_path,file_name):

    dssp_file = os.path.join(file_path, file_name[:-4] + ".dssp")
    aa_prob_file = os.path.join(file_path, "aaProbData.txt")
    prob_dir = os.path.join(file_path, "probs")
    seqs_dir = os.path.join(file_path, "seqs")

    if os.path.isfile(dssp_file):
        os.remove(dssp_file)
    
    if os.path.isfile(aa_prob_file):
        os.remove(aa_prob_file)

    if os.path.exists(prob_dir):
        shutil.rmtree(prob_dir)

    if os.path.exists(seqs_dir):
        shutil.rmtree(seqs_dir)   

def assessAngle(phi, psi, angleType):
    
    #Append residues to respective normal or outlier dictionary for scatterplot
    #Write residue ifnormation to ramaOutput txt file
    if RAMA_PREF_VALUES[angleType][int(math.degrees(psi)) + 180][int(math.degrees(phi)) + 180] < \
        RAMA_PREFERENCES[angleType]["bounds"][1]:                                
            angAssessment = "Questionable"                                                                        
    #if residue falls within partially preferred bounds
    elif RAMA_PREF_VALUES[angleType][int(math.degrees(psi)) + 180][int(math.degrees(phi)) + 180] > \
        RAMA_PREFERENCES[angleType]["bounds"][1] and RAMA_PREF_VALUES[angleType][int(math.degrees(psi)) + 180][int(math.degrees(phi)) + 180] < \
            RAMA_PREFERENCES[angleType]["bounds"][2]:                               
        angAssessment = "Acceptable"                                                                        
                                        
        #if residue falls within the bounds of preferrable 
    elif RAMA_PREF_VALUES[angleType][int(math.degrees(psi)) + 180][int(math.degrees(phi)) + 180] > \
        RAMA_PREFERENCES[angleType]["bounds"][2]:                              
        angAssessment = "Preferable"                                                                        
    
    return angAssessment

def color_map(row):
     
    # Otherwise color according to proteinMPNN criteria
    if float(row['ProteinMPNN <br> Probability']) < ACCEPTABLE_CUTOFF_MPNN:
        return [''] * len(row)
    
    elif 'Pro Angle' in row:
        # If residue is proline color row light grey
        if row['Residue'] == 'PRO':
            return [f'background-color:' + PRO_COLOR for col in row]
        elif ACCEPTABLE_CUTOFF_MPNN <= float(row['ProteinMPNN <br> Probability']) <= OPTIMAL_CUTOFF_MPNN and row['Pro Angle'] != 'Questionable':
            return [f'background-color:' + '#ffe599' for col in row]
        elif float(row['ProteinMPNN <br> Probability']) > OPTIMAL_CUTOFF_MPNN and row['Pro Angle'] != 'Questionable' and row['Pro Angle'] != 'No_Angle':
            return [f'background-color:' + OPTIMAL_COLOR for col in row]
        else:
            return [''] * len(row)

    elif 'Gly Angle' in row:
        if row['Residue'] == 'GLY':
            return [f'background-color:' + PRO_COLOR for col in row]
        if ACCEPTABLE_CUTOFF_MPNN <= float(row['ProteinMPNN <br> Probability']) <= OPTIMAL_CUTOFF_MPNN and row['Gly Angle'] != 'Questionable':
            return [f'background-color:' + '#ffe599' for col in row]
        elif float(row['ProteinMPNN <br> Probability']) > OPTIMAL_CUTOFF_MPNN and row['Gly Angle'] != 'Questionable' and row['Gly Angle'] != 'No_Angle':
            return [f'background-color:' + OPTIMAL_COLOR for col in row]
        else:
            return [''] * len(row)

def output_to_HTML(file_path, job_name, results_file_type):
    
    #Whether full results or partial results is converted depends on results_file_type
    if results_file_type == "all":
        data_path = os.path.join(file_path, job_name)
    elif results_file_type == "top":
        data_path = os.path.join(file_path, job_name)

    rama_list = pd.read_csv(data_path, sep = '\t', header = 0)
    rama_list['Phi'] = rama_list['Phi'].astype(str)
    rama_list['Psi'] = rama_list['Psi'].astype(str)
    rama_list['ProteinMPNN_Prob'] = rama_list['ProteinMPNN_Prob'].astype(str)
    
    if 'DDG_Pred' in rama_list.columns:
        rama_list['DDG_Pred'] = rama_list['DDG_Pred'].astype(str) 
    if 'Pre-Pro_Angle' in rama_list.columns:
        rama_list.rename(columns = {'Pre-Pro_Angle': 'Pre-Pro <br> Angle'}, inplace = True)
    if 'Pro_Angle' in rama_list.columns:
        rama_list.rename(columns = {'Pro_Angle': 'Pro Angle'}, inplace = True)
    if 'Gly_Angle' in rama_list.columns:
        rama_list.rename(columns = {'Gly_Angle': 'Gly Angle'}, inplace = True)
    rama_list.rename(columns = {'ProteinMPNN_Prob': 'ProteinMPNN <br> Probability'}, inplace = True)
    rama_list.rename(columns = {'Sec_Struct': 'Sec Struct'}, inplace = True)
    
    #Format HTML table
    rama_list_styled = rama_list.style.apply(color_map, axis = 1).hide(axis = 'index')
    rama_list = rama_list_styled.to_html(index = False, justify = "center", dtable_id = 'results_table')
    rama_list = rama_list.replace('<table', '<table style="text-align:center; font-family: sans-serif;"')   
    
    if results_file_type == "all":
        html_path = os.path.join(file_path, "rama_results.html")
    
    elif results_file_type == "top":
        html_path = os.path.join(file_path, "rama_results_top.html")
        
    with open(html_path, 'w') as html_output:   
        html_output.write(rama_list)

def retrievePDB(pdb_code, file_path):
       
    pdb_code = pdb_code.lower()
    file_name = pdb_code + ".cif"
    getPDB = PDB.PDBList(pdb = file_path, obsolete_pdb = file_path)
    downloaded_file = getPDB.retrieve_pdb_file(pdir = file_path, pdb_code = pdb_code, file_format = "mmCif")
    
    if not os.path.isfile(downloaded_file):
        return("retrieveFailed")
        
    newName = os.path.join(file_path, file_name)
    os.rename(downloaded_file, newName)    
    
    return file_name

def reorder_pdb(pdb_file_name, file_path, MPNN_scores, run_gly):
   
    pdb_file = os.path.join(file_path, pdb_file_name[:-4] + '.pdb')
    ordered_pdb = os.path.join(file_path, 'sorted_struct.pdb')
    
    with open(pdb_file, 'r') as infile:
        atom_lines = [line for line in infile if line.startswith('ATOM') or line.startswith('HETATM')]

    # Sort atom lines based on chain ID and residue number
    sorted_atom_lines = sorted(atom_lines, key=lambda line: (line[21], int(line[22:26])))

    # Write the reordered ATOM lines to the output PDB file
    with open(ordered_pdb, 'w') as output_file:
        cryst_record = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          "
        output_file.write("HEADER" + "\n")
        output_file.write(cryst_record + "\n")
        output_file.writelines(sorted_atom_lines)
    
    apply_b_factor(file_path, MPNN_scores, run_gly)

def apply_b_factor(file_path, MPNN_scores, run_gly):

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('sorted_struct', os.path.join(file_path, 'sorted_struct.pdb'))
    curr_res = 0    
    model = structure[0]
    for chain in model:
        for residue in chain:
            if PDB.is_aa(residue):
                if residue.get_resname() == "PRO" and run_gly != "True":
                    curr_MPNN = float(101)
                elif residue.get_resname() == "GLY" and run_gly == "True":
                    curr_MPNN = float(101)    
                else:
                    try:
                        curr_MPNN = float(MPNN_scores[curr_res]) * 100 
                    except IndexError:
                        error_file = os.path.join(file_path,'error.txt')
                        with open(error_file,'w') as error:
                            error.write("ProteinMPNN Error")
                for atom in residue:
                    if MPNN_scores[curr_res] == '?':
                        atom.set_bfactor(0)
                    else:
                        atom.set_bfactor(curr_MPNN)
                curr_res += 1
            else:
                for atom in residue:
                    atom.set_bfactor(0)

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(os.path.join(file_path, 'sorted_struct.pdb'))

def get_secondary_struct(pdb_file_name,res_info_dict):
    
    sec_codes = {'E':'B_Strand', 'H':'A_Helix', 'T':"Turn", 'G':'3-10_Helix', 'I':'Pi Helix',
                'B':'B-Bridge', 'S':'Bend', 'P':'Poly_Proline', '-':'No_Sec', ' ':'No_Sec', 'C':'No_Sec'}

    sec_dict = {}
    dssp_file_name = pdb_file_name[:-4] + ".dssp"
    pdb_file_name = pdb_file_name[:-4] + '.pdb'
    subprocess.run(["/www/cgi-bin/ProScan/bin/mkdssp", pdb_file_name, dssp_file_name])     
    
    if not os.path.isfile(dssp_file_name):
        return sec_list #If dssp fails to parse pdb/cif return empty list
    
    with open(dssp_file_name, 'r') as dssp_file:
    
        #Skip over lines until residue information starts    
        for line in dssp_file:       
            if line.startswith("  #"):
                break
        
        #Append residue secondary codes to secList and skip chain breaks (!)
        for curr_res in dssp_file:              
            if curr_res[13] == "!":        
                continue 
            
            curr_id = curr_res[6:10] 
            curr_chain = curr_res[11] 
            curr_sec = sec_codes[curr_res[16]]
            curr_key = str(curr_chain).strip() + str(curr_id).strip()
            sec_dict[curr_key] = curr_sec                          
            
    for key in res_info_dict:
        if key in sec_dict:
            res_info_dict[key].append(sec_dict[key])
        else:
            res_info_dict[key].append('?')

    return res_info_dict

def run_proteinMPNN(file_name, file_path, run_gly):
  
    # Get directory of proteinMPNN helper script and main script
    run_script_dir = '/www/cgi-bin/ProScan/ProteinMPNN'
    help_script_dir = os.path.join(run_script_dir, 'helper_scripts')
    output_path = file_path + "/parsed_pdbs.jsonl"
    
    # Run proteinMPNN
    subprocess.call(['/www/cgi-bin/miniconda3/envs/wsgi/bin/python', 'parse_multiple_chains.py', '--input_path',file_path
    , '--output_path', output_path, '--seq_path', file_path], cwd=help_script_dir)
    subprocess.call(['/www/cgi-bin/miniconda3/envs/wsgi/bin/python', 'protein_mpnn_run.py', '--jsonl_path', output_path, 
    '--out_folder', file_path, '--num_seq_per_target', '1', '--sampling_temp', '0.2', '--seed', '37', '--batch_size', '1', '--save_probs', '1'], cwd=run_script_dir)   
    prob_list = parse_aa_probs(file_name,file_path, run_gly)
    return prob_list

def parse_aa_probs(pdb_file_name, file_path, run_gly):

    prob_list = []
    all_probs = []
    
    MPNN_seq = ''  #Read in sequence parsed out by MPNN 
    with open(os.path.join(file_path, 'MPNN_seq.txt'), 'r') as parsed_seq:      
        MPNN_seq = parsed_seq.readline()
    
    data_dir = os.path.join(file_path, "probs")
    data_file = os.path.join(data_dir, pdb_file_name[:-4] + ".npz")  
         
    if os.path.isfile(data_file) == 0:
        return prob_list #If MPNN fails to parse file return empty list
    
    all_data = np.load(data_file)
    prob_data = all_data["log_probs"]
    np.savetxt(os.path.join(file_path, "aaProbData.txt"), np.exp(prob_data).mean(0), fmt="%.5f")    
    
    with open(os.path.join(file_path, "aaProbData.txt"), "r") as prob_file:
        #Read in proline probabilities
        for line in prob_file:
            curr_line = line.split(" ")
            all_probs.append(curr_line)
            if run_gly == 'True':
                prob_list += (line[41:48]),
            else:
                prob_list += (line[97:104]),
    
    #Ignore probabilities of gap residues
    try:
        filtered_prob_list = [prob_list[i] for i in range(len(MPNN_seq)) if MPNN_seq[i] != '-']
        return filtered_prob_list
    except IndexError:
        error_file = os.path.join(file_path,'error.txt')
        with open(error_file,'w') as error:
            error.write("ProteinMPNN Error")

#If user requests chains to be exclude, removes specified chains from pdb
def remove_chains(pdb_file_name, file_path):

    pdb_file = os.path.join(file_path, pdb_file_name)
    with open(os.path.join(file_path, 'chain_info.txt'), 'r') as chain_file:
        chain_list = chain_file.readline().strip()
        chain_list_type = chain_file.readline().strip()
      
    if chain_list == 'EMPTY':
        chain_list = []
    else:
        chain_list = chain_list.split(',')

    if pdb_file_name[-4:] == ".pdb":
        
        structure = PDB.PDBParser().get_structure("input_structure", "%s" % pdb_file)
          
        if len(chain_list) != 0: #If user did not put any chains, avoid detaching any
            chains = list(structure.get_chains())
            if chain_list_type == "chainIgnoreList":    
                for chain in chains:
                    if chain.get_id() in chain_list:
                        structure[0].detach_child(chain.get_id())
        
            elif chain_list_type == "chainOnlyList":
                pdb_chains = []
                for chain in chains:
                    pdb_chains.append(chain.get_id())
                check_chain_input = bool(set(chain_list) & set(pdb_chains)) #Make sure user specified a chain that actually exist so entire pdb/cif is not deleted
                if check_chain_input == True:
                    for chain in chains:
                        if chain.get_id() not in chain_list:
                            structure[0].detach_child(chain.get_id())
        
        io = PDB.PDBIO()      
        with open(os.path.join(file_path,"tempStruct.pdb"),"w") as f:
            
            #Write random cryst1 record line required for dssp functionality
            cryst_record = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          "
            f.write("HEADER" + "\n")
            f.write(cryst_record + "\n")
            io.set_structure(structure)
            io.save(f,preserve_atom_numbering=True)

        os.replace(os.path.join(file_path, "tempStruct.pdb"), pdb_file)
        return None
    
    elif pdb_file_name[-4:] == ".cif":
    
        structure = PDB.MMCIFParser().get_structure("input_structure", "%s" % pdb_file)       
        
        chains = list(structure.get_chains())    
        
        if len(chain_list) != 0: 
            if chain_list_type == "chainIgnoreList":    
                for chain in chains:
                    if chain.get_id() in chain_list:
                        structure[0].detach_child(chain.get_id())
        
            elif chain_list_type == "chainOnlyList":
                pdb_chains = []
                for chain in chains:
                    pdb_chains.append(chain.get_id())
                check_chain_input = bool(set(chain_list) & set(pdb_chains)) #Make sure user specified a chain that actually exist so entire pdb/cif is not deleted
                if check_chain_input == True:
                    for chain in chains:
                        if chain.get_id() not in chain_list:
                            structure[0].detach_child(chain.get_id())

        io = PDB.PDBIO()
        with open(os.path.join(file_path,"tempStruct.pdb"),"w") as f:
          
           #Write random cryst1 record line required for dssp functionality
           cryst_record = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          "
           f.write("HEADER" + "\n")
           f.write(cryst_record + "\n")
           io.set_structure(structure)
           io.save(f,preserve_atom_numbering=True)
        shutil.copy(pdb_file,os.path.join(file_path, 'originalPDB.cif'))
        os.replace(os.path.join(file_path, "tempStruct.pdb"), os.path.join(file_path, pdb_file_name[:-4] + ".pdb"))
    
    return None

def write_high_scoring(file_path, job_name, file_name_top):

    data_path = os.path.join(file_path, job_name)
    filtered_path=os.path.join(file_path, file_name_top)
    rama_df = pd.read_csv(data_path, sep='\t', header=0)
    
    rama_df['ProteinMPNN_Prob'] = rama_df['ProteinMPNN_Prob'].astype(float)

    if 'Pro_Angle' in rama_df.columns:  
        best_hits_df = rama_df[(rama_df['ProteinMPNN_Prob'] >= ACCEPTABLE_CUTOFF_MPNN) & (rama_df['Pro_Angle'] != 'No_Angle') & (rama_df['Pro_Angle'] != 'Questionable') & (rama_df['Residue'] != 'PRO')]
        best_hits_df.to_csv(filtered_path, sep='\t', index=False)

    elif 'Gly_Angle' in rama_df.columns:
        best_hits_df = rama_df[(rama_df['ProteinMPNN_Prob'] >= ACCEPTABLE_CUTOFF_MPNN) & (rama_df['Gly_Angle'] != 'No_Angle') & (rama_df['Gly_Angle'] != 'Questionable') & (rama_df['Residue'] != 'GLY')]
        best_hits_df.to_csv(filtered_path, sep='\t', index=False)
    
    return None

def create_rama_plot(file_path, file_name):

    PLOTLY_PREF_COLOR = '#D2042D'
    PLOTLY_ACCEPT_COLOR = '#FAFA33'
    WT_res = []; res_nums = []; phi_angles = []; psi_angles = []; MPNN_scores = []; chains = []; primary_angles = []

    with open(os.path.join(file_path,file_name),'r') as results_file:
          results_file.readline()
          for line in results_file:   
            line_split = line.split('\t')
            if float(line_split[3]) != 0 and float(line_split[4]) != 0:
                WT_res.append(line_split[0])
                res_nums.append(line_split[1])
                chains.append(line_split[2])
                phi_angles.append(float(line_split[3]))
                psi_angles.append(float(line_split[4]))
                MPNN_scores.append(float(line_split[6]))
                primary_angles.append(line_split[7])
    #Store info that will displayed on hover in pandas dict
    all_results_df = pd.DataFrame({'res_num':res_nums, 'WT_res':WT_res, 'phi_angles':phi_angles, 'psi_angles':psi_angles, 'MPNN_scores':MPNN_scores, 'Chains':chains, 'primary_angles':primary_angles})
    
    non_pro_df = all_results_df[((all_results_df['WT_res'] != 'PRO') & (all_results_df['MPNN_scores'] < ACCEPTABLE_CUTOFF_MPNN)) | ((all_results_df['MPNN_scores'] > OPTIMAL_CUTOFF_MPNN) & (all_results_df['primary_angles'] == 'Questionable'))]
    medium_only_df = all_results_df[(all_results_df['WT_res'] != 'PRO') & (all_results_df['MPNN_scores'] >= ACCEPTABLE_CUTOFF_MPNN) & (all_results_df['MPNN_scores'] < OPTIMAL_CUTOFF_MPNN)]
    high_only_df = all_results_df[(all_results_df['WT_res'] != 'PRO') & (all_results_df['MPNN_scores'] >= OPTIMAL_CUTOFF_MPNN)  & (all_results_df['primary_angles'] != 'Questionable')]
    only_pro_df = all_results_df[all_results_df['WT_res'] == 'PRO']

    #Create seperate scatterplot traces for WT prolines, low scoring, medium scoring, and high scoring residues. Set information displayed on hover    
    scatter_plot_nopro = go.Scatter(x = non_pro_df['phi_angles'], y = non_pro_df['psi_angles'],
                        mode = 'markers', showlegend = False, marker = dict(color = "#abb2b9", size = 7, line = dict(color = 'black', width = 1)),
                        customdata = non_pro_df[['res_num', 'WT_res', 'MPNN_scores', 'Chains']],
                        hovertemplate = 'Residue Number: %{customdata[0]}<br>Chain ID: %{customdata[3]}<br>WT Residue: %{customdata[1]}<br>MPNN Score: %{customdata[2]}<extra></extra>',
                        hoverlabel = dict(bgcolor = "white", font_color = "black", bordercolor = "#000000"))

    scatter_plot_pro = go.Scatter(x = only_pro_df['phi_angles'], y = only_pro_df['psi_angles'],
                        mode = 'markers', showlegend = False, marker = dict(color = PRO_GRAPH_COLOR, size = 7, line = dict(color = 'black', width = 1)),
                        customdata = only_pro_df[['res_num', 'WT_res', 'MPNN_scores', 'Chains']],
                        hovertemplate = 'Residue Number: %{customdata[0]}<br>Chain ID: %{customdata[3]}<br>WT Residue: %{customdata[1]}<br>MPNN Score: %{customdata[2]}<extra></extra>',
                        hoverlabel = dict(bgcolor = "white", font_color = "black", bordercolor = PRO_GRAPH_COLOR))

    scatter_plot_medium = go.Scatter(x = medium_only_df['phi_angles'], y = medium_only_df['psi_angles'],
                        mode = 'markers', showlegend = False, marker = dict(color = ACCEPTABLE_COLOR, size = 7, line = dict(color = 'black', width = 1)),
                        customdata = medium_only_df[['res_num', 'WT_res', 'MPNN_scores', 'Chains']],
                        hovertemplate = 'Residue Number: %{customdata[0]}<br>Chain ID: %{customdata[3]}<br>WT Residue: %{customdata[1]}<br>MPNN Score: %{customdata[2]}<extra></extra>',
                        hoverlabel = dict(bgcolor = "white", font_color = "black", bordercolor = ACCEPTABLE_COLOR))

    scatter_plot_high = go.Scatter(x = high_only_df['phi_angles'], y = high_only_df['psi_angles'],
                        mode = 'markers', showlegend = False, marker = dict(color = OPTIMAL_COLOR, size = 7, line = dict(color = 'black', width = 1)),
                        customdata = high_only_df[['res_num', 'WT_res', 'MPNN_scores', 'Chains']],
                        hovertemplate = 'Residue Number: %{customdata[0]}<br>Chain ID: %{customdata[3]}<br>WT Residue: %{customdata[1]}<br>MPNN Score: %{customdata[2]}<extra></extra>',
                        hoverlabel = dict(bgcolor = "white", font_color = "black", bordercolor = OPTIMAL_COLOR))
 
    #Create graph quadrants
    shapes = [dict(type='line', x0=-180, x1=180, y0=0, y1=0, line=dict(dash='dash')),
    dict(type='line', x0=0, x1=0, y0=-180, y1=180, line=dict(dash='dash'))]

    #Create figure, fix axis range, apply quadrants, set background to transparent    
    rama_fig = go.Figure(data=[scatter_plot_nopro,scatter_plot_pro,scatter_plot_medium,scatter_plot_high],
                         layout=dict(yaxis_range=[-180, 180], xaxis_range=[-180, 180],
                         xaxis=dict(showgrid=False,fixedrange=True,dtick=60), yaxis=dict(showgrid=False,fixedrange=True,dtick=60),
                         shapes=shapes, plot_bgcolor='white', margin=dict(l=40, r=40, t=50, b=40),xaxis_title='Phi',yaxis_title='Psi'))
    
    pro_image = '/www/ProScan/pro_image_plt_greyscale.png'
    pre_image = '/www/ProScan/pre_image_plt_greyscale.png'
    gly_image = '/www/ProScan/gly_image_plt_greyscale.png'
    pro_plot_image = base64.b64encode(open(pro_image, 'rb').read())
    pre_plot_image = base64.b64encode(open(pre_image, 'rb').read())
    gly_plot_image = base64.b64encode(open(gly_image, 'rb').read())

    rama_fig.add_layout_image(
    dict(source='data:image/png;base64,{}'.format(pro_plot_image.decode()),xref="paper",yref="paper",
        x=0,y=1,sizex=1,sizey=1,sizing="contain",opacity=1,layer="below"))
        
    rama_fig.add_layout_image(
    dict(source='data:image/png;base64,{}'.format(pre_plot_image.decode()),xref="paper",yref="paper",
        x=0,y=1,sizex=1,sizey=1,sizing="contain",opacity=1,visible=False,layer="below"))
    
    rama_fig.add_layout_image(
    dict(source='data:image/png;base64,{}'.format(gly_plot_image.decode()),xref="paper",yref="paper",
        x=0,y=1,sizex=1,sizey=1,sizing="contain",opacity=1,visible=False,layer="below"))

    scatter_options = [{'label': 'Show All Residues', 'method': 'update', 'args': [{'visible': [True, True, True, True]}]},
                       {'label': 'Hide WT Prolines', 'method': 'update', 'args': [{'visible': [True, False,True, True]}]},
                       {'label': 'Hide Low Scoring', 'method':'update','args':[{'visible':[False,False,True,True]}]}]
    
    # Add dropdown menu to switch between contour plots
    scatter_button = go.layout.Updatemenu(
    type='dropdown',
    showactive=True,
    buttons=scatter_options,
    x=0.25,
    xanchor='left',
    y=1.15,
    yanchor='top',
    )

    contour_buttons = [
    dict(method="relayout",
        args=[{"images[2].visible": False, "images[1].visible": False, "images[0].visible":True},
              {"title": "Proline"}],
        label="Proline"
    ),
    dict(
        method="relayout",
        args=[{"images[2].visible": False, "images[1].visible": True, "images[0].visible":False},
              {"title": "Pre-Proline"}],
        label="Pre-Proline"
    ),
    dict(
        method="relayout",
        args=[{"images[2].visible": True, "images[1].visible": False, "images[0].visible":False},
              {"title": "Glycine"}],
        label="Glycine"
    )
    ]

    updatemenu = go.layout.Updatemenu(type="dropdown",showactive=True,buttons=contour_buttons,
                                    x=-0.1,xanchor='left',y=1.15,yanchor='top',)

    rama_fig.update_layout(updatemenus=[updatemenu,scatter_button])    
    
    json_data = pio.to_json(rama_fig)   
    json_path = os.path.join(file_path,'rama_plot.json')   
    with open(json_path,'w') as rama_plot:
        rama_plot.write(json_data)
       
    return None

def create_rama_plot_gly(file_path, file_name):

    PLOTLY_PREF_COLOR = '#D2042D'
    PLOTLY_ACCEPT_COLOR = '#FAFA33'
    WT_res = []; res_nums = []; phi_angles = []; psi_angles = []; MPNN_scores = []; chains = []; primary_angles = []

    with open(os.path.join(file_path,file_name),'r') as results_file:
          results_file.readline()
          for line in results_file:   
            line_split = line.split('\t')
            if float(line_split[3]) != 0 and float(line_split[4]) != 0:
                WT_res.append(line_split[0])
                res_nums.append(line_split[1])
                chains.append(line_split[2])
                phi_angles.append(float(line_split[3]))
                psi_angles.append(float(line_split[4]))
                MPNN_scores.append(float(line_split[6]))
                primary_angles.append(line_split[7])
    #Store info that will displayed on hover in pandas dict
    all_results_df = pd.DataFrame({'res_num':res_nums, 'WT_res':WT_res, 'phi_angles':phi_angles, 'psi_angles':psi_angles, 'MPNN_scores':MPNN_scores, 'Chains':chains, 'primary_angles':primary_angles})
    
    non_gly_df = all_results_df[((all_results_df['WT_res'] != 'GLY') & (all_results_df['MPNN_scores'] < ACCEPTABLE_CUTOFF_MPNN)) | ((all_results_df['MPNN_scores'] > OPTIMAL_CUTOFF_MPNN) & (all_results_df['primary_angles'] == 'Questionable'))]
    medium_only_df = all_results_df[(all_results_df['WT_res'] != 'GLY') & (all_results_df['MPNN_scores'] >= ACCEPTABLE_CUTOFF_MPNN) & (all_results_df['MPNN_scores'] < OPTIMAL_CUTOFF_MPNN)]
    high_only_df = all_results_df[(all_results_df['WT_res'] != 'GLY') & (all_results_df['MPNN_scores'] >= OPTIMAL_CUTOFF_MPNN)  & (all_results_df['primary_angles'] != 'Questionable')]
    only_gly_df = all_results_df[all_results_df['WT_res'] == 'GLY']

    #Create seperate scatterplot traces for WT prolines, low scoring, medium scoring, and high scoring residues. Set information displayed on hover    
    scatter_plot_nopro = go.Scatter(x = non_gly_df['phi_angles'], y = non_gly_df['psi_angles'],
                        mode = 'markers', showlegend = False, marker = dict(color = "#abb2b9", size = 7, line = dict(color = 'black', width = 1)),
                        customdata = non_gly_df[['res_num', 'WT_res', 'MPNN_scores', 'Chains']],
                        hovertemplate = 'Residue Number: %{customdata[0]}<br>Chain ID: %{customdata[3]}<br>WT Residue: %{customdata[1]}<br>MPNN Score: %{customdata[2]}<extra></extra>',
                        hoverlabel = dict(bgcolor = "white", font_color = "black", bordercolor = "#000000"))

    scatter_plot_pro = go.Scatter(x = only_gly_df['phi_angles'], y = only_gly_df['psi_angles'],
                        mode = 'markers', showlegend = False, marker = dict(color = PRO_GRAPH_COLOR, size = 7, line = dict(color = 'black', width = 1)),
                        customdata = only_gly_df[['res_num', 'WT_res', 'MPNN_scores', 'Chains']],
                        hovertemplate = 'Residue Number: %{customdata[0]}<br>Chain ID: %{customdata[3]}<br>WT Residue: %{customdata[1]}<br>MPNN Score: %{customdata[2]}<extra></extra>',
                        hoverlabel = dict(bgcolor = "white", font_color = "black", bordercolor = PRO_GRAPH_COLOR))

    scatter_plot_medium = go.Scatter(x = medium_only_df['phi_angles'], y = medium_only_df['psi_angles'],
                        mode = 'markers', showlegend = False, marker = dict(color = ACCEPTABLE_COLOR, size = 7, line = dict(color = 'black', width = 1)),
                        customdata = medium_only_df[['res_num', 'WT_res', 'MPNN_scores', 'Chains']],
                        hovertemplate = 'Residue Number: %{customdata[0]}<br>Chain ID: %{customdata[3]}<br>WT Residue: %{customdata[1]}<br>MPNN Score: %{customdata[2]}<extra></extra>',
                        hoverlabel = dict(bgcolor = "white", font_color = "black", bordercolor = ACCEPTABLE_COLOR))

    scatter_plot_high = go.Scatter(x = high_only_df['phi_angles'], y = high_only_df['psi_angles'],
                        mode = 'markers', showlegend = False, marker = dict(color = OPTIMAL_COLOR, size = 7, line = dict(color = 'black', width = 1)),
                        customdata = high_only_df[['res_num', 'WT_res', 'MPNN_scores', 'Chains']],
                        hovertemplate = 'Residue Number: %{customdata[0]}<br>Chain ID: %{customdata[3]}<br>WT Residue: %{customdata[1]}<br>MPNN Score: %{customdata[2]}<extra></extra>',
                        hoverlabel = dict(bgcolor = "white", font_color = "black", bordercolor = OPTIMAL_COLOR))
 
    #Create graph quadrants
    shapes = [dict(type='line', x0=-180, x1=180, y0=0, y1=0, line=dict(dash='dash')),
    dict(type='line', x0=0, x1=0, y0=-180, y1=180, line=dict(dash='dash'))]

    #Create figure, fix axis range, apply quadrants, set background to transparent    
    rama_fig = go.Figure(data=[scatter_plot_nopro,scatter_plot_pro,scatter_plot_medium,scatter_plot_high],
                         layout=dict(yaxis_range=[-180, 180], xaxis_range=[-180, 180],
                         xaxis=dict(showgrid=False,fixedrange=True,dtick=60), yaxis=dict(showgrid=False,fixedrange=True,dtick=60),
                         shapes=shapes, plot_bgcolor='white', margin=dict(l=40, r=40, t=50, b=40),xaxis_title='Phi',yaxis_title='Psi'))
    
    pro_image = '/www/ProScan/pro_image_plt_greyscale.png'
    pre_image = '/www/ProScan/pre_image_plt_greyscale.png'
    gly_image = '/www/ProScan/gly_image_plt_greyscale.png'
    pro_plot_image = base64.b64encode(open(pro_image, 'rb').read())
    pre_plot_image = base64.b64encode(open(pre_image, 'rb').read())
    gly_plot_image = base64.b64encode(open(gly_image, 'rb').read())

    rama_fig.add_layout_image(
    dict(source='data:image/png;base64,{}'.format(gly_plot_image.decode()),xref="paper",yref="paper",
        x=0,y=1,sizex=1,sizey=1,sizing="fill",opacity=1,layer="below"))

    rama_fig.add_layout_image(
    dict(source='data:image/png;base64,{}'.format(pro_plot_image.decode()),xref="paper",yref="paper",
        x=0,y=1,sizex=1,sizey=1,sizing="fill",opacity=1,visible=False,layer="below"))
        
    rama_fig.add_layout_image(
    dict(source='data:image/png;base64,{}'.format(pre_plot_image.decode()),xref="paper",yref="paper",
        x=0,y=1,sizex=1,sizey=1,sizing="fill",opacity=1,visible=False,layer="below"))
    
    scatter_options = [{'label': 'Show All', 'method': 'update', 'args': [{'visible': [True, True, True, True]}]},
                       {'label': 'Hide WT Gly', 'method': 'update', 'args': [{'visible': [True, False,True, True]}]},
                       {'label': 'Hide Low Scoring', 'method':'update','args':[{'visible':[False,False,True,True]}]}]
    
    # Add dropdown menu to switch between contour plots
    scatter_button = go.layout.Updatemenu(
    type='dropdown',
    showactive=True,
    buttons=scatter_options,
    x=0.25,
    xanchor='left',
    y=1.15,
    yanchor='top',
    )

    contour_buttons = [
    dict(
        method="relayout",
        args=[{"images[0].visible": True, "images[1].visible": False, "images[2].visible":False},
              {"title": "Glycine"}],
        label="Glycine"
    ),
    dict(
        method="relayout",
        args=[{"images[0].visible": False, "images[1].visible": True, "images[2].visible":False},
              {"title": "Proline"}],
        label="Proline"
    ),
    dict(
        method="relayout",
        args=[{"images[0].visible": False, "images[1].visible": False, "images[2].visible":True},
              {"title": "Pre-Proline"}],
        label="Pre-Proline"
    )
    ]

    updatemenu = go.layout.Updatemenu(type="dropdown",showactive=True,buttons=contour_buttons,
                                    x=-0.1,xanchor='left',y=1.15,yanchor='top',)

    rama_fig.update_layout(updatemenus=[updatemenu,scatter_button],autosize=True)    
    
    json_data = pio.to_json(rama_fig)   
    json_path = os.path.join(file_path,'rama_plot.json')   
    with open(json_path,'w') as rama_plot:
        rama_plot.write(json_data)
       
    return None

def get_hbonds(file_path, file_name, run_gly):

    pdb_file = os.path.join(file_path,file_name[:-4]+'.pdb')
    disulfide_lines = []
    ### See if there is better way to do this ###
    curr_dirr = os.getcwd()
    os.chdir(file_path)
    hbplus_cmd = ["/www/cgi-bin/ProScan/hbplus/hbplus",pdb_file] #Go into hbplus directory and execute f/ there
    process = subprocess.Popen(hbplus_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
   
    #For some reason hbplus will not output disulfide bonds to .hb2 file. So capture them from output
    for line in iter(process.stdout.readline, ''):
    
        line = line.strip()
        if line.startswith("Disulphide"):       
            disulfide_lines.append(line)
    
    ### See if there is better way to do this ###
    os.chdir(curr_dirr) #Switch back to ProScan root directory 

    disulfide_list = parse_disulfides(disulfide_lines)
    filter_hbonds(file_path, file_name, disulfide_list, run_gly)
    return None

def get_glycans(file_path, file_name):

    model_file = os.path.join(file_path,'originalPDB.cif')
    glycan_list = []
    with open(model_file,'r') as cif_file:

        for line in cif_file:
            if line[:6] == 'covale' and 'NAG' in line:
                if line[68:71] == 'ASN': 
                    glycan_list.append([line[66],str(int(line[72:76].zfill(4))),'Glycan_NAG'])
    
    return glycan_list

def parse_disulfides(disulfide_lines):

    disulfide_list = []

    for line in disulfide_lines:
        disulfide_list.append([line[25],line[26:30],'Disulfide_' + line[44:49]])
        disulfide_list.append([line[44],line[45:49],'Disulfide_' + line[25:30]])

    return disulfide_list

def filter_hbonds(file_path, file_name, disulfide_list, run_gly):
    
    if file_name[-4:] == '.cif':
        glycan_list = get_glycans(file_path, file_name)
    
    hbond_file_name = file_name[:-4] + '.hb2'
    hbond_file = os.path.join(file_path, hbond_file_name)
    hbond_filtered_file = os.path.join(file_path, 'effected_hbonds.txt')
    hbond_list = []

    with open(hbond_file,'r') as hbond_predictions:

        for header_lines in range(8): ##Skip header lines
            next(hbond_predictions)
        
        #Ignore solvent H-bonds (!= HOH). If side-chain (SC) or N-atom main chain (MC) h-bonds are interrupted, append to list
        for res_line in hbond_predictions:
            if res_line[33] == 'S' and res_line[21:24] != 'HOH' and res_line[0:5] != res_line[14:19] and res_line[33] != 'H' and res_line[34] != 'H': 
                hbond_list.append([res_line[0], res_line[1:5],'H-bond_' + res_line[14] + res_line[15:19]])
            if res_line[34] == 'S' and res_line[6:9] != 'HOH' and res_line[0:5] != res_line[14:19] and res_line[33] != 'H' and res_line[34] != 'H':
                hbond_list.append([res_line[14], res_line[15:19],'H-bond_' + res_line[0] + res_line[1:5]])
            if res_line[33] == 'M' and res_line[10] == 'N' and res_line[21:24] != 'HOH' and res_line[0:5] != res_line[14:19] and res_line[33] != 'H' and res_line[34] != 'H':
                hbond_list.append([res_line[0], res_line[1:5], 'H-bond_'+ res_line[14] + res_line[15:19]])
            else:
                continue

    if file_name[-4:] == '.cif' and len(glycan_list) != 0:
        for glycan in glycan_list:
            hbond_list.append(glycan)
    
    for dis_bond in disulfide_list: #Disulfides and h-bonds are stored as a list [Chain ID, Res Num, Bond information]
        hbond_list.append(dis_bond)

    with open(hbond_filtered_file,'w') as hbond_predictions:
        for residue in hbond_list: #Write bonds to a text dock to be parsed out
            res_line = '\t'.join(map(str, residue))
            hbond_predictions.write(res_line + '\n')
        
    return None

def append_disrupted_hbonds(file_path, run_data_dict):

    hbond_filtered = os.path.join(file_path, 'effected_hbonds.txt')
    hbond_dict = {}

    with open(hbond_filtered,'r') as hbond_predictions:

        for res_line in hbond_predictions:
            res_line_split = res_line.split('\t')
            res_dict_id = res_line_split[0] + str(int(res_line_split[1])) #Individual residues (ID'ed by chain + Res_num) are keys of dictionary. Each key holds a list with their bond info
            if res_dict_id in hbond_dict:
                if res_line_split[2].strip() not in hbond_dict[res_dict_id]:
                    hbond_dict[res_dict_id].append(res_line_split[2].strip())
            else:
                hbond_dict[res_dict_id] = []
                hbond_dict[res_dict_id].append(res_line_split[2].strip())
    
    for key in run_data_dict:
              
        if key in hbond_dict:
            if len(hbond_dict[key]) == 1:     
               run_data_dict[key].append(hbond_dict[key][0])
            else:
                note_string = " ".join(hbond_dict[key])
                run_data_dict[key].append(note_string)
        else:
            run_data_dict[key].extend([str('-')])
        
    return run_data_dict


def clean_runs():
    
    run_directory = os.listdir('/www/ProScan/runs')

    for run in run_directory:
        entry_path = os.path.join('/www/ProScan/runs', run)
    
        if os.path.isdir(entry_path):
            try:
                shutil.rmtree(entry_path)
                print(f"Deleted directory and its contents: {entry_path}")
            except Exception as e:
                print(f"Error deleting directory {entry_path}: {e}")


def make_mut_file(pdb_file_name, file_path, run_gly):

    aa_code_dict = {'ALA': 'A','ARG': 'R','ASN': 'N','ASP': 'D','CYS': 'C','GLN': 'Q',
            'GLU': 'E','GLY': 'G','HIS': 'H','ILE': 'I','LEU': 'L','LYS': 'K','MET': 'M',
    'PHE': 'F','PRO': 'P','SER': 'S','THR': 'T','TRP': 'W','TYR': 'Y','VAL': 'V',}

    parser = PDB.PDBParser()
    structure = parser.get_structure('structure', os.path.join(file_path,'sorted_struct.pdb'))
    residue_info_list = []
    seen_residues = set()
    for model in structure:
        for chain in model:
            chain_id = chain.id
            for residue in chain:
                residue_id = residue.id[1]
                residue_type = residue.get_resname()
                residue_info = [chain_id, residue_id, residue_type]
                if tuple(residue_info) not in seen_residues:
                    if residue_type in aa_code_dict:
                        residue_info_list.append(residue_info)
                        seen_residues.add(tuple(residue_info))
    
    num_residues = len(residue_info_list)
    base_file_name = 'sorted_struct.muts.txt'
    mut_file_path = os.path.join(file_path,base_file_name)
    
    with open(mut_file_path,'w') as mut_file:
        
        mut_file.write('START' + '\n')
        if num_residues < 10:
            mut_file.write('   ' + str(num_residues) + '\n')
        elif num_residues >= 10 and num_residues < 1000:
            mut_file.write('  ' + str(num_residues) + '\n')
        else:
            mut_file.write(' ' + str(num_residues) + '\n')
        
        if run_gly == 'True':
            mut_res = 'G'
        else:
            mut_res = 'P'

        for residue in residue_info_list:
            mut_file.write('  ' + str(1) + '\n')
            if int(residue[1]) < 10:
                mut_line = 'MUTATIONS   ' + str(residue[1]) + '  ' + residue[0] + '  ' + aa_code_dict[residue[2]] + '  ' + mut_res + '\n'  
            elif int(residue[1]) >= 10 and int(residue[1]) < 100:
                mut_line = 'MUTATIONS  ' + str(residue[1]) + '  ' + residue[0] + '  ' + aa_code_dict[residue[2]] + '  ' + mut_res + '\n'  
            elif int(residue[1]) >= 100 and int(residue[1]) < 1000:
                mut_line = 'MUTATIONS ' + str(residue[1]) + '  ' + residue[0] + '  ' + aa_code_dict[residue[2]] + '  ' + mut_res + '\n'
            else:
                mut_line = 'MUTATIONS' + str(residue[1]) + '  ' + residue[0] + '  ' + aa_code_dict[residue[2]] + '  ' + mut_res + '\n'

            mut_file.write(mut_line)

    mut_track_file = os.path.join(file_path, 'mut_track.txt')

    with open (mut_track_file,'a') as mut_tracker:

        mut_tracker.write(str(num_residues) + '\n')
    
def execute_ros(pdb_file_name,file_path,run_gly):
    
    home_path = os.path.abspath(os.getcwd())
    os.chdir(file_path)
    make_mut_file(pdb_file_name,file_path, run_gly)
   
    
    path_dest = os.path.join(file_path,'paths.txt')
    shutil.copyfile('/www/cgi-bin/ProScan/paths.txt',path_dest)
    
    ros_bin = '/www/cgi-bin/ProScan/rosetta2.3/rosetta++/rosetta.gcc34'
    code = 'aa'
    pdb = 'sorted_struct'
    out = pdb + '.ddg.ros.out'
    mut = pdb + '.muts.txt'
    seed = str(12)
    ros_cmd = [
    ros_bin,
    code,
    pdb,
    "_",
    "-interface",
    "-intout",
    out,
    "-monomeric_protein",
    "-ignore_unrecognized_res",
    "-safety_check",
    "-skip_missing_residues",
    "-mutlist",
    mut,
    "-constant_seed",
    "-jran",
    seed,
    "-yap",
    "-s",
    pdb
    ]

    ros_cmd_str = " ".join(ros_cmd)
    ddg_parse_script = '/www/cgi-bin/ProScan/get_ddg.pl'
    with subprocess.Popen(ros_cmd_str, stdout=subprocess.PIPE, shell=True, text=True, bufsize=1, universal_newlines=True) as process: #Clear buffer in short intervals so number of mutants is updated more frequently
        for line in process.stdout:
            if "examining" in line:
                # Append to file if needed
                with open("mut_track.txt", "a") as file:
                    file.write(line)
    ddg_output = subprocess.run(['perl', ddg_parse_script, 'sorted_struct.ddg.ros.out'], capture_output=True, text=True)
    ros_ddg_output_file = os.path.join(file_path,'ros_ddg.txt')

    with open (ros_ddg_output_file,'w') as ros_output:
        ros_output.write(ddg_output.stdout)

    os.chdir(home_path)
    
def parse_ros(file_path,res_info_dict):

    ddg_file = os.path.join(file_path,"ros_ddg.txt")
    ddg_vals = {}

    with open(ddg_file,'r') as ddg_values:

        for line in ddg_values:
            line_split = line.split(',')
            ddg_res = line_split[1].strip()
            ddg_chain = line_split[3].strip()
            ddg_value = line_split[4].strip()
            ddg_key = str(ddg_chain) + str(ddg_res)
            ddg_vals[ddg_key] = [ddg_value]
    
    for key in res_info_dict:
        if key in ddg_vals:
            res_info_dict[key].extend(ddg_vals[key])
        else:
            res_info_dict[key].extend([0])

    return res_info_dict