from flask import Flask, render_template, request, send_file, url_for, redirect, jsonify
from werkzeug.utils import secure_filename
import os
import proscan_utils
import random
import string
import datetime
import json
from multiprocessing import Process
app = Flask(__name__)
app.secret_key = 'dnfaoiduf98au34'
app.config['UPLOAD_EXTENSIONS'] = ['.pdb','.cif']

RUNDIR_PATH = os.path.join(app.root_path, 'runs')
EXAMPLE_PATH = "/www/ProScan/static/examples/4z0x.cif_CUD_240312_125507"
EXAMPLE_DIR = "4z0x.cif_CUD_240312_125507"

def get_unique_name(base_name):
    
    random_tag = '_' + ''.join(random.SystemRandom().choice(
        string.ascii_uppercase + string.digits) for _ in range(3))
    datetime_tag = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
    uniquedirname = "_".join([("".join([base_name, random_tag])), datetime_tag])
    return uniquedirname
 
def create_and_cd_to_unique_dir(base_name):
       
    unique_job_name = get_unique_name(base_name)
    unique_dir_name = os.path.join(RUNDIR_PATH, unique_job_name)
    os.makedirs(unique_dir_name, exist_ok=True)

    return unique_job_name

#Save uploaded pdb and calculate results with calc_ramachandran
def run_proscan(job_name):

    file_name = job_name[:-18]
    file_name_all = job_name[:-18] + "_ProScan_results_all.txt"
    file_name_top = job_name[:-18] + '_ProScan_results_top.txt'
    
    run_path = os.path.join(RUNDIR_PATH, job_name)
    proscan_utils.remove_chains(file_name, run_path)
    proscan_utils.calc_ramachandran(file_name, run_path)
    proscan_utils.write_high_scoring(run_path, file_name_all, file_name_top)   

    return None

#Create tables for results page and plotly graph
def create_graphs(job_name):

    file_name_all = job_name[:-18] + "_ProScan_results_all.txt"
    file_name_top = job_name[:-18] + '_ProScan_results_top.txt'
    run_path = os.path.join(RUNDIR_PATH, job_name)
    
    user_options_list = get_user_option_inputs(job_name)
    if user_options_list[3] == "True":
        proscan_utils.create_rama_plot_gly(run_path, file_name_all)
    else:
        proscan_utils.create_rama_plot(run_path, file_name_all)    
    proscan_utils.output_to_HTML(run_path, file_name_all, "all")
    proscan_utils.output_to_HTML(run_path, file_name_top, "top")
    return None

#User input is written to a textfile to be accessed later in the run
def write_chain_list_file(fileName, chain_list, run_path, chain_list_type, ros_check, gly_check):

    
    chain_info_file = os.path.join(run_path, 'chain_info.txt')
    chain_check = chain_list.split(',')
    chain_info = open(chain_info_file,'w')
    
    if chain_check == ['']: 
        chain_info.write('EMPTY\n')
        chain_info.write(chain_list_type + '\n')
    else:
        chain_info.write(chain_list + '\n')
        chain_info.write(chain_list_type + '\n')
    
    if ros_check == 'off':
        chain_info.write('True' + '\n')
    else:
        chain_info.write('False' + '\n')
    
    if gly_check == 'on':
        chain_info.write('True')
    else:
        chain_info.write('False')

    chain_info.close()

    return None

def get_user_option_inputs(job_name):

    run_path = os.path.join(RUNDIR_PATH, job_name)
    user_options = os.path.join(run_path,'chain_info.txt')
    user_options_list = []

    with open (user_options,'r') as options:
        for line in options:
            user_options_list.append(line.strip())
    
    return user_options_list

#Render initial page. File upload handled by index.html
@app.route('/')  
def home():  
    
    return render_template("index.html") 

@app.route('/upload',methods = ['POST'])
def upload():

    idPDB = request.form["idPDB"]
    chain_list = request.form["chainList"]
    user_file = request.files['file']
    chain_list_type = request.form['chainType']
    ros_check = request.form.get('rosCheckbox','off')
    gly_check = request.form.get('glyCheckbox','off')

    #If user uploads a file and puts in PDB give them an error
    if user_file and idPDB:              
        error_message = "Both submission fields filled. PDB code and file upload detected. Please upload file OR enter PDB code"            
        return render_template("error.html", errorMessage = error_message)   
    #If user does not upload a file give them another error      
    elif not user_file and not idPDB:            
        error_message = "No pdb uploaded. Please enter a valid pdb code or upload pdb file"            
        return render_template("error.html", errorMessage = error_message)
    #If enters pdb code, retrieve pdb file from rcsb and save it as cif
    elif idPDB and not user_file:
        idPDB = idPDB.lower()
        job_name = create_and_cd_to_unique_dir(idPDB + '.cif')
        run_path = os.path.join(RUNDIR_PATH, job_name)          
        fileName = proscan_utils.retrievePDB(idPDB, run_path)              
        
        if fileName == "retrieveFailed":
            error_message = "PDB file could not be fetched. Please check PDB ID and try again"            
            return render_template("error.html", errorMessage = error_message)          
    else:         
        fileName = user_file.filename  
        if fileName[-4:] != '.cif' and fileName[-4:] != '.pdb': #Check that file type is allowed
            error_message = "Unrecognized file format uploaded. Please upload structure file in pdb or mmcif/pdbx format"            
            return render_template("error.html", errorMessage = error_message)             
        
        job_name = create_and_cd_to_unique_dir(fileName)
        run_path = os.path.join(RUNDIR_PATH, job_name)
        user_file.save(os.path.join(run_path, secure_filename(user_file.filename)))             
    
    write_chain_list_file(fileName,chain_list, run_path, chain_list_type, ros_check, gly_check)
    return redirect(url_for('scanning_structure', jobName = job_name))

@app.route('/run_example')
def run_example():

    example_path =  EXAMPLE_PATH
    jobName = EXAMPLE_DIR
    json_path = os.path.join(example_path, "rama_plot.json")  
    with open(json_path,'r') as plotly_json:
        plot_data = json.load(plotly_json)
    return render_template("results.html", plotData = plot_data, jobName = jobName)

@app.route('/scanning_structure/<jobName>')
def scanning_structure(jobName):
    
    proscan_process = Process(target = run_proscan, args = (jobName,))
    proscan_process.start()
    return render_template('running_scan.html', jobName = jobName)

@app.route('/check_status/<jobName>')
def check_status(jobName):

    file_name_top = jobName[:-18] + '_ProScan_results_top.txt'
    run_path = os.path.join(RUNDIR_PATH, jobName)
    top_file = os.path.join(run_path, file_name_top)
    mut_file = os.path.join(run_path,'sorted_struct.muts.txt')
    error_file = os.path.join(run_path,'error.txt')
    
    if os.path.exists(top_file):
        status = 'C'
    elif os.path.exists(error_file):
        status = 'E'
    elif os.path.exists(mut_file): 
        progress_file = os.path.join(run_path,'mut_track.txt')
        if os.path.exists(progress_file) and os.stat(progress_file).st_size != 0:
            with open(progress_file, 'rb') as progress:
                lines = progress.readlines() 
                if len(lines) > 3:
                    total_muts = lines[0].decode('utf-8').strip()
                    second_to_last_line = lines[-2].decode('utf-8').strip()  # Get the second-to-last line
                    last_mut = second_to_last_line.split()[-1]  # Extract the mut number of that line
                    status = 'L'
                    return jsonify({'status': status, 'lastMut': str(last_mut), 'totalMut': str(total_muts)})
        status = 'G'
    else:
        status = 'R'

    return jsonify({'status': status})

@app.route('/loading_graphs/<jobName>')
def loading_graphs(jobName):
    return render_template('loading_graphs.html', jobName = jobName)

@app.route('/graphs_loading/<jobName>')
def graphs_loading(jobName):

    create_graphs(jobName)
    return render_template('loading_complete.html', jobName = jobName)

@app.route('/loading_complete/<jobName>')
def loading_complete(jobName):

    return redirect(url_for('results', jobName = jobName))

@app.route('/results/<jobName>',methods = ['POST','GET'])  
def results(jobName):
    
    run_path = os.path.join(RUNDIR_PATH, jobName)
    json_path = os.path.join(run_path, "rama_plot.json")  
    with open(json_path, 'r') as plotly_json:
        plot_data = json.load(plotly_json)    
    return render_template("results.html", plotData = plot_data, jobName = jobName)
   
#Allows user to download full results file
@app.route('/download/<jobName>', methods = ['GET','POST'])
def download(jobName):
    
    file_name = jobName[:-18] + "_ProScan_results_all.txt"
    if jobName == EXAMPLE_DIR:
        run_path = EXAMPLE_PATH
    else:
        run_path = os.path.join(RUNDIR_PATH, jobName)
    file_path = os.path.join(run_path, file_name)
    return send_file(file_path, as_attachment = True)

#Allows user to download best scoring results file
@app.route('/download_top/<jobName>',methods = ['GET','POST'])
def download_top(jobName):

    if jobName == EXAMPLE_DIR:
        run_path = EXAMPLE_PATH
    else:
        run_path = os.path.join(RUNDIR_PATH, jobName)
    file_path = os.path.join(run_path,jobName[:-18] + '_ProScan_results_top.txt') 
    return send_file(file_path, as_attachment = True)
         
#Return to first page and clear the contents of the session
@app.route('/goHome',methods = ['POST'])
def go_home():
      
    return redirect(url_for('home'))
   
@app.route('/fetch_table/<jobName>')
def fetch_table(jobName):

    if jobName == EXAMPLE_DIR:
        run_path = EXAMPLE_PATH
    else:
        run_path = os.path.join(RUNDIR_PATH, jobName)
    HTML_path = os.path.join(run_path, "rama_results.html")
    
    with open(HTML_path) as f:
        html = f.read()    
    return jsonify({'table': html})

@app.route('/fetch_table_top/<jobName>')
def fetch_table_top(jobName):

    if jobName == EXAMPLE_DIR:
        run_path = EXAMPLE_PATH
    else:
        run_path = os.path.join(RUNDIR_PATH, jobName)

    HTML_path = os.path.join(run_path, "rama_results_top.html")
    
    with open(HTML_path) as f:
        html = f.read()    
    return jsonify({'table': html})

@app.route('/get_pdb/<jobName>')
def get_pdb(jobName):
    
    if jobName == EXAMPLE_DIR:
        run_path = EXAMPLE_PATH
    else:
        run_path = os.path.join(RUNDIR_PATH, jobName)
    pdb_path = os.path.join(run_path, 'sorted_struct.pdb')   
    return send_file(pdb_path, as_attachment = False)

@app.route('/about')
def about():
    return render_template("about.html")
    
@app.route('/help')
def help():
    return render_template('help.html')

@app.route('/error')
def error():
    
    error_message = 'ProteinMPNN Parsing Error. Please check the format of your input and resubmit'
    return render_template('error.html', errorMessage = error_message)

if __name__ == '__main__':  
    app.run(debug=True)


