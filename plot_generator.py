def adjust_rama_z_vals(x_data,y_data,z_data,rama_type):

    if rama_type == 'PRO':
        
        #Replace absolute path later...
        with open('/www/ProScan/data/pref_proline.data', 'r') as pref_file:
            phi_psi_vals = pref_file.readlines()
        
            for curr_coord_set in phi_psi_vals:
            
                if curr_coord_set[0] == '#':
                    continue
                else:
                    coords = curr_coord_set.split()
                    x_data.append(float(coords[0]))
                    y_data.append(float(coords[1]))
                    z_data.append(float(coords[2]))
 
    elif rama_type == 'PREPRO':
        
        #Replace absolute path later...
        with open('/www/ProScan/data/pref_preproline.data', 'r') as pref_file:
            phi_psi_vals = pref_file.readlines()
        
            for curr_coord_set in phi_psi_vals:
            
                if curr_coord_set[0] == '#':
                    continue
                else:
                    coords = curr_coord_set.split()
                    x_data.append(float(coords[0]))
                    y_data.append(float(coords[1]))
                    z_data.append(float(coords[2]))
                    
    return x_data,y_data,z_data

#Create_rama_plot imports in contour plots so new ones do not have to be generated each time
#If for some reason the rama plots need to be regenerated use the below function
def generate_fresh_contour_plots(file_path):

     #Add x, y, and z values of all phi-psi coord pairs
    x_data = []; y_data = []; z_data = []
    x_grid,y_grid,z_grid = adjust_rama_z_vals(x_data,y_data,z_data,'PRO')
    
    x_data_pre = []; y_data_pre = []; z_data_pre = []
    x_grid_pre,y_grid_pre,z_grid_pre = adjust_rama_z_vals(x_data_pre,y_data_pre,z_data_pre,'PREPRO')
    
    #Define contour graph resolution and color intervals
    contour_levels = dict(start=0, end=1, size=0.01,showlines=False)
    contour_colors_pro = [[0, '#FFFFFF'], [0.001,'#FFFFFF'],[0.002, '#7FFF8C'],
    [0.02, '#6cd487'],[1, '#6cd487']]
    
    contour_colors_pre = [[0, '#FFFFFF'], [0.001,'#FFFFFF'],[0.002, '#FFEDA0'],
    [0.02, '#FEB24C'],[1, '#FEB24C']]

    
    #Initialize contour plots
    contour_plot_pro = go.Contour(z=z_grid, x=x_grid, y=y_grid,showscale=False,
                              contours_coloring='fill', contours=contour_levels,
                              colorscale=contour_colors_pro,hoverinfo='none',connectgaps=False)
    
    contour_plot_pre = go.Contour(z=z_grid_pre, x=x_grid_pre, y=y_grid_pre,showscale=False,
                              contours_coloring='fill', contours=contour_levels,
                              colorscale=contour_colors_pre,hoverinfo='none',connectgaps=False)
    
    contour_pro_fig = go.Figure(data=[contour_plot_pro])
    contour_pre_fig = go.Figure(data=[contour_plot_pre])

    
    pro_path = os.path.join(file_path,"pro_contour.json")
    pre_path = os.path.join(file_path,"pre_contour.json")
    #GENERATED PLOTS WILL BE SAVED TO RUNS FOLDER BY DEFAULT#
    contour_pro_fig.write_json(pro_path)
    contour_pre_fig.write_json(pre_path)
