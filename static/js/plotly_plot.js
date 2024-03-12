function grab_plotly(graphData) {
    document.addEventListener("DOMContentLoaded", function () {
        // Assuming graphData is in the correct format for Plotly
        Plotly.newPlot('pro_plot', graphData.data, graphData.layout,{displayModeBar:false});
    });
}


var submitButton = document.getElementById('submitButton');
submitButton.addEventListener('click', function () {
    applyFilter();
});

var submitButton = document.getElementById('clearButton');
submitButton.addEventListener('click', function () {
    clearFilter();
});

window.addEventListener('resize', function() {
    Plotly.Plots.resize(document.getElementById('pro_plot'));
});

//This function loads in the chain/residue from inputRes input of results.html
//This will parse out the chain/residue pair and check if its present in the data of the plotly plot object grabbed by div (graphDiv)
//If res/chain pair are in the plot data, it will create a new trace with just that datapoint. It will hide all other traces, so only this new trace is displayed 
function applyFilter() {
    
    var graphDiv = document.getElementById('pro_plot')
    var inputRes = document.getElementById('inputRes');
    var inputValue = inputRes.value;

    //Check the user actually inputted something in inputRes
    if (inputValue.trim() !== '') {
        var chain = inputValue.charAt(0);
        var id = inputValue.slice(1);
    }
    
    // Hide all points
    var style_update = {
        opacity: 0.0,
    };
    
    // Hide all hover annotations 
    var layout_update = {
        hovermode: false
    };

    var resFound = 0; //Checks if residue is found in any trace
    
    //Loop through all four current traces
    for (var i = 0; i < graphDiv.data.length; i++) {
        
        var xCoordinates = graphDiv.data[i].x;
        var yCoordinates = graphDiv.data[i].y;
        for (var j = 0; j < xCoordinates.length; j++) { //Loop through the member of each trace
            
            var currRes = graphDiv.data[i].customdata[j][0];
            var currChain = graphDiv.data[i].customdata[j][3];
            if (currRes === id && currChain === chain) { //If residue is found by matching up parsed id/chain to customdata array
                
                resFound = 1;
                var currX = xCoordinates[j];
                var currY = yCoordinates[j];
                var currResType = graphDiv.data[i].customdata[j][1];
                var currMPNNScore = graphDiv.data[i].customdata[j][2];
                //Create new trace with just the found residue/chain datapoint
                var scatterPlotSelection = {x: [currX],y: [currY],mode: 'markers',showlegend: false,
                    marker: {color: 'blue',size: 8,line: {color: 'black',width: 2} },
                    customdata: [currRes, currResType, currMPNNScore, currChain]
                };
            }
        }
    }

    //If the residue has been found
    if (resFound === 1){
        //If another residue is already filtered out, clear it off first and reset the plots
        if (graphDiv.data.length > 4){
            clearFilter();
            }
        //Then apply filters, first hiding all points and their annotations. Then add the new layer with the single point   
        Plotly.restyle(graphDiv, style_update);
        Plotly.relayout(graphDiv,layout_update);
        Plotly.addTraces(graphDiv, scatterPlotSelection);      
    }    
}

//Resets the plotly view so all original datapoints and annotations are shown.
function clearFilter() {
    
    var graphDiv = document.getElementById('pro_plot')
    
    //Check to ensure that there is a fifth trace (Ensures user only clears custom filter trace)
    if (graphDiv.data.length > 4){
    Plotly.deleteTraces(graphDiv, [-1]); //Remove most recent trace, which will be the custom filter trace
    var update = {
        opacity: 1.0 
    };

    var layout_update = {
        hovermode: true
    };
    
    Plotly.restyle(graphDiv, update); //Add back visibility of all points and their annotations
    Plotly.relayout(graphDiv,layout_update);
    }
}

