var component;
document.addEventListener("DOMContentLoaded", function () {
 
  fetch(`/get_pdb/${jobName}`)
    .then(response => {
      if (!response.ok) {
        throw new Error('Failed to fetch PDB file');
      }
        return response.blob(); 
    })
    .then(data => {
      visualizePDBWithBFactor(data);
    })
    .catch(error => {
      console.log("Error fetching PDB file:", error);
    });
});

function visualizePDBWithBFactor(pdbData) {
  
  var stage = new NGL.Stage("mut_viewport");
  stage.setParameters({ backgroundColor: "white" });
    
  // handle window resizing
  window.addEventListener("resize", function(event) {
    stage.handleResize();
  }, false );

  
  stage.loadFile(pdbData, { ext: 'pdb' }).then(function(o){
    
    o.addRepresentation("cartoon")
    o.autoView();
  }); 
}

function recenterView(stage) {
    stage.autoView();
}


