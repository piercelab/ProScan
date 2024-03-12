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
  
  var BAD_COLOR = "0xD3D3D3"; var ACCEPT_COLOR = "0xFAFA33"; var PRO_COLOR = "0x6F8FAF"; var GOOD_COLOR = "0x00bb00"
  var BAD_BOUNDS = 3; var ACCEPT_BOUNDS = 20; var PRO_BFACTOR = 101; 

  var stage = new NGL.Stage("ngl_viewer");
  stage.setParameters({ backgroundColor: "white" });
    
  // handle window resizing
  window.addEventListener("resize", function(event) {
    stage.handleResize();
  }, false );

  var highBFactorResidues = []; 

  //Access all the atoms. The b-factors of all residues are based on ProteinMPNN probs, so recolor based on bfactor
  function createColormakerScheme() {  
    var schemeId = NGL.ColormakerRegistry.addScheme(function(params) {
      this.atomColor = function(atom) {
        if (atom.bfactor < BAD_BOUNDS) {
            return BAD_COLOR;
        } else if (atom.bfactor >= BAD_BOUNDS && atom.bfactor < ACCEPT_BOUNDS) {
              highBFactorResidues.push(atom.resno);
            return ACCEPT_COLOR;
        } else if (atom.bfactor == PRO_BFACTOR) { 
            return PRO_COLOR;
        } else {
            highBFactorResidues.push(atom.resno);
            return GOOD_COLOR;
        }
      };
    });
      return schemeId;
  }
   
  var recenterButton = document.getElementById('recenterButton');
  recenterButton.addEventListener('click', function () {
      recenterView(stage);
  });

  stage.loadFile(pdbData, { ext: 'pdb' }).then(function(o){
      
    schemeId = createColormakerScheme()
    o.addRepresentation("cartoon", {color: schemeId})
    var highBFactorSelection = highBFactorResidues.join(" + ");
      
    if (highBFactorResidues.length !== 0) {
      o.addRepresentation("ball+stick", { sele: highBFactorSelection, color: schemeId });  
    }
    stage.autoView();

    var rep = 0;
    var submitButton = document.getElementById('submitButton');

    submitButton.addEventListener('click', function () {
      rep = goToResidue(o,rep);
    }); 

    var clearButton = document.getElementById('clearButton');
    clearButton.addEventListener('click', function () {
      rep = resetSeleView(o,rep);   
      var inputRes = document.getElementById('inputRes');
      inputRes.value = ""
    }); 
  }); 
}

function recenterView(stage) {
    stage.autoView();
}

function goToResidue(o,rep) {
  
  var inputRes = document.getElementById('inputRes');
  var inputValue = inputRes.value;

  if (inputValue.trim() !== '') {
    var chain = inputValue.charAt(0);
    var id = inputValue.slice(1);
    chain = ":" + chain;
    curr_sele = chain + " and " + id
    
    var atomCheck = o.structure.getView(new NGL.Selection(curr_sele))
    var numAtom = atomCheck.atomCount
    
    if (numAtom !== 0){
      if (rep !== 0){
        rep.setVisibility(false);
      }
     
      rep = o.addRepresentation("ball+stick", {sele: curr_sele, color:schemeId});
      o.autoView(curr_sele);
      return rep; 
    }
  }
  return 0;
}

function resetSeleView(o,rep) {
  
  if (rep !== 0){
    rep.setVisibility(false);
    o.autoView();
    selectorActive = 0;
    rep = 0;
    return rep;
  }
  
  
  o.autoView();
  return rep;
}


