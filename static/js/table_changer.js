var previousScript = null;

function tableChanger(scriptSrc, callback) {
  // Get rid of existing table (why is this so annoying to accomplish)
  var existingTable = document.getElementById('table_display');
  if (existingTable) {
    existingTable.innerHTML = '';
  }

  // Remove the currently loaded changer script at the head of the HTML page
  if (previousScript) {
    document.head.removeChild(previousScript);
  }

  // Create a changer script element
  var script = document.createElement('script');
  script.src = scriptSrc;

  // Set up a callback for when the script is loaded
  script.onload = function () {
    // Update the reference to the current script
    previousScript = script;

    // Call the callback function
    if (callback && typeof callback === 'function') {
      callback();
    }
  };

  // Append the new changer script to the head
  document.head.appendChild(script);
}

//Handle hiding and showing of searched residues
function handleTableOperations() {
  var tableDisplay = document.getElementById('table_display');
  var inputRes = document.getElementById('inputRes');
  var submitButton = document.getElementById('submitButton');
  var clearButton = document.getElementById('clearButton');

  submitButton.addEventListener('click', function () {
    var inputValue = inputRes.value;
    //Make sure that the user actually inputs something in inputRes
    if (inputValue.trim() !== '') { 
      var chain = inputValue.charAt(0);
      var id = inputValue.slice(1);
      var table = tableDisplay.querySelector('table');
      
      //If table exists (it should) gp through all rows, lining up res and chain ids with corresponding row indicies
      if (table) { 
        var rows = table.getElementsByTagName('tr');
        
        var matchingRow = Array.from(table.getElementsByTagName('tr')).find(function (row) {
          return row.cells[1].textContent === id && row.cells[2].textContent === chain;
        });

        //If the chain/id are in the table, add 'hidden-row class to all rows'. Styling of this class sets row visibility to false
        if (matchingRow){ 
          for (var i = 1; i < rows.length; i++) {
            rows[i].classList.add('hidden-row');
          }
        }
        
        if (matchingRow) {
          var errorMessage = document.getElementById("res_not_found");
          errorMessage.style.visibility = "hidden";
          matchingRow.classList.remove('hidden-row'); //Unhide the matching row
          console.log('Matching row found and shown.');
        } else {
          var errorMessage = document.getElementById("res_not_found");
          errorMessage.style.visibility = "visible"; 
          console.log('No matching row found.');
        }
      }
    } else {
      console.log('Input is empty. Please enter a value.');
    }
  });

  //This function removes hidden-row class from all rows, making them all visible again
  clearButton.addEventListener('click', function () {
    var errorMessage = document.getElementById("res_not_found");
    errorMessage.style.visibility = "hidden"; 
    
    var table = tableDisplay.querySelector('table');
    if (table) {
      var rows = table.getElementsByTagName('tr');
      for (var i = 1; i < rows.length; i++) {
        rows[i].classList.remove('hidden-row');
      }
    }
  });
}

document.addEventListener('DOMContentLoaded', function () {
  // Initially load all results
  tableChanger(document.getElementById('display_all').dataset.scriptPath, handleTableOperations);

  // Load the full results on display_all button click
  document.getElementById('display_all').addEventListener('click', function () {
    tableChanger(this.dataset.scriptPath, handleTableOperations);
  });

  // Load the top results on display_top button click
  document.getElementById('display_top').addEventListener('click', function () {
    tableChanger(this.dataset.scriptPath, handleTableOperations);
  });
});
