var container = document.getElementById("table-container");
var button = document.getElementById("expand_table");


button.addEventListener("click", tableExpander);

function tableExpander() {
    if (container.style.height === "450px") {
        // If the current height is 450px, remove the height and change button text to "Collapse"
        container.style.height = "auto";
        button.textContent = "Collapse Table";
      } else {
        // If the current height is not 450px, set it to 500px and change button text to "Expand"
        container.style.height = "450px";
        button.textContent = "Expand Table";
      }
}
