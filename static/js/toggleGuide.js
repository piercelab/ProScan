
var guideDiv = document.getElementById("quick_guide_text");
var guideToggle = document.getElementById("guideToggle");
guideToggle.addEventListener("click", function() {
    if (guideDiv.style.display === "none" || guideDiv.style.display === "") {
        guideDiv.style.display = "flex";
      guideToggle.textContent = "Hide Results Quick Guide";
    } else {
        guideDiv.style.display = "none";
      guideToggle.textContent = "Show Results Quick Guide";
    }
  });
