document.getElementById("view_tester").addEventListener("click", function() {

    var url = '/mut_viewer/' + jobName; 
    var newTab = window.open(url, "_blank");
    newTab.focus();
});
