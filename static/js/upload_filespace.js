var fileDropzone = document.getElementById('file-dropzone');
var fileInput = document.getElementById('file-input');
var fileInfo = document.getElementById('file-info');
var fileName = document.getElementById('file-name');

fileDropzone.addEventListener('click', function() {
  fileInput.click();
});

fileDropzone.addEventListener('dragover', function(event) {
  event.preventDefault();
  fileDropzone.style.border = '2px dashed #888';
});

fileDropzone.addEventListener('dragleave', function() {
  fileDropzone.style.border = '2px dashed #ccc';
});

fileDropzone.addEventListener('drop', function(event) {
  event.preventDefault();
  fileDropzone.style.border = '2px dashed #ccc';
  fileInput.files = event.dataTransfer.files;
  displayFileInfo(fileInput.files[0]);
});

fileInput.addEventListener('change', function() {
  displayFileInfo(fileInput.files[0]);
});

function displayFileInfo(file) {
  fileInfo.style.display = 'block';
  fileName.textContent = 'Selected File: ' + file.name;
}
