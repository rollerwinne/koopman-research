function mvtemp(path)
sourcePath='./temp';
targetPath=['/Users/mmc/WorkSpace/bupt/GraduateThesis/images/',path];
fileList = dir(sourcePath);  
  
for k = 3:length(fileList)  
    movefile([sourcePath,'/',fileList(k).name],targetPath);
end 
end