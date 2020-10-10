function mvtemp
sourcePath='./temp';
targetPath='/Users/mmc/WorkSpace/bupt/KoopmanThesis/images/new';
targetPath='./recycle';
%targetPath=['/Users/mmc/WorkSpace/bupt/GraduateThesis/images/',path];
fileList = dir(sourcePath);  
  
for k = 3:length(fileList)  
    movefile([sourcePath,'/',fileList(k).name],targetPath);
end 
end