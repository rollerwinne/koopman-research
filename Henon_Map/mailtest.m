clear
tic;timestart=char(datetime('now'));

%% these are program
pause(1);

%% copy these code to send an E-mail
%% send an E-mail to me
[a,b]=qqmail2me(timestart,mfilename('fullpath'),[]); %程序开始时间、文件名、附件