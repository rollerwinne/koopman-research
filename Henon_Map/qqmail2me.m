function [subject,content]=qqmail2me(timestart,filename,attachments)
timeuse=toc;
if timeuse>3600
    timestr=[num2str(timeuse/3600,4),' hours'];
elseif timeuse>60
    timestr=[num2str(timeuse/60,4),' minites'];
else
    timestr=[num2str(timeuse),' seconds'];
end
timeend=char(datetime('now'));
[~,name,~]=fileparts(filename);
IPaddress=char(java.net.InetAddress.getLocalHost);
subject=['''',name,'.m'' has been finished'];
content={['Time used: ',timestr,'.'];
    ['Start time: ',timestart,'.'];
    ['Finish time: ',timeend,'.'];
    ['Path: ''',filename,'.m''.']
    ['Address: ',IPaddress,'.'];
    ['----------'];
    ['This is an auto E-mail from MATLAB.']};

MailAddress = '502358401@qq.com';
password = 'rygppwbxdikjcbeb';
receiver='13269482397@163.com';
setpref('Internet','E_mail',MailAddress);
setpref('Internet','SMTP_Server','smtp.qq.com');
setpref('Internet','SMTP_Username',MailAddress);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
% props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
% props.setProperty('mail.smtp.socketFactory.port','465');
sendmail(receiver,subject,content,attachments);
end