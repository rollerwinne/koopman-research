%% Question before Summer vacation
% 随机梯度算法并不随机 https://www.cnblogs.com/lliuye/p/9451903.html
% 
% 约束条件写到argmin中
% 

%% 2019.7.7
% 1.鲁棒性决定分界点(加噪声)
% 2.用不同的精度(本征函数的个数)去逐步精确本征函数的划分

%% Summer vacation task
% 随机梯度算法 
% 写报告for SVG
% 论文结构

%% 2019.8.26 Lecture for Domenico
% 1.报告中描述基函数不够清晰(slides需要改进)
% 2.为什么要取这样(高斯 or else)的基函数,为什么这么取基函数就是对的(正交的投影更有效,但是也不一定正交,函数空间尽可能完备)
% 3.本征函数的极小值点为什么反映的边界点,极大值点能反映出什么性质(interval 中间 ,有可能跟周期轨道有关)
% 4.这种做法具有鲁棒性吗,加了噪声会不会有好的划分

%% 2019.8.31 Discuss with Lan
% 1.本征函数加了噪声
% 2.加了噪声还能不能通过argmin来找到"结构简单"的本征函数

% 噪声的分辨率和基函数的分辨率
% robustness 基函数 噪声 系统不一样