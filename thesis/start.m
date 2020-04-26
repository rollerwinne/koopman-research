function start(~)
clear;clc;close all;
if nargin<1
    addpath(genpath(pwd));
else
    rmpath(genpath(pwd));
end
end