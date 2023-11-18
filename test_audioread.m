clc; clear all;

% from MATLAB example code
load handel.mat
filename = 'test_handel.wav';
audiowrite(filename,y,Fs);
clear y Fs
[y, Fs] = audioread("test_handel.wav");
sound(y,Fs);