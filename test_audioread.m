clc; clear all;

% from MATLAB example code
[y, Fs] = audioread("test_440_pure.mp3");
sound(y,Fs);