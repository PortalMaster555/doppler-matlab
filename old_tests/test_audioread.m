clc; clear all;

% from MATLAB example code
[y, Fs] = audioread("test_440_48khz.wav");
sound(y,Fs);