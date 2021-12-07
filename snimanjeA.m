%zadatak 2 

clc 
clear all
close all

recObj = audiorecorder;
disp('Start speaking');
t = 2;
recordblocking(recObj,t);
disp('End of Recording');

x = getaudiodata(recObj);
play(recObj);
%%
Fs = 8000;
audiowrite('nbA.wav', x, Fs);

