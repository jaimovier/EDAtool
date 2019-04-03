%% Example script for EDAtool
% by Javier Jaimovich (2012)

clc; help EDAtool.m

%% Load example EDA

load('EDAtoolExamples');

choice = menu('Select EDA signal',...
    cellfun(@(y) y.name,Examples,'UniformOutput',false));

disp(Examples{choice})

EDA = Examples{choice}.signal;
SR = Examples{choice}.SR;
range = Examples{choice}.range;

%% Run function

%PARAMETERS
Threshold = 4;
ArtefactWindow = 0.2;
ResampleOption = 0;
InterOption = 1;
NaNSize = 2.5;

% Make parameter vector
par = [range ArtefactWindow Threshold NaNSize];

help EDAtool.m
[EDA_f, EDA_x, EDA_P, EDA_T, Q] = ...
    EDAtool(EDA, SR, 1, InterOption, ResampleOption, par);
