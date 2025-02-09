% Demo of response spectra computation for given ground motion acceleration
clear clc close all;

% example ground-motion waveform from the 2014 M6.0 South Napa Eq.
load('./data/gacc.mat'); 

% dt = sampling interval in seconds
dt = 0.005;

% xi = ratio of critical damping
xi = 0.05;  

% sPeriod = spectral period vector
sPeriod = [0.01,0.02,0.022,0.025,0.029,0.03,0.032,0.035,0.036,...
    0.04,0.042,0.044,0.045,0.046,0.048,0.05,0.055,0.06,0.065,0.067,0.07,...
    0.075,0.08,0.085,0.09,0.095,0.1,0.11,0.12,0.125,0.13,0.133,0.14,0.15,...
    0.16,0.17,0.18,0.19,0.2,0.22,0.24,0.25,0.26,0.28,0.29,0.3,0.32,0.34,...
    0.35,0.36,0.38,0.4,0.42,0.44,0.45,0.46,0.48,0.5,0.55,0.6,0.65,0.667,...
    0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,...
    2,2.2,2.4,2.5,2.6,2.8,3,3.2,3.4,3.5,3.6,3.8,4,4.2,4.4,4.6,4.8,5,7.5,10];

[PSA, PSV, SD, SA, SV, OUT] = responseSpectra(xi, sPeriod, gacc, dt);

addpath ('libs'); 

% plot PSA, PSV and SD spectrum 
plotSpectra(1,sPeriod,PSA,PSV,SD,'./plots/','responseSpectra'); 

