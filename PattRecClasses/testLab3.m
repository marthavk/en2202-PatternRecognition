% script for testing backward algorithm
clc;
clear all;
q=[1; 0];
A=[0.9 0.1 0; 0 0.9 0.1];
b1 = GaussD('Mean', 0, 'StDev', 1);
b2 = GaussD('Mean', 3, 'StDev', 2);
B=[b1, b2];
h=HMM;
mc = MarkovChain;
mc.InitialProb = q;
mc.TransitionProb = A;
h.StateGen = mc;

%% construct pX (pX(j,t)= P( X(t)= observed x(t) | S(t)= j ))
% from arbitrary observation sequence x
x = [-0.2  2.6 1.3];
pX = B.prob(x);

% declare test c vector
c = [1 0.1625 0.8266 0.0581];

%% test backward algorithm
betaHat = backward(mc, pX, c)