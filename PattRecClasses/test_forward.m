clc;
clear;
q=[1; 0];
A=[0.9 0.1 0; 0 0.9 0.1];
b1 = GaussD('Mean', 0, 'StDev', 1);
b2 = GaussD('Mean', 3, 'StDev', 2);
B=[b1, b2];
h=HMM;
mc = MarkovChain;
mc.InitialProb = q;
mc.TransitionProb = A;
h.OutputDistr = B;
h.StateGen = mc;
x = [-0.2 2.6 1.3];

%% Test Forward Algorithm
pX = B.prob(x);
[alfaHat, c]=forward(mc,pX)

alfaHat_correct = [1 0.3847 0.4189; 0 0.6153 0.5811];
c_correct = [1 0.1625 0.8266 0.0581];
error_alfaHat = sum(sum(alfaHat_correct-alfaHat))
error_c = sum(c_correct-c)


%% Test LogProb
[logP] = logprob(h, x);
logP_correct = -9.1877;
error_logP = logP_correct - logP


