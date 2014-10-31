function [alfaHat, c]=forward_fixme(mc,pX)
%calculates state and observation probabilities for one single data sequence,
%using the forward algorithm, for a given single MarkovChain object,
%to be used when the MarkovChain is included in a HMM object.
%
%Input:
%mc= MarkovChain object
%pX= matrix with state-conditional likelihood values,
%   without considering the Markov depencence between sequence samples.
%	pX(j,t)= myScale(t)* P( X(t)= observed x(t) | S(t)= j ); j=1..N; t=1..T
%NOTE: pX may be arbitrarily scaled, as defined externally,
%   i.e., pX may not be a properly normalized probability density or mass.
%
%NOTE: If the HMM has Finite Duration, it is assumed to have reached the end
%after the last data element in the given sequence, i.e. S(T+1)=END=N+1.
%
%Result:
%alfaHat=matrix with normalized state probabilities, given the observations:
%	alfaHat(j,t)=P[S(t)=j|x(1)...x(t), HMM]; t=1..T
%c=row vector with observation probabilities, given the HMM:
%	c(t)=P[x(t) | x(1)...x(t-1),HMM]; t=1..T
%	c(1)*c(2)*..c(t)=P[x(1)..x(t)| HMM]
%   If the HMM has Finite Duration, the last element includes
%   the probability that the HMM ended at exactly the given sequence length, i.e.
%   c(T+1)= P( S(T+1)=N+1| x(1)...x(T-1), x(T)  )
%Thus, for an infinite-duration HMM:
%   length(c)=T
%   prod(c)=P( x(1)..x(T) )
%and, for a finite-duration HMM:
%   length(c)=T+1
%   prod(c)= P( x(1)..x(T), S(T+1)=END )
%
%NOTE: IF pX was scaled externally, the values in c are 
%   correspondingly scaled versions of the true probabilities.
%
%--------------------------------------------------------
%Code Authors: Course TAs
%Last Edited: Gerasimos Markatos, Martha Vlachou-Konchylaki 25/10/2014
%--------------------------------------------------------
%%

% T represent the number of Time instances
T=size(pX,2);       

% N represents the number of states                    
N = length(mc.InitialProb);     

% Parameters of the HMM model
q = mc.InitialProb;             
A = mc.TransitionProb;          
B = pX;                                                                                                                        

% Scaling Factors
c = zeros(1,T);  

[rows,columns] = size(A);   
% set parameter FINITE to 1 if the model is FINITE
FINITE = 0;
if(rows ~= columns)     
    FINITE = 1;
    c = [c 0]; 
end

alfaHat = zeros(N,T); 
alfaTemp = zeros(N,T); 

%% INITIALIZATION
alfaTemp(:,1) = q.*B(:,1);   
c(1) = sum(alfaTemp(1,:),2); 
alfaHat(:,1) = alfaTemp(:,1)/c(1);  


%%  FORWARD STEP
for t=2:T
    for j=1:N
        alfaTemp(j,t) = B(j,t)*sum(alfaHat(:,t-1).*A(:,j)); 
    end
    c(t) = sum(alfaTemp(:,t));  
    alfaHat(:,t) = alfaTemp(:,t)/c(t); 
end

%% TERMINATION
if(FINITE)    
    c(end) = sum(alfaHat(:,end).*A(:,end));
end
