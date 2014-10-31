function S=rand(mc,T)
%S=rand(mc,T) returns a random state sequence from given MarkovChain object.
%
%Input:
%mc=    a single MarkovChain object
%T= scalar defining maximum length of desired state sequence.
%   An infinite-duration MarkovChain always generates sequence of length=T
%   A finite-duration MarkovChain may return shorter sequence,
%   if END state was reached before T samples.
%
%Result:
%S= integer row vector with random state sequence,
%   NOT INCLUDING the END state,
%   even if encountered within T samples
%If mc has INFINITE duration,
%   length(S) == T
%If mc has FINITE duration,
%   length(S) <= T
%
%---------------------------------------------
%Code Authors:
%---------------------------------------------

S=zeros(1,T);%space for resulting row vector
nS=mc.nStates;

init = mc.InitialProb;
A = mc.TransitionProb;
t = 1;
pD = DiscreteD(init);

if (mc.finiteDuration==1)
    endState = nS+1;
    lastState = rand(pD,1);
    if (lastState~=endState)
            S(t)=lastState;

        for t=2:T  
            pD = DiscreteD(A(lastState,:));
            lastState=rand(pD,1);
            if (lastState~=endState)
                S(t)=lastState;
            else
                break;
            end    
        end

    end
    
else
    lastState = rand(pD,1);
    S(t) = lastState;
    for t=2:T
        pD = DiscreteD(A(lastState,:));
        lastState=rand(pD,1);
        S(t)=lastState;
    end
end

% remove zero elements (if chain is finite and length(S) < T)
S( :, ~any(S,1) ) = [];

  

