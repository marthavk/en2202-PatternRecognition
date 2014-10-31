function [alfaHat, c]=forward_fixed(mc,pX)
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
%--------------------------------------------------------
%%

% general comments and changes: 
% * we eliminated for loops when possible
% * we preallocated matrices, setting the sizes, when possible
% * we corrected logic errors
% * we corrected dynamic semantic errors 
% * we erased everything that had to do with Z.

% T represent the number of Time instances
T=size(pX,2);       % CORRECTED: pX is a matrix with state-conditional 
                    %likelihood values over time

% N represents the number of states                    
N = length(mc.InitialProb);     % CORRECTED : we changed numberOfStates into N and added a comment since it will be used a lot through the following lines
q = mc.InitialProb;             % CORRECTED : Use of [] not recommended. 
A = mc.TransitionProb;          % OK
B = pX;                         % TO CHECK : ERASE B=pX because pX does not correspond 
                                % to the matrix B of an HMM. ???
                                % pX will be referred with its name from
                                % now on

c = zeros(1,T);  % CORRECTED: c=zeros(numberOfStates) returns a square matrix. We need size (1,numberOfStates)

[rows,columns] = size(A);   %   OK
% for better code reading we set a parameter of FINITE to 1 if the model is
% FINITE
FINITE = 0;
% here it is checked if the model is finite or infinite                       
if(rows ~= columns)     
    FINITE = 1;
    c = [c 0];              % Corrected: c is the one to have one more 
                            % element if the model is finite, not q.
                            % Also, c is a row, not a column vector.            
end

alfaHat = zeros(N,T);              % TO CHECK: Set the size of this matrix as well - good programming practice
alfaTemp = zeros(N,T);     % CORRECTED: The size of the matrix should be dynamic.
% we used only one alfaTemp, there was no reason to use two. In the
% initialization step we can only set the first column of the alphaTemp
% matrix

%% INITIALIZATION
alfaTemp(:,1) = q.*B(:,1);      % CORRECTED: set only the first column 
                                % of initAlfaTemp. Also, the foor loop 
                                % can be avoided      


c(1) = sum(alfaTemp(1,:),2);        % CORRECTED

% we can change that and avoid the for loop
alfaHat(:,1) = alfaTemp(:,1)/c(1);  


%%  FORWARD STEP
for t=2:T
    % alfaTemp = [];            no need to declare alfaTemp again
    for j=1:N
        alfaTemp(j,t) = B(j,t)*sum(alfaHat(:,t-1).*A(:,j)); % changed the 
        % multiplication. we don't want matrix multiplication, we want 
        % element multiplication for each column and then the sum.
    end
    
    % we can avoid the previous for loop here as well but then the code
    % will be not easily readble. 
    
    c(t) = sum(alfaTemp(:,t));  % that was totally wrong. first, c(t) - 
    % not c(2) - must be set in each iteration. Then the sum should be of 
    % the sum of the right column of alfaTemp, the one that corresponds to 
    % the current t.
    
    alfaHat(:,t) = alfaTemp(:,t)/c(t); % changed: the first part of the 
    % line should be alfaHat not alfaTemp. Also, the time pointer
    % should be added and the division should be made by the respective
    % c (i.e. c(t)) not c(1)
    % also, we can avoid the for loop here
    
    % alfaHat = [alfaHat alfaTemp']; unecessary and wrong
end

%% TERMINATION
% [rows,columns] = size(A); 
% instead of checking the rows and columns, we check the FINITE parameter
if(FINITE)
    % c(max(rows,columns)) = 0.0581; % wtf is that again???? 
    c(end) = sum(alfaHat(:,end).*A(:,end));
end
