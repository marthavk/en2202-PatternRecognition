function R=rand(pD,nData)
%R=rand(pD,nData) returns random scalars drawn from given Discrete Distribution.
%
%Input:
%pD=    DiscreteD object
%nData= scalar defining number of wanted random data elements
%
%Result:
%R= row vector with integer random data drawn from the DiscreteD object pD
%   (size(R)= [1, nData]
%
%----------------------------------------------------
%   Sampling from a discrete distribution: 
% Created by Dahua Lin, On Oct 27, 2008
% Modified by Martha Vlachou-Konchylaki, Gerasimos Markatos, On Sep 25,
% 2014
%----------------------------------------------------
K = numel(pD);
p=pD.ProbMass';
if K>1
    error('Method works only for a single DiscreteD object');
end;



%*** Insert your own code here and remove the following error message 




% construct the bins

edges = [0, cumsum(p)];
s = edges(end);
if abs(s - 1) > eps
    edges = edges * (1 / s);
end

% draw bins

rv = rand(1, nData);
c = histc(rv, edges);
ce = c(end);
c = c(1:end-1);
c(end) = c(end) + ce;

% extract samples

xv = find(c);

if numel(xv) == nData  % each value is sampled at most once
    x = xv;
else                % some values are sampled more than once
    xc = c(xv);
    d = zeros(1, nData);
    dv = [xv(1), diff(xv)];
    dp = [1, 1 + cumsum(xc(1:end-1))];
    d(dp) = dv;
    x = cumsum(d);
end

% randomly permute the sample's order
R = x(randperm(nData));

end


