function  S = SignificantDiff( x1,s1,x2,s2 )
%This function test to see if x1 is significantly better than x2. X1, x2 
%are the mean of measurements ,S1 and s2 are the respective standard deviation.
% Significant difference is observed if the respective confidence interval 
%have no intersection. 
%I think that the squre of 100 in formula is the number of experiments.
if (x1-0.196*s1>x2-0.196*s2)
    S=1;
elseif (x1-0.196*s1<x2-0.196*s2)
    S=-1;
else
    S=0;
end
