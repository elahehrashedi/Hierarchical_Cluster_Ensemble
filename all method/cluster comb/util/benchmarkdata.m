function [X, group, maxg, namestr]=benchmarkdata(data_idx,PreProcess)
% load the ith dataset
% normalize it between [0,1]
%'banana'	'breast_cancer'	'diabetis'	'flare_solar'	'german'	'heart'	
%'image'	'ringnorm'	'splice'	'thyroid'	'titanic'	'twonorm'	'waveform'
if nargin < 2
    PreProcess=1;   %preprocess the data to [0 1] interval
end
load benchmarks;
namestr    =benchmarks{1,data_idx}
xpression  =[namestr,'.x'];
X          =eval(xpression);

xpression  = [namestr,'.t'];
target     =eval(xpression);

%data set must be of moderate size to be used in our ensembles
if (size(X,1)>1000)
   X=X(1:1000,:); 
   target=target(1:1000,:);
end

[group, gnames] = grp2idx(target);
maxg = max(group);

if (PreProcess==1)
    % means of zero and standard deviations of 1.
    p=X';
    [pn,meanp,stdp] = prestd(p);
    X=pn';
    gscatter(X(:,1),X(:,2),group);
    %each col of X must be one feature
end