function [X, group, maxg, namestr]=benchmarkdata_mock(data_idx,PreProcess)
% load the ith dataset
% normalize it between [0,1]
if nargin < 2
    PreProcess=1;   %preprocess the data to [0 1] interval
end
data_names={'fourty','long1','sizes5','spiral','square1', 'square4', 'twenty'};
load mock;
namestr=data_names{data_idx}
xpression  =namestr;
X          =eval(xpression);
target     =X(:,3);
X=X(:,1:2);

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