% function P = ConvertEnsemble2myFormat(cls)
%
% DESCRIPTION 
% CreateClusterEnsemble and CreateClusterer have a structure syntax but
% StrehlFormat uses cls = [cl1; cl2;...]note that StrehlFormat can not
% handle measurement level clustering such as Fuzzy c-means clustering
% Description
% P the structure containing clusterer result
% cls clustering results every on in one row
function P= ConvertEnsemble2myFormat(cls)
for i=1:size(cls,1)
    %P(1,i).partition=cls(i,:)';
    P.partition{1,i}=cls(i,:)';
end
P.type='abstract';