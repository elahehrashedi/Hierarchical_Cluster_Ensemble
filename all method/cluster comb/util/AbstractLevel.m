function d=AbstractLevel(mu,c)
%convert measurement level to abstract level for example 2=abstractlevel(0 1 0 0 0)
%p: the clustering(classifier) output in measuement level
%d: the clustering(classifier) output in abstract level
%c: number of cluster or class
if (size(mu,2)==1) %if it is not in measurementlevel
    error('may be the partition is currently in abstractlevel');
end
[Y,d] = max(mu,[],2);
