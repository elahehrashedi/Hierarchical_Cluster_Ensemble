function d=MeasurementLevel(p,c)
%if we use none fuzzy clustring we must convert abstract level to vector for example 2=0 1 0 0 0
%p: the clustering(classifier) output in abstract level
%d: the clustering(classifier) output in measurement level
%c: number of cluster or class
if (size(p,2)~=1) %if it is not in abstract level
    error('may be the partition is currently in measurementlevel');
end
d=zeros(size(p,1),c);   %creating d of pstar
for i=1:c
    a=zeros(1,c);
    a(1,i)=1;
    d(find(p==i),:)=repmat(a,size(find(p==i),1),1);
end