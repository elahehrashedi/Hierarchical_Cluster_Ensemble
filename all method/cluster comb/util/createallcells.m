function A=createallcells(k,Dim)
%function that defines cells mf for a Dim dimensional problem and number them
%Dim number of feature
%k   number of mf for each feature
m=(1:k)';
for D=2:Dim
    res=[];
    for i=1:k
        tmp=[i*ones(size(m,1),1) m];
        res=[res;tmp];
    end
    m=res;
end
A=m;