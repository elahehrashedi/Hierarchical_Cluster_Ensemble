function desmat=DendrogramDescriptor(Z,type)
% generates a descriptor matrix of the hierarchical, binary cluster tree represented by Z.
% Z is an m-1-by-3 matrix, generated by the linkage function, or randomagglomeration 
% where m is the number of objects in the original data set
if nargin < 1
    Z=[1,2,2;
        4,5,3;
        8,3,4;
        6,7,5;
        10,9,7;
        12,11,10;]
    type='cd';
    dendrogram(Z);
end

switch lower(type)
    case 'pmd'
        desmat=pmd(Z);
    case 'cmd'
        desmat=cmd(Z);
    case 'smd'
        desmat=smd(Z);
    case 'cd'
        desmat=cd(Z);
    case 'pd'
        %         desmat=pd(Z);
        error('not written yet');
    case 'my'    
        desmat=my(Z);
    otherwise
        error('unknown DendrogramDescriptor');
end

%---------------------------------------------------------------------------------------------------
function pmdmatrix=pmd(Z)
% generates a pmd (partition membership divergence) matrix of the hierarchical, binary cluster tree 
% represented by Z. Th enumber of partition(level) that jk are toghether are counted and is
% subtracted from m possible partition to calculate the # of level that jk are not toghether
m=size(Z,1)+1;
SimMatrix=zeros(m,m);
for i=1:m
    newp= cluster(Z,'maxclust',i); %create the ith partition with i cluster in it
    if (max(newp)~=i)
        warning('the number of cluster in step i is not i');
    end
    for j=1:i
        ii=find (newp==j);
        for k=1:size(ii,1)
            SimMatrix(ii(k,1),ii)=SimMatrix(ii(k,1),ii)+1;
        end
    end
end
pmdmatrix=ones(m,m)*m-SimMatrix;
%---------------------------------------------------------------------------------------------------
function cmdmatrix=cmd(Z)
% generates a cmd (cluster membership divergence) matrix of the hierarchical, binary cluster tree 
% represented by Z. The clusters(nodes) are created and the number of item in each of them is
% reperesented as diffrence for the newly combined members. because the cluster size is repeatedly
% increased if we compute the min between the former cmdmatrix and the new calculated one the items
% that are previously joined are not changed. 
m = size(Z,1)+1;
cmdmatrix=ones(m,m)*inf;
for i=1:1:m
    cluster{i}.X=[i];     
    cmdmatrix(i,i)=1;     %every item has a cluster only one itm
end
for k=1:1:m-1
    i=Z(k,1); j=Z(k,2);
    X=[cluster{i}.X, cluster{j}.X];
    cluster{m+k}.X=X;
    len=size(X,2);       %for every node find the cluster size
    for i=1:len
        cmdmatrix(X(1,i),X)=min(cmdmatrix(X(1,i),X),len); %set the dissimilarity according to min cluster size having jk together
    end
end

%---------------------------------------------------------------------------------------------------
function smdmatrix=smd(Z)
% generates a smd (subtree membership divergence) matrix of the hierarchical, binary cluster tree 
% represented by Z. we create the subtrees (nodes) and count in which of them the jk are combined
% toghether and subtract it from the number of possible node to compute the # of node in which they
% are not toghether.
m = size(Z,1)+1;
SimMatrix=zeros(m,m);
for i=1:1:m
    cluster{i}.X=[i];     
end
for k=1:1:m-1
    i=Z(k,1); j=Z(k,2);
    X=[cluster{i}.X, cluster{j}.X];
    cluster{m+k}.X=X;
    len=size(X,2);       %for every node find the cluster size
    %for every node check to see if jk are together
    for i=1:len
        SimMatrix(X(1,i),X)=SimMatrix(X(1,i),X)+1;
    end
end
smdmatrix=ones(m,m)*(m-1)-SimMatrix;
%---------------------------------------------------------------------------------------------------
function cdmatrix=cd(Z)
% generates a cd (cophenetic diffrence) matrix of the hierarchical, binary cluster tree 
% represented by Z. we create the nodes in the tree and set the diss according to the min level
% which jk are combined
m = size(Z,1)+1;
cdmatrix=ones(m,m)*inf;
for i=1:1:m
    cluster{i}.X=[i];     
    cdmatrix(i,i)=0;     %every comine itself at level 0
end
for k=1:1:m-1
    i=Z(k,1); j=Z(k,2);
    X=[cluster{i}.X, cluster{j}.X];
    cluster{m+k}.X=X;
    len=size(X,2);       %for every node find the cluster size
    for i=1:len
        cdmatrix(X(1,i),X)=min(cdmatrix(X(1,i),X),Z(k,3)); %set the dissimilarity according to min cluster level that combine the item jk
    end
end

%---------------------------------------------------------------------------------------------------
function mymatrix=my(Z)
% generates a cd (cophenetic diffrence) matrix of the hierarchical, binary cluster tree 
% represented by Z.
m = size(Z,1)+1;
mymatrix=zeros(m,m);
for i=1:1:m
    cluster{i}.X=[i];     
    cluster{i}.id=0;
end
for k=1:1:m-1
    i=Z(k,1); j=Z(k,2);
    cluster{m+k}.X=[cluster{i}.X, cluster{j}.X];     
    id=max(cluster{i}.id, cluster{j}.id)+1;
    cluster{m+k}.id=id;
    mymatrix(cluster{i}.X,cluster{j}.X)=id;
    mymatrix(cluster{j}.X,cluster{i}.X)=id;
end