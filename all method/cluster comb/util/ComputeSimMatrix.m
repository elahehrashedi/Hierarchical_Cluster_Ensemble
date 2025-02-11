function [SimMatrix,DiffMatrix]=ComputeSimMatrix(P,L,c)
Errorcheck=1;
switch lower(P.type)
    case 'hierarchical'
        for i=1:1:L
            Z=P.partition{1,i}.dend;
            subsampleidx=P.partition{1,i}.subsampleidx;
            P.partition{1,i}=NaN*zeros(P.N,1);
            P.partition{1,i}(subsampleidx,:)=cluster(Z,'maxclust',c,'criterion','distance'); 
            if (max(P.partition{1,i})~=c) error('error converting dendrogram to cluster');end
        end
    case 'abstract'
        %nop
    case 'measurement'
%        if (size(P(1,1).partition,2)~=1)   % if the input is in measurementlevel
            % convert it to abstractlevel
            disp('this function can work with abstract level only');
            disp('converting to abstract...');
            for i=1:1:L
                P.partition{1,i}=AbstractLevel(P.partition{1,i},c);
            end
%        end
    case 'strehl'
        %nop
    otherwise
end
if (~exist('c')||(Errorcheck==1))
    cls = ConvertEnsemble2StrehlFormat(P,L);
end
if ~exist('c'),
 c = max(max(cls));  
end;

%my SimMatrix calculation only works when the number of cluster in the
%clusterer is a predetermined value c so in cases that they where not equal
%the strehl version is used
% in the direct optimization is also when the number of cluster differ may
% cause division by zero

%N=size(P(1,1).partition,1); %the number of training samples
N=P.N;                       %the number of training samples
SimMatrix=zeros(N,N);
for i=1:L
    %newp=P(1,i).partition;
    newp=P.partition{1,i};
    for j=1:c
        ii=find (newp==j);
        for k=1:size(ii,1)
            SimMatrix(ii(k,1),ii)=SimMatrix(ii(k,1),ii)+1;
        end
    end
end
SimMatrix=SimMatrix/L;
SimMatrix=checks(SimMatrix);
DiffMatrix=1-SimMatrix;

if (Errorcheck)    %error cheking with Strehl code 
    clbs = clstoclbs(cls);
    s = clbs' * clbs;
    s = checks(s./size(cls,1));
    if (length(find(s~=SimMatrix))>0)
        warning(' error in calculating simmatrix Strehl and my sim matrix are not the same');
        SimMatrix=s;
        DiffMatrix=1-SimMatrix;
    end
end
