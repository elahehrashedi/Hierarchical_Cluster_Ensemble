function P=CreateClusterer(X,L,Type,c)
% Partitions the points in the data matrix
% X into c clusters( if c is determied first otherwise randomely specify it) L times and reurn the ensemble result
% if type=random then the X value will not be used (only its dimension will be used ) and L random clustering will be returned

if nargin < 1
    load kmeansdata;
    size(X);
    X=X(1:10,:);
    c=3;
    L=2;              %the number of learner in ensemble
    Type='all linkage+subsample';       %out put of clustring will be at abstract level(only the cluster number)
end
if  isequal(Type,'abs all linkage')
    L=7;             %set the number of clusterer equal to the number of diffrent possible hierarchical methods
end
if  isequal(Type,'nosubsample all linkage')
    L=10;             %set the number of clusterer equal to the number of diffrent possible hierarchical methods and distance metrics
end
%eli added
if  isequal(Type,'nosubsample single linkage')
    L=7;             %set the number of clusterer equal to the number of diffrent possible hierarchical methods and distance metrics
end

np=size(X,1);        %the number of input pattern

if nargin < 4        %if the number of cluster is not determined specify it randomly in range a,b
    rand('state',sum(100*clock));
    a = 2; b = 22;        % the range for cluster numbers
    c = a + (b-a) * rand(1,L);
    c=round(c);
else                % if the number is specified covert it to matrix format
    c=c*ones(1,L);
end

switch lower(Type)
    case 'k-means'
        disp('k-means+random initialization clusterers')
        %P = repmat(struct('partition',zeros(size(X,1),1)), 1, L); %The ensemble P
        for i=1:1:L
            IDX = kmeans(X,c(1,i),'emptyaction','singleton');
            %P(1,i).partition=IDX;
            P.partition{1,i}=IDX;
            P.type='abstract';
            P.N=np;
        end

    case 'k-means+subsample'
        disp('k-means+random initialization+subsample clusterers')
        %P = repmat(struct('partition',NaN*zeros(size(X,1),1)), 1, L); %The ensemble P
        r=size(X,1);
        for i=1:1:L
            perm = randperm(r);
            subsampleidx=perm(1:round(0.8*r));% use randomly selected 80% data to create cluster
            % the cluster of th eremaning data will set to NaN
            % this causee their cluster be different to each other in
            % direct optimization for example getting mean cause the
            % resulting pstar to be nan to
            % we can also run a kmeans then using its centers and a diffrent kmeans with one iteration to
            % force all instances have a clustering
            [IDX,Center]= kmeans(X(subsampleidx,:),c(1,i),'emptyaction','singleton');
            %P(1,i).partition(subsampleidx,:)=IDX;
            IDX = kmeans(X,c(1,i),'emptyaction','singleton','maxiter',10,'start',Center);
            %P(1,i).partition=IDX;
            P.partition{1,i}=IDX;
            P.type='abstract';
            P.N=np;
        end
    case 'random'
        disp('random clusterers')
        %P = repmat(struct('partition',zeros(size(X,1),1)), 1, L); %The ensemble P
        for i=1:1:L
            a = 1; b = c(1,i);        % the range for cluster numbers
            IDX = a + (b-a) * rand(size(X,1),1);
            IDX=round(IDX);
            %P(1,i).partition=IDX;
            P.partition{1,i}=IDX;
            P.type='abstract';
            P.N=np;
        end
    case 'fcm'
        disp('Fuzzy c-means measurement ensemble')
        c=c(1,1);   %in this case c can not be different between clusterer
        %P = repmat(struct('partition',zeros(size(X,1),c)), 1, L); %The ensemble P
        for i=1:1:L
            [center, U, obj_fcn] = fcm(X, c,[NaN,NaN,NaN,0]);     %set the display off
            %P(1,i).partition=U';  %membership
            P.partition{1,i}=U';  %membership
            P.type='measurement';
            P.N=np;
        end
    case 'abs sl'
        disp('abstract single linkage+subsample clusterers')
        P = linkage_cluster(X,L,c,'single');

    case 'abs avgl'
        disp('abstract average linkage+subsample clusterers')
        P = linkage_cluster(X,L,c,'average');

    case 'abs cl'
        disp('abstract complete linkage+subsample clusterers')
        P = linkage_cluster(X,L,c,'complete');

    case 'abs weighted'
        disp('abstract weighted linkage+subsample clusterers')
        P = linkage_cluster(X,L,c,'weighted');

    case 'abs centroid'
        disp('abstract centroid linkage+subsample clusterers')
        P = linkage_cluster(X,L,c,'centroid');

    case 'abs median'
        disp('abstract median linkage+subsample clusterers')
        P = linkage_cluster(X,L,c,'median');

    case 'abs ward'
        disp('abstract ward linkage+subsample clusterers')
        P = linkage_cluster(X,L,c,'ward');

    case 'abs all linkage'
        disp('abstract all linkage+subsample clusterers')
        P = linkage_cluster(X,L,c,'all');

    case 'h sl'
        disp('hierarchical single linkage+subsample clusterers')
        P = linkage_dendro(X,L,c,'single');

    case 'h avgl'
        disp('hierarchical average linkage+subsample clusterers')
        P = linkage_dendro(X,L,c,'average');

    case 'h cl'
        disp('hierarchical complete linkage+subsample clusterers')
        P = linkage_dendro(X,L,c,'complete');

    case 'h weighted'
        disp('hierarchical weighted linkage+subsample clusterers')
        P = linkage_dendro(X,L,c,'weighted');

    case 'h centroid'
        disp('hierarchical centroid linkage+subsample clusterers')
        P = linkage_dendro(X,L,c,'centroid');

    case 'h median'
        disp('hierarchical median linkage+subsample clusterers')
        P = linkage_dendro(X,L,c,'median');

    case 'h ward'
        disp('hierarchical ward linkage+subsample clusterers')
        P = linkage_dendro(X,L,c,'ward');

    case 'h all linkage'
        disp('hierarchical all linkage+subsample clusterers')
        P = linkage_dendro(X,L,c,'all');

    case 'nosubsample all linkage'
        disp('hierarchical all linkage+with out subsampling clusterers')
        P = nosubsample_linkage_dendro(X);
    %eli added
    case 'nosubsample single linkage'
        disp('hierarchical single linkage+with out subsampling clusterers')
        P = nosubsample_linkage_dendro_2(X);        

    case 'random dendrogram'
        disp('random dendrograms')
        for i=1:1:L
            P.partition{1,i}.dend=RandomAgglomeration(np);
            P.type='hierarchical';
            P.partition{1,i}.subsampleidx=1:np;
            P.N=np;
        end


    otherwise
        error('Unknown type');
end


%%%%%%%%%%%%%%%%%%%%%%%%%  linkage clustering with subsample%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = linkage_cluster(X,L,c,method)
if isequal(method,'all')            %create an ensemble using different hierarchical methods
    %P = repmat(struct('partition',NaN*zeros(size(X,1),1)), 1, L); %The ensemble P
    r=size(X,1);
    methods={'single','complete','average','weighted','centroid', 'median','ward'}
    for i=1:1:L
        perm = randperm(r);
        subsampleidx=perm(1:round(0.8*r));% use randomly selected 80% data to create cluster
        % the cluster of th eremaning data will set to NaN
        % this causee their cluster be different to each other
        IDX = clusterdata(X(subsampleidx,:),'linkage',methods(mod(i,7)+1),'maxclust',c(1,i));
        P.partition{1,i}=NaN*zeros(size(X,1));
        P.partition{1,i}(subsampleidx,:)=IDX;
        P.type='abstract';
        P.N=r;
    end
else
    %P = repmat(struct('partition',NaN*zeros(size(X,1),1)), 1, L); %The ensemble P
    r=size(X,1);
    for i=1:1:L
        perm = randperm(r);
        subsampleidx=perm(1:round(0.8*r));% use randomly selected 80% data to create cluster
        % the cluster of th eremaning data will set to NaN
        % this causee their cluster be different to each other
        IDX = clusterdata(X(subsampleidx,:),'linkage',method,'maxclust',c(1,i));
        %P(1,i).partition(subsampleidx,:)=IDX;
        P.partition{1,i}=NaN*zeros(size(X,1));
        P.partition{1,i}(subsampleidx,:)=IDX;
        P.type='abstract';
        P.N=r;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%  linkage dendrogram with subsample%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   T = CLUSTERDATA(X,CUTOFF) is the same as
%      Y = pdist(X, 'euclid');
%      Z = linkage(Y, 'single');
%      T = cluster(Z, 'cutoff', CUTOFF);
%
function P = linkage_dendro(X,L,c,method)
if isequal(method,'all')            %create an ensemble using different hierarchical methods
    %P = repmat(struct('partition',NaN*zeros(size(X,1),1)), 1, L); %The ensemble P
    r=size(X,1);
    methods={'single','complete','average','weighted','centroid', 'median','ward'};
    for i=1:1:L
        perm = randperm(r);
        subsampleidx=perm(1:round(0.8*r));% use randomly selected 80% data to create cluster
        % the cluster of th eremaning data will set to NaN
        % this causee their cluster be different to each other
        %IDX = clusterdata(X(subsampleidx,:),'linkage',methods(mod(i,7)+1),'maxclust',c(1,i));
        Y = pdist(X(subsampleidx,:), 'euclid');
        Z = linkage(Y, methods(mod(i,7)+1));

        P.partition{1,i}.dend=Z;
        P.partition{1,i}.subsampleidx=subsampleidx;
        P.type='hierarchical';
        P.N=r;
    end
else
    %P = repmat(struct('partition',NaN*zeros(size(X,1),1)), 1, L); %The ensemble P
    r=size(X,1);
    for i=1:1:L
        perm = randperm(r);
        subsampleidx=perm(1:round(0.8*r));% use randomly selected 80% data to create cluster
        % the cluster of th eremaning data will set to NaN
        % this causee their cluster be different to each other
        Y = pdist(X(subsampleidx,:), 'euclid');
        Z = linkage(Y, method);
        P.partition{1,i}.dend=Z;
        P.partition{1,i}.subsampleidx=subsampleidx;
        P.type='hierarchical';
        P.N=r;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%  linkage dendrogram without subsample%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = nosubsample_linkage_dendro(X)
r=size(X,1);
methods={'single','complete','average','weighted','ward','centroid', 'median'};
dis_metrics={'euclidean','cityblock','seuclidean','mahalanobis','minkowski',...
    'cosine','correlation','hamming','jaccard','chebychev'};
k=0;
for i=1:7
    for j=1:1
        k=k+1;
        Y = pdist(X, dis_metrics{j});
        Z = linkage(Y, methods{i});
        %p.part1,1= single , part1,2=complete
        P.partition{1,k}.dend=Z;
        P.partition{1,k}.subsampleidx=1:r;
        P.type='hierarchical';
        P.N=r;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%% single linkage dendrogram without subsample%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = nosubsample_linkage_dendro_2(X)
r=size(X,1);
methods={'single'};
dis_metrics={'euclidean','cityblock','seuclidean','mahalanobis','minkowski',...
    'cosine','correlation','hamming','jaccard','chebychev'};
k=0;
for i=1:1
    for j=1:10
        k=k+1;
        Y = pdist(X, dis_metrics{j});
        Z = linkage(Y, methods{i});
        % part: single,elu + single,city , single,
        P.partition{1,k}.dend=Z;
        P.partition{1,k}.subsampleidx=1:r;
        P.type='hierarchical';
        P.N=r;
    end
end