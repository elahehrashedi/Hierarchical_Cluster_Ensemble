function [P,subsampleidx]=Create_Clusterer2(X,L,Type,c,w,notrand,percent)
% Partitions the points in the data matrix
% X into c clusters( if c is determied first otherwise randomely specify it) L times and reurn the ensemble result
% if type=random then the X value will not be used (only its dimension will be used ) and L random clustering will be returned

subsampleidx=[];
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
    L=10;             %set the number of clusterer equal to the number of diffrent possible hierarchical methods and distance metrics
end
%eli added
if  isequal(Type,'nosubsample average linkage')
    L=10;             %set the number of clusterer equal to the number of diffrent possible hierarchical methods and distance metrics
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
            
            %==================================
            [subsampleidx,w] = sampler (r,w,notrand,percent);
            %perm = randperm(r);
            %subsampleidx=perm(1:round(0.8*r));
            %==================================
            
            % use randomly selected 80% data to create cluster
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
        [P,subsampleidx] = linkage_cluster(X,L,c,'single',w,notrand,percent);

    case 'abs avgl'
        disp('abstract average linkage+subsample clusterers')
        [P,subsampleidx] = linkage_cluster(X,L,c,'average',w,notrand,percent);

    case 'abs cl'
        disp('abstract complete linkage+subsample clusterers')
        [P,subsampleidx] = linkage_cluster(X,L,c,'complete',w,notrand,percent);

    case 'abs weighted'
        disp('abstract weighted linkage+subsample clusterers')
        [P,subsampleidx] = linkage_cluster(X,L,c,'weighted',w,notrand,percent);

    case 'abs centroid'
        disp('abstract centroid linkage+subsample clusterers')
        [P,subsampleidx] = linkage_cluster(X,L,c,'centroid',w,notrand,percent);

    case 'abs median'
        disp('abstract median linkage+subsample clusterers')
        [P,subsampleidx] = linkage_cluster(X,L,c,'median',w,notrand,percent);

    case 'abs ward'
        disp('abstract ward linkage+subsample clusterers')
        [P,subsampleidx] = linkage_cluster(X,L,c,'ward',w,notrand,percent);

    case 'abs all linkage'
        disp('abstract all linkage+subsample clusterers')
        [P,subsampleidx] = linkage_cluster(X,L,c,'all',w,notrand,percent);

    case 'h sl'
        disp('hierarchical single linkage+subsample clusterers')
        [P,subsampleidx] = linkage_dendro(X,L,c,'single',w,notrand,percent);

    case 'h avgl'
        disp('hierarchical average linkage+subsample clusterers')
        [P,subsampleidx] = linkage_dendro(X,L,c,'average',w,notrand,percent);

    case 'h cl'
        disp('hierarchical complete linkage+subsample clusterers')
        [P,subsampleidx] = linkage_dendro(X,L,c,'complete',w,notrand,percent);

    case 'h weighted'
        disp('hierarchical weighted linkage+subsample clusterers')
        [P,subsampleidx] = linkage_dendro(X,L,c,'weighted',w,notrand,percent);

    case 'h centroid'
        disp('hierarchical centroid linkage+subsample clusterers')
        [P,subsampleidx] = linkage_dendro(X,L,c,'centroid',w,notrand,percent);

    case 'h median'
        disp('hierarchical median linkage+subsample clusterers')
        [P,subsampleidx] = linkage_dendro(X,L,c,'median',w,notrand,percent);

    case 'h ward'
        disp('hierarchical ward linkage+subsample clusterers')
        [P,subsampleidx] = linkage_dendro(X,L,c,'ward',w,notrand,percent);

    case 'h all linkage'
        disp('hierarchical all linkage+subsample clusterers')
        [P,subsampleidx] = linkage_dendro(X,L,c,'all',w,notrand,percent);

    case 'nosubsample all linkage'
        disp('hierarchical all linkage+with out subsampling clusterers')
        P = nosubsample_linkage_dendro(X);
    %eli added
    case 'nosubsample single linkage'
        disp('hierarchical single linkage+with out subsampling clusterers')
        P = nosubsample_linkage_dendro_2(X,1);       
    %eli added
    case 'nosubsample average linkage'
        disp('hierarchical average linkage+with out subsampling clusterers')
        P = nosubsample_linkage_dendro_2(X,2); 
    %eli added
    case 'nosubsample complete linkage'
        disp('hierarchical complete linkage+with out subsampling clusterers')
        P = nosubsample_linkage_dendro_2(X,3); 
    %eli added
    case 'nosubsample weighted linkage'
        disp('hierarchical weighted linkage+with out subsampling clusterers')
        P = nosubsample_linkage_dendro_2(X,4);
            %eli added
    case 'nosubsample centroid linkage'
        disp('hierarchical centroid linkage+with out subsampling clusterers')
        P = nosubsample_linkage_dendro_2(X,5);
            %eli added
    case 'nosubsample median linkage'
        disp('hierarchical median linkage+with out subsampling clusterers')
        P = nosubsample_linkage_dendro_2(X,6);
            %eli added
    case 'nosubsample ward linkage'
        disp('hierarchical ward linkage+with out subsampling clusterers')
        P = nosubsample_linkage_dendro_2(X,7);
        
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

w1=w;

end

%%%%%%%%%%%%%%%%%%%%%%%%%  linkage clustering with subsample%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,subsampleidx] = linkage_cluster(X,L,c,method,w,notrand,percent)


if isequal(method,'all')            %create an ensemble using different hierarchical methods
    %P = repmat(struct('partition',NaN*zeros(size(X,1),1)), 1, L); %The ensemble P
    r=size(X,1);
    methods={'single','complete','average','weighted','centroid', 'median','ward'};
    for i=1:1:L
        
        %==================================
        [subsampleidx,w] = sampler (r,w,notrand,percent);
        %perm = randperm(r);
        %subsampleidx=perm(1:round(0.8*r));
        %==================================
        
        % use randomly selected 80% data to create cluster
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
        
        %==================================
        [subsampleidx,w] = sampler (r,w,notrand,percent);
        %perm = randperm(r);
        %subsampleidx=perm(1:round(0.8*r));
        %==================================
            
        % use randomly selected 80% data to create cluster
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%  linkage dendrogram with subsample%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   T = CLUSTERDATA(X,CUTOFF) is the same as
%      Y = pdist(X, 'euclid');
%      Z = linkage(Y, 'single');
%      T = cluster(Z, 'cutoff', CUTOFF);
%
function [P,subsampleidx] = linkage_dendro(X,L,c,method,w,notrand,percent)


if isequal(method,'all')            %create an ensemble using different hierarchical methods
    %P = repmat(struct('partition',NaN*zeros(size(X,1),1)), 1, L); %The ensemble P
    r=size(X,1);
    methods={'single','complete','average','weighted','ward'};%'centroid',, 'median'
    for i=1:1:L
        %==================================
        [subsampleidx,w] = sampler (r,w,notrand,percent);
        %perm = randperm(r);
        %subsampleidx=perm(1:round(0.8*r));
        %==================================

        % use randomly selected 80% data to create cluster
        % the cluster of th eremaning data will set to NaN
        % this causee their cluster be different to each other
        %IDX = clusterdata(X(subsampleidx,:),'linkage',methods(mod(i,7)+1),'maxclust',c(1,i));
        Y = pdist(X(subsampleidx,:), 'euclid');
        
        %ela added random
        %Z = linkage(Y, methods(mod(i,7)+1));
        n = 5;
        f = ceil(n.*rand(100,1));        
        Z = linkage(Y, methods(f(1)));
        display(methods(f(1)));

        P.partition{1,i}.dend=Z;
        P.partition{1,i}.subsampleidx=subsampleidx;
        P.type='hierarchical';
        P.N=r;
    end
else
    %P = repmat(struct('partition',NaN*zeros(size(X,1),1)), 1, L); %The ensemble P
    r=size(X,1);
    for i=1:1:L

        %==================================
        [subsampleidx,w] = sampler (r,w,notrand,percent);
        %perm = randperm(r);
        %subsampleidx=perm(1:round(0.8*r));
        %==================================
            
        % use randomly selected 80% data to create cluster
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
end
%%%%%%%%%%%%%%%%%%%%%%%%% single linkage dendrogram without subsample%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = nosubsample_linkage_dendro_2(X,methodnum)


r=size(X,1);
methods={'single','average','complete','weighted','centroid', 'median','ward'};
dis_metrics={'euclidean','cityblock','seuclidean','mahalanobis','minkowski',...
    'cosine','correlation','hamming','jaccard','chebychev'};
k=0;
for i=1:1
    for j=methodnum:methodnum%0 %create just one auclidean
        k=k+1;
        Y = pdist(X, dis_metrics{i});
        Z = linkage(Y, methods{j});
        % part: single,elu + single,city , single,
        P.partition{1,k}.dend=Z;
        P.partition{1,k}.subsampleidx=1:r;
        P.type='hierarchical';
        P.N=r;
    end
end
end

% %%%%%%%%%%%%%%%%%%%%%%%%% sampler %%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [subsampleidx,w] = sampler (r,w,notrand,percent)
% 
% b= percent;%selection
% a=1-b;
% 
% %w is not changed
% %just send out subsample
% %upper weighted patterns have lower chance to be selected
% 
% if (notrand<1) %a*r or 1
%         display('random selection');
%         w=w/sum(w);
%         perm = randperm(r);
%         subsampleidx=perm(1:round(b*r));
% else
% %             display('random selection');
% %         w=w/sum(w);
% %         perm = randperm(r);
% %         subsampleidx=perm(1:round(b*r));
%         
%         display('weighted selection');
%         minimum = min(w);
%         if (minimum<0)
%             added=zeros(1,r);
%             added(1,:) = -minimum; 
%             w=w+added;
%         end
% 
% 
%         w = w/sum(w);
%         
%         %revers w 
%         %ww = ones(1,r) - w ;
%         ww=ones(1,r)+w;
%         
%         n=1;
%         for i=1:r
%             if(ww(1,i)~=0)
%                 for j=1:ww(1,i)*100
%                     p(n)=i;
%                     n=n+1;
%                 end
%             else
%                 p(n)=i;
%                 n=n+1;               
%             end
%         end
%         sublenght=round(b*r);
%         i=1;
%         while i<sublenght
%             j=randperm(size(p,2));
%             sample=p(1,j(1));
%             subsampleidx(1,i)=sample;
%             %remove a particular value from a vector
%             index = find(p == sample );
%             p(index) = [];
%             i = i+1;
%         end
% end

%%%%%%%%%%%%%%%%%%%%%%%%% sampler %%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is a new sampler
function [subsample,w] = sampler (r,w,notrand,percent)
b= percent;%selection
a=1-b;
n=size(w,2);
subsample=[];

%w is not changed
%just send out subsample
%upper weighted patterns have lower chance to be selected
    if (notrand==0) %a*r or 1
        %display('random selection');
        w=w/sum(w);
        perm = randperm(r);
        subsample=perm(1:round(b*r));
    elseif notrand==1    
        %display('weighted selection');
        minimum = min(w);
        wsum=w/sum(w);
        if (minimum<0)
            w2=w-minimum;
            w2 = w2/sum(w2);
            wsum=w2;
        end

        for i=2:n
            wsum(i)=wsum(i-1)+ wsum(i);   
        end
        
        allsamples=combntns(1:n,1);%integers between 1 to num of
        
        %&&&&&&&&&&&&&&&&&&&&&&&&&&
        %&&&&&&&&&&&&&&&&&&&&&&&&&&
        for i=1:b*r
            %display(i);
            randnumber=rand();
            flag=0;%sample not found
            for j=size(wsum,2):-1:1
                if (flag==0) && (randnumber>wsum(j))
                    %display('2sample');display(allsamples(j+1));
                    %j+1 is selected
                    %add sample to subsamples
                    subsample=[subsample,allsamples(j+1)];
                    %change wsum from sample j+1, sub w(j)
                    newlenght=size(wsum,2);
                    if j+1 < newlenght
                        wsum(j+1:newlenght)=wsum(j+1:newlenght)-w(allsamples(j+1));
                    end
                    %delete sample from wsum
                    wsum(j+1)=[];
                    wsum = wsum/max(wsum);
                    %delete sample from allsamples
                    allsamples(j+1)=[]; 
                    %break the loop
                    flag=1;%sample found
                end
            end
            if ((flag==0) && (j == 1))
                %display('1sample');display(allsamples(1));                   
                %1 is selected
                %add sample to subsamples
                subsample=[subsample,allsamples(1)];
                %change wsum from sample j+1, sub w(j)
                newlenght=size(wsum,2);
                if 1 < newlenght
                    wsum(1:newlenght)=wsum(1:newlenght)-w(allsamples(1));
                end
                %delete sample from wsum
                wsum(1)=[];
                wsum = wsum/max(wsum);
                %delete sample from allsamples
                allsamples(1)=[]; 
            end
        end
        %&&&&&&&&&&&&&&&&&&&&&&&&&&
        %&&&&&&&&&&&&&&&&&&&&&&&&&&
    end
end


