function CreateandSaveClusterer(L,niter,ClustererTypes,descriptors,Data_Names,normalize,fillnum)
% a newer version (1.5) which has a flag for normalization
if nargin == 0
    clc; clear all;
    L=10;
    niter=1;%00; % number of iteration
	ClustererTypes={'h sl','h avgl','h cl','h weighted','h centroid','h median','h ward','h all linkage','nosubsample all linkage'};
    descriptors={'pmd','cmd','smd','cd','my'};
    Data_Names= {'breast_cancer','flare_solar','titanic','ionosphere',...
     'wpbc','image_segmentation','liver_disorders','wine',...
     'laryngeal1','laryngeal3','weaning','contractions',...
        'diabetis','german','balance_scale'...
           'vehicle','laryngeal2' };
end
if nargin < 6
   % normalize is used to enforce that the dendrogram descriptor rows make a probability function. 
   % if we use information theoretic combination method this normalization is nessecary.
   % in simple hierarchical comb. we did not use this normalization and the
   % normalization is performed in combination algorithm. in fuzzy
   % consensus we decide to normalize the descriptor here.
   % in simple hier.com we fill unspecified entry of dendrograms with zero
   %and in information theoretic we filled them with ones 
   %the default is information theoretic
   %information theoretic 
           normalize=1; fillnum=1;
   %simple hierarchical comb
%             normalize=0 ;fillnum=0;
   %fuzzy comb 
   %       normalize=2 ;fillnum=1;
end

for data_idx=1:length(Data_Names)
    data_name=Data_Names{data_idx};
    %X : data
    %group : cluster 
    %maxg : num of groups
    %N : the number of input pattern
    [X, group, maxg]=LoadBenchmarkdata(data_name,1,1);
    c=maxg;
    for j=1:niter     % number of iteration
        
        for nomofctype=1:length(ClustererTypes)  
            CType=ClustererTypes{nomofctype};
            P=CreateClusterer(X,L,CType,c);

            %filename=sprintf('C:\\Clusterers\\%s_%d.mat',data_name,j);
            filename=sprintf('F:\\eli\\cluster comb\\Clusterers\\%s_%d.mat',data_name,j);

            % ela: commented
            
            
            %==================================================================
            N=P.N;  
            for descriptor_id=1:length(descriptors) 
                %l boode felan 7
                for i=1:1:L
                    %correcteddendro=correct(DendrogramDescriptor(P.partition{1,i},'pmd')); 
                    %make the dendrogram N*N and with the same order as original data
                    if (fillnum==1)
                        correcteddendro=ones(N,N);
                    else
                        correcteddendro=zeros(N,N);
                    end

                    subsampleidx=P.partition{1,i}.subsampleidx;
                    first=DendrogramDescriptor(P.partition{1,i}.dend, descriptors{descriptor_id});



                % ####eq.1  this is very critical expression inserted in this version
                first=first/max(first(:))*.95;
                
                
                    for jj=1:1:size(subsampleidx,2)
                        correcteddendro(subsampleidx(1,jj),subsampleidx)=first(jj,:);
                    end
                    
                    for jj=1:size(correcteddendro,1),
                        correcteddendro(jj,jj) = 0;
                    end

                    if normalize==1
                        correcteddendro=correcteddendro./repmat(sum(correcteddendro,2),1,N);
                    end
                    if normalize==2
                        correcteddendro=correcteddendro/max(correcteddendro(:));
                    end

                    %descriptors are 'pmd','cmd','smd','cd','my'
                    eval(sprintf('P.%s{1,i}=correcteddendro;',descriptors{descriptor_id}));
                end
            end
            %==================================================================
%             Y=squareform(P.pmd{1,i});         %is similar to Y = pdist(X);
%             Zstar = linkage(Y,'single');
%             acc1=HierCombMethodAccuracy(Zstar,group,X,'cophenet');
% 
%             
%             acc3=HierCombMethodAccuracy(P.partition{1,i}.dend,group,X,'cophenet');
            %==================================================================           
            
            save(filename, 'P','L','c','X','group');
        end;
    end
end

