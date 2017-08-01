function [subsampleidx,return_matrix,N] = Create_Clusterer1(X, group, c ,ClustererType,descriptor,w,notrand,percent,fillnum,normalize)

subsampleidx=[];

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

        L=1;
        i=1;
        j=1;
        
        [P,subsampleidx]=Create_Clusterer2(X,L,ClustererType,c,w,notrand,percent);                       
        N=P.N;  
        %correcteddendro=correct(DendrogramDescriptor(P.partition{1,i},'pmd')); 
        %make the dendrogram N*N and with the same order as original data
        if (fillnum==1)
            correcteddendro=ones(N,N);
        elseif (fillnum==0)
            correcteddendro=zeros(N,N);
        elseif (fillnum==NaN)
            correcteddendro=zeros(N,N)/0;
        end

        subsampleidx=P.partition{1,i}.subsampleidx;
        %first=DendrogramDescriptor(P.partition{1,i}.dend, descriptors{descriptor_id});
        first=DendrogramDescriptor(P.partition{1,i}.dend, descriptor);


        % ####eq.1  this is very critical expression inserted in this version
        first=first/max(first(:))*.95;


        for jj=1:1:size(subsampleidx,2)
            correcteddendro(subsampleidx(1,jj),subsampleidx)=first(jj,:);
        end

        for jj=1:size(correcteddendro,1),
            correcteddendro(jj,jj) = 0;
        end

        %how to normalize matrix
        correcteddendro = normalfunc(correcteddendro,normalize,N);
        return_matrix=correcteddendro;
end