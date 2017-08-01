function [ T,Zstar ] = hEnsembleCluster( X,L,c,ClustererType,ConsensusParam )
%This function creates an ensemble of hierarchical base clusterings and then combine
%them to a final h. clustering using hierarchical clustering combination
%methods and finally drive a flat clustering out of the resulted tree with
%the c cluster in it.
if nargin < 4
    ClustererTypes={'h sl','h avgl','h cl','h weighted','h centroid','h median','h ward','h all linkage','nosubsample all linkage'};
    ClustererType=ClustererTypes{9};
    %ConsensusParams={hiff,ds,linkagemethod,descriptor,matCombtype,alpha}
    %%%%===============================Information Theoretic
    ConsensusParam{1,1}= 4;                                   %hiff
    ConsensusParam{1,2}= 0;                                   %ds
    ConsensusParam{1,3}= 'average';                            %linkagemethod
    ConsensusParam{1,4}= 'cd';                                %descriptor
    ConsensusParam{1,5}= 'min';                               %matCombtype
    ConsensusParam{1,6}= 0;                                   %alpha
    %%%%===============================Fuzzy
%     ConsensusParam{1,1}= 3;                                   %hiff
%     ConsensusParam{1,2}= 0;                                  %ds
%     ConsensusParam{1,3}= 'single';                           %linkagemethod
%     ConsensusParam{1,4}= 'cd';             %descriptor
%     ConsensusParam{1,5}= '';                                  %matCombtype
%     ConsensusParam{1,6}= 0; %             %alpha
    
end
P=createEnsemble(X,L,ClustererType,c);
Zstar=allHierarchicalConsensus(P,L,c,ConsensusParam{:});
T = cluster(Zstar,'maxclust',c);
end


function P=createEnsemble(X,L,ClustererType,c)
%This function creates an ensemble of hierarchical clusterings and also
%creates their desciptors,notice that we only used the cd desc, changing to
%include other descriptors is easy.

% normalize is used to enforce that the dendrogram descriptor rows make a probability function.
% if we use information theoretic combination method this normalization is nessecary.
% in simple hierarchical comb. we did not use this normalization and the
% normalization is performed in combination algorithm. in fuzzy
% consensus we decide to normalize the descriptor here.
% in simple hier.com we fill unspecified entry of dendrograms with zero
%and in information theoretic we filled them with ones
% information theoretic
%         normalize=1 ;fillnum=1
%simple hierarchical comb
%         normalize=0 ;fillnum=0
%fuzzy comb
          normalize=2 ;fillnum=1;

%descriptors={'pmd','cmd','smd','cd','my'};
descriptors={'cd'};

P=CreateClusterer(X,L,ClustererType,c);
N=P.N;
%create dedescriptors for the ensemble
for descriptor_id=1:length(descriptors)
    for i=1:1:L
        % make the dendrogram N*N and with the same order as original data
        if (fillnum==1)
            correcteddendro=ones(N,N);
        else
            correcteddendro=zeros(N,N);
        end
        subsampleidx=P.partition{1,i}.subsampleidx;
        first=DendrogramDescriptor(P.partition{1,i}.dend, descriptors{descriptor_id});

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
        eval(sprintf('P.%s{1,i}=correcteddendro;',descriptors{descriptor_id}));
    end
end
end