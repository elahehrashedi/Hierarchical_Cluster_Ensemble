%function Z=HierarchicalConsensus(P,L,c,hiff,ds,linkagemethod,descriptor,combtype,alpha)
% this function create one hierarchical clustering from the clustering ensemble
% using the coverted method specified in kuncheva paper entitled "Experimental
% comparison of cluster ensemble methods"
% Description
% P             %the clustering ensemble P.partition{1,i} is the ith hierarchical clustering
% L=10;         %the number of learner in ensemble 
% c             %the number of cluster
% ds            %the DiffMatrix is assumed as a similarity or data (1=data,0=sim)
% hf            %the hierarchical/flat method is used for creating DiffMatrix 
% linkagemethod %the linkagemethod which is used in creating final clustering
% descriptor    %is used in hierarchical detrmination of SimMatrix
% only the methods can be used that produce hierarchical clustering as output
function [Zstar,iter]=allHierarchicalConsensus(P,L,c,hiff,ds,linkagemethod,descriptor,matCombtype,alpha)


iter=0; % the default value. if a combination method returns number of iteration it needed this value will be overwritten
if ~exist('c'),
    c = max(max(cls));
end;
%hierarchical,information theoretic,fuzzy,flat,Genetic
if (hiff==1)        
    %'simpleHierarchical' 
    % the two followig methods are the same but the
    %second one uses the saved clusterers to speed it up
   %[SimMatrix,DiffMatrix]=ComputeHierarchicalSimMatrix(P,L,c,descriptor);
   [SimMatrix,DiffMatrix]=ComputeSimpleHierarchicalSimMatrix(P,L,c,descriptor);
   
elseif (hiff==2)    
    %'Informationtheoretic'
   [SimMatrix,DiffMatrix]=ComputeInformationtheoreticSimMatrix(P,L,c,descriptor,matCombtype,alpha);
elseif (hiff==3)    
    %'FuzzyConsensus' 
   [SimMatrix,DiffMatrix,iter]=FuzzyConsensus(P,L,c,descriptor);
elseif (hiff==4)    
    %'Flat' 
   [SimMatrix,DiffMatrix]=ComputeSimMatrix(P,L,c);  
elseif (hiff==5)    
    %'Genetic' 
   %[SimMatrix,DiffMatrix]=ComputeGeneticSimMatrix(P,L,c,descriptor,matCombtype,alpha);   
	[SimMatrix,DiffMatrix]=ComputeGeneticSimMatrix(P,L,c,X,group,descriptor,matCombtype,alpha); 
 
elseif (hiff==6)    
   %'AggregateOper'
   [SimMatrix,DiffMatrix]=ComputeAggregateOperSimMatrix(P,L,c,descriptor,matCombtype,alpha);     
end

% note that in the case of fuzzy combination the output diffmatrix is
% ultrametric so any hierarchical method can retrive its topology without
% any difference. another point is that it is not rational to interpret it as
% a data

if (ds==1)
    %here SimMatrix(M) is interpreted as data
    Y = pdist(SimMatrix);             %is similar to Y = pdist(X);
else
    %assume DiffMatrix(1-M) is dissimilarity to make DiffMatrix suitable for
    %linkage we must take upper triangular matrix
    Y=squareform(DiffMatrix);         %is similar to Y = pdist(X);
end
Zstar = linkage(Y,linkagemethod);

% %knwon types
% 'single linkage on m(similarity)'
%     (P,L,c,hf=0,ds=0,linkagemethod='single',descriptor='')
% 'single linkage on m(data)'
%     (P,L,c,hf=0,ds=1,linkagemethod='single',descriptor='')
% 'mean linkage on m(similarity)'
%     (P,L,c,hf=0,ds=0,linkagemethod='average',descriptor='')
% 'mean linkage on m(data)'
%     (P,L,c,hf=0,ds=1,linkagemethod='average',descriptor='')
% 'hierarchical single linkage on m(similarity)'
%     (P,L,c,hf=1,ds=0,linkagemethod='single',descriptor='pmd')
% 'hierarchical single linkage on m(data)'
%     (P,L,c,hf=1,ds=1,linkagemethod='single',descriptor='pmd')
% 'hierarchical mean linkage on m(similarity)'
%     (P,L,c,hf=1,ds=0,linkagemethod='single',descriptor='pmd')
% 'hierarchical mean linkage on m(data)'
%     (P,L,c,hf=1,ds=1,linkagemethod='single',descriptor='pmd')
