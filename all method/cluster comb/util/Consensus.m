%function pstar=Consensus(P,L,Type,c)
% this function calculate one partitioning from the clustering ensemble
% using the method specified in kuncheva paper entitled "Experimental
% comparison of cluster ensemble methods"
% Description
% L=100;    %the number of learner in ensemble 
% P         %the clustering ensemble P.partition{1,i} is the ith partition
% c         %the number of cluster
% Type      %the consensus function used

function pstar=Consensus(P,L,Type,c)

if (~exist('c')||isequal(Type,'hgpa'))
    cls = ConvertEnsemble2StrehlFormat(P,L);
end
if ~exist('c'),
    c = max(max(cls));
end;

switch lower(Type)
    case 'single linkage on m(similarity)'
        %assume DiffMatrix(1-M) is dissimilarity to make DiffMatrix suitable for
        %linkage we must take upper triangular matrix
        disp('single linkage on M(similarity)');
        [SimMatrix,DiffMatrix]=ComputeSimMatrix(P,L,c);
        Y=squareform(DiffMatrix);           %is similar to Y = pdist(X);
        Z = linkage(Y,'single');
        pstar = cluster(Z,'maxclust',c);

    case 'single linkage on m(data)'
        %here SimMatrix(M) is interpreted as data
        disp('single linkage on M(data)');
        [SimMatrix,DiffMatrix]=ComputeSimMatrix(P,L,c);
        Y = pdist(SimMatrix);             %is similar to Y = pdist(X);
        Z = linkage(Y,'single');
        pstar = cluster(Z,'maxclust',c);
        
    case 'mean linkage on m(data)'
        %here SimMatrix(M) is interpreted as data
        disp('mean linkage on M(data)');
        [SimMatrix,DiffMatrix]=ComputeSimMatrix(P,L,c);
        Y = pdist(SimMatrix);             %is similar to Y = pdist(X);
        Z = linkage(Y,'average');
        pstar = cluster(Z,'maxclust',c);

    case 'k-means on m(data)'
        %here SimMatrix(M) is interpreted as data
        %in the consensus step we can perform several kmeans which one was
        %beter is returned
        disp('k-means on M(data)');
        [SimMatrix,DiffMatrix]=ComputeSimMatrix(P,L,c);
        pstar= kmeans(SimMatrix,c,'emptyaction','singleton','replicates',5);

    case 'cspa on m(similarity)'
        %here SimMatrix(M) is interpreted as similarity graph
        disp('cspa on m(similarity)');
        [SimMatrix,DiffMatrix]=ComputeSimMatrix(P,L,c);
        pstar = metis(SimMatrix,c);

    case 'hgpa'
        disp('HGPA')
        pstar = hgpa(cls,c);
        
    case 'direct optimization(book)'
        pstar=DirectOptimization(P,L,c,'book heuristic');
    case 'direct optimization(my)'
        pstar=DirectOptimization(P,L,c,'my heuristic');
    case 'direct optimization(opt)'
        pstar=DirectOptimization(P,L,c,'opt');
%---------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------
    case 'direct simple voting(book)'
        pstar=Directvoting(P,L,c,'book heuristic','simple voting');
    case 'direct simple voting(my)'
        pstar=Directvoting(P,L,c,'my heuristic','simple voting');
    case 'direct simple voting(opt)'
        pstar=Directvoting(P,L,c,'opt','simple voting');

    case 'direct wighted voting(book)'
        pstar=Directvoting(P,L,c,'book heuristic','wighted voting');
    case 'direct wighted voting(my)'
        pstar=Directvoting(P,L,c,'my heuristic','wighted voting');
    case 'direct wighted voting(opt)'
        pstar=Directvoting(P,L,c,'opt','wighted voting');

    case 'direct selected voting(book)'
        pstar=Directvoting(P,L,c,'book heuristic','selected voting');
    case 'direct selected voting(my)'
        pstar=Directvoting(P,L,c,'my heuristic','selected voting');
    case 'direct selected voting(opt)'
        pstar=Directvoting(P,L,c,'opt','selected voting');

    case 'direct selected weighted voting(book)'
        pstar=Directvoting(P,L,c,'book heuristic','selected weighted voting');
    case 'direct selected weighted voting(my)'
        pstar=Directvoting(P,L,c,'my heuristic','selected weighted voting');
    case 'direct selected weighted voting(opt)'
        pstar=Directvoting(P,L,c,'opt','selected weighted voting');
%---------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------
    %###################Divide and conquere    
    %book and my are local ones using ghosh heuristic to find relabling
    %opt finds the optimal relabing 
    case 'dc incremental_2d(book)'
        pstar=DivisiveCombination(P,L,c,'book heuristic','Incremental_2d');
    case 'dc incremental_2d(my)'
        pstar=DivisiveCombination(P,L,c,'my heuristic','Incremental_2d');
    case 'dc incremental_2d(opt)'
        pstar=DivisiveCombination(P,L,c,'opt','Incremental_2d');

    case 'dc incremental_3d(book)'
        pstar=DivisiveCombination(P,L,c,'book heuristic','Incremental_3d');
    case 'dc incremental_3d(my)'
        pstar=DivisiveCombination(P,L,c,'my heuristic','Incremental_3d');
    case 'dc incremental_3d(opt)'
        pstar=DivisiveCombination(P,L,c,'opt','Incremental_3d');

    case 'dc hierarchical_2d(book)'
        pstar=DivisiveCombination(P,L,c,'book heuristic','hierarchical_2d');
    case 'dc hierarchical_2d(my)'
        pstar=DivisiveCombination(P,L,c,'my heuristic','hierarchical_2d');
    case 'dc hierarchical_2d(opt)'
        pstar=DivisiveCombination(P,L,c,'opt','hierarchical_2d');

    case 'dc hierarchical_3d(book)'
        pstar=DivisiveCombination(P,L,c,'book heuristic','hierarchical_3d');
    case 'dc hierarchical_3d(my)'
        pstar=DivisiveCombination(P,L,c,'my heuristic','hierarchical_3d');
    case 'dc hierarchical_3d(opt)'
        pstar=DivisiveCombination(P,L,c,'opt','hierarchical_3d');
    
    case 'dc incremental_2d_2pass(book)'
        pstar=DivisiveCombination(P,L,c,'book heuristic','incremental_2d_2pass');
    case 'dc incremental_2d_2pass(my)'
        pstar=DivisiveCombination(P,L,c,'my heuristic','incremental_2d_2pass');
    case 'dc incremental_2d_2pass(opt)'
        pstar=DivisiveCombination(P,L,c,'opt','incremental_2d_2pass');
    
        %###################Divide and conquere
%---------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------
    case 'direct optimization 3d(book)'
        pstar=DirectOptimization_3d(P,L,c,'book heuristic');
    case 'direct optimization 3d(my)'
        pstar=DirectOptimization_3d(P,L,c,'my heuristic');
    case 'direct optimization 3d(opt)'
        pstar=DirectOptimization_3d(P,L,c,'opt');
    
    otherwise
        error('Unknown type');
end

