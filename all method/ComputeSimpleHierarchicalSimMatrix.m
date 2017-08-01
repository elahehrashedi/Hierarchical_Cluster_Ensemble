function [SimMatrix,DiffMatrix]=ComputeSimpleHierarchicalSimMatrix(P,L,c,descriptor)
%This function is the same as ComputeHierarchicalSimMatrix but it uses the
%saved clusterers to speed it up.
if nargin < 4
    descriptor='pmd';
end
switch lower(P.type)
    case 'hierarchical'
        N=P.N;
        DiffMatrix=zeros(N,N);
        for i=1:1:L
            %DiffMatrix=DiffMatrix+correcteddendro;
            eval(sprintf('DiffMatrix=DiffMatrix+P.%s{1,i};',descriptor));
        end
        DiffMatrix=DiffMatrix/max(DiffMatrix(:));
        SimMatrix=1-DiffMatrix;
        % to check and make diagonal equal to 0 and 1
        SimMatrix=checks(SimMatrix);
        DiffMatrix=1-SimMatrix;
    case 'abstract'
        error('iput clustering must be hierarchical');
    case 'measurement'
        error('iput clustering must be hierarchical');
    case 'strehl'
        error('iput clustering must be hierarchical');
    otherwise
end
