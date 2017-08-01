function [ Pens ] = comb( DT,matCombtype,alpha)
% In the experiments using renyi alpha=0 lines 28 and 30 are inserted now
% and in experiments on all data sets these lines 
% notice that alpha in the renyi is not realy alpha in the derived formula
% i.e. 1-alpha is assumed alpha
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here
if nargin < 3
    matCombtype='renyi';
    alpha=20;
end

[L,c]=size(DT);
switch lower(matCombtype)
    case 'max'
        Pens=max(DT,[],1);
    case 'min'
        Pens=min(DT,[],1);
    case 'average'
        Pens=nanmean(DT,1);
    case 'product'
        Pens=prod(DT);
    case 'rankavg' %Average of Rank     
        [a,I]=sort(DT,2); % replace the original DT with Rank of them 
        DT=I;
        Pens=mean(DT,1);
    case 'renyi'
        if (alpha==0) % renyi with alpha==0 is equivalent to product
            DT=DT*c;  % makes the values greater so their product not be zero(due to round of error)
            Pens=prod(DT);
            Pens=Pens/sum(Pens);
        else
            DTalpha=DT.^alpha;
            s=nansum(DTalpha,1);
            Pens=s.^(1/alpha);
        end
    case 'jsd'
        for m=1:c
            p=[1,sum(DT(:,m))/(1-2^L)];
            for i=2:L
                idx=creatediffpermut(L,i);
                pp=0;
                for j=1:size(idx,1)
                    pp=pp+prod(DT(idx(j,:),m));
                end
                p=[p,pp/(1-2^L)];
            end;
            r = roots(p);
            id=testroot(r,DT(:,m),L);
            if (id==-1) error('*****no real roots between 0 and 1 found*****');end
            Pens(1,m)=r(id);
        end
        % test the roots
        %         for tt=1:c
        %             prod(Pens(1,tt)+DT(:,tt))-Pens(1,tt)^L*2^L;
        %         end
    otherwise
        error('unknown comb type');
end
%Pens

function id=testroot(r,X,L)
% this function test that the found root is a real root of our eq. and also
% is a real number between 0 and 1
id=-1;
for ii=1:length(r)
    if (prod(r(ii,1)+X)-r(ii,1)^L*2^L<0.001 && (0 <= r(ii,1) & r(ii,1) <= 1)&& imag(r(ii,1))==0) 
        id=ii; 
        if (r(ii,1)~=0) return; end
    end;
end