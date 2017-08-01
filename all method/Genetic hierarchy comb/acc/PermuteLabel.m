function newp=PermuteLabel(pstar,p,c,heuristic,inout)
%permute the labels of p to find the labeling that maximise the similarity between pstar and p
%c is the number of cluster
% inout determines the input output type a=abstract m=measurement
%creating confusion matrix
switch lower(inout)
    case 'a'
        cnf=zeros(c,c);
        A=createallcells(c,2);
        for i=1:size(A,1)
            cnf(A(i,1),A(i,2))=size(find(pstar ==A(i,1)& p==A(i,2)),1);
        end
        cnf;
    case 'm'     % the 2 clusterings are in measurement level
        cnf=pstar'*p;
    otherwise
        disp('Unknown inout type.');
end

switch lower(heuristic)
    case 'book heuristic'
        %no preprocessing
        for j=1:c
            max1=max(cnf(:));
            [ii,jj]=find(cnf==max1);
            ii=ii(1,1);                %if we have two identical digit consider the first one
            jj=jj(1,1);
            %newp(find(p==jj),:)=ii;
            relables(1,ii)=jj;
            cnf(:,jj)=-1*ones(1,c);
            cnf(ii,:)=-1*ones(c,1);
        end

    case 'my heuristic'
        k=cnf;   
        %ktemp=repmat(max(k,[],1),c,1).*repmat(max(k,[],2),1,c);
        ktemp=repmat(sum(k,1),c,1).*repmat(sum(k,2),1,c);
        zerocount=length(find(ktemp==0));
        if (zerocount>0) 
                zeroind=find(ktemp==0);
                ktemp(zeroind)=1;
                save('error_ktempzero');
           %     error('some elements in ktemp are zero');
                warning('some elements in ktemp are zero');
        end
        ktemp=k./ktemp;
        cnf=ktemp; 
        for j=1:c
            max1=max(cnf(:));
            [ii,jj]=find(cnf==max1);
            ii=ii(1,1);                %if we have two identical digit consider the first one
            jj=jj(1,1);
            %newp(find(p==jj),:)=ii;
            relables(1,ii)=jj;
            cnf(:,jj)=-1*ones(1,c);
            cnf(ii,:)=-1*ones(c,1);
        end
    case 'opt'
        % This version differs with PermuteLabel in that it find the optimum relabling
        %permute the labels of p to find the labeling that maximise the similarity between pstar and p
        %c is the number of cluster
        A=createallpermutation(c);
        relables=[];
        bestres=0;
        for ii=1:length(A)
            res=0;
            for cc=1:c
                res=res+cnf(cc,A(ii,cc));
            end
            if res>bestres
                relables=A(ii,:);
                bestres=res;
            end
        end
%         for ii=1:c
%             jj=relables(1,ii);
%             newp(find(p==jj),:)=ii;
%         end
    otherwise
        error('Unknown type');
end

%constructing output 
switch lower(inout)
    case 'a'
        newp=zeros(size(p));
        for ii=1:c
            jj=relables(1,ii);
            newp(find(p==jj),:)=ii;
        end
    case 'm'     % the 2 clusterings are in measurement level
        newmup=zeros(size(p,1),c);
        for ii=1:c
            jj=relables(1,ii);
            newmup(:,ii)=p(:,jj);
        end
        newp=newmup;
    otherwise
        disp('Unknown inout type.');
end
