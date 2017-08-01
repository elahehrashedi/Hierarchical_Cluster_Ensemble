function acc=CombMethodAccuracy(pstar,c,group)
%This function evaluate how much the clustring is similar to ground truth 
%information. group is the class label of every point in benchmark, pstar 
%is the label extracted from ensemble
maxg=max(group);
target=zeros(size(group));
for i=1:c                 % for every cluster find the coresponding label
    inx=find(pstar==i);
    temp=zeros(1,maxg);   % count every class in this cluster the cluster label will be the one with the highest count
    for j=1:maxg
        temp(1,j)= size(find(group(inx,1)==j),1);
    end
    [cc,NewTargat]=max(temp,[],2);
    target(inx,1)=NewTargat*ones(size(inx));
end
err=size(find(group~=target),1)/length(group);
acc=1-err;