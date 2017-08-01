clc
clear all;close all;
% Data_Names= {'long1','sizes5','spiral','square1', 'square4', ...
%     'circle','longspiral','smile','triangle','helal','longsquare','spiralsquare'};

Data_Names= {'breast_cancer'};

   linkagemethods={'single','average','weighted','complete','ward'};
   L=10;
   for data_idx=1:length(Data_Names)
       data_name=Data_Names{data_idx};
       [X, group, maxg]=LoadBenchmarkdata(data_name,1,1);
       lx=length(X);
       if lx>300
           stp=round(lx/300);
           idx=1:stp:lx;
           X=X(idx,:);
           group=group(idx,:);
       end
       c=maxg;
%        Y = pdist(X,'euclid');
%        Z = linkage(Y,'single');
%        T = cluster(Z,'maxclust',c);

       T=hEnsembleCluster( X,L,c);

       T=PermuteLabel(group,T,c,'book heuristic','a');
       acc(data_idx)=CombMethodAccuracy(T,c,group);
       
       figure;

%       subplot(1,2,1),
%       gscatter(X(:,1),X(:,2),group,'','+o*.xsd^v<>ph');
%       title('original');
%       daspect([1 1 1])
%       subplot(1,2,2),
        gscatter(X(:,1),X(:,2),T,'','+o*.xsd^v<>ph');
%        title('clustered');
        daspect([1 1 1])
       
   end
   mean(acc)