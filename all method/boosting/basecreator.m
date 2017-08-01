function [Zstar,X] = basecreator()
%yasamin612001@yahoo.com
clc; 
clear all;


ClustererTypes={'nosubsample single linkage','nosubsample average linkage',...
    'nosubsample complete linkage','nosubsample weighted linkage','nosubsample centroid linkage',...
    'nosubsample median linkage','nosubsample ward linkage'};
Data_Names={'wine','contractions','titanic','laryngeal2','vehicle',...
    'laryngeal1','laryngeal3','breast_cancer','flare_solar'};
descriptors={'cd'};%,'cmd','smd','pmd','my'};
linkagemethods={'average','single','complete','ward','centroid', 'median'};

for xx=1:length(Data_Names)
    [X, group, c]= LoadBenchmarkdata (Data_Names{xx},1,1);
    w=[];%weights
        for k=1:length(linkagemethods)
            for i=1:length(ClustererTypes)        
                for j=1:length(descriptors)

                    cnum = 1;
                    linkagemethod = linkagemethods(k) ;

                    [P(cnum),subsampleidx,created_matrix] = B_CreateClusterer1(X, group, c,ClustererTypes{i},descriptors{j},w,0,100,0);
                    N=P.N;
                    w=ones(1,N)/N;
                    % get its accuracy
                    %eval(sprintf('R{cnum}.correcteddendro=P.%s{1,i};',descriptor));
                    eval(sprintf('R{cnum}.correcteddendro=P(cnum).%s{1,1};',descriptors{j}));
                    %difference matrix
                    BaseY=squareform(R{cnum}.correcteddendro); 
                    Zstar = linkage(BaseY,linkagemethod);
                    
                    %prev accuracy 3 places changes
                    baseacc=HierCombMethodAccuracy(Zstar,group,X,'cophenet');
                    baseacc%=cophenet(Zstar,Y);
                    
                end
                end
            end
        end
   
end
   