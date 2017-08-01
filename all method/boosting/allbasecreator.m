function allbasecreator()
ClustererTypes={'nosubsample single linkage','nosubsample average linkage',...
    'nosubsample complete linkage','nosubsample weighted linkage','nosubsample centroid linkage',...
    'nosubsample median linkage','nosubsample ward linkage'};
% Data_Names={'contractions','thyroid','laryngeal2','heart',...
% 'breast_cancer','laryngeal3','liver_disorders',...
% 'weaning','diabetis',...
% 'pima_indians_diabetes','wpbc',...
% 'splice',...
% 'ringnorm'};
Data_Names={'respiratory','thyroid_aeberhard'};
linkagemethods={'weighted','centroid','average','single','complete','ward','median'};%, 
w=[];

for xx=1:length(Data_Names)
    [X, group, c]= LoadBenchmarkdata (Data_Names{xx},1,1);
    for yy=1: length(linkagemethods)
        Y = pdist(X, 'euclid');
        Z = linkage(Y, linkagemethods{yy});
        b=cophenet(Z,Y);
        eval ( sprintf('baseacc_%s(yy)=b;',Data_Names{xx}));
    end
end

  display('salam');      
       
end