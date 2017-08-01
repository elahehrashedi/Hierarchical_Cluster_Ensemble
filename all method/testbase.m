%evaluate the base clusteres
clc; clear all;
% ela: all are the same
L=10;
niter=1; % number of iteration

Data_Names=  {'breast_cancer','flare_solar','titanic','laryngeal1','laryngeal2','contractions','wpbc','vehicle','wine'};%{};%, };
%Data_Names= {'diabetis','german','banana','splice','vehicle'}
%         'diabetis','german','balance_scale'...
%            'vehicle','laryngeal2' };
%Data_Names= {'contractions'};
%Data_Names= {'heart_switzerland','glass','page_blocks','soybean_large','voice_9'};

for data_idx=1:length(Data_Names)
    data_name=Data_Names{data_idx};
    k=0;
    for j=1:niter     % number of iteration
        filename=sprintf('clusterers\\%s_%d.mat',data_name,j);
        load(filename, 'P','L','c','X','group');
        % ela comment
        CP=complete_dend(P,L);
        %L
        for i=1:7
            k=k+1;
            %ela comment
            %acc(data_idx,k)=HierCombMethodAccuracy(CP.partition{1,i}.dend,group,X,'cophenet');
            acc(data_idx,k)=HierCombMethodAccuracy(P.partition{1,i}.dend,group,X,'cophenet');
            
        end
   
            % SL with all descriptors
            %p.pmd1,1 p.cmd1.1
%             acc(data_idx,1)=HierCombMethodAccuracy(P.pmd{1,1}.dend,group,X,'cophenet');
%             acc(data_idx,2)=HierCombMethodAccuracy(P.cmd{1,1}.dend,group,X,'cophenet');
%             acc(data_idx,3)=HierCombMethodAccuracy(P.smd{1,1}.dend,group,X,'cophenet');
%             acc(data_idx,4)=HierCombMethodAccuracy(P.cd{1,1}.dend,group,X,'cophenet');
%             acc(data_idx,5)=HierCombMethodAccuracy(P.my{1,1}.dend,group,X,'cophenet');
            
%             acc(data_idx,1)=HierCombMethodAccuracy(P.partition{1,1}.dend,group,X,'cophenet');


    end
end

m=mean(acc,2)';
st=std(acc,0,2)';
mn=min(acc,[],2)';
mx=max(acc,[],2)';
filename='base_hierarchical_large_no_sub';
save(filename);
