function [Zstar] = Boost()

clc; 
clear all;
L=1;
dorefine=1;

%breast_cancer 4;contractions 1;flare_solar 0.1;titanic 0.0;laryngeal1 0.1,laryngeal2 1,
%'laryngeal3' 1; wine 1;'heart_hungarian' 1;
%Data_Names= {'breast_cancer','contractions','flare_solar','titanic','laryngeal1','laryngeal2','laryngeal3','titanic','wine'};
%Data_Names= {'diabetis','splice','image_segmentation','iris','voice_3',...
% 'heart_cleveland','german','heart','thyroid','heart_switzerland','ionosphere','wpbc','waveform21','balance_scale',...
% 'glass','page_blocks','pima_indians_diabetes','soybean_large',...
% 'vehicle','waveform_noise','liver_disorders','weaning','soybean_small',...
% 'ringnorm','twonorm','waveform','twenty','fourty',...
% 'banana'};
% Data_Names={'splice'};%2900
% Data_Names={'ringnorm'};%7400
%Data_Names={'titanic'};



%Data_Names={'contractions','thyroid','laryngeal2','heart'};%98+215+691
% Data_Names={'breast_cancer','laryngeal3','liver_disorders'};%266+353+345
% Data_Names={'weaning','diabetis'};%302+768
Data_Names={'wpbc'};%768+198
% Data_Names={'respiratory','thyroid_aeberhard'};
%Data_Names={'respiratory','thyroid_aeberhard'};


    for xx=1:length(Data_Names)
        cnum = 0 ; % number of clusters which are combined
        data_name=Data_Names{xx};
        [X, group, c]= LoadBenchmarkdata (data_name,1,1);
        pdistx=pdist(X);
        spdistx=squareform(pdist(X));
        
        %[Zstar,combacc,terror]=boost_combine(X, L, c,group,Data_Names{xx});
        boost_combine(X, L, c,group,Data_Names{xx});
%         if dorefine==1
%             [ZZ,Y,c1,c2]=refine(Zstar,spdistx);
%             Zstar=ZZ;
%         end
%         selacc=cophenet(Zstar,pdist(X));
%         eval(sprintf('combacc_%s=combacc;',data_name));
%         eval(sprintf('terror_%s=terror;',data_name));
%         eval(sprintf('selcophenetic_%s=selacc;',data_name)); 
%         filename=sprintf('%s.mat',data_name);
%         save(filename ,'selacc','combacc', 'terror')
    end

    
end
