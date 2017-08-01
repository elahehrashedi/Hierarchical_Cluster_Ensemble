% function testdatasets()
clear all;
clc;
Data_Names= {'fourty','long1','sizes5','spiral','square1', 'square4', 'twenty',...      % disp('MOCK datasets')
      'banana','breast_cancer','diabetis','flare_solar','german','heart',...     %disp('Benchmarks datasets')
      'ringnorm','splice','thyroid','titanic','twonorm','waveform',...
      'heart_cleveland','heart_hungarian','heart_va','heart_switzerland','ionosphere',...%disp('UCI ML datasets')
      'thyroid_aeberhard','wpbc','waveform21','balance_scale','glass','image_segmentation',...
      'iris','liver_disorders','page_blocks','pima_indians_diabetes','soybean_large',...
      'soybean_small','vehicle','waveform_noise','wine',...
      'laryngeal1','laryngeal2','laryngeal3','respiratory','voice_3','voice_9','weaning','contractions'}; %disp('real data kuncheva') 
for data_idx=1:length(Data_Names)
    data_name=Data_Names{data_idx};
    [X, group, maxg]=LoadBenchmarkdata(data_name,1,0);
    aa{data_idx,1}=sprintf('%s',data_name);
    bb{data_idx,1}=sprintf('%d',size(X,1));
    cc{data_idx,1}=sprintf('%d',size(X,2));
    dd{data_idx,1}=sprintf('%d',maxg);
    
end


