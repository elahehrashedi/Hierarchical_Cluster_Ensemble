function [X, group, maxg]=LoadBenchmarkdata(data_name,PreProcess,cut)
% find the dataset containing the data_name
% normalize it between [0,1] 
if nargin < 2
    PreProcess=1;   %preprocess the data to [0 1] interval
    cut=1;
end

switch lower(data_name)
    case {'circle','longspiral','smile','triangle','helal','longsquare','spiralsquare'}
        disp('hand created datasets')
        load hand_created;
        xpression  =data_name;
        X          =eval(xpression);
        target     =X(:,3);
        X=X(:,1:2);

    case {'fourty','long1','sizes5','spiral','square1', 'square4', 'twenty'}
       disp('MOCK datasets')
       load mock;
       xpression  =data_name;
       X          =eval(xpression);
       target     =X(:,3);
       X=X(:,1:2);

    case {'banana','breast_cancer','diabetis','flare_solar','german','heart',...
            'ringnorm','splice','thyroid','titanic','twonorm','waveform'}
        disp('Benchmarks datasets')
        load benchmarks;
        namestr    =data_name;
        xpression  =[namestr,'.x'];
        X          =eval(xpression);

        xpression  = [namestr,'.t'];
        target     =eval(xpression);

    case {'heart_cleveland','heart_hungarian','heart_va','heart_switzerland','ionosphere',...
            'thyroid_aeberhard','wpbc','waveform21','balance_scale','glass','image_segmentation',...
            'iris','liver_disorders','page_blocks','pima_indians_diabetes','soybean_large',...
            'soybean_small','vehicle','waveform_noise','wine'}
        disp('UCI ML datasets')
        load UCI_ML;
        xpression  =[data_name,'.x'];
        X          =eval(xpression);

        xpression  = [data_name,'.t'];
        target     = eval(xpression);

    case {'laryngeal1','laryngeal2','laryngeal3','respiratory','voice_3','voice_9','weaning','contractions'}
        disp('real data kuncheva')
        load realdata_kuncheva;
        xpression  =[data_name,'.x'];
        X          =eval(xpression);

        xpression  = [data_name,'.t'];
        target     =eval(xpression);
    
    otherwise
        disp('Unknown method.')
end

%data set must be of moderate size to be used in our ensembles
if ((size(X,1)>1000)&&(cut==1))
    X=X(1:1000,:);
    target=target(1:1000,:);
end

[group, gnames] = grp2idx(target);
maxg = max(group);

% fill missing value with the average of that column
if length(find(isnan(X)))>0
    idx=find(isnan(X));
    [i,j]=find(isnan(X));
    m= nanmean(X,1);
    for k=1:length(idx)
        X(idx(k))=m(1,j(k));        
    end
end

if (PreProcess==1)
    % means of zero and standard deviations of 1.
    p=X';
    [pn,meanp,stdp] = prestd(p);
    X=pn';
    %gscatter(X(:,1),X(:,2),group);
    %each col of X must be one feature
end

