%function EnsembleComparison()
clc; clear all;
L=10;
niter=10; % number of iteration
ClustererTypes={'k-means','k-means+subsample','single linkage+subsample',...
    'average linkage+subsample','complete linkage+subsample','weighted linkage+subsample',...
    'centroid linkage+subsample','median linkage+subsample','ward linkage+subsample',...
    'all linkage+subsample','fcm','random'};

ConsensusTypes={ 'single linkage on m(similarity)','single linkage on m(data)',...
    'mean linkage on m(data)','k-means on m(data)','cspa on m(similarity)','hgpa',...
    'direct optimization(book)','direct optimization(my)',...
    'direct simple voting(book)','direct simple voting(my)','direct wighted voting(book)',...
    'direct wighted voting(my)','direct selected voting(book)','direct selected voting(my)',...
    'direct selected weighted voting(book)','direct selected weighted voting(my)',...
    'direct optimization 3d(book)','direct optimization 3d(my)' };

ClustererType=ClustererTypes{2};

for data_idx=1:7
    %data_idx=7;
    %'banana'	'breast_cancer'	'diabetis'	'flare_solar'	'german'	'heart'	'image'	'ringnorm'
    %'splice'	'thyroid'	'titanic'	'twonorm'	'waveform'
    %[X, group, maxg, namestr]=benchmarkdata(data_idx,1);

    %'fourty','long1','sizes5','spiral','square1', 'square4', 'twenty'
    
    [X, group, maxg, namestr]=benchmarkdata_mock(data_idx,1);
    c=maxg;
    for j=1:niter     % number of iteration
        %P=CreateClusterer(X,L,ClustererType); if the cluster number is not
        %defined the cnf computation of
        %direct optimization will be
        %incorrect
        P=CreateClusterer(X,L,ClustererType,c);
        for i = 7:length(ConsensusTypes);
            ConsensusType=ConsensusTypes{i};
            pstar=Consensus(P,L,ConsensusType,c);
            err(i,j)=CombMethodAccuracy(pstar,c,group);
            disp(['clusterensemble: ' ConsensusType ' at ' num2str(err(i,j))]);
        end;
        disp(['++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++']);
    end;
    err=err(7:12,:);
    [qual, best] = min(err)
    sum(err,2)/niter
    filename=sprintf('%s_L%d',namestr,L);
    save(filename);
end
