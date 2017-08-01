%this file is desined to compare the result of fuzzy clustering with other
%clustering combination method namely the EA and simple combination method.
%thte result is reported in fuzzy journal.
%#############
warning off
%#############

clc; clear all;
L=50;
niter=1; % number of iterationXXXXXXXXXXXXXXXXX
ClustererTypes={'h sl','h avgl','h cl','h weighted','h centroid','h median','h ward','h all linkage'};
ClustererType=ClustererTypes{1};


% Data_Names= {'breast_cancer','flare_solar','titanic','diabetis','german','banana','splice',...
%     'ionosphere','wpbc','image_segmentation','liver_disorders','wine','vehicle',...
%     'heart_switzerland','glass','page_blocks','soybean_large','voice_9',... 
%     'laryngeal1','laryngeal2','laryngeal3','weaning','contractions'  };
Data_Names= {'titanic'};

CombMethod={'simpleHierarchical','Informationtheoretic','FuzzyConsensus','flat'};
linkagemethods={'single','average'};
descriptors={'cd','pmd','cmd','smd','my'};

% %==========================================================================
% %'SimpleConsensus'
methodno=0;
    methodno=methodno+1;
    ConsensusParams{methodno,1}= 1;                                   %hiff
    ConsensusParams{methodno,2}= 0;                                   %ds
    ConsensusParams{methodno,3}= 'single';                            %linkagemethod
    ConsensusParams{methodno,4}= 'cd';                                %descriptor
    ConsensusParams{methodno,5}= '';                                  %matCombtype
    ConsensusParams{methodno,6}= 0;                                   %alpha
% %'FlatConsensus'
    methodno=methodno+1;
    ConsensusParams{methodno,1}= 4;                                   %hiff
    ConsensusParams{methodno,2}= 0;                                   %ds
    ConsensusParams{methodno,3}= 'single';                            %linkagemethod
    ConsensusParams{methodno,4}= 'cd';                                %descriptor
    ConsensusParams{methodno,5}= '';                                  %matCombtype
    ConsensusParams{methodno,6}= 0;                                   %alpha
    
% %========================================================================
%CreateandSaveClusterer2(L,niter,ClustererType,descriptors,Data_Names);
for data_idx=1:length(Data_Names)
    data_name=Data_Names{data_idx};
    for j=1:niter     % number of iteration
        %'F:\\MATLAB\\R2008a\\work\\clusterers\\%s_%d.mat'
        filename=sprintf('F:\\eli\\cluster comb\\Clusterers\\%s_%d.mat',data_name,j);
        load(filename, 'P','L','c','X','group');

        for i = 1:size(ConsensusParams,1)
            varargin=ConsensusParams(i,:);
             [Zstar,iter]=allHierarchicalConsensus(P,L,c,varargin{:});
             acc(i,data_idx)=HierCombMethodAccuracy(Zstar,group,X,'cophenet');
            %saving star to mesure their suitability later
            %F:\\MATLAB\\R2008a\\work\\fuzzy\\journal large data\\other\\
            Zstarfilename=sprintf('res\\other_%s_con%d_iter%d.mat',data_name,i,j);
            save(Zstarfilename, 'Zstar','iter');
        end;
        disp(['++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++']);
    end;
end
%F:\\MATLAB\\R2008a\\work\\fuzzy\\journal large data\\other
filename='res\\other_hierarchical2';
save(filename);