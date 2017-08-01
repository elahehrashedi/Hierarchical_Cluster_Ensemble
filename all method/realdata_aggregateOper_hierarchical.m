%function realdata_hierarchical()

%#############
warning off
%#############

clc; clear all;
L=10;
niter=10; % number of iteration
ClustererTypes={'h sl','h avgl','h cl','h weighted','h centroid','h median','h ward','h all linkage'};
ClustererType=ClustererTypes{1};


Data_Names= {'breast_cancer','flare_solar','titanic','ionosphere',...
     'wpbc','image_segmentation','liver_disorders','wine',...
     'laryngeal1','laryngeal3','weaning','contractions'};
%         'diabetis','german','balance_scale'...
%            'vehicle','laryngeal2' };


CombMethod={'simpleHierarchical','Informationtheoretic','FuzzyConsensus','flat','Genetic','AggregateOper'};
linkagemethods={'single'};
descriptors={'cd','pmd','cmd','smd','my'};

% %==========================================================================
% %'AggregateOper'
methodno=0;
for descriptor=1:5
    for agg_type=1:18
        switch agg_type
            case {1,8}
                % (agg_type==1) //dombi s-norm
                % (agg_type==8) //dombi t-norm
                param=0:.1:1;
                param=[param, 2, 4,8];
                for pp=1:length(param)
                    methodno=methodno+1;
                    ConsensusParams{methodno,1}= 6;                                   %hiff
                    ConsensusParams{methodno,2}= 0;                                   %ds
                    ConsensusParams{methodno,3}= 'single';                            %linkagemethod
                    ConsensusParams{methodno,4}= descriptors{descriptor};             %descriptor
                    ConsensusParams{methodno,5}= agg_type;                            %matCombtype
                    ConsensusParams{methodno,6}= param(1,pp);                         %alpha
                end
            case {2,9}
                % (agg_type==2) //dubois prade s-norm
                % (agg_type==9) //dubois prade t-norm
                param=0:.1:1;
                for pp=1:length(param)
                    methodno=methodno+1;
                    ConsensusParams{methodno,1}= 6;                                   %hiff
                    ConsensusParams{methodno,2}= 0;                                   %ds
                    ConsensusParams{methodno,3}= 'single';                            %linkagemethod
                    ConsensusParams{methodno,4}= descriptors{descriptor};             %descriptor
                    ConsensusParams{methodno,5}= agg_type;                            %matCombtype
                    ConsensusParams{methodno,6}= param(1,pp);                         %alpha
                end
            case{3,10}
                % (agg_type==3) //Yager s-norm
                % (agg_type==10) //Yager t-norm
                param=1:.1:1.9;
                param=[param,2:10];
                for pp=1:length(param)
                    methodno=methodno+1;
                    ConsensusParams{methodno,1}= 6;                                   %hiff
                    ConsensusParams{methodno,2}= 0;                                   %ds
                    ConsensusParams{methodno,3}= 'single';                            %linkagemethod
                    ConsensusParams{methodno,4}= descriptors{descriptor};             %descriptor
                    ConsensusParams{methodno,5}= agg_type;                            %matCombtype
                    ConsensusParams{methodno,6}= param(1,pp);                         %alpha
                end
            case {4,5,6,7,11,12,13,14}
                % (agg_type==4) //drasticSum s-norm
                % (agg_type==5) //EinsteinSum s-norm
                % (agg_type==6) //algebraicSum s-norm
                % (agg_type==7) //max s-norm

                % (agg_type==11) //drasticProd t-norm
                % (agg_type==12) //EinsteinSum s-norm
                % (agg_type==13) //algebraicSum s-norm
                % (agg_type==14) //min t-norm
                % //================================averaging operator======================
                methodno=methodno+1;
                ConsensusParams{methodno,1}= 6;                                   %hiff
                ConsensusParams{methodno,2}= 0;                                   %ds
                ConsensusParams{methodno,3}= 'single';                            %linkagemethod
                ConsensusParams{methodno,4}= descriptors{descriptor};             %descriptor
                ConsensusParams{methodno,5}= agg_type;                            %matCombtype
                ConsensusParams{methodno,6}= 0;                                   %alpha
                
            case {15,17,18}
                % (agg_type==15) //averaging max-min operator
                % (agg_type==17) //averaging FuzzyAnd operator
                % (agg_type==18) //averaging FuzzyOr operator
                param=0:0.1:1;
                for pp=1:length(param)
                    methodno=methodno+1;
                    ConsensusParams{methodno,1}= 6;                                   %hiff
                    ConsensusParams{methodno,2}= 0;                                   %ds
                    ConsensusParams{methodno,3}= 'single';                            %linkagemethod
                    ConsensusParams{methodno,4}= descriptors{descriptor};             %descriptor
                    ConsensusParams{methodno,5}= agg_type;                            %matCombtype
                    ConsensusParams{methodno,6}= param(1,pp);                         %alpha
                end
               
            case {16}
                % (agg_type==16) //averaging Generalized mean operator
                param=-10:1:-1;
                param=[param, 1:1:10;];
                for pp=1:length(param)
                    methodno=methodno+1;
                    ConsensusParams{methodno,1}= 6;                                   %hiff
                    ConsensusParams{methodno,2}= 0;                                   %ds
                    ConsensusParams{methodno,3}= 'single';                            %linkagemethod
                    ConsensusParams{methodno,4}= descriptors{descriptor};             %descriptor
                    ConsensusParams{methodno,5}= agg_type;                            %matCombtype
                    ConsensusParams{methodno,6}= param(1,pp);                         %alpha
                end

        end
    end
end
% %========================================================================
%CreateandSaveClusterer2(L,niter,ClustererType,descriptors,Data_Names);
for data_idx=1:length(Data_Names)
    data_name=Data_Names{data_idx};
    for j=1:niter     % number of iteration
        filename=sprintf('E:\\fuzzy\\clusterers\\%s_%d.mat',data_name,j);
        load(filename, 'P','L','c','X','group');

        for i = 1:size(ConsensusParams,1)
            varargin=ConsensusParams(i,:);
            Zstar=allHierarchicalConsensus(P,L,c,varargin{:});

   
            disp([data_name '(' num2str(j) ')' ' clusterensemble: ' num2str(i) ]);

            Zstarfilename=sprintf('E:\\fuzzy\\aggregateOper result\\zstar\\%s_con%d_iter%d.mat',data_name,i,j);
            save(Zstarfilename, 'Zstar');
        end;
        disp(['++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++']);
    end;
    
end
