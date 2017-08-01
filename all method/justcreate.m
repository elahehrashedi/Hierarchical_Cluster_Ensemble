function justcreate( )
clc; clear all;
L=50;
niter=1;%0; % number of iteration

%ClustererTypes={'nosubsample single linkage'};
ClustererTypes={'h sl'};
%ClustererTypes={'h sl','h avgl','h cl','h weighted','h centroid','h median','h ward','h all linkage','nosubsample all linkage'};

%Data_Names=  {'breast_cancer','flare_solar','titanic','laryngeal1','laryngeal2','contractions','wpbc','vehicle','wine'};%{};%, };
Data_Names=  {'titanic'};
%Data_Names= {'wpbc'};%,'','laryngeal2','','wpbc','vehicle'};
%     Data_Names= {'breast_cancer','flare_solar','titanic','ionosphere',...
%      'wpbc','image_segmentation','liver_disorders','wine',...
%      'laryngeal1','laryngeal3','weaning','contractions',...
%         'diabetis','german','balance_scale'...
%            'vehicle','laryngeal2' }

%CombMethod={'genetic'};

%linkagemethods={'single','complete','average','weighted','ward','centroid', 'median'};
descriptors={'cd','pmd','cmd','smd','my'};
%matCombtypes={'product','jsd','max','min','average','rankavg'}; % renye is inserted by hand

%alphas ={-32 -24	-16	-8	-4	-2	-1	0 1	2	4	8	16	24	32};

% methodno=0;
% 
%  for ds=0:1
%     for linkagemethod=1:5
%         for descriptor=1:5
%             for alpha=1:15
%                 methodno=methodno+1;
%                 ConsensusParams{methodno,1}= 2;                                   %hiff
%                 ConsensusParams{methodno,2}= ds;                                  %ds
%                 ConsensusParams{methodno,3}= linkagemethods{linkagemethod};       %linkagemethod
%                 ConsensusParams{methodno,4}= descriptors{descriptor};             %descriptor
%                 ConsensusParams{methodno,5}= 'renyi';                             %matCombtype
%                 ConsensusParams{methodno,6}= alphas{alpha};                       %alpha
%             end
%         end
%     end
% end


%for numoftypes=1:length(ClustererTypes)
%    CreateandSaveClusterer(L,niter,ClustererTypes{numoftypes},descriptors,Data_Names);
%     disp(numoftypes);
% end

CreateandSaveClusterer(L,niter,ClustererTypes,descriptors,Data_Names);

end