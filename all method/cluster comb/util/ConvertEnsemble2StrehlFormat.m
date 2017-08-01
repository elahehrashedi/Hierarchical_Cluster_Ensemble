% function cls = ConvertEnsemble2StrehlFormat(P,L)
%
% DESCRIPTION 
% CreateClusterEnsemble and CreateClusterer have a structure syntax but
% StrehlFormat uses cls = [cl1; cl2;...]note that StrehlFormat can not
% handle measurement level clustering such as Fuzzy c-means clustering
% Description
% L the number of clusterer in ensemble
% P the structure containing clusterer result
function cls = ConvertEnsemble2StrehlFormat(P,L)
if strcmp(P.type,'measurement') % if the input is in measurementlevel convert it to abstractlevel
%if (size(P(1,1).partition,2)~=1)  
    disp('Strehl Format can work with abstract level only');
    disp('converting to abstract...');
    for i=1:1:L
        %P(1,i).partition=AbstractLevel(P(1,i).partition,c);
        P.partition{1,i}=AbstractLevel(P.partition{1,i},c);
        
    end
end
cls=[];
for i=1:L
    %newp=P(1,i).partition;
    newp=P.partition{1,i};
    cls=[cls;newp'];
end
%P.type='Strehl';
