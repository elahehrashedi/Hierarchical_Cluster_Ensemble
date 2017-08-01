function [Zstar_out,combacc,terror]=boost_combine(X, L, c , group,dataname)
%function boost_combine(X, L, c , group,dataname)
global sCC pdistx spdistx coeficient fillnum normalize show;

ClustererTypes={'h all linkage'};%,'h all linkage''h avgl''h sl','','h cl','h weighted','h median','h ward'};%,,'nosubsample all linkage'};%
CTypeInFileName={'hcentroid','hsl','havgl','hcl','hweighted','hmedian','hward'};
descriptors={'cd'};%,'cmd','smd','pmd','my'};
linkagemethods={'complete'};%,'single','average','ward','centroid', 'median','centroid'};
matCombtypes={'min'};%,'max','average'};
alphas ={-1,1,2};

% %best parameters 'h all linkage',cd,average,min and function 4
% ClustererTypes={'h centroid'};
% CTypeInFileName={'hcentroidl'};
% descriptors={'cd'};
% linkagemethods={'average'};
% matCombtypes={'min'};%,'jsd','product','rankavg'};
% alphas ={-1,1,2};

sCC=[];
L=1;
niter=30;%100; % number of iteration
notrand=1;%1 weighted selection  0 bagg
func = 3;%  1 : wout not based on win ,2: wout base on win , 3: paper, 4:book
percent=0.4; %percent of sample selection
fillnum=1;
if_refine_base=1; %1:refine z , 0: do not   refine base matrices
if_refine_comb=1; %1:refine z , 0: do not   refine combination
normalize = 2;
bestCC=[];
worstCC=[];
%save_W=[];
%save_CC=[];
w=[];%weights
subsampleidx=[];
combmatrix=[];
created_matrix=[];
 
pdistx=pdist(X);
spdistx=squareform(pdist(X));


    for l=1:length(matCombtypes)
        for k=1:length(linkagemethods)
            for i=1:length(ClustererTypes)
                for j=1:length(descriptors)

                    cnum = 1;
                    ConsensusParams{3}= ClustererTypes{i} ;
                    ConsensusParams{4}= descriptors{j} ;  
                    linkagemethod = linkagemethods(k) ;
                    ConsensusParams{5}= matCombtypes(l) ;

                    if strcmp(char(matCombtypes(l)),'min')                        
                        fillnum=1;
                    elseif strcmp(char(matCombtypes(l)),'max')
                        fillnum=0;
                    elseif strcmp(char(matCombtypes(l)),'average')
                        fillnum=1;
                    end

%                     % first create clustering with all samples    
%                     [subsampleidx,matrix,N] = Create_Clusterer1(X, group, c,'nosubsample average linkage',ConsensusParams{4},w,0,100,fillnum,normalize);
% 
%                     BaseY=squareform(matrix); 
%                     Zstar = linkage(BaseY,linkagemethod);
% 
% %                     if if_refine_base==1
% %                         [ZZ,Y,c1,c2]=refine(Zstar,spdistx);
% %                         Zstar=ZZ;
% %                         ZD=DendrogramDescriptor(ZZ, 'cd');
% %                         ZD=normalfunc (ZD,normalize,N) ;
% %                         matrix=ZD;
% %                     end
%                     baseacc=HierCombMethodAccuracy(Zstar,group,X,'cophenet');
%                     baseterror= Terror(Zstar,group,c);
                    N=size(X,1);                    
                    baseacc=1;
                    
                    
                    %combination dendrogram
                    if strcmp(char(matCombtypes(l)),'min')
                        combmatrix=ones(N,N);
                    else
                        combmatrix= zeros(N,N);
                    end
                    itemselectioncounter=zeros(1,N);
                    coeficient=zeros(N,N);
                    
                    % weight of samples
                    w=ones(1,N)/N;
                    acc(1,:)=0;
                    %%save
                    %save_W(1,:)=ones(1,N)/N;%w_out;
                    %save_CC(1,:)=zeros(1,N)/N;
                    %save_itemselected(1,:)=zeros(1,N);
                    bestCC=zeros(1,N);
                    worstCC=ones(1,N);
                    CCsum=zeros(1,N);
                    e=zeros(1,N);

                    %%%%                                       
                    % change one of inputs
                    for input=1:niter

                        cnum = cnum+1;% in first loop it is 2

                        %Create Clusterer with subsample
                        [subsampleidx,matrix,N] = Create_Clusterer1(X, group, c,ConsensusParams{3},ConsensusParams{4},w,notrand,percent,fillnum,normalize);

                        %count how many times an item selected
                        %%%%saving
                        %save_itemselected(cnum,:)=counts (zeros(1,N),subsampleidx,N);
                        %%%%%%%%%%
                        itemselectioncounter = counts (itemselectioncounter,subsampleidx,N);

                        Y=squareform(matrix);
                        Zstar = linkage(Y,linkagemethod);

                        if if_refine_base==1
                            [ZZ,Y,c1,c2]=refine(Zstar,spdistx);
                            if(c2>c1)
                                Zstar=ZZ;
                                ZD=DendrogramDescriptor(ZZ, 'cd');
                                ZD=normalfunc (ZD,normalize,N) ;
                                matrix=ZD;
                            end
                        end                       
                        acc(cnum,j)=cophenet(Zstar,pdistx);

                        % change weights befor
                        %[w_out,cc,do_combine,bestCC,worstCC,CCsum,e_out]=weight_patterns(w,matrix,X,func,save_CC(cnum-1,:),save_itemselected(cnum,:),bestCC,worstCC,CCsum,cnum);
                        [w_out,cc,do_combine,bestCC,worstCC,CCsum,e_out]=weight_patterns(w,matrix,X,func,zeros(1,N),zeros(1,N),bestCC,worstCC,CCsum,cnum);
                        %%%save
                        %save_W(cnum,:)=w_out;%w_out-w
                        %save_CC(cnum,:)=cc;
                        e(cnum)=e_out;
                        %w=w_out;change place

                        if do_combine
                            % combine dedrograms
                            sume=sum(e)-e_out;
                            combmatrix2=combine2(combmatrix,matrix,ConsensusParams{5},N,cnum,sume,e_out);
                            Y=squareform(combmatrix2); 
                            Zstar = linkage(Y,linkagemethod);

                            if if_refine_comb==1
                                [ZZ,Y,c1,c2]=refine(Zstar,spdistx );
                                if c2>c1
                                    Zstar=ZZ;
                                    ZD=DendrogramDescriptor(ZZ, 'cd');
                                    ZD=normalfunc (ZD,normalize,N) ;
                                    combmatrix2=ZD;
                                end
                            end
                            combacc(cnum,j)=cophenet(Zstar,pdistx);
                            terror(cnum,j)= Terror(Zstar,group,c);
                            if combacc(cnum,j)<combacc(cnum-1,j)
                                combacc(cnum,j)=combacc(cnum-1,j);
                                %combmatrix not changed
                                %w not changed
                                w=w_out;                                
                            else
                                combmatrix=combmatrix2;
                                %w changed
                                w=w_out;
                            end

                        else
                            %combacc(cnum,j)=0;
                            combacc(cnum,j)=combacc(cnum-1,j);
                        end
                        
                        if(cnum==2)%save max Zstar
                            Zstar_out= Zstar;
                        else
%                             if (cophenet(Zstar_out,pdistx)<combacc(cnum,j))
%                                 Zstar_out= Zstar;
%                             end
                            if (Terror(Zstar_out,group,c)<terror(cnum,j))
                                Zstar_out= Zstar;
                            end
                        end
                    end
                    clear R; clear P;                      

 
%                      data=char(dataname);cluster=char(CTypeInFileName{i});desc=char(descriptors{j});mat=char(matCombtypes(l));ll=char(linkagemethods(k));    
% %                     %filename=sprintf('res\\b%s_%s_%s_%s_%s',data,cluster,desc,mat,ll);
%                     filename=sprintf('z%s_%s_%s_%s_%s',data,cluster,desc,mat,ll);
% %                     %save(filename ,'acc','combacc', 'baseacc','meanACC','StdACC','LAST_combacc');
% %                     %save(filename ,'acc','combacc', 'baseacc','save_W','save_CC','save_itemselected','itemselectioncounter');
%                     save(filename ,'acc','combacc', 'baseacc')
% %                     clear combacc; clear acc;
                end
            end
        end
    end        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%w weights
function [w_out,CC,do_combine,bestCC,worstCC,CCsum,loge] = weight_patterns(w_in,created_matrix,X,func,CC_prev,itemselected,bestCC,worstCC,CCsum,cnum)

          global sCC pdistx spdistx show;
            cnum=cnum-1;
          % the base matriix to be compared with
          compared_matrix=spdistx;
          w_out=[];
          %clustering quality
          CC=[];
          CQ=[];
          do_combine=1;
          e = 0 ;
          psudoloss = 0 ;
          
          %for all patternd calculate weight
          N=size(X,1);
          for i=1:N
              %weight   
              %compared_matrix = the ordinary Euclidean distance between the ith and jth observations.    
              %createdmatrix=cophenetic_matrix = the dendrogrammatic distance between the model points Ti and Tj. 
              %This distance is the height of the node at which these two points are first joined together
              
              %mohasebe mojaddad
              Abar =  sum(compared_matrix(i,:))/N;
              Bbar = sum (created_matrix(i,:))/N;
              a=Zigma(compared_matrix(i,:),Abar,created_matrix(i,:),Bbar,N); 
              b=(Zigma(compared_matrix(i,:),Abar,compared_matrix(i,:),Abar,N)*Zigma(created_matrix(i,:),Bbar,created_matrix(i,:),Bbar,N));
              cci=abs(a/b.^0.5);
              

              %cq2=corrcoef(compared_matrix(i,:),created_matrix(i,:));
              %correlation near to one : good clustering So
               
              if b==0   CC(i)=0;
              else      CC(i)=abs(cci);     end              
              %if itemselected(i)==1     sCC=[sCC;CC_prev(i),CC(i)];     end
              %the higher the value of CQti the worse the clustering quality of instance xi.\
          end

           for i=1:N
%               if CC(i)>bestCC(i)
%                 bestCC(i)=CC(i);
%               end
%               if CC(i)<worstCC(i)
%                 worstCC(i)=CC(i);
%               end
              CCsum(i)=CCsum(i)+CC(i);              
          end
          
          for i=1:N
              %CQ(i)=1-(bestCC(i)-worstCC(i));
              %CQ(i)=1-CC(i);
              CQ(i)=1 - (CCsum(i)/cnum); %remember that maxCQ = 1
          end
          
%           for i=1:N
%               %CQ(i)=1-(bestCC(i)-worstCC(i));
%               CQ_prev(i)=1-CC_prev(i);
%               %CQ_prev(i)=1- ((CQ(i)*cnum)-CC(i))*(cnum-1);
%           end

        %==================================================================
        for i=1:N
            e = e + w_in(i)*(1-CC(i));
        end
                
        % reverse weights so that patterns taht have higher correlation get lower weight
        %correlation between -1 to 1
        %==================================================================
         if func==1
              %me1
              for i=1:N
                  % better ones get lower weight
                  if (CQ(i)>0)
                    %w_out(i)= 1 - CQ(i) ;
                    w_out(i)= CQ(i);
                  else
                    %w_out(i)= 1 + CQ(i) ;
                    display('negative weight');
                    w_out(i)= CQ(i);
                  end
              end  
        %==================================================================              
          elseif func == 2          
              %me 2
              for i=1:N
                  % better ones get lower weight ; higher value of CQti the worse the clustering quality of instance xi
                  if (CQ(i)>0)
                    %w_out(i)= w_in(i) - CQ(i) ;
                    w_out(i)= w_in(i) * CQ(i) ;%where CQ=1 -mean(CCsum);
                  else
                    w_out(i)= w_in(i) * CQ(i) ;
                  end
              end 
        %==================================================================              
        %paper
         elseif func == 3
            for i=1:N
                %loss for selecteds
                %psudoloss = psudoloss + w_in(i)* itemselected(i) * (CC(i)-CC_prev(i)) ;
                %psudoloss = psudoloss +(CC(i)-CC_prev(i)) ;
                %psudoloss = psudoloss +(   1 - w_in(i)*(CC(i)-bestCC(i))   );
                psudoloss = psudoloss + w_in(i)* CQ(i); % where CQ=1 -mean(CCsum);
            end
            psudoloss = psudoloss/2;
            Beta = (1-psudoloss)/psudoloss;
          
            if Beta<0.5
                 %do_combine=0;
                 display('hight value of loss in 3 ');
                 w_out=ones(1,N)/N;
            else 
              for i=1:N
                w_out(i)= (w_in(i)* Beta.^CQ(i)) ;
              end
            w_out= w_out / sum(w_out) ;
            end
        %==================================================================            
         elseif func== 4              
             if(e==0 || e>=0.9)
                w_out=ones(1,N)/N;
                do_combine=0;
                display('hight value of loss in 4 '); 
             else
                 beta=e/(1-e);
                 for i=1:N
                    w_out(i)= (w_in(i)* beta.^CC(i)) ;                
                 end
                w_out= w_out / sum(w_out) ;               
             end         
         end
          %==========================================================================
          % Zed is a normalization constant such that W is a distribution
          w_out= w_out / sum(w_out) ;
          %loge=abs(log((1-e)/e));%true
          loge=1-e;
          %loge=Beta; 
          %show=[show;1-e,1/Beta];
            
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function combmatrix = combine2(combmatrix,correcteddendro,matCombtype,N,cnum,sum_e,last_e)

global coeficient fillnum normalize;

        type= char(matCombtype);
        alpha=2;
        if strcmp(char(matCombtype),'waverage')
            for i=2:N
                for j=1:i-1
                    if(correcteddendro(i,j) ~= fillnum)
                        m1=coeficient(i,j);
                        m2=m1+last_e;
                        combmatrix(i,j) = (correcteddendro(i,j)*last_e + combmatrix(i,j)*m1)/m2;
                        combmatrix(j,i)=combmatrix(i,j);
                        coeficient(i,j)=m2;
                        coeficient(j,i)=m2;
                    end
                end
            end
        elseif strcmp(char(matCombtype),'average')
        	for i=2:N
                for j=1:i-1
                    %if(correcteddendro(i,j) ~= fillnum)
                        m1=sum_e;
                        m2=sum_e+last_e;
                        combmatrix(i,j) = (correcteddendro(i,j)*last_e + combmatrix(i,j)*m1)/m2;
                        combmatrix(j,i)=combmatrix(i,j);
                    %end
                end
            end                                    
        else      
            for j=1:N
                DT=[];
                if strcmp(char(matCombtype),'average1')
                    if cnum==2 % first matrix treaed as NAN
                        combmatrix(j,:) = correcteddendro(j,:);
                    else
                        combmatrix(j,:) = (correcteddendro(j,:) + combmatrix(j,:)*(cnum-2))/(cnum-1);
                    end
                else
                    DT=[DT; correcteddendro(j,:)];
                    DT=[DT; combmatrix(j,:)];
                    combmatrix(j,:) = comb( DT ,type,alpha);
                end
            end
        end
        %combmatrix=normalfunc (combmatrix,normalize,N);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z=Zigma(X,Xbar,Y,Ybar,N)
    z=0;
    for i=1:N
        z=z+(X(i)-Xbar)*(Y(i)-Ybar);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function itemselectioncounter = counts (itemselectioncounter,subsampleidx,N)
    for i=1:size(subsampleidx,2)
       itemselectioncounter(subsampleidx(i)) =  itemselectioncounter(subsampleidx(i)) +1;        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function terror = Terror(Z,group,c)

       T = cluster(Z,'maxclust',c);
       T=PermuteLabel(group,T,c,'book heuristic','a');
       terror=CombMethodAccuracy(T,c,group);
       
end


