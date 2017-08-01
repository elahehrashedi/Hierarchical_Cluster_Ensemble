function [Z,Y,c1,c2]=refine(Z,Y)
% generates a cmd (cluster membership divergence) matrix of the hierarchical, binary cluster tree 
% represented by Z. The clusters(nodes) are created and the number of item in each of them is
% reperesented as diffrence for the newly combined members. because the cluster size is repeatedly
% increased if we compute the min between the former cmdmatrix and the new calculated one the items
% that are previously joined are not changed.

c1=cophenet(Z,squareform(Y));
if nargin < 1
    %     Z=[1,2,2;
    %         4,5,3;
    %         8,3,4;
    %         6,7,5;
    %         10,9,7;
    %         12,11,10;]
    %     type='cd';
    %     dendrogram(Z);
        clc; clear all;
        %X = [rand(300,3); rand(200,3)+2; rand(200,3)+3];
        load matlab;
        Y = pdist(X);
    %   X= [4,4;
    %     8 4;
    %     15 8;
    %     24 4;
    %     24 12];
    %     Y = pdist(X);

    %     load fisheriris
    %     Y = pdist(meas);
        Z = linkage(Y,'single');
        %figure;dendrogram(Z);
        c = cophenet(Z,Y);
        Y = squareform(Y);
end

m = size(Z,1)+1;
% cmdmatrix=ones(m,m)*inf;
for i=1:1:m
    cluster{i}.P=[i];     
    cmdmatrix(i,i)=1;     %every item has a cluster only one itm
end
for k=1:1:m-1
    i=Z(k,1); j=Z(k,2);
    P=[cluster{i}.P,cluster{j}.P];
    cluster{m+k}.P=P;
    len=size(P,2);       %for every node find the cluster size
    %=============================================================
%     S=0;
%     for i=1:len
%         S=S+sum(Y(P(1,i),P)); %set the dissimilarity according to min cluster size having jk together
%     end
%     Z(k,3)=S/(2*len);
    %===========================================================
    
    %===============================================
    S=0;
    for ii=1:length(cluster{i}.P)
        for jj=1:length(cluster{j}.P)
            S=S+Y(cluster{i}.P(1,ii),cluster{j}.P(1,jj));
        end
    end
    Z(k,3)=S/(length(cluster{i}.P)*length(cluster{j}.P));
    if k>1
        if Z(k,3)<Z(k-1,3)
            Z(k,3)=Z(k-1,3)+0.00000001*Z(k-1,3);
            %The above formula guarantees the structure of the tree but may
            %not be the optimal one
            %Z(k-1,3)=Z(k,3)-0.00000001*Z(k,3);
        end
    end
    %=================================================
end
Y = squareform(Y);
c2 = cophenet(Z,Y);

%figure;dendrogram(Z);

end

