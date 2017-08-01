function res=NMI(p1,p2,c)
%this function calculate the normalized mutual information between two cluster
%c is the number of cluster
% 
N=size(p1,1);           %total number of training samples
%new method
if (size(p1,2)==1) %if the partition are in abstract level
    mp1=MeasurementLevel(p1,c);
else
    mp1=p1;
end
if (size(p2,2)==1) %if the partition are in abstract level
    mp2=MeasurementLevel(p2,c);
else
    mp2=p2;
end
cnf=mp1'*mp2;
%cnf=[1 1 0 ;0 0 2 ;1 1 0 ;0 2 0 ]; for test only (this work correcly)

[crow,ccol]=size(cnf);
Ni=sum(cnf,2); %sum of all columns for row i
Nj=sum(cnf,1); %sum of all rows for column j
Ntemp=repmat(Ni,1,ccol).*repmat(Nj,crow,1);
%we try to compute res=-2*sum(aa)/(s1+s2)
aa=cnf.*log((cnf.*N)./Ntemp);
aa(find(isnan(aa)))=0;
bb=-2*sum(sum(aa));
s1=sum(Ni.*log(Ni/N));
s2=sum(Nj.*log(Nj/N));
res=bb/(s1+s2);