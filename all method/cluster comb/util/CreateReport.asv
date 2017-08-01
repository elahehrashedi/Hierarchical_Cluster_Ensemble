% function a=CreateReport(meanf,Stdf)
% This function create a table with the format meanf+Stdf
meanf=meanACC;
Stdf=StdACC;
if (size(meanf)~=size(Stdf))
    eror('the size of mean and std must be the same');
end
for i=1:size(meanf,1)
    for j=1:size(meanf,2)
        aa{i,j}=sprintf('%4.4f+%4.4f',meanf(i,j),Stdf(i,j));
    end
end
 bb=char(aa(:,1)); %columnn 1
 bb=char(aa(:,2)); %columnn 2
 