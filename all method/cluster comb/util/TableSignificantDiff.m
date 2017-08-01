function  [b,w,s,bi,wi,ti]=TableSignificantDiff( AR,SAR )
%This function compares different method to see if one is significantly better than the other. 
% Significant difference is observed if the respective confidence interval have no intersection. 
%I think that the squre of 100 in formula is the number of experiments.
%AR is the performance evaluation matrix in which each row corresponds to a data set and each
%coloumn is one method 
% AR(i,j) and SAR(i,j) are the mean and standard deviation of the measurement in 100 diferent run
% the number of datasets where x has been significantly better than denoted b(x,y)
% the number of datasets where x has been significantly worse than y denoted w(x,y)
% the number of datasets where the diff is not significant denoted s(x,y)
% an index of total performance denoted t
[r,c]=size(AR);
s=zeros(c,c);
w=s;
b=s;
for x=1:c           % one method
    for y=1:c       % another method
        if (x~=y)
            for k=1:r               % the number of datasets
                if (AR(k,x)-0.196*SAR(k,x)>AR(k,y)-0.196*SAR(k,y))
                    b(x,y)=b(x,y)+1;
                elseif (AR(k,x)-0.196*SAR(k,x)<AR(k,y)-0.196*SAR(k,y))
                    w(x,y)=w(x,y)+1;
%                 else
%                     ss(x,y)=ss(x,y)+1;
                end
            end
        end
    end
end
s=r-b-w;
ti=sum((b-w),2);
bi=sum((b),2);
wi=sum((w),2);