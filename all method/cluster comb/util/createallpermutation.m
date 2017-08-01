function A=createallpermutation(Dim)
%function that defines permutation for a Dim dimensional problem and number them
%Dim number of feature
%k   number of mf for each feature
% for  dim4 
A=[];
if (Dim==7)
    for i=1:Dim
        for j=1:Dim
            for k=1:Dim
                for l=1:Dim
                    for m=1:Dim
                        for n=1:Dim
                            for o=1:Dim
                                if ((i~=j)&&(i~=k)&&(i~=l)&&(i~=m)&&(i~=n)&&(i~=o)&&(j~=k)&&(j~=l)&&...
                                        (j~=m)&&(j~=n)&&(j~=o)&&(k~=l)&&(k~=m)&&(k~=n)&&(k~=o)...
                                        &&(l~=m)&&(l~=n)&&(l~=o)&&(m~=n)&&(m~=o)&&(n~=o))
                                    A=[A;i,j,k,l,m,n,o];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if (Dim==6)
    for i=1:Dim
        for j=1:Dim
            for k=1:Dim
                for l=1:Dim
                    for m=1:Dim
                        for n=1:Dim
                            if ((i~=j)&&(i~=k)&&(i~=l)&&(i~=m)&&(i~=n)&&(j~=k)&&(j~=l)&&...
                                    (j~=m)&&(j~=n)&&(k~=l)&&(k~=m)&&(k~=n)&&(l~=m)&&(l~=n)&&(m~=n))
                                A=[A;i,j,k,l,m,n];
                            end
                        end
                    end
                end
            end
        end
    end
end
if (Dim==5)
    for i=1:Dim
        for j=1:Dim
            for k=1:Dim
                for l=1:Dim
                    for m=1:Dim
                        if ((i~=j)&&(i~=k)&&(i~=l)&&(i~=m)&&(j~=k)&&(j~=l)&&(j~=m)&&(k~=l)&&(k~=m)&&(l~=m))
                            A=[A;i,j,k,l,m];
                        end
                    end
                end
            end
        end
    end
end
if (Dim==4)
    for i=1:Dim
        for j=1:Dim
            for k=1:Dim
                for l=1:Dim
                    if ((i~=j)&&(i~=k)&&(i~=l)&&(j~=k)&&(j~=l)&&(k~=l))
                        A=[A;i,j,k,l];
                    end
                end
            end
        end
    end
end
if (Dim==3)
    for i=1:Dim
        for j=1:Dim
            for k=1:Dim
                if ((i~=j)&&(i~=k)&&(j~=k))
                    A=[A;i,j,k];
                end
            end
        end
    end
end
if (Dim==2)
    for i=1:Dim
        for j=1:Dim
                if (i~=j)
                    A=[A;i,j];
                end          
        end
    end
end