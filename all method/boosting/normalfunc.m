function matrix = normalfunc (matrix,normalize,N)                    
    if normalize==1
        matrix=matrix./repmat(sum(matrix,2),1,N);
    end
    if normalize==2
        matrix=matrix/max(matrix(:));
    end
    if normalize==3
        m=mean(matrix(:));
        v=var(matrix(:));
        v=v.^(1/2);
        correcteddendro=(matrix-m)/v;
    end
    if normalize==4
    end
end