[X, group, c]= LoadBenchmarkdata ('thyroid_aeberhard',1,1);
[subsampleidx,matrix,N] = Create_Clusterer1(X, group, c,'nosubsample centroid linkage','cd',[],0,100,1,2);
                    BaseY=squareform(matrix); 
                    Zstar = linkage(BaseY,'average');
                            pdistx=pdist(X);
                    c1=cophenet(Zstar,pdistx);
                    Zstar2=linkage(BaseY,'weighted');
                    c2=cophenet(Zstar,pdistx);
                    c1
                    c2
                    c1