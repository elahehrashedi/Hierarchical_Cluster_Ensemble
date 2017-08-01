function res=CombMethodNMI(pstar,P,L,c)
% s=0;
% for i=1:L
%     %s=s+NMI(pstar,P(1,i).partition,c);
%     s=s+NMI(pstar,P.partition{1,i},c);
% 
% end
% res=s/L;

%The folowing code do the same as above but with error checking
cls = ConvertEnsemble2StrehlFormat(P,L);
res = ceevalmutual(cls,pstar');