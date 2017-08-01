%from bi wi and ti create sorted table with descending ti
[i,j]=sort(tiACC,'descend')
%Method	b	w	t
colmethod=ConsensusTypes(j);
colb=biACC(j);
colw=wiACC(j)
colt=i;

[i,j]=sort(tiNMI,'descend')
%Method	b	w	t
colmethod=ConsensusTypes(j);
colb=biNMI(j);
colw=wiNMI(j)
colt=i;