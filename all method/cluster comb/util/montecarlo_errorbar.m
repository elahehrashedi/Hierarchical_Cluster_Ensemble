
plottools;
hold on;
 
color={'blue','green','red','cyan','magenta','yellow','black'};
%XScale, YScale, ZScale       {linear} | log 
%Axis scaling. Linear or logarithmic scaling for the respective axis. See also loglog, semilogx, and semilogy.
axes1 = axes(...
  'XMinorTick','on',...
  'XScale','linear');
xlim(axes1,[1 7]);
hold(axes1,'all');
for i=1:7
errorbar (1:7,meanNMI(:,i),StdNMI(:,i), 'DisplayName', 'i','LineWidth',2,'Color',color{mod(i,7)}); figure(gcf) 
end