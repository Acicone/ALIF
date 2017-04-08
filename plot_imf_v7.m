function[h]=plot_imf_v7(imf1,imf2,T,legText1,legText2)


[m,n]=size(imf1);
if nargin<3, T=1:n; end

if isempty(T), T=1:n; end

TextSize=20;

h=figure;
for i=1:m
    subplot(m,1,i); 
    I2=plot(T,imf2(i,:),'r','LineWidth',2);
     hold on
    I1=plot(T,imf1(i,:),'k','LineWidth',2);
    set(gca,'fontsize', TextSize);         
end

% if nargin >= 3  && not(isempty(titLeg))
%     set(gcf,'NextPlot','add');
%     axes;
%     tt= title(titLeg,'FontSize',TextSize);
%     set(gca,'Visible','off');
%     set(tt,'Visible','on');
% end

if nargin >= 4 && not(isempty(legText1)) && not(isempty(legText2))
hh=legend([I1,I2],{legText1,legText2},'Location','Best');
set(hh,'Interpreter','latex')
end

end