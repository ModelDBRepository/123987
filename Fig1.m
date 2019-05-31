%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.
%Fig.1

clear all
%global vpos1 vpos2 vpos3 vpos4
%global Tsim 

figure(1)
set(gcf,'position',[100 100 520 420])
if exist('chain5loop3N&R.jpg')
    subplot('position',[0.03 0.78 0.45 0.16]);
    fig1=imread('chain5loop3N&R.jpg');
    imagesc(fig1);
    axis off;
    text(380,-30,['A'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','normal','FontSize',16,'FontAngle','normal','Color',[0 0 0]) 
    text(1245,-30,['B'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','normal','FontSize',16,'FontAngle','normal','Color',[0 0 0])   
    text(-30,30,['a'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','normal','FontSize',16,'FontAngle','normal','Color',[0 0 0]) 
end    
    
vpos1=[0.07 0.4 0.4 0.28];
vpos2=[0.07 0.09 0.4 0.25];
Tsim=40;
Chain5loop3N

Izhinit
vpos3=[0.55 0.4 0.4 0.28];
vpos4=[0.55 0.09 0.4 0.25];
ChainloopIzhi
text(-42720,3,['b'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','normal','FontSize',16,'FontAngle','normal','Color',[0 0 0])
text(-42720,0.2,['c'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','normal','FontSize',16,'FontAngle','normal','Color',[0 0 0]) 
set(gcf,'Color',[1 1 1],'position',[100 100 520 420]);