%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.
%Fig.4

clear all

global  beta gamma maxAc Em0  Eml1 delta 
epsilon=0;
gamma=0;
graphic=1;
delta=0;
beta=1.05 ;
maxAc=8;
Em0=0.1;
Tsim=400;

fig4params
figure
set(gcf,'position', [100 100 520 420]);
subplot('position', [0.05 0.78 0.35 0.1])

plot(actot,'k','LineWidth', 0.2);
set(gca,'XLim', [0 Tsim],'YLim',[0  8]);
set(gca,'XTick',[0: Tsim/2: Tsim])
set(gca,'YTick',[0: 4: 8])
set(gca,'XTickLabel',[])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(15,7,['OS'],'FontWeight','bold','FontSize',9,'FontAngle','italic')
text(490,0,['a'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','normal','FontSize',16,'FontAngle','normal','Color',[0 0 0])
text(60,10,['e = 0'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','symbol','FontWeight','bold','FontSize',10,'FontAngle','italic','Color',[0 0 0])
text(160,10,['g = 0'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','symbol','FontWeight','bold','FontSize',10,'FontAngle','italic','Color',[0 0 0])
text(260,10,['d = 0'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','symbol','FontWeight','bold','FontSize',10,'FontAngle','italic','Color',[0 0 0])
text(160,14,['A'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','normal','FontSize',16,'FontAngle','normal','Color',[0 0 0])
box off
subplot('position', [0.05 0.64 0.35 0.1])
plot(Eml1,'k')
set(gca,'XLim', [0 Tsim],'YLim',[0  8]);
set(gca,'XTick',[0: Tsim/2: Tsim])
set(gca,'XTickLabel',[])
set(gca,'YTick',[0: 4: 8])

set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(15,7,['Em'],'FontWeight','bold','FontSize',9,'FontAngle','italic')
text(55,6,['l1'],'FontWeight','bold','FontSize',8,'FontAngle','italic')
box off
set(gcf,'Color',[1 1 1])
text(450,5,['b=1.05'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','symbol','FontWeight','bold','FontSize',10,'FontAngle','italic','Color',[0 0 0])
text(450,2,['MaxAc=8'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','ariall','FontWeight','normal','FontSize',9,'FontAngle','italic','Color',[0 0 0])


graphic=1; 
beta=1.02 ; 
maxAc=8;
Tsim=400;
fig4params

subplot('position', [0.05 0.50 0.35 0.1])
plot(actot,'k','LineWidth', 0.2);
set(gca,'XLim', [0 Tsim],'YLim',[0  8]);
set(gca,'XTick',[0: Tsim/2: Tsim])
set(gca,'XTickLabel',[])
set(gca,'YTick',[0: 4: 8])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(15,7,['OS'],'FontWeight','bold','FontSize',9,'FontAngle','italic')
box off
text(490,0,['b'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
  'FontName','arial','FontWeight','normal','FontSize',16,'FontAngle','normal','Color',[0 0 0])
subplot('position', [0.05 0.36 0.35 0.1])
plot(Eml1,'k')
set(gca,'XLim', [0 Tsim],'YLim',[0  8]);
set(gca,'XTick',[0: Tsim/2: Tsim])
set(gca,'XTickLabel',[])
set(gca,'YTick',[0: 4: 8])

set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(15,7,['Em'],'FontWeight','bold','FontSize',9,'FontAngle','italic')
text(55,6,['l1'],'FontWeight','bold','FontSize',8,'FontAngle','italic')
text(450,5,['b=1.02'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','symbol','FontWeight','bold','FontSize',10,'FontAngle','italic','Color',[0 0 0])
text(450,2,['MaxAc=8'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','ariall','FontWeight','normal','FontSize',9,'FontAngle','italic','Color',[0 0 0])

        
box off
set(gcf,'Color',[1 1 1])

                
graphic=1;
beta=1.05 ; 
maxAc=4;
Tsim=400;
fig4params


subplot('position', [0.05 0.22 0.35 0.1])
plot(actot,'k','LineWidth', 0.2);
set(gca,'XLim', [0 Tsim],'YLim',[0  8]);
set(gca,'XTick',[0: Tsim/2: Tsim])
set(gca,'XTickLabel',[])
set(gca,'YTick',[0: 4: 8])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(15,7,['OS'],'FontWeight','bold','FontSize',9,'FontAngle','italic')
box off
text(490,0,['c'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
      'FontName','arial','FontWeight','normal','FontSize',16,'FontAngle','normal','Color',[0 0 0])
subplot('position', [0.05 0.08 0.35 0.1])
plot(Eml1,'k')
set(gca,'XLim', [0 Tsim],'YLim',[0  8]);
set(gca,'XTick',[0: Tsim/2: Tsim])
set(gca,'XTickLabel',[0:20:40])
set(gca,'YTick',[0: 4: 8])

set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(15,7,['Em'],'FontWeight','bold','FontSize',9,'FontAngle','italic')
text(55,6,['l1'],'FontWeight','bold','FontSize',8,'FontAngle','italic')
text(250,-4,['time (s)'],'FontWeight','bold','FontSize',9,'FontAngle','normal')
box off
set(gcf,'Color',[1 1 1])
text(450,5,['b=1.05'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','symbol','FontWeight','bold','FontSize',10,'FontAngle','italic','Color',[0 0 0])
text(450,2,['MaxAc=4'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','ariall','FontWeight','normal','FontSize',9,'FontAngle','italic','Color',[0 0 0])



epsilon=0.15;
gamma=0.2;
delta=6;
graphic=1;
beta=1.05 ;
maxAc=8;
Tsim=400;

fig4params
subplot('position', [0.6 0.78 0.35 0.1])
plot(actot,'k','LineWidth', 0.2);
set(gca,'XLim', [0 Tsim],'YLim',[0 8]);
set(gca,'XTick',[0: Tsim/2: Tsim])
set(gca,'XTickLabel',[])
set(gca,'YTick',[0: 4: 8])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(15,7,['OS'],'FontWeight','bold','FontSize',9,'FontAngle','italic')
text(160,14,['B'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','normal','FontSize',16,'FontAngle','normal','Color',[0 0 0])
box off

text(60,10,['e = 0.15'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','symbol','FontWeight','bold','FontSize',10,'FontAngle','italic','Color',[0 0 0])
text(180,10,['g = 0.2'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','symbol','FontWeight','bold','FontSize',10,'FontAngle','italic','Color',[0 0 0])
text(300,10,['d = 6'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','symbol','FontWeight','bold','FontSize',10,'FontAngle','italic','Color',[0 0 0])
subplot('position', [0.6 0.64 0.35 0.1])
plot(Eml1,'k')
set(gca,'XLim', [0 Tsim],'YLim',[0  8]);
set(gca,'XTick',[0: Tsim/2: Tsim])
set(gca,'XTickLabel',[])
set(gca,'YTick',[0: 4: 8])

set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(15,7,['Em'],'FontWeight','bold','FontSize',9,'FontAngle','italic')
text(55,6,['l1'],'FontWeight','bold','FontSize',8,'FontAngle','italic')
box off
set(gcf,'Color',[1 1 1])
 
graphic=1;
beta=1.02 ; 
maxAc=8;

Tsim=400;
fig4params

subplot('position', [0.6 0.50 0.35 0.1])
plot(actot,'k','LineWidth', 0.2);
set(gca,'XLim', [0 Tsim],'YLim',[0  8]);
set(gca,'XTick',[0: Tsim/2: Tsim])
set(gca,'XTickLabel',[])
set(gca,'YTick',[0: 4: 8])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9 )
text(15,7,['OS'],'FontWeight','bold','FontSize',9,'FontAngle','italic')
box off
   
subplot('position', [0.6 0.36 0.35 0.1])
plot(Eml1,'k')
set(gca,'XLim', [0 Tsim],'YLim',[0  8]);
set(gca,'XTick',[0: Tsim/2: Tsim])
set(gca,'XTickLabel',[])
set(gca,'YTick',[0: 4: 8])

set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(15,7,['Em'],'FontWeight','bold','FontSize',9,'FontAngle','italic')
text(55,6,['l1'],'FontWeight','bold','FontSize',8,'FontAngle','italic')


box off
set(gcf,'Color',[1 1 1])


delta=6;
graphic=1;
beta=1.05 ; 
maxAc=4;
Tsim=400;
fig4params


subplot('position', [0.6 0.22 0.35 0.1])
plot(actot,'k','LineWidth', 0.2);
set(gca,'XLim', [0 Tsim],'YLim',[0  8]);
set(gca,'XTick',[0: Tsim/2: Tsim])
set(gca,'YTick',[0: 4: 8])
set(gca,'XTickLabel',[])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(15,7,['OS'],'FontWeight','bold','FontSize',9,'FontAngle','italic')
box off

subplot('position', [0.6 0.08 0.35 0.1])
plot(Eml1,'k')
set(gca,'XLim', [0 Tsim],'YLim',[0  8]);
set(gca,'XTick',[0: Tsim/2: Tsim])
set(gca,'YTick',[0: 4: 8])
set(gca,'XTickLabel',[0:20:40])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(15,7,['Em'],'FontWeight','bold','FontSize',9,'FontAngle','italic')
text(55,6,['l1'],'FontWeight','bold','FontSize',8,'FontAngle','italic')
text(250,-4,['time (s)'],'FontWeight','bold','FontSize',9,'FontAngle','normal')

box off
set(gcf,'Color',[1 1 1],'position',[100 100 520 420]);


