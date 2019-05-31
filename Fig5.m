%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.
%Fig5

clear all
global beta gamma maxAc delta actot
global Tsim epsilon 

epsilon=0;
gamma=0;
delta=0;
beta=1.05 ;
maxAc=8;
Tsim=12000;
scaledelta=(2:2:10);
scalebeta=(1.01:0.02:1.1);
durbeta=zeros(size(scalebeta'));
interbeta=zeros(size(scalebeta'));
toutbeta=zeros(size(scalebeta'));
durdelta=zeros(size(scaledelta'));
interdelta=zeros(size(scaledelta'));
toutdelta=zeros(size(scaledelta'));
compteurbeta=0;
compteurdelta=0;
for beta=scalebeta
    compteurbeta=compteurbeta+1;
    fig5params
    y=interp(actot,10);
    [a,b]=detectionblsimmcp(y,100,2,30);
    toutbeta(compteurbeta)=beta;
    durbeta(compteurbeta)=a;
    interbeta(compteurbeta)=b;
end
beta=1.05;

for delta=scaledelta
    compteurdelta=compteurdelta+1;
    fig5params
    y=interp(actot,10);
    [a,b]=detectionblsimmcp(y,100,2,30);
    toutdelta(compteurdelta)=delta;
    durdelta(compteurdelta)=a;
    interdelta(compteurdelta)=b;
end


% Create figure
figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],...
    'Color',[1 1 1]);

% Create textbox
annotation(figure1,'textbox','String',{'A'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.015 0.97 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'B'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.495 0.97 0.004 0.001]);

%figure(1)
subplot(1,2,2)
plot(toutbeta, 60./interbeta,'ko', toutbeta,durbeta,'k*')
hold on
subplot(1,2,1)
hold on
subplot2 = subplot(1,2,1);
pl1=plot(toutdelta, 60./interdelta, 'ko');
pl2=plot(toutdelta, durdelta, 'k*');
legend1 = legend(subplot2,'show');
set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
   'Position',[0.14 0.45 0.39 0.14],'FontSize',8,...
    'FontName','Arial');
set(pl1,'Marker','o','DisplayName','L episode frequency [1/min]');
set(pl2,'Marker','*','DisplayName','L episode duration [s]' );

hold on
    
Xbeta=[ones(length(toutbeta),1) toutbeta];
tetabeta=Xbeta\(60./interbeta);
Xaffbeta=[ones(2,1) toutbeta([1,end],1)];
yaffbeta=Xaffbeta*tetabeta;
subplot(1,2,2)
plot(toutbeta([1,end],1), yaffbeta,'k:')
text(1.025,17,[' MaxAc = 8 '],...  
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','normal','FontSize',10,'FontAngle','italic')
text(1.05,-1.6,[' b '],...  
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','symbol','FontWeight','normal','FontSize',12,'FontAngle','italic')
set(gca,'Xlim',[1 1.1])
set(gca,'Ylim',[0 20])
box off

Xdelta=[ones(length(toutdelta),1) toutdelta];
tetadelta=Xdelta\(durdelta);
Xaffdelta=[ones(2,1) toutdelta([1,end],1)];
yaffdelta=Xaffdelta*tetadelta;
subplot(1,2,1)
plot(toutdelta([1,end],1), yaffdelta,'k:')
set(gca,'Xlim',[0 10])
set(gca,'Ylim',[0 10])
text(3.5,8.5,[' b = 1.05 '],...  
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','symbol','FontWeight','normal','FontSize',10,'FontAngle','italic')
text(5,-0.8,[' MaxAc  '],...  
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','normal','FontSize',12,'FontAngle','italic')

box off