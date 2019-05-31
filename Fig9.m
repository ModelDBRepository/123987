%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.
%Fig9

clear all
global beta gamma maxAc Em0 cumul Ac Eml1 delta

P =99.8778;                                             %Tonic interspike time
erreurP =0.0833;
Delex =1.6000;                                          %Transmission delay
errdelex =5.6843e-014;
dt=0.125;
Idep=4.8;
nNeuron = 29 ;
nCurrent = nNeuron;
tEnd =15000;                                            %maximum simulation time
time = dt:dt:tEnd ;
pars=[0.02      0.2     -65     8       14 ];           %1-tonic spiking
nType=1*ones(1,nNeuron);
a=pars(nType,1)';
b=pars(nType,2)';
c=pars(nType,3)';
d=pars(nType,4)';
Is=0*ones(nNeuron,1)';
SynRatex = 1*ones(1,nNeuron);% 3 msec
SynRatin = 0.026*ones(1,nNeuron);
SynDecayex = 1 - SynRatex*dt;
SynDecayin = 1 - SynRatin*dt;
Gmaxsynex=50;
Gmaxsynin=30;
groupebc=[5 11 13 18 24 26];                    %Bc and Bd groups defined in Fig7.m
groupebd=[7 9 15 20 22 28];
seedmodul=1;
rand('seed',seedmodul);
modul=1;        
epsilon=0;              
beta=1.00013;
gamma=0.002;              
maxAc=8;
delta=0;
Em0=0.1;
Ac=0;
Eml1=0;
    
tindex = 1;                 %pointer for output arrays
Istate = 0;                 %state variable used for synaptic currents decay
vmod=zeros(1,nNeuron);
vmod(1)=1;
I = Is ;
v = c;                          % reset voltage;
s = zeros(1,nNeuron);           % spike occurrence
u = zeros(1,nNeuron);           % recovery variable


maxI=50;
minI=-100;
%W is the connection matrix used in Fig7
w = ...
[0          -1           0           0           0          -1           0          -1           0          -1           0          -1           0          -1           0          -1           0           0          -1           0          -1           0          -1           0          -0           0          -0        0.15          -1
0          -1           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0       -0.15           0           0           0           0           0           0           0           0           0           0
0           0           0           0           0           0           0           0           0       -0.15           0       -0.15           0           0           0       -0.15           0           0           0           0           0           0           0           0           0           0           0           0           0
1           0           0           0           0          -1           0           0           0           0           0       -0.15        0.15           0           0           0           0           0           0           0           0           0           0           0           0           0       -0.15           0           0
1           0           0           1           0          -1           0          -0           0       -0.15        0.15           0           0           0           0           0           0           0          -1           0           0           0           0           0           0           0           0           0           0
1           0           0           0           1           0           0           0           0           0           0           0           0           0           0           0           0           0           0        0.15           0           0       -0.15           0           0           0           0           0           0
1           0           0        0.15           1           0        0.15          -1           0          -1           0           0           0           0           0           0           0           0           0           0          -1        0.15           0           0           0           0           0           0           0
0           0           0           0           0           0           1           0           0           0           0           0           0       -0.15           0           0           0           0           0           0           0           0           0        0.15           0           0           0           0           0
1           0           0           0           0           0           1           0           0          -1        0.15          -1           0           0           0           0           0           0           0           0           0           0          -1        0.15           0           0       -0.15           0           0
0       -0.15           0           0           0           0           0           0           1       -0.15           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0
1           0           0           0           0           0           0           0           1           0           0          -1           0          -1           0           0           0           0           0           0           0           0           0           0          -1           0           0           0           0
1           0           0           0           0           0           0           0           0           0           1       -0.15           0           0           0           0           0           0           0           0           0           0       -0.15           0           0           0           0           0           0
1           0           0           0       0.15           0           0           0           0           0           1           0           0          -1           0          -1           0           0           0           0           0           0           0           0       -0.15           0          -1           0           0
1           0           0           0           0           0           0           0        0.15           0           0           0           1           0           0           0           0           0           0           0           0           0           0           0           0        0.15       -0.15           0           0
1       -0.15           0           0           0           0           0           0           0           0        0.15           0           1           0           0          -1           0           0           0        0.15           0           0           0           0           0           0           0           0          -1
0       -0.15           0           0        0.15           0           0           0           0           0           0           0           0           0           1           0           0           0           0           0           0           0           0           0           0           0           0        0.15           0
1           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0          -1           0           0           0           0           0           0           0       -0.15           0           0
1       -0.15           0           0           0          -1           0           0           0           0           0           0           0           0           0           0           0           0          -0           0          -1           0           0           0           0           0           0           0           0
1           0           0           0           0           0           0           0           0       -0.15           0           0           0           0           0           0           0           1           0           0           0        0.15           0           0           0           0           0           0           0
1           0           0           0           0           0           0          -1           0           0           0           0           0           0           0           0           0           1           0           0          -1           0          -1           0           0           0           0           0           0
0           0           0           0           0           0           0           0           0           0           0           0        0.15           0           0           0           0        0.15           0           1           0           0           0           0           0           0           0           0           0
1           0           0           0           0           0           0           0           0          -1           0           0           0           0           0           0           0           0           0           1           0           0          -1           0          -1           0           0        0.15           0
0           0           0           0           0           0        0.15           0           0           0           0           0        0.15           0           0           0           0           0           0           0           0           1           0           0           0           0           0           0           0
1           0           0           0           0           0           0           0           0           0        0.15          -1           0       -0.15           0           0           0           0           0           0           0           1           0        0.15          -1           0          -1           0           0
1           0           0           0           0           0           0           0           0           0           0       -0.15           0           0           0           0           0        0.15           0           0           0           0           0           1           0           0           0           0           0
1           0           0           0           0           0           0           0           0           0           0           0           0          -1           0           0           0           0           0           0       -0.15           0           0           1           0           0          -1           0          -0
1           0           0           0           0           0           0           0           0           0           0           0        0.15           0           0           0           0           0           0           0       -0.15           0           0           0           0           0           0           0           0
1           0           0           0           0           0           0           0           0           0           0           0           0           0           0          -1           0           0           0           0           0           0           0           0           0           1           0           0          -1
0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0       -0.15           0           0           1           0] ;
w((w==0.15))=0.007;
w((w==-0.15))=-0.01;
SynStrength =Gmaxsynex* abs(w)'.*(w'>0) - Gmaxsynin* abs(w)'.*(w'<0); 
groupex=find(sum(w)>=0);
groupin=find(sum(w)<0);
type=zeros(1,nNeuron);
type(groupex)=1;
type(groupin)=-1;

groupe=zeros(1,nNeuron);
groupe(groupebc)=1;
groupe(groupebd)=2;
groupe(1)=-1;

R=zeros(1,nNeuron);    
R([2 4 17])=1;
CurrentStrength = diag(R);
CurrentInput = zeros(length(time), nNeuron);
size(CurrentInput)
CurrentInput(:,[2 4 17]) =Idep* ones(size(CurrentInput(:,[2 4 17])));

v = c;                                                  % reset voltage;
s = zeros(1,nNeuron);                                   % spike occurance
u = zeros(1,nNeuron);                                       % revocery variable
Vout = zeros(length(time),nNeuron) ;
Iout = zeros(length(time),nNeuron) ;
Init=zeros(nNeuron,1);
type0=type((Init==1));
groupe0=groupe((Init==1));
if ~isempty(groupe0)>0
    firings=[zeros(sum(Init),1),find(Init==1),type0,groupe0];
else
    firings=zeros(0,4);  
end
cumul=0;
hbar = waitbar(0,'En cours');
for t=time
    format short
    tempsavante=t-P+(2)*Delex;            
    tempsavanti=t-P+(3)*Delex;

    randI=randn(1,nNeuron);
    indquie=find(firings(:,1)<=tempsavante & firings(:,1)>(tempsavante-dt) & firings(:,3)==1);
    indquii=find(firings(:,1)<=tempsavanti & firings(:,1)>(tempsavanti-dt) & firings(:,3)==-1);
    qui=firings([indquie;indquii],2);
    s(qui)=1;
    Istate =  s*SynStrength+Istate .* SynDecayex.*(Istate>0) +Istate .* SynDecayin.*(Istate<0); 
    I = Istate +  Is + CurrentInput(tindex,:)*CurrentStrength +vmod*Eml1(end)*(modul==1)+epsilon*randI;
    v = v + dt * (0.04*v.^2+5*v+140-u+I);
    u = u + dt * a .* (b.*v-u);
    Vout(tindex,:) = v ;
    Iout(tindex,:) = I ;
    s = zeros(1,nNeuron) ;
    fired=find(v>=30); 
    lfire=sum(fired==1);
    modulation(lfire)
    if ~isempty(fired)
        firings=[firings; t+0*fired', fired',(type(fired))',(groupe(fired))'];
        v(fired) = c(fired) ;
        u(fired) = u(fired)+d(fired) ;
    end;
    tindex = tindex+1;
    waitbar(t/tEnd,hbar);
end 


% Create figure
figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],...
    'Color',[1 1 1]);

% Create textbox
annotation(figure1,'textbox','String',{'A'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.005 0.97 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'B'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.005 0.65 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'C'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.005 0.40 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'D'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.005 0.16 0.004 0.001]);


% Create textbox
annotation(figure1,'textbox','String',{'excitatory neurone activities'},...
    'FontWeight','bold',...
    'FontSize',9,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.09 0.95 0.4 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'Bc'},...
    'FontWeight','bold',...
    'FontAngle','italic',...
    'FontSize',8,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.30 0.16 0.004 0.001]);
% Create textbox
annotation(figure1,'textbox','String',{'FS'},...
    'FontWeight','bold',...
    'FontAngle','italic',...
    'FontSize',9,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.27 0.18 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'Bd'},...
    'FontWeight','bold',...
    'FontAngle','italic',...
    'Color',[0.5 0.5 0.5],...
    'FontSize',8,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.73 0.16 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'FS'},...
    'FontWeight','bold',...
    'Color',[0.5 0.5 0.5],...
    'FontAngle','italic',...
    'FontSize',9,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.7 0.18 0.004 0.001]);

figure(1)
subplot('Position',[0.08 0.65 0.90 0.25]);
plot(firings(find(firings(:,4)==-1),1), firings(find(firings(:,4)==-1),2),'diamond','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0], 'MarkerSize', 1 )
hold on
plot(firings(find(firings(:,4)==1),1), firings(find(firings(:,4)==1),2),'diamond','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0], 'MarkerSize', 1 )
plot(firings(find(firings(:,4)==2),1), firings(find(firings(:,4)==2),2),'diamond','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerSize', 1)
set(gca,'XLim', [0 tEnd]);
set(gca,'YLim',[0 nNeuron+1]);
set(gca,'YDir','reverse')
%set(gca,'XTick',[0: Tsim/3: Tsim])
set(gca,'XTick',[])
%set(gca,'XLim', [10000 Tsim]);
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
% text(11000,-5,['Neuron activities'],...
%     'VerticalAlignment','middle',...
%     'HorizontalAlignment','left',...     
%     'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','normal','Color',[0 0 0])
box off


    subplot('Position',[0.08 0.41 0.90 0.22]);
    plot(time,Vout(:,1),'LineWidth',0.8,'Color',[0 0 0], 'LineStyle','-')
    hold on
    plot(time,0.5*Iout(:,1)-150,'LineWidth',0.2,'Color',[0 0 0], 'LineStyle','-')
    plot(time,Vout(:,7)-450,'LineWidth',0.8,'Color',[0 0 0], 'LineStyle','-')
    plot(time,0.5*Iout(:,7)-600,'LineWidth',0.2,'Color',[0 0 0], 'LineStyle','-')
    plot(time,Vout(:,26)-800,'LineWidth',0.8,'Color',[0.5 0.5 0.5], 'LineStyle','-')
    plot(time,0.5*Iout(:,26)-950,'LineWidth',0.2,'Color',[0.5 0.5 0.5], 'LineStyle','-')
    set(gca,'ylim', [-1020 100])
    box off
    set(gca,'XTick',[])
    set(gca,'YTick',[])
set(gcf,'Color',[1 1 1])
text(1000,0,['l1 '],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','italic','Color',[0 0 0])
text(1000,-320,['Bc7 '],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','italic','Color',[0 0 0])
text(1000,-750,['Bd26 '],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','italic','Color',[0.5 0.5 0.5])
box off
minV=-55;
Voutf=(Vout>=minV).*Vout +(Vout<minV)*minV;
vtotbc=sum(Voutf(:,groupebc),2);
vtotbd=sum(Voutf(:,groupebd),2);
subplot('Position',[0.08 0.25 0.90 0.13]);
hold on
plot(vtotbc+50,'LineWidth',0.5,'Color',[0 0 0], 'LineStyle','-');
plot(vtotbd-350,'LineWidth',0.5,'Color',[0.5 0.5 0.5], 'LineStyle','-');
set(gca,'XTickLabel',[])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
set(gcf,'Color',[1 1 1])
set(gca,'ylim', [-750 200])
box off


set(gca,'XTick',[])
set(gca,'YTick',[])

set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(1000,100,['sum Bc '],...
     'VerticalAlignment','middle',...
     'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','italic','Color',[0 0 0])
text(1000,-320,['sum Bd '],...
     'VerticalAlignment','middle',...
     'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','italic','Color',[0.5 0.5 0.5])
box off

window1=round(P/dt);          

N=length(vtotbc);
fe=1000/dt;
Dt=1/fe;
tmax=Dt*N;                                      % secondes
t = 0:Dt:tmax-Dt;                               % secondes
filtvtotbc=filtfilt(ones(1,window1)/window1,1,vtotbc);
filtvtotbc=filtvtotbc-mean(filtvtotbc);
filtvtotbd=filtfilt(ones(1,window1)/window1,1,vtotbd);
filtvtotbd=filtvtotbd-mean(filtvtotbd);
subplot('Position',[0.08 0.08 0.90 0.15]);
plot(t(1000:end-1000),filtvtotbc(1000:end-1000),'LineWidth',0.5,'Color',[0 0 0], 'LineStyle','-')
hold on
plot(t(1000:end-1000),filtvtotbd(1000:end-1000),'LineWidth',0.5,'Color',[0.5 0.5 0.5], 'LineStyle','-')

set(gca,'YLim',[-1 10])

set(gca,'FontName','arial','FontWeight','bold','FontSize',9)

box off
xlabel('time (s)');
close(hbar)