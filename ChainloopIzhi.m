%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.
%Models from  Izhikevich E : http://nsi.edu/users/izhikevich/publications/spikes.htm
%Modified by Bruce Land -- BioNB 222, March 2008 

format short
nNeuron = 11 ;
nCurrent = nNeuron;
tEnd =4000;                         %maximum simulation time
time = dt:dt:tEnd ;
pars=[0.02      0.2     -65     8       14 ];    %1-tonic spiking
nType = 1*ones(1,nNeuron) ;
a=pars(nType,1)';
b=pars(nType,2)';
c=pars(nType,3)';
d=pars(nType,4)';
Is=0*ones(nNeuron,1)';
Gmaxsynex=50;
Gmaxsynin=30;
Idep=4.8;
SynRatex = 1*ones(1,nNeuron);                   % 3 msec
SynRatin = 0.026*ones(1,nNeuron);            
SynDecayex = 1 - SynRatex*dt;
SynDecayin = 1 - SynRatin*dt;
tindex = 1;                                     %pointer for output arrays
Istate = 0;                                     %state variable used for synaptic currents decay
vmod=zeros(1,nNeuron);
vmod(1)=1;
I = Is ;
v = c;                                          % reset voltage;
s = zeros(1,nNeuron);                           % spike occurence
u = zeros(1,nNeuron);                           % recovery variable
w = ...
[ 0  0 -1  0  0  0  0  0  0  0  0
1  0 -1  0 -1  0  0  0  0  0  0
0  1  0  0  0  0  0  0  0  0  0
0  1  0  0 -1  0 -1  0  0  0  0
0  0  0  1  0  0  0  0  0  0  0
0  0  0  1  0  0 -1  0 -1  0  0
0  0  0  0  0  1  0  0  0  0  0
0  0  0  0  0  1  0  0 -1  0 -1
0  0  0  0  0  0  0  1  0  0  0
0  0  0  0  0  0  0  1  0  0 -1
0  0  0  0  0  0  0  0  0  1  0] ;
SynStrength =Gmaxsynex* (w'>0) - Gmaxsynin* (w'<0); ...
R=zeros(1,nNeuron);    
R(1)=1;
CurrentStrength = diag(R);
CurrentInput = zeros(length(time), nNeuron);
CurrentInput(:,1) =Idep* ones(size(CurrentInput(:,1)));

v = c;                                              % reset voltage;
s = zeros(1,nNeuron);                               % spike occurrence
u = zeros(1,nNeuron);                               % recovery variable
Vout = zeros(length(time),nNeuron) ;
Iout = zeros(length(time),nNeuron) ;
Init=zeros(nNeuron,1);
groupex=find(sum(w)>=0);
groupin=find(sum(w)<0);
type(groupex)=1;
type(groupin)=-1;
type0=type((Init==1));
maxI=50;
if ~isempty(type0)
    firings=[zeros(sum(Init),1),find(Init==1),type0];
else
    firings=zeros(0,3);  
end
for t=time
  
    tempsavante=t-P+(2)*Delex;            
    tempsavanti=t-P+(3)*Delex;
    indquie=find(firings(:,1)<=tempsavante & firings(:,1)>(tempsavante-dt) & firings(:,3)==1);
    indquii=find(firings(:,1)<=tempsavanti & firings(:,1)>(tempsavanti-dt) & firings(:,3)==-1);
    qui=firings([indquie;indquii],2);
    s(qui)=1;
    Istate =  s*SynStrength+Istate .* SynDecayex.*(Istate>0) +Istate .* SynDecayin.*(Istate<0); 
    I = Istate +  Is + CurrentInput(tindex,:)*CurrentStrength;
    v = v + dt * (0.04*v.^2+5*v+140-u+I);
    u = u + dt * a .* (b.*v-u);
    Vout(tindex,:) = v ;
    Iout(tindex,:) = I ;
    s = zeros(1,nNeuron) ;
    fired=find(v>=30);                      % indices of cells that spiked
    lfire=sum(fired==1);
    if ~isempty(fired)
        firings=[firings; t+0*fired', fired',(type(fired))];
        v(fired) = c(fired) ;
        u(fired) = u(fired)+d(fired) ;
       
    end;
    tindex = tindex+1;
end                                         % main time step for-loop


figure(1)
subplot('Position',[0.55 0.80 0.40 0.13]);
plot(time,Vout(:,2),'k','LineWidth',0.8)
set(gca,'ylim', [-220 100])
hold on
plot(time,Iout(:,2)-150,'k','LineWidth',0.2)
set(gca,'YTickLabel',{'','-100',' 0', '100'});
set(gca,'XLim',[0 2000]);
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
set(gca,'Xtick',[0 1000 2000]);
set(gca,'XTickLabel',[0 :tEnd/4000:2*tEnd/4000])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
box off

subplot('Position',vpos3);
plot(firings(find(firings(:,3)>=0 ),1), firings(find(firings(:,3)>=0  ),2),'diamond','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0], 'MarkerSize', 2 )
hold on
plot(firings(find(firings(:,3)<0 ),1), firings(find(firings(:,3)<0 ),2),'diamond','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerSize', 1)
set(gca,'XLim', [0 4000]);
set(gca,'YLim',[0 nNeuron+1]);
set(gca,'YDir','reverse')
set(gca,'XTick',[0000:1000: 4000])
set(gca,'XTickLabel',[[0 :tEnd/4000:tEnd/1000]])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)

box off

subplot('Position',vpos4);
minV=-55;
Voutf=(Vout>=minV).*Vout +(Vout<minV)*minV;
vtot=sum(Voutf(:,groupex),2);
window1=round(P/dt);              %time steps
N=length(vtot);
fe=1000/dt;
Dt=1/fe;
tmax=Dt*N;                      % time, s
t = 0:Dt:tmax-Dt;               % time, s
filtvtot=filtfilt(ones(1,window1)/window1,1,vtot);
filtvtot=filtvtot-mean(filtvtot);
plot(filtvtot(1000:31000),'k')
set(gca,'YLim',[-1 1])
set(gca,'YTick',[-1:1:1])
set(gca,'XLim',[0000 32000])
set(gca,'XTick',[0000 8000 16000 24000 32000])
set(gca,'XTickLabel',[[0:4]])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9) 
box off
xlabel('time (s)');
text(800,1,['filtered sum of spikes'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','normal','Color',[0 0 0]) 
