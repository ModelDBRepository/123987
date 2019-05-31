%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.
%Fig.6

clear all

global beta gamma maxAc delta  Em0 cumul  Ac Eml1
  
matdep =[ 0  0  -1   
          1  0  -1   
          0  1  0 ];                      
Rdep=[1 0 0 ]; 
matL=[0 -1 0 ; 0 -1 0; 0 0 0];
X=6;
Y=2;
mattoub=boucladj(matdep,Rdep', X, [4,2],1, [2,5], -1);  
Spre=mattoub(:,1:end-1);
Rpre=mattoub(:,end);
Spara=bouclpara(Spre,Rpre, Y); 
R=Spara(:,end);
S=Spara(:,1:end-1);
R=R';
vectv=zeros(size(S,1),1);
vectv([1:Y*(1+X*2)])=1;
vecth=zeros(1, size(S,2));
inhibsl=[];
for y=1:Y
    inhibsl=[inhibsl, [0:X-1]*2+3+(y-1)*((X*2)+1)];
end
inhibsl;
vecth(inhibsl)=-1;
excitateurs=[1,[1:X]*2];
for y=2:Y
    excitateursy=[1+(y-1)*((X*2)+1),[1:X]*2+(y-1)*((X*2)+1)];
    excitateurs=[excitateurs;excitateursy];
end
excitateurs;
taillantiphase=2;
groupebc=[ ];
groupebd=[ ];
for i=1:size(excitateurs,2)
    if mod(floor((i-1)/taillantiphase),2)==0
        groupebc=[groupebc,excitateurs(1,i)];
    else
        groupebd=[groupebd,excitateurs(1,i)];
    end
end
groupebc=groupebc(2:end);

for y=2:Y
    groupebcy=[ ];
    groupebdy=[ ];
    for i=1:size(excitateurs,2)
        if mod(floor((i-1)/taillantiphase),2)==0
            groupebcy=[groupebcy,excitateurs(y,i)];
        else
            groupebdy=[groupebdy,excitateurs(y,i)];
        end
    end
    groupebc=[groupebc,groupebcy(2:end)];
    groupebd=[groupebd,groupebdy];
end
groupebdin=groupebd+1;                  
vectv(groupebd)=1;                      
vectv(groupebdin)=0;           
S=completelop(matL,vectv,vecth,S,1);
R=[0 1 0,R];
inhibsl=inhibsl+3;
excitateurs=excitateurs+3;
groupebc=groupebc+3;
groupebd=groupebd+3;
groupebdin=groupebdin+3;
excitateurs=excitateurs(:,2:end);
excitateurs=reshape(excitateurs',1,size(excitateurs,1)*size(excitateurs,2));
v=1:Y;
for i=1:X
    for j=1:Y
        for k=1:Y
            S(excitateurs((k-1)*X+i),inhibsl((j-1)*X+i))=-1;
        end
    end
end

groupex=find(sum(S)>=0);
groupin=find(sum(S)<0); 
Ne=length(groupex);
Ni=length(groupin);
Ne=length(groupex);     Ni=length(groupin);
retards=ones(Ne+Ni);
th=[0.5*ones(Ne,1); 	0.5*ones(Ni,1)];

modul=1;        
epsilon=0.15;              
beta=1.05 ;            
gamma=0.2 ;             
delta=6;
maxAc=8;
Em0=0.1;

%***************************
vseuiL=[1;zeros(Ne+Ni-1,1)];   
Tsim=200;
seedscalevol=1;
randn('seed',seedscalevol) ;    
rand('seed',seedscalevol);

Init=zeros(Ne+Ni,1);
type=(sign(sum(S)))';
groupe=zeros(size(S,1),1);
groupe(1)=3; 
groupe(groupebc)=1;
groupe(groupebd)=2;
nirings=[zeros(sum(Init),1),find(Init==1),type((Init==1)),groupe((Init==1))];
hstock=zeros(Ne+Ni,Tsim);
hstock(:,1)=Init;
matact=zeros(Ne+Ni,Tsim);
matact(:,1)=Init;
Ac=0;
Eml1=0;
retL=1;   
%*****************************

delaiaff=0;
finaff=Tsim;
Tomax=max(max(retards));
cumul=0;

for t=1:Tsim
    feu=[];
    matfiltre=zeros(Ne+Ni);
    h=R';
    if ~isempty(nirings)
       for ret=1:max(max(retards))
           fincre=zeros(Ne+Ni);                           
           inter=find(nirings(:,1)==t-ret);               
           feu=[feu; ret+0*inter,inter];                  
           finter=[[nirings(inter,2)]', Ne+Ni+1];          
           fincre(1:Ne+Ni,finter)=1;                       
           fincre=fincre(:,1:Ne+Ni);                       
           matinter=(retards==ret);                        
           matfiltre=matfiltre+((matinter==fincre)&(matinter>0));   
       end
       tfire=nirings(feu(:,2),1);                          
       firedt=nirings(feu(:,2),2);                         
       st=sort(firedt); 
       
       if ~isempty(st)>0                                     
           dst=diff(st);                                   
           y=st(dst>0);
           firedtone=[y;st(end)];                          
           Sfiltre=S.*matfiltre;                           
           h=h+sum(Sfiltre(:,firedtone),2);  
       end
       inter1=find(nirings(:,1)==t-retL);   
         finter1=[nirings(inter1,2)];
         finter2=[nirings(inter1,1)];
         l1fire=(sum(finter1==1)>0);  
         lfire=(l1fire>0);
         modulation(lfire);
    end
    randh=randn(Ne+Ni,1);

    hent=h+Eml1(end)*(modul==1)*vseuiL+epsilon*randh-th;                                                        
    hstock(:,t+1)=hent;
    n=(hent>0);                                                
    fired=find(n>0);
    matact(:,t+1)=n; 
    
    if ~isempty(fired);
        nirings=[nirings; t+0*fired,fired,type(fired),groupe(fired)];                
    end
 
end

matactbc=matact(groupebc,:);
matactbd=matact(groupebd,:);
matactl=matact(1,:);
actot=sum(matact(groupex(2:end),:));
actotin=sum(matact(groupin,:));
     indexc=find(nirings(:,4)==1 & nirings(:,1));
     indexd=find(nirings(:,4)==2 & nirings(:,1));
     indexl1=find(nirings(:,4)==3 & nirings(:,1)) ;
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
annotation(figure1,'textbox','String',{'a'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.03 0.85 0.004 0.001]);
% Create textbox
annotation(figure1,'textbox','String',{'RG'},'FontName','Arial',...
    'FontAngle','italic',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.63 0.975 0.004 0.001]);

% Create rectangle
annotation(figure1,'rectangle','LineStyle',':','FaceColor','flat',...
    'Position',[0.61 0.84 0.26 0.14]);

% Create textbox
annotation(figure1,'textbox','String',{'PG'},'FontName','Arial',...
    'FontAngle','italic',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.75 0.81 0.004 0.001]);
% Create textbox
annotation(figure1,'textbox','String',{'Bc'},'FontName','Arial',...
    'FontAngle','italic',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.63 0.70 0.004 0.001]);
% Create textbox
annotation(figure1,'textbox','String',{'Bd'},'FontName','Arial',...
    'FontAngle','italic',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.8 0.70 0.004 0.001],...
    'Color',[0.5 0.5 0.5]);

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
    'Position',[0.005 0.35 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'D'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.005 0.13 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'excitatory neurone activities'},...
    'FontWeight','bold',...
    'FontSize',9,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.1 0.68 0.4 0.001]);

% Create textbox
annotation(figure1,'textbox',...
    'String',{'number of active excitatory Bc and Bd neurones'},...
    'FontWeight','bold',...
    'FontSize',9,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.1 0.42 0.6 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'Bc'},'FontWeight','bold',...
    'FontSize',8,...
    'FontName','Arial',...
    'FontAngle','italic',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.33 0.35 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'OS'},'FontWeight','bold',...
    'FontSize',9,...
    'FontName','Arial',...
    'FontAngle','italic',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.3 0.37 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'Bd'},'FontWeight','bold',...
    'FontSize',8,...
    'FontName','Arial',...
    'FontAngle','italic',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.73 0.35 0.004 0.001],...
    'Color',[0.5 0.5 0.5]);

% Create textbox
annotation(figure1,'textbox','String',{'OS'},'FontWeight','bold',...
    'FontSize',9,...
    'FontName','Arial',...
    'FontAngle','italic',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.7 0.37 0.004 0.001],...
    'Color',[0.5 0.5 0.5]);

% Create textbox
annotation(figure1,'textbox','String',{'Em'},'FontWeight','bold',...
    'FontSize',9,...
    'FontName','Arial',...
    'FontAngle','italic',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.1 0.15 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'l1'},'FontWeight','bold',...
    'FontSize',8,...
    'FontName','Arial',...
    'FontAngle','italic',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.14 0.13 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'c'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.56 0.85 0.004 0.001]);


% Create textbox
annotation(figure1,'textbox','String',{'b'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.4 0.85 0.004 0.001]);
     
set(gcf,'Color',[1 1 1])
if exist('fig6newBcBd.jpg')
subplot('Position',[0.61 0.68 0.30 0.30]);
    fig1=imread('fig6newBcBd.jpg');
    imagesc(fig1);
    axis off;
end  

subplot('Position',[0.12 0.72 0.25 0.25]);
imagesc(S)
notremap=[0 0 0;1 1 1;0.5 0.5 0.5];
colormap(flipud(notremap))
set(gcf,'Color',[1 1 1])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
box on

subplot('Position',[0.48 0.72 0.02 0.25]);
imagesc(R',[-1 1])
set(gcf,'Color',[1 1 1])
set(gca,'XTick',[])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)

box on
tdep=50;
subplot('Position',[0.08 0.43 0.85 0.2]);
plot(nirings(indexc,1), nirings(indexc,2),'diamond','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0], 'MarkerSize', 1 )
hold on
plot(nirings(indexl1,1), nirings(indexl1,2),'diamond','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0], 'MarkerSize', 1)
plot(nirings(indexd,1), nirings(indexd,2),'o','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerSize', 1)
set(gca,'XLim', [tdep Tsim],'YLim',[0 Ne+Ni+1]);
set(gca,'YDir','reverse')

set(gca,'XTick',[])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
box off

subplot('Position',[0.08 0.18 0.85 0.22]);
actotbc=sum(matact([1 groupebc],:));
  actotbd=sum(matact(groupebd,:));
     plot(actotbc,'LineWidth',0.5,'Color',[0 0 0], 'LineStyle','-');
     hold on
     plot(actotbd,'LineWidth',0.5,'Color',[0.5 0.5 0.5], 'LineStyle','-');
     hold off
     set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
     set(gca,'XLim',[tdep Tsim]);
     set(gca,'XTick',[]);

box off

subplot('Position',[0.08 0.08 0.85 0.04]);
plot(Eml1,'k')
set(gca,'XLim', [tdep Tsim],'YLim',[min(Eml1)  10]);
set(gca,'YTick',[0 10])
set(gca,'XTickLabel',[0 5 10 15])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
xlabel('time (s)');


box off

set(gcf,'Color',[1 1 1],'position',[100 100 520 420]);