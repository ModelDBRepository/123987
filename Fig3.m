%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.
%Fig.3

clear all
global beta gamma cumul maxAc Em0 Ac Eml1 delta

format long
matdep =[ 0  0  -1   
          1  0  -1   
          0  1  0 ];                      
Rdep=[1 0 0 ]; 
matL=[0 -1 0 ; 0 -1 0; 0 0 0];
X=5;
Y=1;
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

S=completelop(matL,vectv,vecth,S,1);
R=[0 1 0,R];
inhibsl=inhibsl+3;
excitateurs=excitateurs+3;
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
Ne=length(groupex);     Ni=length(groupin);
retards=ones(Ne+Ni);
th=[0.5*ones(Ne,1); 	0.5*ones(Ni,1)];

modul=1;        
epsilon=0;              
beta=1.05 ;             
gamma=0.00;            
maxAc=6;
delta=0;
Em0=0.1;
            

vseuiL=[1;zeros(Ne+Ni-1,1)];    
seedscalevol=3760;
randn('seed',seedscalevol)     
rand('seed',seedscalevol)

Tsim=150;
Init=zeros(Ne+Ni,1);
type=(sign(sum(S)))';
nirings=[zeros(sum(Init),1),find(Init==1),type((Init==1))];
hstock=zeros(Ne+Ni,Tsim);
hstock(:,1)=Init;
matact=zeros(Ne+Ni,Tsim);
matact(:,1)=Init;
Ac=0;
Eml1=0;
retL=1; 


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
    
       if ~isempty(st)                                     
           dst=diff(st);                                   
           y=st(dst>0);
           firedtone=[y;st(end)];                         
           Sfiltre=S.*matfiltre;                           
           h=h+sum(Sfiltre(:,firedtone),2);  
       end
       inter1=find(nirings(:,1)==t-retL);   
       finter1=[nirings(inter1,2)];
       finter2=[nirings(inter1,1)];
       lfire=(sum(finter1==1)>0);   
       modulation(lfire);
       
    end
    randh=randn(Ne+Ni,1);
    hent=h+Eml1(end)*(modul==1)*vseuiL+epsilon*randh-th;                                                            
    hstock(:,t+1)=hent;
    n=(hent>0);                                               
    fired=find(n>0);
    matact(:,t+1)=n;
    if ~isempty(fired);
        nirings=[nirings; t+0*fired,fired,type(fired)]; 
    end
    
end

actot=sum(matact(groupex,:));
actotin=sum(matact(groupin,:));


   % Create figure
figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],...
    'Color',[1 1 1]);

% Create textbox
annotation(figure1,'textbox','String',{'A'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.003 0.97 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'a'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.03 0.85 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'B'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.003 0.65 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'C'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.003 0.4 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'D'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.005 0.22 0.004 0.001]);


% Create textbox
annotation(figure1,'textbox','String',{'b'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.39 0.85 0.004 0.001]);

% Create textbox
annotation(figure1,'textbox','String',{'c'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.52 0.85 0.004 0.001]); 
if exist('model1onechain5loop3N&R.jpg')
subplot('Position',[0.58 0.70 0.38 0.25]);
    fig1=imread('model1onechain5loop3N&R.jpg');
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

subplot('Position',[0.08 0.45 0.85 0.2]);
plot(nirings(find(nirings(:,3)>=0),1), nirings(find(nirings(:,3)>=0),2),'diamond','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0], 'MarkerSize', 1.5 )
hold on
plot(nirings(find(nirings(:,3)<0),1), nirings(find(nirings(:,3)<0),2),'o','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerSize', 1)
set(gca,'XLim', [0 Tsim],'YLim',[0 Ne+Ni+1]);
set(gca,'YDir','reverse')
set(gca,'XTick',[0: Tsim/3: Tsim])
set(gca,'XTick',[])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(4,0,['neurone activities'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','normal','Color',[0 0 0])

box off
line([82 90],[15,15],'LineWidth',4,'Color',[0 0 0], 'LineStyle','-')
text(84,17,['Tp'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','italic','Color',[0 0 0])

subplot('Position',[0.08 0.25 0.85 0.15])
plot([0:Tsim],actot,'k');
set(gca,'XLim', [0 Tsim],'YLim',[min(actot)  max(actot)+1]);
set(gca,'XTick',[0: Tsim/3: Tsim])
set(gca,'XTick',[])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(4,max(actot)+1,['number of active excitatory neurones'],...
    'VerticalAlignment','middle',...
   'HorizontalAlignment','left',...
   'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','normal','Color',[0 0 0])
text(30,max(actot)-2,['OS'],...
    'VerticalAlignment','middle',...
   'HorizontalAlignment','left',...
   'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','italic','Color',[0 0 0])
box off

subplot('Position',[0.08 0.08 0.85 0.12]);
plot(Eml1,'k')
set(gca,'XLim', [0 Tsim],'YLim',[min(Eml1)  max(Eml1)]);
set(gca,'YTick',[0: max(Eml1)/2:max(Eml1)])
set(gca,'YTick',[])
set(gca,'XTick',[0: Tsim/3: Tsim])
set(gca,'XTickLabel',[0 5 10 15])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
xlabel('time (s)');
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(4,max(Eml1),['Em'],'FontWeight','bold','FontSize',9,'FontAngle','italic')
text(11,max(Eml1)-0.5,['l1'],'FontWeight','bold','FontSize',8,'FontAngle','normal')
box off
set(gcf,'Color',[1 1 1],'position',[100 100 520 420]);


    
    
    
    
        