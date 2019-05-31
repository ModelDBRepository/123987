%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.
%Fig.7
clear all

global beta gamma maxAc delta Em0 cumul  Ac  Eml1

germe=18   ;
facteur=0.15;
nombrerreur= 60;
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

S1=matmodif(S,nombrerreur,facteur,germe);
B=nnz(S1);
groupex=find(sum(S1)>=0);
groupin=find(sum(S1)<0); 
Ne=length(groupex);
Ni=length(groupin);
Ne=length(groupex);     Ni=length(groupin);
retards=ones(Ne+Ni);
th=[0.5*ones(Ne,1); 	0.5*ones(Ni,1)];

modul=1;        
epsilon=0.15;              
beta=1.07 ;            
gamma=0.2 ;             
delta=6;
maxAc=6;
Em0=0.1;


vseuiL=[1;zeros(Ne+Ni-1,1)];   

seedscalevol=1;
randn('seed',seedscalevol)     
rand('seed',seedscalevol)

%initialisations 
Init=zeros(Ne+Ni,1);
type=(sign(sum(S1)))';
groupe=zeros(size(S1,1),1);
groupe(1)=3; 
Tsim=150;            
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
       
       if ~isempty(st)                                   
           dst=diff(st);                                   
           y=st(dst>0);
           firedtone=[y;st(end)];                          
           Sfiltre=S1.*matfiltre;                           
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
     indexc1=find(nirings(:,4)==1 & nirings(:,1));
     indexd1=find(nirings(:,4)==2 & nirings(:,1));
     indexl1=find(nirings(:,4)==3 & nirings(:,1)) ;

figure
set(gcf,'Color',[1 1 1])    

set(gca,'XTickLabel',[])
subplot('position',[0.12 0.68 0.28 0.28]);
imagesc(S);
 set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
notremap=[1 0 0;1 1 1;0 0 1];
colormap(flipud(notremap))
notremap=[0 0 0;1 1 1;0.5 0.5 0.5];
colormap(flipud(notremap))
subplot('position',[0.55 0.68 0.28 0.28]);

imagesc(S1, [-1 1] );
 set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
notremap=[0.5 0.5 0.5; repmat([0.85 0.85 0.85],8,1) ;1 1 1;repmat([0.85 0.85 0.85],8,1) ; 0 0 0];
colormap(notremap)
text(-55,1,['A'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','normal','FontSize',16,'FontAngle','normal','Color',[0 0 0])
text(-52,10,['a'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','normal','FontSize',16,'FontAngle','normal','Color',[0 0 0])
text(-8,10,['b'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','normal','FontSize',16,'FontAngle','normal','Color',[0 0 0])
subplot('position',[0.08 0.28 0.85 0.30])
plot(nirings(indexl1,1),nirings(indexl1,2),'diamond','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0], 'MarkerSize', 1 )
 hold on
plot(nirings(indexc1,1),nirings(indexc1,2),'diamond','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0], 'MarkerSize', 1 )
plot(nirings(indexd1,1),nirings(indexd1,2),'diamond','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerSize', 1 )
hold off
     set(gca,'YLim',[0 Ne+Ni+1]);
     set(gca,'YDir','reverse')
     set(gca,'XLim',[0 Tsim]);
     set(gca,'XTick',[0: Tsim/4: Tsim])
     set(gca,'XTickLabel',[])
     set(gca,'FontName','arial','FontWeight','bold','FontSize',10)
 text(-12,5,['B'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','normal','FontSize',16,'FontAngle','normal','Color',[0 0 0])
text(5,-2,['excitatory neurone activities'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','normal','Color',[0 0 0])    
     actot=sum(matact(groupex(1:end),:));
     actotbc1=sum(matact([1 groupebc],:));
     actotbd1=sum(matact(groupebd,:));
    
    box off
subplot('position',[0.08 0.09 0.85 0.14])
     plot(actotbc1,'LineWidth',0.5,'Color',[0 0 0], 'LineStyle','-');
     hold on
     plot(actotbd1,'LineWidth',0.5,'Color',[0.5 0.5 0.5], 'LineStyle','-');
     hold off
     set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
set(gca,'XLim', [0 Tsim],'YLim',[0  10]);
set(gca,'XTickLabel',[0 5 10 15])
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
xlabel('time (s)');
 box off
     text(37,4.5,['Bc'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','bold','FontSize',8,'FontAngle','italic','Color',[0 0 0])
text(77,4.5,['Bd'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','bold','FontSize',8,'FontAngle','italic','Color',[0.5 0.5 0.5])

text(5,11,['number of active Bc and Bd neurones'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','normal','Color',[0 0 0]) 
text(30,5,['OS'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','italic','Color',[0 0 0]) 
text(70,5,['OS'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','italic','Color',[0.5 0.5 0.5]) 
 text(-12,7,['C'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...     
    'FontName','arial','FontWeight','normal','FontSize',16,'FontAngle','normal','Color',[0 0 0])
box off
set(gcf,'Color',[1 1 1],'position',[100 100 520 420]);
