%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.

%global nirings seedscalevol Init
%global vpos1 vpos2 Tsim

matdep =[ 0  0   -1  
          1  0   -1   
          0  1  0 ];                %single loop
Rdep=[1 0 0]; 
mattoub=Boucladj(matdep,Rdep', 5, [4,2],1, [2,5], -1);  
S=mattoub(:,1:end-1);
R=mattoub(:,end);
R=R';
groupex=find(sum(S)>0);             %group of excitatory neurons
groupin=find(sum(S)<0);             %group of inhibitory neurons

Ne=length(groupex);     Ni=length(groupin);
th=[0.5*ones(Ne,1); 	0.5*ones(Ni,1)];
retards=ones(Ne+Ni);                %delay matrix

seedscalevol=3760;
randn('seed',seedscalevol)     
rand('seed',seedscalevol)
Init=zeros(Ne+Ni,1);
type=(sign(sum(S)))';
nirings=[zeros(sum(Init),1),find(Init==1),type((Init==1))]; % neuron states
hstock=zeros(Ni+Ne,Tsim+1);        %potential matrix
hstock(:,1)=Init;
matact=zeros(Ni+Ne,Tsim+1);         %binary activity matrix
matact(:,1)=Init;
epsilon=0;                          %dynamical noise parameter


Tomax=max(max(retards));
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
                             
    end
    randh=randn(Ne+Ni,1);
    hent=h+epsilon*randh-th;                                                                                 
    hstock(:,t+1)=hent;
    n=(hent>0);                                                
    fired=find(n>0);
    matact(:,t+1)=n;                                      
    if ~isempty(fired);
        nirings=[nirings; t+0*fired,fired,type(fired)];                 
    end
end

actot=sum(matact(groupex,:));
subplot('position',vpos1)
set(gcf,'Color',[1 1 1])
plot(nirings((nirings(:,3)>=0),1), nirings((nirings(:,3)>=0),2),'diamond','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0], 'MarkerSize', 2 )
hold on
plot(nirings((nirings(:,3)<0),1), nirings((nirings(:,3)<0),2),'o','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerSize',1 )
set(gca,'XLim', [0 Tsim],'YLim',[0 Ne+Ni+1]);
set(gca,'YDir','reverse')
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
text(10,-2,['MCP model   <-    neurone activities    ->   Izhikevich model'],...
    'VerticalAlignment','middle',...
   'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','normal','Color',[0 0 0])
    box off
subplot('position',vpos2)
plot([0:Tsim],actot,'k');
set(gca,'XLim', [0 Tsim],'YLim',[1 4]);
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
xlabel('time steps');
text(2,4,['number of active excitatory neurones     '],...
    'VerticalAlignment','middle',...
   'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','normal','Color',[0 0 0])
box off


    
    
        