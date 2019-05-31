%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.
%Fig.2

clear all

global  beta gamma delta Em0 cumul Ac maxAc Eml1 
S = [0    -1     
     0    -1 ];
R = [0    1  ];
groupex=find(sum(S)>=0);
groupin=find(sum(S)<0); 
Ne=length(groupex);     Ni=length(groupin);
retards=ones(Ne+Ni);
th=[0.5*ones(Ne,1); 	0.5*ones(Ni,1)];

modul=1;            %l1 automodulation if modul=1      
epsilon=0.06;       %dynamical noise
beta=1.05 ;          
gamma=0;            %modulation noise
maxAc=6;            
delta=0;            %noise on maxAc     
Em0=0.1;              
vseuiL=[1];
seedscalevol=5;randn('seed',seedscalevol) ;    rand('seed',seedscalevol);
Tsim=150;
Init=zeros(Ne+Ni,1);
type=(sign(sum(S)))';
nirings=[zeros(sum(Init),1),find(Init==1),type((Init==1))];
hstock=zeros(Ni+Ne,Tsim);
hstock(:,1)=Init;
matact=zeros(Ni+Ne,Tsim);
matact(:,1)=Init;
Ac=0;                   %l1 spike counter
Eml1=0;                 %l1 automodulation initialization 
retL=1;                 %l1 postsynaptic delay
%*****************************
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

         if ~isempty(feu)                                           
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
         else
             h=R';
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
    n=(hent>=0);                                                
    fired=find(n>0);
    matact(:,t+1)=n; 
  if ~isempty(fired)
        nirings=[nirings; t+0*fired,fired,type(fired)];                 
  end
end



% Create figure
figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],...
    'Color',[1 1 1]);
% Create rectangle
annotation(figure1,'rectangle','FaceColor',[0 0 0],...
    'Position',[0.2786 0.7476 0.04 0.04]);
% Create line
annotation(figure1,'line',[0.2186 0.2286],[0.8976 0.8776],'LineWidth',0.6);
% Create line
annotation(figure1,'line',[0.2286 0.2186],[0.8776 0.8576],'LineWidth',0.6);
% Create line
annotation(figure1,'line',[0.1886 0.2186],[0.8976 0.8976],'LineWidth',0.6);
% Create line
annotation(figure1,'line',[0.2186 0.1886],[0.8576 0.8576],'LineWidth',0.6);
% Create line
annotation(figure1,'line',[0.1886 0.1786],[0.8576 0.8776],'LineWidth',0.6);
% Create line
annotation(figure1,'line',[0.1786 0.1886],[0.8776 0.8976],'LineWidth',0.6);
% Create arrow
annotation(figure1,'arrow',[0.2286 0.3336],[0.8776 0.8776],'HeadLength',7,...
    'HeadWidth',7,...
    'LineStyle','--');
% Create arrow
annotation(figure1,'arrow',[0.2686 0.2786],[0.7676 0.7676],'HeadLength',7,...
    'HeadWidth',7);
% Create arrow
annotation(figure1,'arrow',[0.2858 0.2686],[0.6676 0.6676],'HeadLength',7,...
    'HeadWidth',7);
% Create rectangle
annotation(figure1,'rectangle','LineStyle','--','FaceColor','flat',...
    'Position',[0.1786 0.6676 0.12 0.1]);
% Create rectangle
annotation(figure1,'rectangle','FaceColor',[1 1 1],...
    'Position',[0.1886 0.6376 0.08 0.06]);
% Create textbox
annotation(figure1,'textbox','String','Mod','FontWeight','bold',...
    'FontSize',9,...
    'FontName','Arial',...
    'FontAngle','italic',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.1886 0.6976 0 0.008]);
% Create arrow
annotation(figure1,'arrow',[0.3386 0.2986],[0.8676 0.7876],'HeadLength',7,...
    'HeadWidth',7,...
    'Color',[0.5 0.5 0.5]);
% Create textbox
annotation(figure1,'textbox','String',{'E'},'FontWeight','bold',...
    'FontSize',9,...
    'FontName','Arial',...
    'FontAngle','italic',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.1886 0.9476 0.004 0.001]);
% Create arrow
annotation(figure1,'arrow',[0.3836 0.3736],[0.9076 0.8976],'HeadLength',7,...
    'HeadWidth',7,...
    'Color',[0.5 0.5 0.5]);
% Create textbox
annotation(figure1,'textbox','String',{'l1'},'FontWeight','bold',...
    'FontName','Arial',...
    'FontAngle','italic',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.3186 0.8176 0.001 0.001]);
% Create textbox
annotation(figure1,'textbox','String',{'l2'},'FontWeight','bold',...
    'FontName','Arial',...
    'FontAngle','italic',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.3586 0.8676 0.001 0.001],...
    'Color',[0.5 0.5 0.5]);
% Create ellipse
annotation(figure1,'ellipse','Position',[0.3286 0.8876 0.06 0.07],...
    'Color',[0.5 0.5 0.5]);
% Create ellipse
annotation(figure1,'ellipse','FaceColor',[0.5 0.5 0.5],...
    'Position',[0.3336 0.8576 0.05 0.05],...
    'Color',[0.5 0.5 0.5]);
% Create textbox
annotation(figure1,'textbox','String',{'A'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.005 0.85 0.004 0.001]);
% Create textbox
annotation(figure1,'textbox','String',{'B'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'EdgeColor',[1 1 1],...
    'Position',[0.005 0.2262 0.004 0.001]);
%figure (1) 
plot(Eml1,'LineWidth',1,'Color',[0 0 0])


set(gca,'XTick',[0: Tsim/3: Tsim])
text(-13,2.5,['Em'],'FontWeight','bold','FontSize',10,'FontAngle','italic')
text(-5,2.4,['l1'],'FontWeight','bold','FontSize',9,'FontAngle','italic')
text(105,3.85,['max'],...  
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','normal','FontSize',12,'FontAngle','italic')
hold on
%set(gca,'YTickLabel',[])
set(gca,'FontName','arial','FontWeight','bold','FontSize',10)
plot(Ac/(2*max(Ac))-1,'LineWidth',0.5,'Color',[0 0 0])
set(gca,'XTick',[0: Tsim/3: Tsim])
set(gca,'XTickLabel',[0 5 10 15]);
set(gca,'YTickLabel',{'0' '6' '0' '0.5' '1' '1.5' '2' ''});
hold on
set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
xlabel('time (s)');
line([10 length(Eml1)],[-0.2 -0.2],'LineWidth',0.5,'Color',[0.5 0.5 0.5], 'LineStyle',':')
plot(nirings((nirings(:,2)==1),1), nirings((nirings(:,2)==1),2)-1.2,'diamond','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0], 'MarkerSize', 2)
set(gca,'XTick',[0: Tsim/3: Tsim])
set(gca,'LineWidth',0.2)

set(gca,'FontName','arial','FontWeight','bold','FontSize',9)
reperexmax=find(Eml1==max(Eml1));
reperexmax=reperexmax(1);
repxmax=find(Ac==max(Ac));
for i=1:length(repxmax)
line(repxmax(i)*ones(1,2),[-1 0.5],'LineWidth',1,'Color',[0 0 0], 'LineStyle','-.')
end
line([0,Tsim],zeros(1,2),'LineWidth',0.2,'Color',[0 0 0], 'LineStyle','-')
line([0,0.00001],[-0.01,-0.499],'LineWidth',20,'Color',[1 1 1], 'LineStyle','-')


text(65,-0.5,...
     '\leftarrowMaxAc',...
     'FontName','arial','FontWeight','normal','FontSize',10,'FontAngle','italic')
 text(5,0.1,...
     '\leftarrow-----Emo',...
     'FontName','arial','FontWeight','normal','FontSize',10,'FontAngle','italic')
text(-13, -0.5,[' Ac'],...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','bold','FontSize',10,'FontAngle','italic' )
text(80,0.7,[' b'],...  
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','symbol','FontWeight','bold','FontSize',12,'FontAngle','italic')
text(-18,-0.2,[' l1 firings '],...  
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontName','arial','FontWeight','bold','FontSize',9,'FontAngle','italic')
box off
set(gcf,'Color',[1 1 1],'position',[100 100 520 420]);


  