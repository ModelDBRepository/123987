%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.

global hstock nirings seedscalevol Init matactl
global Em0 cumul  Eml1  actot
global Tsim epsilon

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
Em0=0.1;
vseuiL=[1;zeros(Ne+Ni-1,1)];
seedscalevol=1;
randn('seed',seedscalevol)     
rand('seed',seedscalevol)

Init=zeros(Ne+Ni,1);
type=(sign(sum(S)))';
nirings=[zeros(sum(Init),1),find(Init==1),type(Init==1)];
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

matactl=matact(1,:);

actot=sum(matact(groupex,:));
actotin=sum(matact(groupin,:));







    
    
    
    
        