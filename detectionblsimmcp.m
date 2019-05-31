%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.

%detection of lung episodes

function [durmoy,intermoy]=detectionblsimmcp(signal,fe,fbase,pourcent)

global amplitude difamplitude
global maxetiq maxetiqcor indiceun dindiceun indebepil
signal=signal(100:end);
dsignal=diff(signal);
N1=length(signal);
maxlocal=signal(1);
tempsmaxlocal=2;
minlocal=signal(1);
tempsminlocal=2;
maxilocaux=maxlocal;
tempsmaxilocaux=1;
minilocaux=minlocal;
tempsminilocaux=1;

dt=1/fe;
signalfnu = fft(signal-mean(signal),N1);
nu1=[-N1/2:N1/2-1]/(N1*dt);

signalfnu2=fftshift(signalfnu);
signalfnu3=abs(signalfnu2(round(end/2)+1:end));
nu1cor=nu1(round(end/2)+1:end);

nu2cor=nu1cor(round(length(nu1cor)*0.5*2/fe):end);                                      
signalfnu4=signalfnu3(round(length(nu1cor)*0.5*2/fe):end);
maxfreq=abs(nu2cor((signalfnu4==max(signalfnu4))));

persignal=1/fbase;
perpointssignal=persignal*fe;
tempsref=perpointssignal/2;
              
compteurtempsmax=0;
compteurtempsmin=0;
for i=2:length(signal)
    if signal(i)>maxlocal
        maxlocal=signal(i);
        tempsmaxlocal=i;
         compteurtempsmax=0;
    else compteurtempsmax=compteurtempsmax+1;
    end
    
     if signal(i)<minlocal
        minlocal=signal(i);
        tempsminlocal=i;
         compteurtempsmin=0;
    else compteurtempsmin=compteurtempsmin+1;
    end
     
     
     if compteurtempsmax>tempsref
         a=dsignal(tempsmaxlocal);
         b=dsignal(tempsmaxlocal-1);
         if   a*b<0                                            
            maxilocaux=[maxilocaux,maxlocal];
            tempsmaxilocaux=[tempsmaxilocaux,tempsmaxlocal];           
       end
         maxlocal=min(signal);
         tempsmaxlocal=i;
         compteurtempsmax=0;   
     end
    
    if compteurtempsmin>tempsref
         a=dsignal(tempsminlocal);
         b=dsignal(tempsminlocal-1);
         if  a*b<0
            minilocaux=[minilocaux,minlocal];
            tempsminilocaux=[tempsminilocaux,tempsminlocal];
         end
         minlocal=max(signal);
         tempsminlocal=i;
         compteurtempsmin=0;
         
    end
end

vmin=[(tempsminilocaux(2:end))',(minilocaux(2:end))'];
vmax=[(tempsmaxilocaux(2:end))',(maxilocaux(2:end))'];


vmininterp=interp1([tempsminilocaux(2:end), length(signal)],[minilocaux(2:end), signal(end)],1:tempsref/2:length(signal));
vmaxinterp=interp1([tempsmaxilocaux(2:end),length(signal)] ,[maxilocaux(2:end), signal(end)],1:tempsref/2:length(signal));
amplitude=vmaxinterp-vmininterp;
difamplitude=diff(amplitude);


[xref,vtot]=smoothist(amplitude,100);
refb=xref((vtot==max(vtot)))+mean(vmin(:,2));
pourcentb=pourcent;                                                                  
amplb=refb+pourcentb*refb/100;
pourcentl=pourcent;                                                                   
ampll=(max(vmax(:,2))-amplb)*(pourcentl/100)+amplb;

indmaxb=find(vmax(:,2)<=amplb);
vmaxb=vmax(indmaxb,:);
indmaxl=find(vmax(:,2)>ampll);
vmaxl=vmax(indmaxl,:);
indmaxnibnil=find(vmax(:,2)>amplb & vmax(:,2)<=ampll);
vmaxnibnil=vmax(indmaxnibnil,:);


prepmaxetiq=[vmaxb(:,1), zeros(size(vmaxb(:,1))); vmaxl(:,1), ones(size(vmaxl(:,1))); vmaxnibnil(:,1), zeros(size(vmaxnibnil(:,1)))];
[maxetiq,ind]=sort(prepmaxetiq(:,1));
maxetiq=[maxetiq,prepmaxetiq(ind,2)];

maxetiqcor=maxetiq(:,2)+circshift(maxetiq(:,2),[1,0])+circshift(maxetiq(:,2),[-1,0]);
maxetiqcor=(maxetiqcor>0);

indiceun=find(maxetiqcor==1) ;                                          
dindiceun=diff(indiceun);
indebepil=find(dindiceun>1);
nombrepil=length(indebepil)+1;

if nombrepil>1  ,
debutepil(1)=maxetiq(indiceun(1),1);
finepil(1)=maxetiq(indiceun(indebepil(1),1));
for i=2:nombrepil-1
    debutepil(i)=maxetiq(indiceun(indebepil(i-1)+1,1));
    finepil(i)=maxetiq(indiceun(indebepil(i),1));
end
debutepil(nombrepil)=maxetiq(indiceun(indebepil(nombrepil-1)+1,1));
finepil(nombrepil)=maxetiq(indiceun(end,1));

durepil=finepil-debutepil;
durmoy=mean(durepil)/fe;                                           %secondes
interepil=debutepil(2:end)-finepil(1:end-1);
intermoy=mean(interepil)/fe;                                        %secondes

vectzones=zeros(debutepil(1),1);
vectzones=[vectzones; ones(finepil(1)-debutepil(1), 1)];
for i=2:length(debutepil)
    vectzones=[vectzones; zeros(debutepil(i)-finepil(i-1), 1)];
    vectzones=[vectzones; ones(finepil(i)-debutepil(i), 1)];
end
vectzones=[vectzones; zeros(length(signal)-finepil(end), 1)];

signalb=signal((vectzones==0));
signall=signal((vectzones==1));

else
    disp('no lung episode');
end







