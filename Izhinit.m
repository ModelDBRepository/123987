%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.
%Models from  Izhikevich E : http://nsi.edu/users/izhikevich/publications/spikes.htm
%Modified by Bruce Land -- BioNB 222, March 2008 

%this program performs the computation of interspike time P and mdeltat
format long
nNeuron = 2 ;
nCurrent = nNeuron;
dt = 0.125 ; %millisecond -- time step  for computation accuracy dt must be (1/2^n)
tEnd =1500; %maximum simulation time
time = dt:dt:tEnd ;
pars=[0.02      0.2     -65     8       14 ];    %1-tonic spiking
nType = 1*ones(1,nNeuron) ;
a=pars(nType,1)';
b=pars(nType,2)';
c=pars(nType,3)';
d=pars(nType,4)';
Is=0*ones(nNeuron,1)';

Idep=4.8;                               % steady depolarizing current input 
Gmaxsynex=50;                            %excitation weight 

depmes=5;                               %number of first spikes suppressed in the computation of P

w = [0   0;1 0] ;
SynStrength =Gmaxsynex*abs(w').* (w'>0); ...
R=zeros(1,nNeuron);
R(1)=1;
CurrentStrength = diag(R);
CurrentInput = zeros(length(time), nNeuron);
CurrentInput(:,1) =Idep* ones(size(CurrentInput(:,1)));
SynRatex = 1*ones(1,nNeuron);                       
SynDecayex = 1 - SynRatex*dt;
v = c;                                              % reset voltage;
s = zeros(1,nNeuron);                               % spike occurance
u = zeros(1,nNeuron);                               % revocery variable
Istate = 0;                                         %state variable used for synaptic currents decay
I = Is ;
tindex = 1;                                         %time pointer for output arrays
                                                    
Vout = zeros(length(time),nNeuron) ;
Init=zeros(nNeuron,1);
firings=[zeros(sum(Init),1),find(Init==1)];
for t=time
    I = Is + CurrentInput(tindex,:)*CurrentStrength ;
    v = v + dt * (0.04*v.^2+5*v+140-u+I);
    u = u + dt * a .* (b.*v-u);
    Vout(tindex,:) = v ;
    Iout(tindex,:) = I ;
    s = zeros(1,nNeuron) ;
    fired=find(v>=30);                                  % indices of cells that spiked
    if ~isempty(fired)
        firings=[firings; t+0*fired', fired'];
        v(fired) = c(fired) ;
        u(fired) = u(fired)+d(fired) ;
        
    end
tindex = tindex+1;                                      %update time index
end
indquand=find(firings(:,2)==1);
indquand=indquand(depmes:end);          %suppression of the first spikes
tempsfire1=firings(indquand,1);
interval1=diff(tempsfire1);
P=mean(interval1);                 %tonic response interspike
erreurP=std(interval1);
tindex = 1;                             %pointer for output arrays
Istate = 0;                             %state variable used for synaptic currents decay
I = Is ;
v = c;                                  % reset voltage;
s = zeros(1,nNeuron);                   % spike occurrence
u = zeros(1,nNeuron);                   % recovery variable

%init plotting variables
Vout = zeros(length(time),nNeuron) ;
Init=zeros(nNeuron,1);
firings=[zeros(sum(Init),1),find(Init==1)];
for t=time
    format short
    tempsavant=t-P;
    indqui=find(firings(:,1)<=tempsavant & firings(:,1)>(tempsavant-dt));
    qui=firings(indqui,2);
    s(qui)=1;
  Istate =  s*SynStrength+Istate .* SynDecayex.*(Istate>0);
  I = Istate +  Is + CurrentInput(tindex,:)*CurrentStrength ;
  v = v + dt * (0.04*v.^2+5*v+140-u+I);
    u = u + dt * a .* (b.*v-u);
    Vout(tindex,:) = v ;
    Iout(tindex,:) = I ;
    s = zeros(1,nNeuron) ;
    fired=find(v>=30);                      % indices of cells that spiked
    if ~isempty(fired)
        firings=[firings; t+0*fired', fired'];
        v(fired) = c(fired) ;
        u(fired) = u(fired)+d(fired) ;
    end;

    tindex = tindex+1;
end 
indquand2=find(firings(:,2)==2);
indquand2=indquand2(depmes:end);
tempsfire2=firings(indquand2,1);
indi=find(Iout(:,2)>=Gmaxsynex);
indi=indi(depmes:end)*dt;
if length(tempsfire2)==length(indi)
        deltat=tempsfire2-indi;
        Delex=mean(deltat);               %mean time delay between pre and post synaptic spikes 
        errdelex=std(deltat);
else
    disp('size problem for deltat')
    Delex=0;
end


