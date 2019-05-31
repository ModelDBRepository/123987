%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.

function Smod=matmodif(S,nombrerreur,facteur,germe)
                     
global memoirej nappelbontir nouvtir
Smod=S; 
memoirej=[];
nappelbontir=0;
conibmod=[];
for k=1:nombrerreur
        nouvtir=0;
        z=bontirageg(size(S,1)*size(S,1),germe);     
        conibmod=[conibmod,z];
        y=floor(z/size(S,1));
        x=rem(z,size(S,1));
        a=y+1-(x==0);
        b=x+(x==0)*size(S,1);
        Smod(a,b)=facteur*(Smod(a,b)==0)*(2*(sum(S(:,b))>=0)-1);
end