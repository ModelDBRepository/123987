%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.

function j=bontirageg(N,germe)

global memoirej nappelbontir nouvtir

if nargin==2
    rand('seed',germe);
end
nappelbontir=nappelbontir-nouvtir+1;

j=zeros(0);
if nappelbontir<N+1
    jessai=ceil(N*rand);
    if sum(memoirej==jessai)==0
        j=jessai;
        memoirej=[memoirej,j];
    else
           nouvtir=1;
           j=bontirageg(N);   
    end
end
    