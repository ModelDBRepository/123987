%Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.

%computation of a smoothed histogram

function [xref,vtot]=smoothist(signal, binmax)


[vref,xref]=hist(signal,binmax);
vtot=vref;
for i=binmax-1:-1:10
    x=min(xref):(max(xref)-min(xref))/i:max(xref);
    v=hist(signal,x);
    vinterp=interp1(x,v,xref);
    vtot=vtot+vinterp;
end
vtot=vtot/(binmax-10);
    