 %Biosystems. 2009 Jul;97(1):35-43.
%Horcholle-Bossavit G, Quenet B.
%Neural model of frog ventilatory rhythmogenesis.


function modulation(lfire)
global beta gamma maxAc Em0 cumul Ac Eml1 delta

       Japp=Eml1(end); 
       cumul=(1-Heavy(cumul-maxAc))*(cumul+lfire);
       terme1=(1-Heavy(cumul-maxAc))*beta*Japp;
       terme2=Heavy(cumul-maxAc)*Em0;
       terme3=gamma*(rand-0.5);  
       
       Japp=terme1+terme2+terme3;
       Japp=(Japp<=0)*(Em0) +Japp*(Japp>0);                        
       Ac=[Ac, cumul];      
       Eml1=[Eml1,Japp];
       if delta>0
            maxAc=maxAc+Heavy(cumul-maxAc)*(delta*(rand-0.5));
       end