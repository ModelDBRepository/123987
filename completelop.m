function matoutou=completel(matL,vectv,vecth,matB,op)

%le module L n'est en contact avec B que par l'excitateur L1
if op==1        %option d'excitation de L1
    matoutou=[matL,[vecth;zeros(size(matL,1)-1,size(matB,2))];[vectv,zeros(size(matB,1),size(matL,2)-1)],matB];
else            %option d'inhibition depuis L2
    matoutou=[matL,[vecth;zeros(size(matL,1)-1,size(matB,2))];[zeros(length(vectv),2),vectv,zeros(size(matB,1),size(matL,2)-3)],matB]
end