function mattoub=Boucladj(matdep, Rdep, nboucles, premf, valf, premb, valb)


%cette fonction permet de construire des matrices de nboucles, la première
%boucle étant définie par matdep, avec les vecteur premf et premb qui sont
%des vecteurs d'une ligne deux colonnes, les premiers maillons de l'association des boucle, maillon forward et
%maillon backward, valf et valb étant respectivmeent les valeurs de ces
%liens. Rdep est un vecteur colonne

matrec=matdep(2:end,2:end);
matinter=matdep;
taille=size(matrec,1);
compr=Rdep;
for i=1:nboucles-1
    matinter=[matinter,zeros(size(matinter,1),size(matrec,2)); zeros(size(matrec,1),size(matinter,2)),matrec];
    compr=[compr;Rdep(2:end)];
end
mattoub=matinter;
mattoub(premf(1),premf(2))=valf;
mattoub(premb(1),premb(2))=valb;
for i=1:nboucles-2
    mattoub(premf(1)+i*taille, premf(2)+i*taille)=valf;
    mattoub(premb(1)+i*taille, premb(2)+i*taille)=valb;
end

mattoub=[mattoub,compr];
