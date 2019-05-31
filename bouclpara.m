function matpara=bouclpara(Spre, Rpre, npara)


%cette fonction permet de construire des matrices de nboucles, la première
%boucle étant définie par matdep, avec les vecteur premf et premb qui sont
%des vecteurs d'une ligne deux colonnes, les premiers maillons de l'association des boucle, maillon forward et
%maillon backward, valf et valb étant respectivmeent les valeurs de ces
%liens. Rdep est un vecteur colonne

matrep=Spre(1:end,1:end);
matinter=matrep;
taille=size(matrep,1);
compr=Rpre;
for i=1:npara-1
    matinter=[matinter,zeros(size(matinter,1),size(matrep,2)); zeros(size(matrep,1),size(matinter,2)),matrep];
    compr=[compr;Rpre(1:end)];
end
matpara=[matinter,compr];
