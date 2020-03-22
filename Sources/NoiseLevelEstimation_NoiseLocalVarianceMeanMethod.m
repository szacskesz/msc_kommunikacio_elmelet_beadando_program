% Minden pixelre veszi azt a 3*3 h�rmas m�trixot, aminek az adott pixel a
% k�z�ppontja.
% Ennek kisz�m�tja az �tlag�t �s a sz�r�s�t (localMean, localVar)
%
% V�g�l �gy az �sszes pixelre megkapott (lok�lis) sz�r�st �tlagolja �s �gy
% megkapjuk hogy mekkora �tlagos sz�r�sa a k�rnyezet�kben lev� pixelekhez
% k�pest.
%
% Az algoritmus az addit�v gauss-os feh�rzaj meg�llap�t�s�ra lehet
% felhaszn�lni, hiszen az fix sz�r�ssal van jelen minden pixelen.
%
% Megjegyz�s: a k�p sz�leinel, ahol a 3*3-as m�trix m�r kil�gna a k�pb�l, a
% kil�g� r�szeket 0-val helyettes�ti (sz�r�s �gy v�ltozatlan marad).
function [nLevels] = NoiseLevelEstimation_NoiseLocalVarianceMeanMethod(img, patchSize)
    nLevels = zeros(1, size(img,3));
    for i = 1:size(img,3)
        channel = img(:,:,i);
        g = im2double(channel);
        
        localMean = filter2(ones([patchSize, patchSize]), g) / prod([patchSize, patchSize]);
        localVar = filter2(ones([patchSize, patchSize]), g.^2) / prod([patchSize, patchSize]) - localMean.^2;
        nLevels(i) = mean2(localVar);
    end
end
