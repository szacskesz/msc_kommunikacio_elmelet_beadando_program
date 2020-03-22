% Minden pixelre veszi azt a 3*3 hármas mátrixot, aminek az adott pixel a
% középpontja.
% Ennek kiszámítja az átlagát és a szórását (localMean, localVar)
%
% Végül így az összes pixelre megkapott (lokális) szórást átlagolja és így
% megkapjuk hogy mekkora átlagos szórása a környezetükben levõ pixelekhez
% képest.
%
% Az algoritmus az additív gauss-os fehérzaj megállapítására lehet
% felhasználni, hiszen az fix szórással van jelen minden pixelen.
%
% Megjegyzés: a kép széleinel, ahol a 3*3-as mátrix már kilógna a képbõl, a
% kilógó részeket 0-val helyettesíti (szórás így változatlan marad).
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
