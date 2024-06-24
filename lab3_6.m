% Zapre vse odprte figure in očisti delovno okolje
close all
clear

    % Pripravi delovno okolje in prebere video.
    % Izračuna in prikaže sredinski okvir ter omogoči izbiro ROI.
    % Ustvari masko za ROI in jo uporabi na izbranem okviru.
    % Za vsak okvir videa izračuna 3D RGB histogram in določi optimalni prag za upragovljanje.
    % Upragovi vsak okvir glede na izračunani prag.
    % Uporabi algoritem poplavljanja (imfill) za zapiranje črnih območij.
    % Prikaže in shrani obdelan video.


    % Za vsak okvir videa se izračuna 3D RGB histogram. Histogram beleži, kolikokrat se posamezne kombinacije rdeče (R), zelene (G) in modre (B) barve pojavijo v okviru.
    % Histogram je razdeljen na "bine" (numBins), ki določajo natančnost barvnega razpona, ki se beleži.
    % Histogram se normalizira tako, da je največja vrednost v histogramu enaka 255. To je potrebno za uskladitev vrednosti barv v razponu, ki se uporablja v digitalni obdelavi slik.
    % Za vsak možen prag (od 0 do 255) se preveri, koliko pikslov bi bilo "napačno označenih" (nepravilno klasificiranih) z uporabo tega praga.
    % Optimalni prag je tisti, pri katerem je število napačno označenih pikslov najmanjše.
    % Ko je optimalni prag določen, se za vsak okvir izvede upragovljanje. To pomeni, da se vsak piksel v okviru primerja s pragom.
    % Če je sivinska vrednost piksla višja ali enaka kot prag, se piksel v izhodnem okviru nastavi na belo (255). Če je nižja, se nastavi na črno (0).



% Preberi video datoteko in izvleči vse okvire
video = VideoReader("ladja.avi"); % Ustvari objekt za branje videa
frames = read(video, [1 Inf]); % Prebere vse okvire videa

% Izberite sredinski okvir videa za analizo
middleFrameIndex = round(video.NumFrames / 2); % Izračuna indeks sredinskega okvira
middleFrame = frames(:, :, :, middleFrameIndex); % Pridobi sredinski okvir

% Prikaz sredinskega okvira in omogoča uporabniku, da izbere 4 točke
figure; imshow(middleFrame); % Prikaže sredinski okvir
title('Kliknite na 4 točke, ki označujejo ROI'); % Dodaj naslov

% Zbira točke od uporabnika in zaokroži vrednosti
[x, y] = ginput(4); % Zbira 4 točke od uporabnika
x = round(x); % Zaokroži x koordinate
y = round(y); % Zaokroži y koordinate

% Ročno ustvari masko za izbrano območje interesa (ROI)
mask = zeros(size(middleFrame, 1), size(middleFrame, 2));
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
mask(ymin:ymax, xmin:xmax) = 1;
% mask(ymin:ymax, xmin:xmax) = 1; nastavi vrednosti znotraj izbranega pravokotnega območja na 1, kar ustvari masko.

% Pripravi ROI iz sredinskega okvira
roi = middleFrame; % Kopira sredinski okvir
% Nastavi vse piksle zunaj ROI na črno
roi(repmat(~mask, [1, 1, 3])) = 0; 

% Pripravi prostor za histogram in izhodno matriko
numBins = 8; % Število intervalov za histogram
histogram = zeros(numBins, numBins, numBins); % Inicializira histogram
video4 = zeros(video.Height, video.Width, 1, video.NumFrames); % Inicializira izhodno matriko

% Obdelava vsakega okvira videa posebej
% Ta zanka obdela vsak okvir videa, izračuna histogram in upragovi vsak okvir glede na izračunani prag.
for frameIndex = 1:video.NumFrames
    currentFrame = frames(:, :, :, frameIndex); % Trenutni okvir za obdelavo

    % Izračun histograma za trenutni okvir
    histogram = zeros(numBins, numBins, numBins); % Resetira histogram
    for i = 1:size(currentFrame, 1)
        for j = 1:size(currentFrame, 2)
            % Pridobi barvne komponente in jih razporedi v histogram
            r = floor(double(currentFrame(i, j, 1)) / (256 / numBins)) + 1;
            g = floor(double(currentFrame(i, j, 2)) / (256 / numBins)) + 1;
            b = floor(double(currentFrame(i, j, 3)) / (256 / numBins)) + 1;
            histogram(r, g, b) = histogram(r, g, b) + 1;
        end
    end

    % Normalizacija histograma
    maxValue = max(histogram, [], 'all'); % Najde največjo vrednost v histogramu
    normalizedHistogram = histogram / maxValue * 255; % Normalizira histogram

    % Določitev optimalnega praga za trenutni okvir
    optimalThreshold = 0;
    minWrongPixels = inf;
    for threshold = 0:255
        wrongPixels = 0; % Štetje napačno označenih pikslov
        % Preverjanje vsakega piksla glede na prag
        for i = 1:size(currentFrame, 1)
            for j = 1:size(currentFrame, 2)
                % Pridobi barvne komponente in izračuna sivinsko vrednost
                r = floor(double(currentFrame(i, j, 1)) / (256 / numBins)) + 1;
                g = floor(double(currentFrame(i, j, 2)) / (256 / numBins)) + 1;
                b = floor(double(currentFrame(i, j, 3)) / (256 / numBins)) + 1;
                grayValue = normalizedHistogram(r, g, b);

                % Preveri, ali je piksel pravilno označen glede na masko
                if mask(i, j) % Piksel znotraj maske
                    if grayValue < threshold
                        wrongPixels = wrongPixels + 1;
                    end
                else % Piksel izven maske
                    if grayValue >= threshold
                        wrongPixels = wrongPixels + 1;
                    end
                end
            end
        end
        % Preveri, če je trenutni prag boljši
        if wrongPixels < minWrongPixels
            minWrongPixels = wrongPixels;
            optimalThreshold = threshold;
        end
    end

    % Upragovljanje trenutnega okvira
    for j = 1:video.Height
        for k = 1:video.Width
            % Pridobi barvne komponente in izračuna sivinsko vrednost
            r = floor(double(currentFrame(j, k, 1)) / (256 / numBins)) + 1;
            g = floor(double(currentFrame(j, k, 2)) / (256 / numBins)) + 1;
            b = floor(double(currentFrame(j, k, 3)) / (256 / numBins)) + 1;
            grayValue = normalizedHistogram(r, g, b);

            % Nastavi piksle glede na prag
            if grayValue >= optimalThreshold
                video4(j, k, 1, frameIndex) = 255; % Bela barva
            else
                video4(j, k, 1, frameIndex) = 0; % Črna barva
            end
        end
    end
end

% Obdelava vsakega okvira videa s funkcijo imfill
for frameIndex = 1:video.NumFrames
    % Uporabi imfill za "poplavljanje" črnih območij v belih
    video4(:, :, 1, frameIndex) = imfill(video4(:, :, 1, frameIndex), 'holes');
end

% Predvajaj obdelan video
implay(video4);

% Zapiši obdelan video v AVI datoteko
v = VideoWriter('6nal.avi', 'Grayscale AVI'); % Ustvari pisalnik videa
open(v); % Odpri datoteko za pisanje
writeVideo(v, uint8(video4 * 255)); % Zapiši video in ga pretvori v 8-bitno obliko
close(v); % Zapri datoteko
