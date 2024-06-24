%%
close all;
clear;

%% 1 in 2
% Preberi sliki
yosemite_RGB = imread('yosemite_meadows.jpg');
pencils_RGB = imread('coloured_pencils.jpg');

% Pretvori sliki v barvna prostora HSV in LAB
yosemite_HSV = rgb2hsv(yosemite_RGB);
pencils_HSV = rgb2hsv(pencils_RGB);
yosemite_LAB = rgb2lab(yosemite_RGB);
pencils_LAB = rgb2lab(pencils_RGB);
%Sliki sta pretvorjeni iz RGB v HSV (Hue, Saturation, Value). HSV je barvni model, ki ločuje barvni ton (hue), nasičenost (saturation) in vrednost svetlosti (value).
%Sliki sta pretvorjeni tudi v LAB barvni prostor. LAB model ločuje svetlost (L) od barvnih informacij (a in b).

% Shrani slike HSV in LAB formatih v polja celic za lažjo iteracijo
images_HSV = {yosemite_HSV, pencils_HSV};
images_LAB = {yosemite_LAB, pencils_LAB};
image_titles = {'yosemite', 'pencils'};

% Oznake komponent za HSV in LAB
component_labels_HSV = ["Hue", "Saturation", "Value"];
component_labels_LAB = ["Lightness", "a", "b"];

% Prikaz posameznih komponent HSV in LAB za vsako sliko
for i = 1:2 
    for j = 1:3
        figure;
        imagesc(images_HSV{i}(:, :, j));
        colormap(gray);
        title([image_titles{i}, ' - ', component_labels_HSV(j)]);

        figure;
        imagesc(images_LAB{i}(:, :, j));
        colormap(gray);
        title([image_titles{i}, ' - ', component_labels_LAB(j)]);
    end
end

% Spremeni in prikaži slike z konstantnimi komponentami v HSV in LAB
for i = 1:2
    % HSV s konstantnimi H, S in V
    for j = 1:3
        temp_img = images_HSV{i};
        temp_img(:, :, j) = mean(temp_img(:, :, j), 'all');
        figure; imshow(hsv2rgb(temp_img));
    end

    % LAB s konstantnimi L, a in b
    for j = 1:3
        temp_img = images_LAB{i};
        temp_img(:, :, j) = mean(temp_img(:, :, j), 'all');
        figure; imshow(lab2rgb(temp_img));
    end
end


%% 3 HSV

% Preberi video in izvleči zadnji okvir
video = VideoReader("ladja.avi");
frames = read(video, [1 Inf]);
lastFrame = frames(:, :, :, video.NumFrames);

% Prikaz zadnjega okvira in omogoča uporabniku, da izbere 4 točke
figure; imshow(lastFrame);
title('Kliknite na 4 točke, ki označujejo ROI');

% Zbiranje točk z ginput
[x, y] = ginput(4); % Zberite 4 točke
x = round(x); % Zaokrožite x koordinate
y = round(y); % Zaokrožite y koordinate

% Ustvari logično masko za izbrano območje interesa (ROI)
roiMask = poly2mask(x, y, size(lastFrame, 1), size(lastFrame, 2));

% Uporabi masko na zadnjem okviru za izolacijo ROI
roi = lastFrame;
roi(repmat(~roiMask, [1, 1, 3])) = 0;

% Pretvori ROI v HSV
roiHSV = rgb2hsv(roi);

% Izračunaj območje HSV za ROI
hRange = [min(roiHSV(:,:,1), [], 'all'), max(roiHSV(:,:,1), [], 'all')];
sRange = [min(roiHSV(:,:,2), [], 'all'), max(roiHSV(:,:,2), [], 'all')];
vRange = [min(roiHSV(:,:,3), [], 'all'), max(roiHSV(:,:,3), [], 'all')];

% Inicializiraj izhodno matriko videa
videoOut = zeros(video.Height, video.Width, video.NumFrames);

% Obdelaj vsak okvir
for i = 1:video.NumFrames
    frameHSV = rgb2hsv(frames(:, :, :, i));
    mask = frameHSV(:,:,1) >= hRange(1) & frameHSV(:,:,1) <= hRange(2) & ...
           frameHSV(:,:,2) >= sRange(1) & frameHSV(:,:,2) <= sRange(2) & ...
           frameHSV(:,:,3) >= vRange(1) & frameHSV(:,:,3) <= vRange(2);
    videoOut(:,:,i) = mask * 255;
end

% Predvajaj obdelan video
implay(videoOut);
% Zapišite grayVideo v AVI datoteko
v = VideoWriter('3nal.avi', 'Grayscale AVI');
open(v);
writeVideo(v, uint8(videoOut * 255)); % Pomnožite z 255, ker mat2gray vrne vrednosti med 0 in 1
close(v);

%% 4

% Izberite okvir za analizo
selectedFrameIndex = 63; % na primer 63. okvir
selectedFrame = frames(:, :, :, selectedFrameIndex);

% Prikaz izbranega okvira
figure; imshow(selectedFrame);
title('Kliknite na 4 točke, ki označujejo ROI');

% Zbiranje točk z ginput
[x, y] = ginput(4); % Zberite 4 točke
x = round(x); % Zaokrožite vrednosti
y = round(y);

% Ustvari logično masko, ki označuje izbrano ROI. Maska je enaka velikosti izbranega okvira.
mask = poly2mask(x, y, size(selectedFrame, 1), size(selectedFrame, 2));

% Uporabite masko na izbranem okviru
roi = selectedFrame;
% Nastavi vse piksle zunaj ROI na 0, kar efektivno "izreže" območje zunaj ROI.
roi(repmat(~mask, [1, 1, 3])) = 0;

% Izračunajte 3D RGB histogram
numBins = 8;
histogram = zeros(numBins, numBins, numBins);
% V dvojni zanki (for i = 1:size(roi, 1) ... for j = 1:size(roi, 2)) prešteje vrednosti barv v ROI in jih razporedi v ustrezne bine histograma.
for i = 1:size(roi, 1)
    for j = 1:size(roi, 2)
        r = floor(double(roi(i, j, 1)) / (256 / numBins)) + 1;
        g = floor(double(roi(i, j, 2)) / (256 / numBins)) + 1;
        b = floor(double(roi(i, j, 3)) / (256 / numBins)) + 1;
        histogram(r, g, b) = histogram(r, g, b) + 1;
    end
end

% Normalizirajte histogram
maxValue = max(histogram, [], 'all'); % Poišče najvišjo vrednost v histogramu.
normalizedHistogram = histogram / maxValue * 255; % Normalizira histogram, tako da ima najvišji bin vrednost 255.

% Inicializirajte izhodno matriko za sivinske slike
video4 = zeros(video.Height, video.Width, 1, video.NumFrames); % Ustvari prazen štiridimenzionalni niz, ki bo shranjeval sivinske slike.

% Obdelajte vsak okvir videa
%  V trojni zanki se za vsak piksel v vsakem okvirju izračuna, v katere bine histograma pade, in uporabi vrednost iz tega bina histograma kot sivinsko vrednost.
for i = 1:video.NumFrames
    for j = 1:video.Height
        for k = 1:video.Width
            r = floor(double(frames(j, k, 1, i)) / (256 / numBins)) + 1;
            g = floor(double(frames(j, k, 2, i)) / (256 / numBins)) + 1;
            b = floor(double(frames(j, k, 3, i)) / (256 / numBins)) + 1;

            % Dodelitev vrednosti iz histograma
            video4(j, k, 1, i) = normalizedHistogram(r, g, b);
        end
    end
end

% Predvajaj obdelan video
implay(video4);

% Zapišite video4 v AVI datoteko
v = VideoWriter('4nal.avi', 'Grayscale AVI');
open(v);
writeVideo(v, uint8(video4 * 255)); % Pomnožite z 255 za zapis v AVI
close(v);



%% 5

% Preberi video
video = VideoReader("ladja.avi");
frames = read(video, [1 Inf]);

% Zadnji okvir
zadnji_frame = frames(:, :, :, video.NumFrames);

% Prikaz zadnjega okvira in omogoča uporabniku, da izbere 4 točke
figure; imshow(zadnji_frame);
title('Kliknite na 4 točke, ki označujejo vogale ROI');
[x, y] = ginput(4); % Zberite 4 točke
x = round(x); % Zaokrožite x koordinate
y = round(y);

% Ustvari logično masko za izbrano območje interesa (ROI)
mask = poly2mask(x, y, size(zadnji_frame, 1), size(zadnji_frame, 2));

% Uporabi masko na zadnjem okviru za izolacijo ROI
vzorec_morja = zadnji_frame;
vzorec_morja(repmat(~mask, [1, 1, 3])) = 0;

% Nato nadaljujte z obdelavo, kot je bilo prvotno načrtovano
st_binov = 4;


histogram_vzorca = zeros(st_binov, st_binov, st_binov);
presek = zeros(st_binov, st_binov, st_binov);

for i = 1:size(vzorec_morja, 1)
    for j = 1:size(vzorec_morja, 2)
        r = floor(double(vzorec_morja(i, j, 1)) / (256 / st_binov)) + 1;
        g = floor(double(vzorec_morja(i, j, 2)) / (256 / st_binov)) + 1;
        b = floor(double(vzorec_morja(i, j, 3)) / (256 / st_binov)) + 1;
        histogram_vzorca(r, g, b) = histogram_vzorca(r, g, b) + 1;
    end
end

% Reset histogram for each 8x8 region
histogram_8x8_podrocja = zeros(st_binov, st_binov, st_binov);

% Obdelaj vsak okvir videa
for i = 1:size(frames, 4)
    for j = 1:8:video.Height - 7
        for k = 1:8:video.Width - 7
            % Trenutni okvir
            trenutni_frame = frames(:, :, :, i);
            
            % Izberi 8x8 področje
            podrocje8x8 = trenutni_frame(j:j+7, k:k+7, :);
            
            % Reset histogram for each 8x8 region
            histogram_8x8_podrocja(:) = 0;
            
            % Izračunaj histogram za 8x8 področje
            for n = 1:8
                for m = 1:8
                    r = floor(double(podrocje8x8(n, m, 1)) / (256 / st_binov)) + 1;
                    g = floor(double(podrocje8x8(n, m, 2)) / (256 / st_binov)) + 1;
                    b = floor(double(podrocje8x8(n, m, 3)) / (256 / st_binov)) + 1;
                    histogram_8x8_podrocja(r, g, b) = histogram_8x8_podrocja(r, g, b) + 1;
                end
            end
            
            % Izračunaj presek histogramov z vzorcem morja
            for a = 1:st_binov
                for b = 1:st_binov
                    for c = 1:st_binov
                        presek(a, b, c) = min(histogram_8x8_podrocja(a, b, c), histogram_vzorca(a, b, c));
                    end
                end
            end
            
            % Vsota preseka za to področje
            vsote((j-1)/8+1, (k-1)/8+1, i) = sum(presek(:));
        end
    end
end

% Normalizacija in priprava izhodnega videa
max_vrednost = max(vsote(:));
norm_vsote = vsote / max_vrednost * 255;
norm_vsote = uint8(norm_vsote);
video2 = zeros(size(norm_vsote, 1), size(norm_vsote, 2), 3, size(norm_vsote, 3), 'uint8');

for i = 1:size(norm_vsote, 3)
    for j = 1:size(norm_vsote, 1)
        for k = 1:size(norm_vsote, 2)
            video2(j, k, 1, i) = norm_vsote(j, k, i);
            video2(j, k, 2, i) = norm_vsote(j, k, i);
            video2(j, k, 3, i) = norm_vsote(j, k, i);
        end
    end
end

% Pretvorba RGB videa v grayscale
grayscale_video = zeros(size(video2, 1), size(video2, 2), 1, size(video2, 4), 'uint8');

for i = 1:size(video2, 4)
    % Uporabite enega od RGB kanalov ali povprečje vseh treh kanalov
    grayscale_video(:, :, 1, i) = 0.299 * video2(:, :, 1, i) + 0.587 * video2(:, :, 2, i) + 0.114 * video2(:, :, 3, i);
end

% Shranjevanje obdelanega videa
v = VideoWriter('5nal.avi', 'Grayscale AVI');
open(v);
writeVideo(v, grayscale_video);
close(v);



