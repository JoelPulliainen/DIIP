clc
clear all
close all


bias1 = imread("_DSC1722.JPG");
bias2 = imread("_DSC1723.JPG");
bias3 = imread("_DSC1724.JPG");
bias4 = imread("_DSC1725.JPG");
bias5 = imread("_DSC1726.JPG");
dark1 = imread("_DSC1762.JPG");
dark2 = imread("_DSC1763.JPG");
dark3 = imread("_DSC1764.JPG");
dark4 = imread("_DSC1765.JPG");
dark5 = imread("_DSC1766.JPG");
flat1 = imread("_DSC1767.JPG");
flat2 = imread("_DSC1768.JPG");
flat3 = imread("_DSC1769.JPG");
flat4 = imread("_DSC1770.JPG");
flat5 = imread("_DSC1771.JPG");


coins = cell(1,12)
[mean_bias,mean_dark,norm_flat] = Imean();
x = linspace(72,83,12);
for i = x
    F = sprintf('_DSC17%d.JPG',i);
    measurement = imread(F);
    %imtool(measurement)
    coins{i-71} = estim_coins1(measurement,mean_bias,mean_dark,norm_flat)
end
% measurement = imread("_DSC1775.JPG");
% coins1 = estim_coins1(measurement,mean_bias,mean_dark,norm_flat)
% 
% % measurement = imread("_DSC1777.JPG");
% % coins2 = estim_coins1(measurement,mean_bias,mean_dark,norm_flat)

% measurement = imread("_DSC1781.JPG");
% coins2 = estim_coins1(measurement,mean_bias,mean_dark,norm_flat)

function coins = estim_coins1(measurement, mu_bias, mu_dark, norm_flat)
    figure("Name","refrence"),imshow(measurement);
    
    I_calibrated = calibration(measurement,mu_bias,mu_dark,norm_flat);
    diff = I_calibrated-measurement;
    figure("Name","Calibrated"),imshow(I_calibrated);

    I_calibrated = db_calibration(I_calibrated);
    [pixel_size_x,pixel_size_y] = size_calc(norm_flat)

    figure("name","segmens32dds"); imshow(I_calibrated);
    

    %I_calibrated = imadjust(I_calibrated,[.2 .2 .2; .5 .6 .8],[]);
    
    I_calibrated = rgb2gray(I_calibrated);
    
    figure("name","segmensd76ds"); imshow(I_calibrated);


    threshold = graythresh(I_calibrated);
    BW = imbinarize(I_calibrated,threshold);
    %BW = imfill(BW,"holes");
    BW = imcomplement(BW); 
    figure("name","segmented image")
    imshow(BW)
    
   
    % Create average filter
    average_filter = fspecial('average',12);
    
    % Filter for sp noise 1 and plot
    BW= imfilter(BW,average_filter);
    BW = imfill(BW,"holes");

    %stats = regionprops('table',BW,'Centroid','Area','Perimeter','MajorAxisLength','MinorAxisLength','EquivDiameter');
    stats = regionprops('table',BW,'all');
    centers = stats.Centroid;
    
    [largest,Idx] = max(stats.Area);
    
    areas = stats.Area;
    peri = stats.Perimeter;
    eqdia = stats.EquivDiameter;
    ans = find(areas>50000);
    centers1 = centers;
    centers = centers(ans,:);
    areas = areas(ans);
    peri = peri(ans);
    eqdia = eqdia(ans);
    areas = areas;
    
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    diameters = diameters(ans);
    
    dia = diameters./pixel_size_x;
    
    figure("name","segmented image with centers marked")
    imshow(BW)
    hold on
    plot(centers(:,1),centers(:,2),'bsquare')
    plot(centers1(:,1),centers1(:,2),'r*')
    hold off
    
    coins = coin_counter(dia);

end
function pred = coin_counter(diameters)
    pred = zeros(1,6);
    for i = 1:length(diameters)
        if diameters(i) >= 25
            pred(1) = pred(1)+1;
        elseif 24.1 > diameters(i) & diameters(i) > 22.8
            pred(2) = pred(2)+1;
        elseif 25 > diameters(i) & diameters(i) > 24.1
            pred(3) = pred(3)+1;
        elseif 22.8 > diameters(i) & diameters(i) > 21.8
            pred(4) = pred(4)+1;
        elseif 20.5 > diameters(i) & diameters(i) > 19
            pred(5) = pred(5)+1;
        else 
            pred(6) = pred(6)+1;
        end
    end
end

% 
function i_cal_no_chec = db_calibration(image)
    mu1 = mean2(image(:,:,1))+20;
    mu2 = mean2(image(:,:,2))+25;
    mu3 = mean2(image(:,:,3))+35;
    for i = 1:1700% 1700
        for j = 2249:size(image,2) % 2249
            if image(i,j,1) < 120 && image(i,j,2) < 120 && image(i,j,3) < 120
                diff_rg = abs(image(i,j,1) - image(i,j,2));
                diff_rb = abs(image(i,j,1) - image(i,j,3));
                diff_gb = abs(image(i,j,2) - image(i,j,3));
                if diff_rg < 20 && diff_rb < 40 && diff_gb < 20
                    image(i,j,1) = mu1;
                    image(i,j,2) = mu2;
                    image(i,j,3) = mu3;
                end
            end
        end
    end
    %imtool(image)
    i_cal_no_chec = image;
    figure("name","Checkerboard")
    imshow(i_cal_no_chec)
end

function i_calibrated = remove_shadows(image)

    for i = 1:size(image,1)
        for j = 1:size(image,2) % 2249
            if image(i,j,1) < 10 && image(i,j,2) < 10 && image(i,j,3) < 10
                diff_rg = abs(image(i,j,1) - image(i,j,2));
                diff_rb = abs(image(i,j,1) - image(i,j,3));
                diff_gb = abs(image(i,j,2) - image(i,j,3));
                if diff_rg < 20 && diff_rb < 40 && diff_gb < 20
                    image(i,j,1) = 255;
                    image(i,j,2) = 200;
                    image(i,j,3) = 150;
                end
            end
        end
    end
    %imtool(image)
    i_calibrated = image;
    figure("name","Checkerboard")
    imshow(i_calibrated)
end

function [pixel_size_x,pixel_size_y] = size_calc(image)
    %Calculate the pixel to mm ratio
    [imagePoints,boardSize] = detectCheckerboardPoints(image);
    dist_y = sqrt(sum((imagePoints(2,:)-imagePoints(1,:)).^2));
    dist_x = sqrt(sum((imagePoints(7,:)-imagePoints(1,:)).^2));
    pixel_size_y = dist_y/12.5;
    pixel_size_x = dist_x/12.5;
    figure("name","Checkerboard")
    imshow(image)
    hold on
    plot(imagePoints(:,1),imagePoints(:,2),'ro')
    hold off

end

function [mean_bias,mean_dark,norm_flat] = Imean()

    bias1 = imread("_DSC1722.JPG");
    bias2 = imread("_DSC1723.JPG");
    bias3 = imread("_DSC1724.JPG");
    bias4 = imread("_DSC1725.JPG");
    bias5 = imread("_DSC1726.JPG");
    dark1 = imread("_DSC1762.JPG");
    dark2 = imread("_DSC1763.JPG");
    dark3 = imread("_DSC1764.JPG");
    dark4 = imread("_DSC1765.JPG");
    dark5 = imread("_DSC1766.JPG");
    flat1 = imread("_DSC1767.JPG");
    flat2 = imread("_DSC1768.JPG");
    flat3 = imread("_DSC1769.JPG");
    flat4 = imread("_DSC1770.JPG");
    flat5 = imread("_DSC1771.JPG");

    mean_bias = (bias1 + bias2 + bias3 + bias4 + bias5)./5;
    mean_dark = ((dark1-mean_bias) + (dark2-mean_bias) + (dark3-mean_bias) + (dark4-mean_bias) + (dark5-mean_bias))./5;
    norm_flat = ((flat1-mean_dark-mean_bias) + (flat2-mean_dark-mean_bias) + (flat3-mean_dark-mean_bias) + (flat4-mean_dark-mean_bias) + (flat5-mean_dark-mean_bias));
end

function calibrated_image = calibration(measurement,bias,dark,flat)

    flatfield = uint8(double(flat)./double(max(flat(:))));
    flat = double(flat)./double(max(flat(:)));
    flat1 = flat(:,:,1);
    mask = imbinarize(flatfield);
    mask1 = mask(:,:,1);
    calibrated_image = double(measurement - bias  - dark ) ./ flat;
    calibrated_image(~mask) = measurement(~mask);
    calibrated_image = uint8(calibrated_image);
    
    figure("name","segmented image poop"),imshow(calibrated_image)
    sad = 23


    %calibrated_image = (measurement - bias - dark) ./ flat;

end