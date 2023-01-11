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
x = linspace(72,83,12);
for i = x
    F = sprintf('_DSC17%d.JPG',i);
    measurement = imread(F);
    %imtool(measurement)
    coins{i-71} = estim_coins1(measurement,bias1,dark1,flat1)
end
% measurement = imread("_DSC1775.JPG");
% coins1 = estim_coins(measurement,mu_bias,mu_dark,norm_flat)
% 
% % measurement = imread("_DSC1777.JPG");
% % coins2 = estim_coins(measurement,mu_bias,mu_dark,norm_flat)

function coins = estim_coins1(measurement, mu_bias, mu_dark, norm_flat)
    figure("Name","refrence"),imshow(measurement);
    
    I_calibrated = calibration(measurement,mu_bias,mu_dark,norm_flat);

    figure("Name","Calibrated"),imshow(I_calibrated);

    [pixel_size_x,pixel_size_y] = size_calc(norm_flat)
    I_calibrated = rgb2gray(I_calibrated);
    threshold = graythresh(I_calibrated);
    BW = imbinarize(I_calibrated,threshold);
    %BW = imfill(BW,"holes");
    BW = imcomplement(BW); 
    figure("name","segmented image")
    imshow(BW)



    
    stats = regionprops('table',BW,'Centroid','Area','Perimeter','MajorAxisLength','MinorAxisLength','EquivDiameter');
    
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
        elseif 24.25 > diameters(i) & diameters(i) > 22.8
            pred(2) = pred(2)+1;
        elseif 25 > diameters(i) & diameters(i) > 24.25
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

function calibrated_image = calibration(measurement,bias,dark,flat)

    flat = flat-bias-dark;

    flat1 = flat(:,:,1);
    flat2 = flat(:,:,2);
    flat3 = flat(:,:,3);
    
    flat1 = flat1(:);
    flat1 = normalize(flat1,'range',[0.5 1.5]);
    flat1 = reshape(flat1, 2592, 3872);
    flat(:,:,1) = uint8(flat1);

    flat2 = flat2(:);
    flat2 = normalize(flat2,'range',[0.5 1.5]);
    flat2 = reshape(flat2, 2592, 3872);
    flat(:,:,2) = uint8(flat2);

    flat3 = flat3(:);
    flat3 = normalize(flat3,'range',[0.5 1.5]);
    flat3 = reshape(flat3, 2592, 3872);
    flat(:,:,3) = uint8(flat3);

    calibrated_image = (measurement - bias - dark) ./ flat;

end