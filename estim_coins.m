function coins = estim_coins(measurement, bias, dark, flat)
% Detect and Count coins in a image
% Using Intensity calibration, geometrical
% calibration, removing checkerboard from the image, grayscaling, and
% filtering and filling holes of a image.
% Input is Measurement image which contains coins to be counted
% Mean bias image, mean dark image and flatfield image.
% Output is s a six-element vector with the first element
% corresponding to the number of 5 cent-coins in the image and the last element corresponding to
% the number of 2 euro-coins in the image
% Additional figures for visualising the process are included as comments

    %figure("Name","refrence"),imshow(measurement);
    
    % Calibrate image intesity
    I_calibrated = calibration(measurement,bias,dark,flat);

    % Calibrate image pixel size using checkerboard
    [pixel_size_x,pixel_size_y] = geo_calibration(flat);

    % Remove checkerboard from the image
    I_calibrated = remove_checkerboard(I_calibrated);

    % Make the image gray
    I_calibrated = rgb2gray(I_calibrated);

    % Plot gray image
    %figure("name","Gray image"); imshow(I_calibrated);

    % Threshold the image, binarize and imcomplement, then plot
    threshold = graythresh(I_calibrated);
    BW = imbinarize(I_calibrated,threshold);
    BW = imcomplement(BW); 

    % figure("name","segmented image")
    % imshow(BW)
   
    % Create average filter
    average_filter = fspecial('average',12);
    
    % Filter image
    BW= imfilter(BW,average_filter);

    % Fill holes in the image
    BW = imfill(BW,"holes");

    % Calculate statistics of the image regions
    %stats = regionprops('table',BW,'Centroid','Area','Perimeter','MajorAxisLength','MinorAxisLength','EquivDiameter');
    stats = regionprops('table',BW,'all');
    centers = stats.Centroid;

    % Find regions that are big enough to be coins
    areas = stats.Area;
    idx = find(areas>50000);
    centers1 = centers;
    centers = centers(idx,:);

    
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    diameters = diameters(idx);
    
    % Calculate size in mm
    dia = diameters./pixel_size_x;
    
    % figure("name","segmented image with centers marked")
    % imshow(BW)
    % hold on
    % plot(centers(:,1),centers(:,2),'bsquare')
    % plot(centers1(:,1),centers1(:,2),'r*')
    % hold off
    
    coins = coin_counter(dia);

end
function pred = coin_counter(diameters)
    % Coint coins based on diameters, Takes input of diameters, outputs an
    % array of the coin counts

    pred = zeros(1,6);
    for i = 1:length(diameters)
        if diameters(i) >= 25
            pred(1) = pred(1)+1;
        elseif 24 > diameters(i) && diameters(i) > 22.8
            pred(2) = pred(2)+1;
        elseif 25 > diameters(i) && diameters(i) > 24
            pred(3) = pred(3)+1;
        elseif 22.8 > diameters(i) && diameters(i) > 21.8
            pred(4) = pred(4)+1;
        elseif 20.5 > diameters(i) && diameters(i) > 19
            pred(5) = pred(5)+1;
        else 
            pred(6) = pred(6)+1;
        end
    end
end


function image_wo_checkerboard = remove_checkerboard(image)
    % Remove checkerboard from the image, Takes input of the image, outputs
    % image with top right corner checkerboard removed

    % Create adjusted mean values for the replacement color
    mu1 = mean2(image(:,:,1))+20;
    mu2 = mean2(image(:,:,2))+25;
    mu3 = mean2(image(:,:,3))+35;

    % Scan over the top right corner
    for i = 1:1700
        for j = 2249:size(image,2)
            if image(i,j,1) < 80 && image(i,j,2) < 80 && image(i,j,3) < 80
                diff_rg = abs(image(i,j,1) - image(i,j,2));
                diff_rb = abs(image(i,j,1) - image(i,j,3));
                diff_gb = abs(image(i,j,2) - image(i,j,3));
                % Check if the RGB difference is small enough -> pixel is
                % black
                if diff_rg < 20 && diff_rb < 40 && diff_gb < 20
                    image(i,j,1) = mu1;
                    image(i,j,2) = mu2;
                    image(i,j,3) = mu3;
                end
            end
        end
    end
    
    image_wo_checkerboard = image;
    % figure("name","Image without checkerboard")
    % imshow(image_wo_checkerboard)
end

function [pixel_size_x,pixel_size_y] = geo_calibration(image)
    %Calculate the pixel to mm ratio with detectCheckerboardPoints. Takes
    %input image to calibrate, outputs pixel size as mm

    [imagePoints,boardSize] = detectCheckerboardPoints(image);
    dist_y = sqrt(sum((imagePoints(2,:)-imagePoints(1,:)).^2));
    dist_x = sqrt(sum((imagePoints(7,:)-imagePoints(1,:)).^2));
    pixel_size_y = dist_y/12.5;
    pixel_size_x = dist_x/12.5;

%     figure("name","Checkerboard points")
%     imshow(image)
%     hold on
%     plot(imagePoints(:,1),imagePoints(:,2),'ro')
%     hold off

end

function calibrated_image = calibration(measurement,bias,dark,flat)
    % Perform intesity calibration for the image using dark bias and
    % flatfield images

    dark = dark-bias;
    flat = flat-bias-dark;
    flatfield = uint8(double(flat)./double(max(flat(:))));
    flat = uint8(double(flat)./double(max(flat(:))));

    % Create mask to not remove the coins on top of the checkerboard
    mask = imbinarize(flatfield);

    calibrated_image = (measurement - bias  - dark ) ./ flat;
    calibrated_image(~mask) = measurement(~mask);
    calibrated_image = uint8(calibrated_image);
    
    % figure("name","Calibrated image"),imshow(calibrated_image)


end