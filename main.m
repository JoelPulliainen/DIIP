clc
clear all
close all

coins = cell(1,12)
[mean_bias,mean_dark,norm_flat] = Imean();
x = linspace(72,83,12);
for i = x
    F = sprintf('_DSC17%d.JPG',i);
    measurement = imread(F);
    coins{i-71} = estim_coins(measurement,mean_bias,mean_dark,norm_flat)
end
% measurement = imread("_DSC1775.JPG");
% coins1 = estim_coins(measurement,mean_bias,mean_dark,norm_flat)

% % measurement = imread("_DSC1777.JPG");
% % coins2 = estim_coins(measurement,mean_bias,mean_dark,norm_flat)

% measurement = imread("_DSC1781.JPG");
% coins2 = estim_coins(measurement,mean_bias,mean_dark,norm_flat)

function [mean_bias,mean_dark,flat] = Imean()

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
    mean_dark = (dark1 + dark2 + dark3 + dark4 + dark5)./5;
    flat = flat1 + flat2 + flat3 + flat4 + flat5;
end