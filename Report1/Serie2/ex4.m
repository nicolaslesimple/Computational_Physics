clear all;
close all;
clc;

% Question 1 : load the image as a red-green-blue (RGB) 3three-dimensional
% array with imread
img_gray = imread('BN.png');
% Question 1 : convert it to two-dimensional grayscale array with 
% rgb2gray and double matlab function
img_gray = double(rgb2gray(img_gray));
img_gray = img_gray - min(min(img_gray));
img_gray = img_gray/max(max(img_gray));
% Question 1 :  plot the image with imshow
figure(1)
imshow(img_gray);
title('Transmission electron microscopy image of an h-BN sample');
filename='./plot/img_gray.eps';
print(gcf,'-depsc',filename)

% Question 2 : Calculate the Fourier transform of the TEM image 
% using the fft2 and fftshift functions of Matlab.
img_fft2 = fft2(img_gray);
img_fftshift = fftshift(img_fft2);
% Question 2 : Visualize the result using imshow
figure(2)
imshow(img_fft2);
title('Fourier transform of TEM image of an h-BN sample');
filename='./plot/fft2.eps';
print(gcf,'-depsc',filename)
figure(3)
imshow(img_fftshift);
title('Fourier transform after shifting of TEM image of an h-BN sample');
filename='./plot/img_fftshift.eps';
print(gcf,'-depsc',filename)
% Question 2 : Make sure the resulting image clearly shows the hexagonal 
% pattern formed by the ensemble of Bragg peaks.
img_perfect = (abs(img_fftshift/3000).^2);
img_perfect(:,263) = 0; % we supress the white line that makes non sense
figure(4)
imshow(img_perfect)
title('Bragg Peaks of FT of TEM image of an h-BN sample');
filename='./plot/bragg_peak.eps';
print(gcf,'-depsc',filename)

% Question 4 : Design a Fourier filter that allows to distinguish 
% the honeycomb lattices of the two layers.
rgb_filter = create_filter();
figure(5)
imshow(rgb_filter);
title('Fourier filter');
filename='./plot/rgb_filter.eps';
print(gcf,'-depsc',filename)
% Question 4 : Apply the color filter to our image and plot it
color_coded_FT = (abs(img_fftshift/100.*rgb_filter).^2);
figure(6)
imshow(color_coded_FT);
title('Application of the color filter on the FT');
filename='./plot/filter_bragg.eps';
print(gcf,'-depsc',filename)

% Question 5 : Perform inverse Fourier transform (using first fftshift, and then ifft2 of MatLab) 
color_coded_FT_shift = fftshift(img_fftshift.*rgb_filter);
color_coded_IFT = 10*abs(ifft2(color_coded_FT_shift));
print(gcf,'-depsc',filename)
figure(7)
imshow(color_coded_IFT);
title('Inverse of the Fourier Transform after application of the color filter');
filename='./plot/output_grayscale.eps';
print(gcf,'-depsc',filename)


function rgb_filter = create_filter()
    % Function allowing the creation of the color filter.
    % This color filter will be apply on the Fourier Transform of th BN
    % image.
    % This filter will allows to distinguish the honeycomb lattices of the two layers.
    
    N = 525; % Number of pixels per dimension on the image
    pts_x_1 = [262,282,283,264,244,243]; % X coordinate of point of the first hexagone
    pts_y_1 = [240,251,274,286,275,253];  % Y coordinate of point of the first hexagone
    pts_x_2 = [270,285,278,256,241,248]; % X coordinate of point of the first hexagone
    pts_y_2 = [241,258,280,285,268,246];  % Y coordinate of point of the first hexagone
    rayon_1 = zeros(6, 1);
    angle_1 = zeros(6, 1);
    angle_2 = zeros(6, 1);
    rayon_2 = zeros(6, 1);
    for i = 1:6
        angle_1(i) = atan2(pts_y_1(i)-N/2 , pts_x_1(i)-N/2);
        angle_2(i) = atan2(pts_y_2(i)-N/2 , pts_x_2(i)-N/2);
        rayon_1(i) = sqrt((pts_y_1(i)-N/2)^2 + (pts_x_1(i)-N/2)^2);
        rayon_2(i) = sqrt((pts_y_2(i)-N/2)^2 + (pts_x_2(i)-N/2)^2);
        
    end
    % Cartesian coordinates 
    [x ,y] = meshgrid(1:N,1:N);
    % Polar coordinate
    phi = atan2(y-N/2 ,x-N/2); 
    hsv = zeros(N,N,3); 
    hsv (: ,: ,1) = 1; hsv (: ,: ,2) = 1; hsv(:,:,3)=0;
    for k=1:N % Iteration on the pixel of the X-axis
        for l=1:N % Iteration on the pixel of the Y-axis
            rayon_tmp = sqrt((N/2 - l)^2 + (N/2 - k)^2);
            for n=1:6 % Iteration on each point of the two hexagones
                % Fill up of the pixel corresponding to the first hexagone
                if (rayon_tmp < rayon_1(n) + 2.5)&&(rayon_tmp > rayon_1(n) - 2.5)&&(mod(phi(k,l)/3/pi,1) > mod(angle_1(n)/3/pi,1) - 0.017)&&(mod(phi(k,l)/3/pi , 1 ) < mod(angle_1(n)/3/pi , 1 ) + 0.017)
                    hsv(k,l,3) = 1;
                % Fill up of the pixel corresponding to the second hexagone
                elseif (rayon_tmp < rayon_2(n) + 2.5)&&(rayon_tmp > rayon_2(n) - 2.5)&&(mod(phi(k,l)/3/pi,1) > mod(angle_2(n)/3/pi,1) - 0.017)&&(mod(phi(k,l)/3/pi,1) < mod(angle_2(n)/3/pi,1) + 0.017)
                    hsv(k,l,1)=119/360; hsv(k,l,2)=0.66;hsv(k,l,3)=1;
                end
            end
        end
    end
    % Convert the HSV filter to RGB and display it
    rgb_filter = hsv2rgb(hsv);
end
