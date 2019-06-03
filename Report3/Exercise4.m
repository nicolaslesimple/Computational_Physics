%% Clear all
clear all;
close all;
clc;

%% Question 1 : Load the photo and apply an svd transformation
figure(1);
imshow('noether.jpg')
image = imread('noether.jpg');
% Transform all the pixels to double 
image_double = double(image);
% Apply the SVD transformation
singval=svd(image_double);

% Plot the sigular values
figure(2);
semilogy(singval,'LineWidth',1,'MarkerSize',3)
set(gca,'FontName','Arial','FontSize',15)
title1=title('Singular Values of the Picture');
set(title1,'FontName','Arial','FontSize',12)
xlabel('Index of the singular value','FontName','Arial','FontSize',16);
ylabel('Singular Values','FontName','Arial','FontSize',16);
grid on;
filename='./plot/singval_log.eps';
print(gcf,'-depsc',filename)

fprintf('The rank of the martix is : ')
disp(rank(image_double));
fprintf('Thus, the matrix has non sigular values.\n')


%% Question 2 
image = imread('noether.jpg');
image_double = double(image);

% Singular Decomposition
[U,S,V]=svd(image_double);

% Select 4, 20 and 100 singular value to represent the original image
[~,N]=size(S);
for i=[4 20 100 220]
    img_matrix(:,:,i)=U*[S(:,1:i) zeros(312,N-i)]*V';
end
img_matrix=uint8(img_matrix(:,:,:));

figure(3);
imshow(image);
hold on;
waitforbuttonpress % Click on the picture to see the one after
imshow(img_matrix(:,:,4))
title('N = 4')
filename='./plot/svd_4.eps';
print(gcf,'-depsc',filename)
waitforbuttonpress
imshow(img_matrix(:,:,20))
title('N = 20')
filename='./plot/svd_20.eps';
print(gcf,'-depsc',filename)
waitforbuttonpress
imshow(img_matrix(:,:,100))
title('N = 100')
filename='./plot/svd_100.eps';
print(gcf,'-depsc',filename)
waitforbuttonpress
imshow(img_matrix(:,:,220))
title('N = 220')
filename='./plot/svd_220.eps';
print(gcf,'-depsc',filename)

%% Additional Study : error in reconstruction in funciton of number of non-zeros columns 
image = imread('noether.jpg');
image_double = double(image);

% SVD compression with the Matlab function
[U,S,V]=svd(image_double);

[~,N]=size(S);
error = zeros(220,1);
for i=[1:1:220]
    img_matrix=U*[S(:,1:i) zeros(312,N-i)]*V';
    error(i)=norm(image_double-img_matrix);
end

figure(4);
% Plot error corresponding to the sigular values
semilogy([1:1:220],error,'LineWidth',1,'MarkerSize',3)
set(gca,'FontName','Arial','FontSize',15)
title1=title('Lower rank approximations');
set(title1,'FontName','Arial','FontSize',12)
xlabel('Number of non-zeros columns','FontName','Arial','FontSize',16);
ylabel('Error with respect to the initial image','FontName','Arial','FontSize',16);
grid on;
filename='./plot/error_log.eps';
print(gcf,'-depsc',filename)