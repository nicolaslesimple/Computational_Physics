clear all;
close all;
clc;

tic
number_of_point = [8,16,32,64, 128, 256, 512,1024];
dft = zeros(8,1);
fft = zeros(8,1);

for i=1:8
    sample = exp(-linspace(-4,4,number_of_point(i)).^2);
    tic
    myfft(sample);
    fft(i,1) = toc;
    tic
    mydft(sample);
    dft(i,1)  = toc;
end

figure(1)
plot(number_of_point, dft, '-b', 'LineWidth', 0.1)
hold on;
plot(number_of_point, fft, '-r', 'LineWidth', 0.1)
title1 = title({'Time needed to compute myfft() and mydft() in function',...
                'of the number of point'});
set(title1,'FontName','Arial','FontSize',12);
leg1 = legend('DFT Time','FFT Time','Location','NorthWest');
set(leg1,'FontName','Arial','FontSize',10)
xlabel('Number of points [-]','FontName','Arial','FontSize',10);
ylabel('Time [s]','FontName','Arial','FontSize',10);
grid on;
hold off 
filename='./plot/time_comparison.eps';
print(gcf,'-depsc',filename)