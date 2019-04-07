clear all;
close all;
clc;

% Question 1 : Import the file in a two-column matrix and merge the real
% and imaginary parts into a single complex vector.
data = importdata('./FID.dat.txt');
data_complex = complex(data(:,1), data(:,2));

% Question 2 : Associate time to each data point, setting the initial 
% time to zero.
delta_t = 83.2*10^-6;
data_size = size(data_complex);
data_size = data_size(1,1);
time_list = (linspace(0,data_size,data_size)'*delta_t);

% Question 3: Calculate the discrete frequencies fn which span the 
% interval [-fc, fc] with fc = 1/2delat_t being the Nyquist frequency.
fc = (1/(2*delta_t));
fn = linspace(-fc,fc,data_size);

% Question 4 : Perform the discrete Fourier transform of the complex FID.
mydft_fid = mydft(data_complex);
fft_fid = fft(data_complex);

% Question 4 : Compare the result of the Fourier transform obtained 
% with the fft function of MatLab.
dif = abs((mydft_fid - fft_fid)./fft_fid);
assert(all(abs(mydft_fid - fft_fid)<1e-6));
fprintf('\tCheck up for mydft passed\n');

% Question 4 : Calculate the chemical shift ?n (in ppm) and plot the
% asked graphs 
f0 = 1067.93;
mu0 = 800.224*10^6;
dirac_n = ((fn - f0)*10^6)/(mu0);

% Question 4 : Calculus of the power spectrum using an existing tool
power_spectrum = abs((mydft_fid).^2/data_size);

% Question 6 : What is the ratio between the number of protons of 
% different types?
ind_first_peak = find(dirac_n(1,:)>3.2 & dirac_n(1,:)<3.5);
ind_second_peak = find(dirac_n(1,:)>4.8 & dirac_n(1,:)<5.2);
big_peak = sqrt(trapz(power_spectrum(ind_first_peak,1)));
small_peak = sqrt(trapz(power_spectrum(ind_second_peak)));
ratio = big_peak/small_peak;
fprintf('\tThe ratio of the integral of the two peaks is :\n');
disp(ratio)

% Plot the real part of the free induction decay as a function of time.
figure(1)
plot(time_list, real(data_complex), '-r', 'LineWidth', 0.1)
title1 = title({'Real part of the free induction',...
             'decay as a function of time'});
set(title1,'FontName','Arial','FontSize',10);
leg1 = legend('R(FID)','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 10);
ylabel('Amplitude [-]','FontName','Arial','FontSize',10);
grid on;
hold off 
filename='./plot/question_2_real_part.eps';
print(gcf,'-depsc',filename)

% Zoom.
figure(2)
plot(time_list, real(data_complex), '-r', 'LineWidth', 0.1)
title1 = title({'Real part of the free induction decay',...
             'as a function of time (Zoom)'});
set(title1,'FontName','Arial','FontSize',12);
leg1 = legend('R(FID)','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
xlabel('Time [s]','FontName','Arial','FontSize',10);
ylabel('Amplitude [-]','FontName','Arial','FontSize',10);
xlim([0,0.04]);
grid on;
hold off 
filename='./plot/question_2_real_part_ZOOM.eps';
print(gcf,'-depsc',filename)

% Plot the imaginary part of the free induction decay as a function of 
% time.
figure(3)
plot(time_list, imag(data_complex), '-b', 'LineWidth', 0.1)
title1  = title({'Imaginary part of the free induction',...
             'decay as a function of time'});
set(title1,'FontName','Arial','FontSize',12)
leg1 = legend('I(FID)','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
xlabel('Time [s]','FontName','Arial','FontSize',10);
ylabel('Amplitude [-]','FontName','Arial','FontSize',10);
grid on;
hold off
filename='./plot/question_2_imag_part.eps';
print(gcf,'-depsc',filename)

%Zoom
figure(4)
plot(time_list, imag(data_complex), '-b', 'LineWidth', 0.1)
title1 = title({'Imaginary part of the free induction',...
              'decay as a function of time'});
set(title1,'FontName','Arial','FontSize',12)
leg1 = legend('I(FID)','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
xlabel('Time[s]','FontName','Arial','FontSize',10);
ylabel('Amplitude [-]','FontName','Arial','FontSize',10);
xlim([0,0.04]);
grid on;
hold off
filename='./plot/question_2_imag_part_ZOOM.eps';
print(gcf,'-depsc',filename)

% Plot the real and imaginary part of the free induction decay as a 
% function of time on the same graph.
figure(5)
plot(time_list, real(data_complex), '-r', 'LineWidth', 0.1)
hold on
plot(time_list, imag(data_complex), '-b', 'LineWidth', 0.1)
leg1 = legend('R(FID)','I(FID)','Location','SouthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Real and Imaginary part of the free induction',...
             'decay as a function of time'});
set(title1,'FontName','Arial','FontSize',12)
xlim([0,0.04]);
xlabel('Time [s]','FontName','Arial','FontSize',10);
ylabel('Amplitude [-]','FontName','Arial','FontSize',10);
grid on;
hold off
filename='./plot/question_2_both_part.eps';
print(gcf,'-depsc',filename)

% Plot of the relative error between the twon implementation of the FFT in
% function of fn
figure(6)
plot(fn, dif, '-b','LineWidth',0.1)
title1 = title({'Relative error between myfft implementation and',...
             'the in build Matlab function fft in function of the frequency:'});
set(title1, 'FontName','Arial','FontSize',12);
leg1 = legend('Relative Error','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
xlabel('Frequency [Hz]', 'FontName','Arial','FontSize',10);
ylabel('Amplitude [-]', 'FontName','Arial','FontSize',10);
grid on;
hold off 
filename='./plot/question_error_comparison.eps';
print(gcf,'-depsc',filename)

% Plot the real part of the discrete Fourier transform of the complex FID 
% as a function of of the chemical shift.
figure(7)
plot(dirac_n, real(mydft_fid),'-r','LineWidth',0.1)
title1 = title({'Real part of the DFT of the complex FID as a',...
             'function of of the chemical shift'});
set(title1,'FontName','Arial','FontSize',12);
leg1 = legend('R(F[FID])','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
xlabel('Chemical shift [ppm]','FontName','Arial','FontSize',10);
ylabel('Amplitude [-]','FontName','Arial','FontSize',10);
grid on;
hold off 
filename='./plot/question_4_real_part.eps';
print(gcf,'-depsc',filename)

% Plot the imaginary part of the discrete Fourier transform of the complex
% FID as a function of of the chemical shift.
figure(8)
plot(dirac_n, imag(mydft_fid),'-b','LineWidth',0.1)
title1 = title({'Imaginary part of the DFT of the complex FID',...
               'as a function of of the chemical shift'});
set(title1,'FontName','Arial','FontSize',12);
leg1 = legend('I(F[FID])','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
xlabel('Chemical shift [ppm]','FontName','Arial','FontSize',10);
ylabel('Amplitude [-]','FontName','Arial','FontSize',10);
xlim([0,6]);
grid on;
hold off 
filename='./plot/question_4_imag_part.eps';
print(gcf,'-depsc',filename)

% Plot both curve in the same plot
figure(9)
plot(dirac_n, real(mydft_fid),'-r','LineWidth',0.1)
hold on
plot(dirac_n, imag(mydft_fid),'-b','LineWidth',0.1)
leg1 = legend('R(F[FID])','I(F[FID])','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'DFT of the complex FID as a function',...
                'of of the chemical shift'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Chemical shift [ppm]','FontName','Arial','FontSize',10);
ylabel('Amplitude [-]','FontName','Arial','FontSize',10);
xlim([0,6]);
grid on;
hold off
filename='./plot/question_4_both_part.eps';
print(gcf,'-depsc',filename)

% Plot the power spectrum of the discrete Fourier transform of the complex FID as a function of of the chemical shift.
figure(10)
plot(dirac_n, power_spectrum,'-b','LineWidth',0.1);
hold on
plot(3.333,2.712*1e9,'rx');
hold on
plot(4.961,1.081*1e9,'rx')
set(gca, 'XDir','reverse')
title1 = title({'Power spectrum (Norm) of the DFT of the complex FID',...
               'as a function of of the chemical shift'});
leg1 = legend('|| F[FID] ||', 'Peaks at 3.33 ppm and 4.96 ppm','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
set(title1,'FontName','Arial','FontSize',12);
xlabel('Chemical shift [ppm]','FontName','Arial','FontSize',10);
ylabel('Amplitude [-]','FontName','Arial','FontSize',10);
xlim([0,6.5]);
grid on;
hold off 
filename='./plot/question_4_power.eps';
print(gcf,'-depsc',filename)


