%%Simple function plot example
%You are recommended to check the help text of each command 
%"help <command>" in the console
clear all;
close all;
clc;

%Generate an indexed figure
figure(1) 

%Some parameters definition
l=0.2;
L=2;

%Define a linear grid and evaluate the function sine over it.
X1=linspace(-1,6,100);
Y1=sin(X1);
%Plot the data series
plot(X1,Y1,'-rx','LineWidth',2,'MarkerSize',6)

%Set general font type and size (it is advisable to do it with reference to the axes tics, since unlike other objects those are hardly modifiable)
set(gca,'FontName','Arial','FontSize',15)

%Title
title1=title('Example plot');
set(title1,'FontName','Arial','FontSize',16)

%Ranges
xlim([-1,6]);
ylim([-2,2]);

%Labels
xlabel('X','FontName','Arial','FontSize',16);
ylabel('Y','FontName','Arial','FontSize',16);

%Add another data series on the same plot

hold on %not to initiate a different plot
%Evaluate a function over an automatic mash
[X2,Y2]=fplot(@(x)sin(2*pi*x/l)*exp(-x/L),[0,5]);
%Plot it...
plot(X2,Y2,'-bo','LineWidth',2,'MarkerSize',4);
hold off %to close the plot

%Legend
leg1=legend('sin(2pi*x/l)*exp(-x/L)  l=0.2  L=2', 'sin(x)','Location','SouthEast');
set(leg1,'FontName','Arial','FontSize',16)

%Size and ratio of the figure [left offset, bottom offset, width, height]
set(gcf,'Position',[0 0 600 500]);

%Export to eps format (a file is created in the current folder)
filename='prova.eps';
print(gcf,'-depsc',filename)

%% You can play with different plotting commands both 2D and 3D
%ex. subplot, semilogx, semilogy, loglog, hist, surf, ezsurf, etc...  

