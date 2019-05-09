clear all;
close all;
clc;
n = 15;
% Initialisation
time = eye(n,1);
time_cycle = eye(n,1);
count_list = eye(n,1);
count_list_cycle = eye(n,1);

% Comparison
for i = 1:n
    rng(i*n)
    name = ['rmg(1:',int2str(n*i),')'];
    try
        [~, m_rl, m_cp] = evalc(name);
    catch e
        [~, m_rl] = evalc(name);
        m_cp = m_rl;
    end
    tic
    [val, count] = eig_j(m_rl);
    time(i) = toc;
    count_list(i) = count;
    tic
    [val, count] = eig_cj(m_rl);
    time_cycle(i) = toc;
    count_list_cycle(i) = count;
end 

x = linspace(5,350,n);
% Plot allowing the time comparison
figure(1)
plot(x,time, '-r', 'LineWidth', 0.1)
hold on 
plot(x,time_cycle, '-b', 'LineWidth', 0.1)
leg1 = legend('Classic method time','Cyclic method time','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Time comparison between the classic and the cyclic method  ',...
             'for a broad range of matrix sizes :'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Size of matrix (n*n)','FontName','Arial','FontSize',10);
ylabel('Time [s]','FontName','Arial','FontSize',10);
grid on;
hold off
filename='./plot/time_comparison.eps';
print(gcf,'-depsc',filename)

% Plot allowing the iteration comparison
figure(2)
plot(x,count_list, '-r', 'LineWidth', 0.1)
hold on 
plot(x,count_list_cycle, '-b', 'LineWidth', 0.1)
leg1 = legend('Classic method number of Jacobi rotations','Cyclic method number of iteration','Location','NorthEast');
set(leg1,'FontName','Arial','FontSize',10)
title1 = title({'Number of Jacobi rotations comparison between the classic',...
     ' and the cyclic method  for a broad range of matrix sizes :'});
set(title1,'FontName','Arial','FontSize',12)
xlabel('Size of matrix (n*n)','FontName','Arial','FontSize',10);
ylabel('Number of Jacobi rotations [-]','FontName','Arial','FontSize',10);
grid on;
hold off
filename='./plot/iter_comparison.eps';
print(gcf,'-depsc',filename)