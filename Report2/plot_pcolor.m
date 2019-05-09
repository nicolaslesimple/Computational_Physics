function plot_pcolor(A, title1, title2, filename1, filename2)
% Plot the wavefunctions which the eigenvector
    figure
    pcolor(A)
    shading interp
    colorbar
    title1_ = title({title1});
    set(title1_,'FontName','Arial','FontSize',12)
    xlabel('Index of columns corresponding to the 3 lowest eigenvalues','FontName','Arial','FontSize',10);
    ylabel('Index corresponding to the position of each point in the vector','FontName','Arial','FontSize',10);
    hold off
    print(gcf,'-depsc',filename1)

    % Plot the wavefunctions which the eigenvector
    figure
    pcolor((A.*A)/norm(A)^2)
    shading interp
    colorbar
    title2_ = title({title2});
    set(title2_,'FontName','Arial','FontSize',12)
    xlabel('Index of columns corresponding to the 3 lowest eigenvalues','FontName','Arial','FontSize',10);
    ylabel('Index corresponding to the position of each point in the vector','FontName','Arial','FontSize',10);
    grid on;
    hold off
    print(gcf,'-depsc',filename2)