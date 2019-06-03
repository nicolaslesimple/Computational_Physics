clear all;
close all;
clc;

%%%%% Non Linear conjugate gradient method %%%%%

% Declaration of the parameters of our algorithm
Number_particles_max=10;
E=zeros(Number_particles_max,1);
x=cell(Number_particles_max,1);
tol = 1e-4;
tol_inside = 1e-12;

% Non linear method to find the configuration to have a min energy
for N=2:Number_particles_max
    
    % First declaration per particles of the several parameters of the
    % simulation
    x{N}=rand(2*N,1); % Random initialisation for the position vector
    r=-grad(x{N}); d=r; alpha=1;
    b=0;
    
    while norm(r)>tol
        while norm(alpha*d)>tol_inside
            hess=hessian(x{N},d); % Hessian calculation
            if hess==0
                x{N}=x{N}-1e-4*grad(x{N})/norm(grad(x{N})); % Special case if the hessian is 0
            end
            % Update of the variable as we are computing an iterative method
            alpha=(-grad(x{N})'*d)/(d'*hess*d);
            x_after=x{N}+alpha*d;
            if energy_total(x_after)>energy_total(x{N}) % If with find a local min, we will search in this direction
                % This conditions allows a faster computation as it
                % optimize the search direction
                x{N}=x{N}-1e-4*grad(x{N})/norm(grad(x{N}));
                b=1;
                break
            else
                x{N}=x_after; % Classic update of x_i
            end
            
        end
        
    % Update of the variable as we are computing an iterative method
    rtmp=-grad(x{N});
    beta=rtmp'*rtmp/(r'*r);
    r=rtmp;
    d=r+beta*d;
    end
    
    % Calculus and storage of the min energy of the system corresponding to
    % the configuration we found with our algorithm
    E(N)=energy_total(x{N});
    
    % For each configuration, we plot the position of the particles
    figure
    title(['Position of each particles in a system composed of N = ' num2str(N) ' particles :'])
    hold on
    plot(x{N}(linspace(1,N*2-1,N)), x{N}(linspace(2,N*2,N)),'o')
    xlim([min(x{N})-2, max(x{N})+2])
    xlabel('x distance')
    ylabel('y distance')
    grid on
    axis equal
    filename=['./plot/' num2str(N) '.eps'];
    print(gcf,'-depsc',filename)
    
    if N == Number_particles_max 
        fprintf('The energy for each configuration (N=1...N_max) is display below :\n')
        disp(E)
    end
    
end
    
% Plot the energy in function of the number of particles in the system
figure
title('Energy of the system in function of the number of particles inside the system :')
hold on
plot(1:Number_particles_max,E)
xlabel('Number of particles N')
ylabel('Total energy E_{tot}')
set(gca,'Yscale','log')
filename='./plot/energy_plot.eps';
print(gcf,'-depsc',filename)




%%%%%%% Gradient Calculus %%%%%%%
function gradient_value=grad(x)
    % This function allows the calculus of the Gradient of the function 
    % using line search and Taylor expansion.
    % ARGS :
    %   - x : Position vector
    % RETURN :
    %   - grad_value : Float corresponding to the value of the gradient 
    %   corresponding to the position vector given as input.
    gradient_value=zeros(size(x));
    step=1e-7;
    for i=1:length(x)
        unity_vector=zeros(size(x));
        unity_vector(i)=1;
        gradient_value(i)=energy_total(x+unity_vector*step)-energy_total(x);
    end
    gradient_value=gradient_value/(step);
end



%%%%%%% Hessian Calculus %%%%%%%
function Hess=hessian(x,d)
    % This function allows the calculus of the Hessian using line search
    % and Taylor expansion.
    % ARGS :
    %   - x : Position vector
    %   - d : Number representing to point on which we want to calculate
    %   the Hessian.
    % RETURN :
    %   - H : Float corresponding to the value of the hessian.
    d=d/norm(d);
    step=1e-4;
    Hess=(energy_total(x+d*step)-2*energy_total(x)+energy_total(x-d*step))/step^2;
end


%%%%%% Total Energy Calculus %%%%%%
function e = energy_total(x)
    % This function allows the calculation of the total energy of the
    % system given a position vector.
    % ARGS :
    %   - x : Position vector
    % RETURN : 
    %   - e : Float representing the total energy of the system
    % Harmonic Energy Calculus
    E_harm = 0;
    for i=1:length(x)   
        E_harm=E_harm+x(i)*x(i);    
    end
    % Coulomb Energy Calculus
    E_coul = 0;
    for i=1:2:length(x)
        for j=1:2:length(x)
            if i~=j
                E_coul=E_coul+1/(norm([x(i:i+1)-x(j:j+1)]));
            end
        end
    end
    % Total Energy Calculus
    e = E_harm + E_coul;
end
