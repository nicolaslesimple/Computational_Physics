clear all;
close all;
clc;

%%%%%%%% First question using N = 7 at maximum %%%%%%%%%

% Declaratio of the parameters of the simulation
gamma=-2.7; 
a=0.4;
nb=200;
N=7;

for N=2:N
    
    % Declaration of each Hamiltonian Matrix
    H=zeros(2*N,2*N);
    l=0;
    
    % Creation of the Hamiltonian Matrix called H
    for k=-pi/a:pi/a/nb:pi/a
        % We first initialized the matrix position corresponding to non
        % null value for the hamiltonian
        for i=1:4:2*N
            i1=i;
            H(i1,i1+1)=1;                % left bound
            H(i1,mod(i1-2-1,2*N)+1)=1;   % lower bound
            H(i1,mod(i1+2-1,2*N)+1)=1;   % upper bound
            i2=i+1;
            H(i2,i2-1)=1;                % right bound
            H(i2,mod(i2+2-1,2*N)+1)=1;   % upper bound
            H(i2,mod(i2-2-1,2*N)+1)=1;   % lower bound
            i3=i+2;
            H(i3,i3+1)=exp(-1i*k*a);     % left bound 
            H(i3,mod(i3+2-1,2*N)+1)=1;   % upper bound
            H(i3,mod(i3-2-1,2*N)+1)=1;   % lower bound
            i4=i+3;
            H(i4,i4-1)=exp(1i*k*a);     % right bound 
            H(i4,mod(i4+2-1,2*N)+1)=1;   % upper bound
            H(i4,mod(i4-2-1,2*N)+1)=1;   % lower bound
        end
        % Increment the counter
        l=l+1;
        % Put some specific values to 0 in the Hamiltonian
        if N~=2
            H(2*N-1,1)=0; 
            H(2*N,2)=0;  
            H(1,2*N-1)=0; 
            H(2,2*N)=0;  
        end
        % Multiply the values corresponding to the definition with gamma
        Hfin=gamma*H(1:2*N,1:2*N);
        % Find column vector E containing the eigenvalues of hamiltonian 
        % and we store it in the matrix E as a column.
        Etemp=eig(Hfin);
        for i=1:2*N
            E(i,l)=Etemp(i); 
        end
     end
    
     % Plot the band structures for each N
     k=-pi/a:pi/a/nb:pi/a;
     figure 
     hold on, grid minor, box on
     for j=1:2*N
         % Plot each columns separtly has each one correspond to another
         % band structure
         plot(k,E(j,1:l))
         cont{j}=['\epsilon_{' num2str(j) '}'];
     end
     title1 = title(['N=' num2str(N)]);
     set(title1,'FontName','Arial','FontSize',12)
     legend(cont)
     xlabel('k vector [nm^{-1}]')
     ylabel('E [eV]')
     grid on;
     filename=num2str(N);
     print(gcf,'-depsc',filename)
end

clear all;
clc;
%%%%%%%% First question using N = 200 at maximum %%%%%%%%%

% Declaratio of the parameters of the simulation
gamma=-2.7; 
a=0.4;
nb=200;
N=200;

for N=2:N
    
    % Declaration of each Hamiltonian Matrix
    H=zeros(2*N,2*N);
    l=0;
    
    % Creation of the Hamiltonian Matrix called H
    for k=-pi/a:pi/a/nb:pi/a 
        for i=1:4:2*N
            % We first initialized the matrix position corresponding to non
            % null value for the hamiltonian
            i1=i;
            H(i1,i1+1)=1;              
            H(i1,mod(i1-2-1,2*N)+1)=1;   
            H(i1,mod(i1+2-1,2*N)+1)=1;   
            i2=i+1;
            H(i2,i2-1)=1;             
            H(i2,mod(i2+2-1,2*N)+1)=1;   
            H(i2,mod(i2-2-1,2*N)+1)=1;   
            i3=i+2;
            H(i3,i3+1)=exp(-1i*k*a);    
            H(i3,mod(i3+2-1,2*N)+1)=1;  
            H(i3,mod(i3-2-1,2*N)+1)=1; 
            i4=i+3;
            H(i4,i4-1)=exp(1i*k*a);   
            H(i4,mod(i4+2-1,2*N)+1)=1;  
            H(i4,mod(i4-2-1,2*N)+1)=1; 
        end
        % Increment the counter
        l=l+1;
        % Put some specific values to 0 in the Hamiltonian
        if N~=2
            H(2*N-1,1)=0; 
            H(2*N,2)=0;  
            H(1,2*N-1)=0; 
            H(2,2*N)=0;  
        end
        % Multiply the values corresponding to the definition with gamma
        Hfin=gamma*H(1:2*N,1:2*N);
        % Find column vector E containing the eigenvalues of hamiltonian 
        % and we store it in the matrix E as a column.
        Etemp=eig(Hfin);
        for i=1:2*N
            E(i,l)=Etemp(i); 
        end
    end
    
    k=-pi/a:pi/a/nb:pi/a;
    % Band Gap assumes the eigenvalues are sorted in ascending order.
    E=sort(E,2);
    % Band gap calculation for each N
    Eg(N-1)=min(min(E(end/2+1:end,:)))-max(max(E(1:end/2,:)));
end

% Select specific values of the Band Gap
Egmax=Eg(Eg>0.01);
N=2:N;
Nmax=N(Eg>0.01);

% Plot the magnitude of the Band Gap
figure(7)
hold on, grid minor, box on
plot(N,Eg,'-..',Nmax,Egmax,'r+')
xlabel('Number of nanoribons N [-]')
ylabel('Energy gap E_g [eV]')
legend('E_g','E_{g\neq0}\equivE_{g,max}')
filename='.\plot\magnitude_band_gap';
print(gcf,'-depsc',filename)