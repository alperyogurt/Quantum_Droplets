%% In this file, one can assign a vector of particle numbers and t_1 = q/n_1 c_1
% values, then calculate the corresponding ground state wavefunction, mu, 
% and energy of the spin-1 droplet. 

% Inputs of rStep, rInterval, quality factors are required. 
% This function provides the figures for the ground state wavefunctions of
% different N and t_1 values that had been executed. 
% 
%%
% Load the look up tables for the I(t) and I'(t)
load('ItLookUpTable.mat');
load('ItDerivativeLookUpTable.mat');

% Assign the N values to be computed (One may choose a single value or a
% vector
NVector= [40];

%Assign the t_1 values to be computed (One may choose a single value or a
% vector
t1Vector = [0.5];

%Determine the radial vector interval and step size
rStep = 0.1;
rInterval = [0 12];

%Determine the time step for the imaginary time evolution
tStep = rStep/1000;

%Determine the number of iterations for the difference equation to check
%whether the resulting wavefunction satisfies the termination criteria
tIteration = 10000;

%Assign the termination criteria for the chemical potential vector
expectedPercentageQuality = 0.1;
expectedMuQuality = 0.01;

%Create the radius vector
r = [rInterval(1):rStep:rInterval(2)];

%Create the matrice that keeps the ground state wavefunctions
phiGroundNQMatrix = zeros(length(NVector),length(t1Vector),length(r));

%Create the matrice that keeps the MF, LHY, KE, mu and total energy values
MFNQ = zeros(length(NVector),length(t1Vector));
LHYNQ = zeros(length(NVector),length(t1Vector));
KENQ = zeros(length(NVector),length(t1Vector));
ENQ = zeros(length(NVector),length(t1Vector));
muNQ = zeros(length(NVector),length(t1Vector));
    
% Sweep the t_1 and N values and save the resulting ground state
% wavefunctions and energy values
for(kk=1:length(t1Vector))
    
    %Assign the executed t_1 value 
    t1 = t1Vector(kk)
    
    %Find the droplet wavefunction and chemical potential
    for ii=1:length(NVector)
        %Keep track of execution time
        tic 
        %Calculate the Ground State Wavefunction
        [phiGroundNQMatrix(ii,kk,:),r1,mu1,varMu1,muPercentageQuality1] = CalculateDropletSpin1GroundStateWavefunction(t1,...
            rStep,rInterval,tStep,tIteration,NVector(ii),...
           expectedPercentageQuality,expectedMuQuality);
       
        muNQ(ii,kk) = mu1;
        
        dummyPhiGround = phiGroundNQMatrix(ii,kk,:);
        

%         % Calculate the energy of the droplet
%         [ENQ(ii,kk), MFNQ(ii,kk), LHYNQ(ii,kk), KENQ(ii,kk)]= dropletSpin1TotalEnergy(t1,r1,rStep,...
%             dummyPhiGround,tVectorExtended,phi1VectorExtended);
        
        %Show values of N, mu varMu to keep track of the process
        NVector(ii)
        mu1
        varMu1
        elapsed_time1 = toc
    end
end



%Demonstrate the resulting wavefunctions for each N values with the all
%possible t_1 values that had been swept
for ii= 1:length(NVector)
    figure(ii);
    for kk= 1:length(t1Vector)
    plot(r,squeeze(phiGroundNQMatrix(ii,kk,:)),'lineWidth',2);
    hold on;
    end
    
    hold off;
    grid on;
    grid minor;
    xlabel('The radial distance $\tilde{r}$','Interpreter','latex');
    ylabel('The droplet wavefunction \psi');
    title(['Droplet Wavefunction $\psi$ vs. Radial Distance $\tilde{r}$ with N= ' num2str(NVector(ii)) ...
        ' for different q or $t_1$: ' num2str(t1Vector(1))...
        ' to ' num2str(t1Vector(length(t1Vector))) ' values'],...
    'Interpreter','latex');
    set(gca,'FontWeight','bold')
    set(gca, 'LineWidth',3);
    
    for(mm=1:length(t1Vector))
        legendCell{mm} = ['t1 = '  num2str(t1Vector(mm)) ' \mu = ' num2str(muNQ(ii,mm))];
%       legendCell{2*ii} = ['TF' num2str(ii)];
    end
    legend(legendCell);
end


%Demonstrate the resulting wavefunctions for each t_1 values with the all
%possible N values that had been swept
for kk= 1:length(t1Vector)
    figure(100+kk);
    for ii= 1:length(NVector)
        plot(r,squeeze(phiGroundNQMatrix(ii,kk,:)),'lineWidth',2);
        hold on;
    end
    
    hold off;
    grid on;
    grid minor;
    xlabel('The radial distance $\tilde{r}$','Interpreter','latex');
    ylabel('The droplet wavefunction \psi');
    title(['Droplet Wavefunction $\psi$ vs. Radial Distance $\tilde{r}$ with $t_1$= ' num2str(t1Vector(kk))...
        ' for different N: ' num2str(NVector(1)) ' to ' num2str(NVector(length(NVector)))  ' values'],...
    'Interpreter','latex');
    set(gca,'FontWeight','bold')
    set(gca, 'LineWidth',3);
    for(mm=1:length(NVector))
        legendCell{mm} = ['N = '  num2str(NVector(mm)) ' \mu = ' num2str(muNQ(mm,kk))];
    end
    legend(legendCell);
end


%Demonstrate the chemical potential vs. t_1 graph for each N values 
for ii= 1:length(NVector)
    figure(200+ii);
    plot(t1Vector,muNQ(ii,:),'lineWidth',2);
    hold off;
    grid on;
    grid minor;
    xlabel('$t_1$','Interpreter','latex');
    ylabel('The chemical potential $\mu$','Interpreter','latex');
    title(['The chemical potential $\mu$ vs. $t_1$ with N= ' num2str(NVector(ii))],...
    'Interpreter','latex');
    set(gca,'FontWeight','bold')
    set(gca, 'LineWidth',3);
end


%Demonstrate the chemical potential vs. N graph for each t_1 values 
for ii= 1:length(t1Vector)
    figure(300+ii);
    plot(NVector,muNQ(:,ii),'lineWidth',2);
    hold off;
    grid on;
    grid minor;
    xlabel('$\tilde{N}$','Interpreter','latex');
    ylabel('The chemical potential $\mu$','Interpreter','latex');
    title(['The chemical potential $\mu$ vs. $\tilde{N}$ with t1= ' num2str(t1Vector(ii))],...
    'Interpreter','latex');
    set(gca,'FontWeight','bold')
    set(gca, 'LineWidth',3);
end



%Save the results of the simulation
save(['DropletNumerical_SPINOR_q_and_N_Sweep_r' num2str(rInterval(1)) 'to' num2str(rInterval(2)) '_rStep_' num2str(rStep)...
    '_tStep_' num2str(tStep) '_N_' num2str(NVector(1))...
    'to' num2str(NVector(length(NVector))) '_t1_' num2str(t1Vector(1)) 'to' num2str(t1Vector(length(t1Vector))) '.mat'], ...
    'muNQ','phiGroundNQMatrix','NVector','t1Vector','r','rStep','tStep','ENQ','MFNQ','LHYNQ','KENQ');