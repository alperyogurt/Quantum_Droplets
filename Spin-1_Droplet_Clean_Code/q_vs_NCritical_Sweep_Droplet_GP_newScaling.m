%% In this file, one can assign a vector of particle numbers and calculate the
% corresponding wavefunction, mu, and energy of the droplet for the inputs
% of rStep, rInterval, quality factors. This function provides the figure
% for the ground state wavefunctions and the corresponding Thomas-Fermi
% wavefunction for the input vector of N values.
% 
%%
load('phi1DerVectors_v1.mat');
load('phi1Vectors_v1.mat');


t1Vector = [0:0.1:1];
%t1Vector = 0;
% Assign the N values to be computed
NDummyOld = 19;
%t1 = q/n_0|c_1|
NCriticalVector = zeros(size(t1Vector));
NResolution = 1;

rStep = 0.1;
rInterval = [0 12];
tStep = rStep/1000;
tIteration = 10000;
expectedPercentageQuality = 1;
expectedMuQuality = 0.1;
r = [rInterval(1):rStep:rInterval(2)];
% phiGroundNQMatrix = zeros(length(NVector),length(t1Vector),length(r));
% phiTFQMatrix = zeros(length(NVector),length(t1Vector),length(r));
% ENQ = zeros(length(NVector),length(t1Vector));
% muNQ = zeros(length(NVector),length(t1Vector));
%     
    
for(kk=1:length(t1Vector))
    t1 = t1Vector(kk)
    complete = 0;
    [dummyPhiGround,r1,mu1,varMu1,muPercentageQuality1] = dropletWavefunctionImagTimeSpin1NewScaling_v2(t1,rStep,rInterval,tStep,tIteration,NDummyOld,...
           expectedPercentageQuality,expectedMuQuality);
    muDummyOld = mu1;
    
    while (complete == 0)
        
        if (abs(muDummyOld) < 0.01)
            NDummyNew = NDummyOld + NResolution;
        end
        if (abs(muDummyOld) > 0.01)
            NDummyNew = NDummyOld - NResolution;
        end

        [dummyPhiGround,r1,mu1,varMu1,muPercentageQuality1] = dropletWavefunctionImagTimeSpin1NewScaling_v2(t1,rStep,rInterval,tStep,tIteration,NDummyNew,...
               expectedPercentageQuality,expectedMuQuality);
        muDummyNew = mu1;

        if(abs(muDummyNew)< 0.01 && abs(muDummyOld)> 0.01)
            NCriticalVector(kk) = NDummyOld;
            complete = 1;
        end
        if(abs(muDummyNew)> 0.01 && abs(muDummyOld)< 0.01)
            NCriticalVector(kk) = NDummyNew;
            complete = 1;
        end
        muDummyOld  = muDummyNew;
        NDummyOld = NDummyNew;
    end
    %Find the droplet wavefunction and chemical potential
        
end

figure(5);
plot(t1Vector,NCriticalVector,'lineWidth',2);
xlabel('$t_1 = \frac{q}{n_1 c_1}$','Interpreter','latex');
    ylabel('Critical $\tilde{N}$ ','Interpreter','latex');
    title(['Critical $\tilde{N}$ vs. $t_1$ '],...
    'Interpreter','latex');
    set(gca,'FontWeight','bold')
    set(gca, 'LineWidth',3);
    grid on;
    grid minor;

% 
% save(['DropletNumerical_SPINOR_q_and_N_Sweep_r' num2str(rInterval(1)) 'to' num2str(rInterval(2)) '_rStep_' num2str(rStep)...
%     '_tStep_' num2str(tStep) '_N_' num2str(NVector(1))...
%     'to' num2str(NVector(length(NVector))) '_t1_' num2str(t1Vector(1)) 'to' num2str(t1Vector(length(t1Vector))) '.mat'], ...
%     'ENQ','muNQ','phiGroundNQMatrix', 'phiTFQMatrix','NVector','t1Vector','r','rStep','tStep');