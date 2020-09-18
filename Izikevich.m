%% The Izhikevich Model

%% 2. Integrate-and-fire neurons

% The work in this section is based on Chapter 5 of Dayan and Abbott. The
% following tutorial is also useful:
% http://neuroscience.ucdavis.edu/goldman/Tutorials_files/Integrate%26Fire.pdf

% Parameters provided
E_L = -65; % resting potential of the model cell
R_m = 90; % total membrane resistance
tau_M = 30; % membrane time constant
V_thresh = -50; % threshold for activation of an action potential
V_reset = -65; % value to which the potential is reset following a spike

a = 0.02; % recovery time scale variable
b = 0.2; % sensitivity variable
d = 2; % after spike reset


% The simulation will last 500ms
dt = 0.01; % length of time step (ms)
simLength = 500/dt; % number of steps in simulation 

% From Dayan and Abbott, p164, 2nA would be a reasonable current pulse
I_e = 2;
I = [zeros(1,0.2*simLength) I_e*ones(1,0.6*simLength) zeros(1,0.2*simLength)];

% We use time steps of 1 ms
timesteps = 0:simLength;

% We want to record the membrane potential, V, for the requisite period,
% assuming that its initial value is equal to the resting potential
V = zeros(1,simLength);
V(1) = E_L;

% First we see how the resting membrane potential would change without
% spiking...
for i=1:simLength
    V(i+1) = E_L + R_m*I(i) + (V(i) - E_L - R_m*I(i))*exp(-dt/tau_M); % update rule is equation 5.9 on D&A p164
end

% ...plotting the results
subThreshVoltageTrace = figure;
subplot(2,1,1)
plot(I,'r')
title('External input to the neuron','FontSize',14)
xlabel('Time (ms)','FontSize',12)
ylabel('Current (nA)','FontSize',12)
axis([0 simLength -0.5 5])
vals = get(gca,'XTick');
set(gca,'XTickLabels',vals*dt) % adjust the scale of the graph to account for the size of the time steps
subplot(2,1,2)
plot(V)
xlim([0 simLength])
vals = get(gca,'XTick'); 
set(gca,'XTickLabels',vals*dt) % adjust the scale of the graph to account for the size of the time steps
title('Change in resting membrane potential (no spiking)','FontSize',14)
xlabel('Time (ms)','FontSize',12)
ylabel('Membrane potential (mV)','FontSize',12)
saveas(gcf, '/Users/JGRaymond/Documents/UniversityTeaching/CCN/LabSolutions/Lab2Report/subThreshVoltage.png', 'png')


% Next we see what happens when we add spiking...
V_spk = zeros(1,simLength);
V_spk(1) = E_L;
V_max = 30; % we use a trick from the ucdavies tutorial to get nicer plots; first set a V_{max} for spikes
V_plot = V_spk; % next make a copy of the vector where we record voltage values
for i=1:simLength
    V_spk(i+1) = E_L + R_m*I(i) + (V_spk(i) - E_L - R_m*I(i))*exp(-dt/tau_M); % update rule is still equation 5.9 on D&A p164
    if (V_spk(i+1) > V_thresh) %cell spiked
        V_spk(i+1) = V_reset; %set voltage back to V_reset
        V_plot(i+1) = V_max; % if we have crossed the threshold value, plot a spike at the chosen maximum membrane potential
    else
        V_plot(i+1) = V_spk(i+1);
    end
end

% ...plotting the results for this as well (note that without the plotting 
% trick the graph would never get above -50 mV)
voltageTraceWithSpks = figure;
subplot(2,1,1)
plot(I,'r')
title('External input to the neuron','FontSize',14)
xlabel('Time (ms)','FontSize',12)
ylabel('Current (nA)','FontSize',12)
axis([0 simLength -0.5 5])
vals = get(gca,'XTick');
set(gca,'XTickLabels',vals*dt)
subplot(2,1,2)
plot(V_plot)
xlim([0 simLength])
vals = get(gca,'XTick');
set(gca,'XTickLabels',vals*dt)
title('Change in resting membrane potential (spiking)','FontSize',14)
xlabel('Time (ms)','FontSize',12)
ylabel('Membrane potential (mV)','FontSize',12)
saveas(gcf, '/Users/JGRaymond/Documents/UniversityTeaching/CCN/LabSolutions/Lab2Report/spikingVoltage.png', 'png')


% We can work out the firing rate using equation 5.11 from p164 of D&A 
% (note that we multiply the results by 1000 to go from mHz to Hz)
r_isi = 1000*(1/(tau_M*log((R_m*I_e + E_L - V_reset)/(R_m*I_e + E_L - V_thresh))));

% Aside: how does this the output of the simulation compare to this? Does 
% the performance of the simulation improve if you reduce the length of the 
% time steps?
r_sim = length(find(V_plot>0))*(10/3); % spike rate during stimulation in simulated data

% For a range of input values, we obtain...
inputs = 0:0.1:2;
rates = zeros(length(inputs));
for i=1:length(inputs)
    if R_m*inputs(i) > V_thresh - E_L % note that there is a condition that must be satisfied for the equation to hold
        rates(i) = 1000*(1/(tau_M.*log((R_m*inputs(i) + E_L - V_reset)/(R_m*inputs(i) + E_L - V_thresh))));
    end
end

% ... which accurately reproduces the graph from p165 of Dayan and Abbott
rateVsCurrent = figure;
plot(inputs,rates,'k')
title('Interspike-interval firing rate as a function of current','FontSize',14)
xlabel('I_{e} (nA)','FontSize',12)
ylabel('r_{ISI} (Hz)','FontSize',12)
saveas(gcf, '/Users/JGRaymond/Documents/UniversityTeaching/CCN/LabSolutions/Lab2Report/ratesVsCurrents.png', 'png')

% We now add spike rate adaptation, which will involve solving the equation
% numerically. For an example of how to do this, see:
% http://people.brandeis.edu/~pmiller/COMP_NEURO/MATLAB/SRA_tref_LIF.m
V_spkadp = zeros(1,simLength);
gsra = zeros(1,simLength);
tau_sra = 100; % time constant for spike rate adaptation from Figure 1(C) 
E_K = -70;
V_spkadp(1) = E_L;
gsra(1) = 0; % initial spike rate adaptation is 0
V_max = 30; % we use a trick from the ucdavies tutorial to get nicer plots; first set a V_{max} for spikes
V_adpplot = V_spkadp; % next make a copy of the vector where we record voltage values
for i=1:simLength
    V_spkadp(i+1) = V_spkadp(i) + dt/tau_M*(E_L-V_spkadp(i)-gsra(i)*(V_spkadp(i)-E_K)+R_m*I(i)); % numerical solution to membrane potential equation 
    if (V_spkadp(i+1) > V_thresh) %cell spiked
        V_spkadp(i+1) = V_reset; %set voltage back to V_reset
        V_adpplot(i+1) = V_max; % if we have crossed the threshold value, plot a spike at the chosen maximum membrane potential
        gsra(i+1) = gsra(i) + 0.06; % whenever the neuron fires a spike we increase the hyperpolarizing current due to spike rate adaptation
    else
        V_adpplot(i+1) = V_spkadp(i+1);
        gsra(i+1) = gsra(i)+dt/tau_sra*-1*gsra(i); % numerical solution to spike rate adaptation equation
    end
end

% We want to check how things change. There are a couple of ways of doing
% this...

% Does spike-rate adaptation slow down the firing rate? It probably
% should...
r_simadp = length(find(V_adpplot>0))*(10/3);

fprintf('\nSpike rate adaptation reduced the firing rate by %4.2f Hz.', r_sim-r_simadp)

% And it should also have an effect on the ISIs...
spikeTimesNoAdap = dt*find(V_plot>0);
spikeTimesAdap = dt*find(V_adpplot>0);

spikeTimesAdapVsNoAdap = figure;
plot(1:length(find(V_adpplot>0))-1,diff(spikeTimesAdap),'r')
hold on
plot(1:length(find(V_plot>0))-1,diff(spikeTimesNoAdap),'b')
title('Evolution of ISIs with and without spike rate adaptation','FontSize',14)
xlabel('Spike Count','FontSize',12)
ylabel('ISI (ms)','FontSize',12)
saveas(gcf, '/Users/JGRaymond/Documents/UniversityTeaching/CCN/LabSolutions/Lab2Report/adaptationISIs.png', 'png')

voltageTraceAdp = figure;
subplot(2,1,1)
plot(I,'r')
title('External input to the neuron','FontSize',14)
xlabel('Time (ms)','FontSize',12)
ylabel('Current (nA)','FontSize',12)
axis([0 simLength -0.5 5])
vals = get(gca,'XTick');
set(gca,'XTickLabels',vals*dt)
subplot(2,1,2)
plot(V_adpplot)
xlim([0 simLength])
vals = get(gca,'XTick');
set(gca,'XTickLabels',vals*dt)
title('Change in resting membrane potential with spike rate adaptation','FontSize',14)
xlabel('Time (ms)','FontSize',12)
ylabel('Membrane potential (mV)','FontSize',12)
saveas(gcf, '/Users/JGRaymond/Documents/UniversityTeaching/CCN/LabSolutions/Lab2Report/adaptationVoltage.png', 'png')