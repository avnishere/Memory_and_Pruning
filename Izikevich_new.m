%Izikevich model

clear all;
a = 0.02;
b = 0.2; 
c = -65;
d = 2;
v = -65;    % Initial value of v
u = b*v;   % Initial value of u

% 1. RegularSpiking <RS>
%a1 = 0.02; b1 = 0.26; 

v_record1 = []; %recording the voltages
v_record2 = [];

%Running the simulation

for t=1:400           % simulation of whatever milliseconds
   
 I = 5; % thalamic input
 
%   if t >= 30
%       I = 10;
%   end
  v_record1 = [v_record1; v]; % appending the record at every iteration
  
  if v>=30 % checking for spikes
    v = c;  
    u = u + d;
  end
  v= v+0.5*(0.04*v.^2+5*v+140-u+I); % voltage update equation
  u= u+ a*(b*v-u);                 
  
end

for t=1:400           % simulation of whatever milliseconds
   
 I = 15; % thalamic input
 
%   if t >= 30
%       I = 10;
%   end
  v_record2 = [v_record2; v]; % appending the record at every iteration
  
  if v>=30 % checking for spikes
    v = c;  
    u = u + d;
  end
  v= v+0.5*(0.04*v.^2+5*v+140-u+I); % voltage update equation
  u= u+ a*(b*v-u);                 
  
end

subplot(1,2,1)
plot(v_record1);
ylabel('Potential (mV)');
xlabel('Time (ms)');
%title('Neuron Response')
ylim([-100,40]);


subplot(1,2,2)
plot(v_record2);
ylabel('Potential (mV)');
xlabel('Time(ms)');
%title('Neuron Response)')
ylim([-100,40]);






