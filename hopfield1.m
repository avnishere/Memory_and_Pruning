%% HOPFIELD     NETWORK

    capacity = 10;

%% Initializing Parameters

%Weight Matrix
weight_mat = zeros(100,100); % preallocation of matrix size
weight_mat = weight_mat - diag(weight_mat); % eliminating self-loops

% Temperature or scaling factor for probabilistic model
T = 4;

% The input is provided as the intial state of the network. 
%% Storage Stage
memories = randi([0,1],capacity,100).*2-1;
weight_matrix = insert_memory(memories);

% construct the network
neuron_network= ones(10,10);

%% Pruning the network

%collect eucledian distance

total_neurons = size(neuron_network,1)*size(neuron_network,2);
euclid = zeros(total_neurons, total_neurons);
p = 1;
[XGrid, YGrid] = meshgrid(1:size(neuron_network,1), 1:size(neuron_network,2));
for i = 1:size(neuron_network,1)
    for j= 1:size(neuron_network,2)
    temp_matrix = sqrt((i - XGrid) .^ 2 + (j - YGrid) .^ 2);
    euclid(p,:) = reshape(temp_matrix,[1, total_neurons]);
    p = p+1;
    
    end
end

% performing pruning

prun_weights  = weight_matrix;
prun_weights( abs(prun_weights) < 0.6*euclid ) = 0;

%% Retreival Stage

% creating the degraded input
test_signal = memories(2,:);

%degrading the pattern
deg = 20;
 for i = 20:20+deg
     if test_signal(i) == 1
         test_signal(i) = -1;
     else
         test_signal(i) = 1;
     end
 end

%transposing for multiplication
test_signal = test_signal';

% updating the nodes  

iter = 1;
hammered= 0;
while hammered<=2 || iter <= 100
    
    for i = 1:100
    
        energies(i) = energy(weight_matrix(i,:), test_signal); % or weight_matrix for nor normal people
        test_signal(i) = classify(energies(i)); % or classify for sign func
    
    end
    hammer = test_signal-memories(2,:)';
    hammered = size(find(hammer~= 0),1);
    iter = iter + 1;
end

hammered

%% Local Functions

% State classification function
function s = classify(energies)

    s = energies;
    s (s>0) = 1;
    s(s<=0) = -1;
 
end

%Energy Computation
function e = energy(weights, states)

    e = weights*states;
    
end

%weights calculation
function weight_mat =  insert_memory(M)

    for i = 1:100
        for j = 1:100
        
            if i ~= j
                weight_mat(i,j) = sum(M(:,i).*M(:,j));
            else 
                weight_mat(i,j) = 0;
            end
        end
    end
    
end

% Probabilistic classification

function s = prob_classify(state_energy)
    T = 4;
    s = 1/(1 + exp((-state_energy/T)));
    s(s>0.5) = 1;
    s(s<=0.5) = -1;

end