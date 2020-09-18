%hopfield
clear all;
close all;
clc;
choose = [1 2 3];
degrade = [0 10 20 30 40 50];
for chosen =1: size(choose,2)
for degrader=1: size(degrade,2)
    
    for runs = 1:500
    
    capacity = 3;

    %Weight Matrix
    weight_mat = rand(100,100); % preallocation of matrix size

    % Storage Stage
    memories = randi([-1,1],capacity,100);
    memories(memories==0) = 1;
    memories;
    weight_matrix = insert_memory(memories);

    % Retreival Stage

    % creating the degraded input
    test_signal = memories(choose(chosen),:);

    %degrading the pattern
    
     for i = 1:degrade(degrader)
         if test_signal(i) == 1
             test_signal(i) = -1;
         else
             test_signal(i) = 1;
         end
     end

    %transposing for multiplication
    test_signal = test_signal';

    % updating the nodes 
stop =0;
        while stop ~= 1
            
            temp = test_signal;
            
            for i = 1:size(weight_matrix,1)
                
                energies(i) = weight_matrix(i,:)*test_signal; % or weight_matrix for nor normal people
                test_signal(i) = classify(energies(i)); % or classify for sign func
                
            end
            
            if test_signal-temp == 0
                stop = 1;
            end

        end

hammer = test_signal - memories(choose(chosen),:)';
ham =find(hammer~= 0);
hammering = size(ham,1);
%end
hammered(degrader,runs,chosen) = mean(hammering);
    end
end
end

for i = 1: size(choose,2)
graph_mat(:,1) = mean(hammered(:,:,chosen),2);
x = degrade;
figure(i)
plot(x,graph_mat(:,1),'g--o');
ylim([-2,5]);
xlabel('Degradation Level in percentages');
ylabel('Mean Hamming Distance');
title('Recall For Memory 1');

end



% Local Functions

% State classification function
function s = classify(energies)

    s = energies;
    s(s>0) = 1;
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