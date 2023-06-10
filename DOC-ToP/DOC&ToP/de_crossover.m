function [child1, child2] = de_crossover(parent1, parent2,parent3, n, lu, pcrossover, de_factor)

parent_1 = parent1(1 : n);
parent_2 = parent2(1 : n);
parent_3 = parent3(1 : n);

% pcrossover=0.9 is the probability of crossover
cross_idx1 = rand(1,n) < pcrossover;
cross_idx1(randi(n)) = 1; % at least one decision component to crossover

child1 = parent_1 + de_factor * (parent_3 - parent_2);
child1(~cross_idx1) = parent_1(~cross_idx1);


cross_idx2 = rand(1,n) < pcrossover;
cross_idx2(randi(n)) = 1; % at least one decision component to crossover
child2 = parent_2 + de_factor * (parent_3 - parent_1);
child2(~cross_idx2) = parent_2(~cross_idx2);

rand_value = rand;
if rand_value < 0.5
    theta = (2*rand_value)^(1/(1+20)) - 1;
else
    theta = 1 - (2 - 2*rand_value)^(1/(1+20));
end
rand_mut = rand(n,1);
muted = rand_mut < 1/n;
y = child1 + theta*(lu(2,:) -lu(1,:));
y(~muted) = child1(~muted);
child1 = y;


rand_value = rand;
if rand_value < 0.5
    theta = (2*rand_value)^(1/(1+20)) - 1;
else
    theta = 1 - (2 - 2*rand_value)^(1/(1+20));
end
rand_mut = rand(n,1);
muted = rand_mut < 1/n;
y = child2+ theta*((lu(2,:) -lu(1,:)));
y(~muted) = child2(~muted);
child2 = y;
% repair the two children

child1 =  max(min(lu(2,:),child1),lu(1,:));
child2 =  max(min(lu(2,:),child2),lu(1,:));
end
