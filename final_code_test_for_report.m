clc
clear

% [bus, branch] = c5;
% [baseMVA, bus, gen, branch, areas, gencost] = c6;
% [baseMVA, bus, gen, branch, areas, gencost] = case14;
% [baseMVA, bus, gen, branch, areas, gencost] = case_ieee30;
% [baseMVA, bus, gen, branch, areas, gencost] = case33bw;
% [bus, branch] = c39;
[baseMVA, bus, gen, branch, areas, gencost] = c57;
% [baseMVA, bus, gen, branch, areas, gencost] = case118;
% [baseMVA, bus, gen, branch, areas, gencost] = case300;

% define size of file
% x = [number of lines/branches]
[line_max,~] = size(branch);
branch(:,4) = .15;

% define number of buses
[size_bus,~] = size(bus);

% save original bus numbers
bus_original = bus(:,1);

% check for parallel lines
% calculate reactance for parallel lines
% adjust size of 

n = 0;

% check for parallel lines
% calculate reactance for parallel lines
% set one parallel line to equivalent parallel reactance
% set other parallel line to zero
for i = 1:line_max
    node1 = branch(i,1);
    node2 = branch(i,2);
    if node1 > 0
        for j = i+1:line_max
            node1_check = branch(j,1);
            node2_check = branch(j,2);
            if node1 == node1_check && node2 == node2_check
                n = n + 1;
                branch(i,4) = 1/(1/branch(i,4)+1/branch(j,4));
                branch(j:line_max-1,:) = branch(j+1:line_max,:);
                branch(line_max-n+1,:) = 0;
            elseif node1 == node2_check && node2 == node1_check
                n = n + 1;
                branch(i,4) = 1/(1/branch(i,4)+1/branch(j,4));
                branch(j:line_max-1,:) = branch(j+1:line_max,:);
                branch(line_max-n+1,:) = 0;
            end
        end
    end
end

% remove extra rows from branch matrix
if n > 0 
    branch(line_max-n+1:line_max,:) = [];
end

% create list of new bus numbers
for i = 1:size_bus
    bus(i,1) = i;
end

[line_max,~] = size(branch);

branch_original(:,1) = branch(:,1);
branch_original(:,2) = branch(:,2);


for i = 1:line_max
    [a,~] = find(bus_original(:,1) == branch(i,1));
    [b,~] = find(bus_original(:,1) == branch(i,2));
    branch(i,1) = a;
    branch(i,2) = b;
end

l_new = zeros(line_max,4);
line_susceptance = zeros(line_max,1);

for k = 1:line_max
    %node i
    l_new(k,1) = branch(k,1);
    
    %node j
    l_new(k,2) = branch(k,2);
    
    %reactance (X) between i and j
    %l_new(k,3) = branch(k,3) + branch(k,4)*i;
    l_new(k,3) = branch(k,4);
    
    %susceptance (B) between i and j
    l_new(k,4) = 1/l_new(k,3);
end

line_susceptance = l_new(:,4);

% calculate Bd (susceptance) matrix
% B = 1 / X 
% susceptance = 1 / reactance
% Bd is a diagonal matrix
% Bd is a L x L matrix
% x = [number of lines]

Bd = zeros(line_max);

for i = 1:line_max
    for j = 1:line_max
        if i == j
            Bd(i,j) = line_susceptance(i);
        else
            Bd(i,j) = 0;
        end
    end
end

%determine number of nodes in grid
% node_max(1) = max(l_new(:,1));
% node_max(2) = max(l_new(:,2));
% node_max = max(node_max);

[node_max, ~] = size(bus);

% give each line a number
% determine number of lines in grid

line = zeros(line_max,3);

for i = 1:line_max
    line(i,1) = i;
    line(i,2) = l_new(i,1);
    line(i,3) = l_new(i,2);
%   line_max = i;
end


%incidence matrix formation [ line x nodes ]
% 1 is for "from node"
% -1 is for "to node"
% 0 is for no connection


% [line_max,y] = size(branch);
A = zeros(line_max,node_max);

for i = 1:line_max
    [a,~] = find(bus(:,1)==line(i,2));
    [b,~] = find(bus(:,1)==line(i,3));
    A(i,a) = 1;
    A(i,b) = -1;
end

% for i = 1:line_max
%     bus(i,1) = i;
% end

% define Pl and Pn
% PTDF = Pl * Pn^(-1)

Pl = Bd*A;
Pn = A'*Bd*A;

max_node_type = 0;

for i = 1:node_max
    if bus(i,2) > max_node_type
        slack_bus = bus(i,1);
        max_node_type = bus(i,2);
    end
end

% define slack bus
% use node type data to determine slack bus
% load nodes = 1
% generator nodes = 2
% slack nodes = 3
% max_node_type = 0;
% for i = 1:node_max
%     if bus(i,2) > max_node_type
%         slack_bus = bus(i,1);
%         max_node_type = bus(i,2);
%     end
% end

% remove one column from Pl
% this column is the slack bus
Pl(:,slack_bus) = 0;

% remove one column and row from Pn
% the column and row are for the slack bus
Pn(:,slack_bus) = 0;
Pn(slack_bus,:) = 0;

% create matrix of zeros
% these matrices will be used to make new Pn and Pl matrices
Pn_new = zeros(node_max-1);
Pl_new = zeros(line_max,node_max-1);

% determine if slack bus is 1, equal to node max or neither
% this if statement will adjust Pn and Pl depending on location of slack bus  
if slack_bus == 1
    Pn_new = Pn(2:node_max,2:node_max);
    
    Pl_new = Pl(:,2:node_max);
elseif slack_bus == node_max
    Pn_new = Pn(1:node_max-1,1:node_max-1);
    
    Pl_new = Pl(:,1:node_max-1);
else
    Pn_new(1:slack_bus-1,slack_bus:node_max-1) = Pn(1:slack_bus-1,slack_bus+1:node_max);
    Pn_new(slack_bus:node_max-1,slack_bus:node_max-1) = Pn(slack_bus+1:node_max,slack_bus+1:node_max);
    Pn_new(slack_bus:node_max-1,1:slack_bus-1) = Pn(slack_bus+1:node_max,1:slack_bus-1);
    Pn_new(1:slack_bus-1,1:slack_bus-1) = Pn(1:slack_bus-1,1:slack_bus-1);
    Pl_new(:,slack_bus:node_max-1) = Pl(:,slack_bus+1:node_max);
    Pl_new(:,1:slack_bus-1) = Pl(:,1:slack_bus-1);
end

% find inverse of Pn
Pn_new = pinv(Pn_new);
% Pn_new_inv = inv(Pn_new);
% Pn_new_pinv = pinv(Pn_new);
% Pn_new_trad = Pn_new^(-1)

% calculate PTDF
ptdf = Pl_new*(Pn_new);

% add reference node to PTDF matrix

ptdf_adj = zeros(line_max,node_max);

if slack_bus == 1
    ptdf_adj(:,2:node_max) = ptdf;
elseif slack_bus == node_max
    ptdf_adj(:,1:node_max-1) = ptdf;
else
    ptdf_adj(:,1:slack_bus-1) = ptdf(:,1:slack_bus-1);
    ptdf_adj(:,slack_bus+1:node_max) = ptdf(:,slack_bus:node_max-1);
end

% % % % disp(ptdf_adj);

% determine generator and load nodes
t = 0;
l = 0;
g = 0;

for i = 1:node_max
    if bus(i,2) == 1
        if bus(i,3) == 0 && bus(i,4) == 0
            t = t + 1;
            tran_list(t) = bus(i,1);
        else
            l = l + 1;
            load_list(l) = bus(i,1);
        end
    else
        g = g + 1;
        gen_list(g) = bus(i,1);
    end
end

% make distribution factor (DF) matrix 
% DF is a [ node x node x line ] matrix
% a_ij = a_ig - a_id

size_gen = size(gen_list');
size_load = size(load_list');
a = zeros(node_max,node_max,line_max);

% % % % % % % why was line_max used

for i = 1:node_max
    for j = 1:node_max
        for l = 1:line_max
            a(i,j,l) = ptdf_adj(l,i) - ptdf_adj(l,j);
        end
    end
end

% calculate the capacity (C) of the system
% matrix is [ node x node ]
% find the maximum PTDF of each gen-load combination
% take the reciprocal of maximum value to find capacity of each gd combo

a = abs(a);
C = zeros(node_max,node_max);
p_max = 100;

for i = 1:node_max
    for j = 1:node_max
        C(i,j) = p_max/max(a(i,j,:));
    end
end

% calculate Ce [ gen x load x node ] matrix

Ce = zeros(size_gen(1),size_load(1),node_max);

for g = 1:size_gen(1)
    for d = 1:size_load(1)
        for e = 1:node_max
            gg = gen_list(g);
            dd = load_list(d);
            if gg == e
                Ce(g,d,e) = 0;
            elseif dd == e
                Ce(g,d,e) = 0;
            else
                Ce(g,d,e) = 1 / (1/C(gg,e) + 1/C(e,dd));
            end
        end
    end
end

% calculate the admittance matrix, then calculate the impedance matrix

adm = zeros(node_max);

for i = 1:node_max
    for j = 1:node_max
        
        %calculate admittances when i ≠ j
        if i ~= j
            for a = 1:line_max
                if l_new(a,1) == i
                    if l_new(a,2) == j
                        adm(i,j) = -l_new(a,4);
                        adm(j,i) = -l_new(a,4);
                    end
                end
            end
        end
        
        %calculate admittances when i = j (diagonal matrix)
        if i == j
            for a = 1:line_max
                if l_new(a,1) == i || l_new(a,2) == i
                    adm(i,i) = l_new(a,4) + adm(i,i);
                end
            end
        end      
    end
end

% removal of last node in order to calculate inverse matrix
% for i = 1:node_max-1
%     for j = 1:node_max-1
%         adjusted_adm(i,j) = adm(i,j);
%     end
% end

% imp = inv(adjusted_adm);

Z = pinv(adm);

% calculate equivalent impedance matrix
% first calculate impedance matrix
% equivalent impedance (Zgd) equals Zgg - 2*Zgd - Zdd

for i = 1:node_max
    for j = 1:node_max
        Zequ(i,j) = Z(i,i) - 2*Z(i,j) + Z(j,j);
    end
end

% calculate Ze [ gen x load x node ] matrix

Ze = zeros(size_gen(1),size_load(1),node_max);

for g = 1:size_gen(1)
    for d = 1:size_load(1)
        for i = 1:node_max
            gg = gen_list(g);
            dd = load_list(d);
            Ze(g,d,i) = 1/(Zequ(gg,i) + Zequ(i,dd));
        end
    end
end

% na_sum(:,1) = 1:node_max;
na_sum(:,1) = bus_original;

na_sum(:,2) = zeros(node_max,1);

for e = 1:node_max
    for g = 1:size_gen(1)
        for d = 1:size_load(1)
            na_sum(e,2) = na_sum(e,2) + Ze(g,d,e) * Ce(g,d,e);
        end
    end
end

% calculate original net ability

grid_net_ability = 0;

for g = 1:size_gen(1)
    for d = 1:size_load(1)
        gg = gen_list(g);
        dd = load_list(d);
        grid_net_ability = grid_net_ability + C(gg,dd) / Zequ(gg,dd);
    end
end

grid_net_ability = grid_net_ability/(size_gen(1)*size_load(1));

display(grid_net_ability)

na_sum(:,2) = na_sum(:,2)/(size_gen(1)*size_load(1)); 

% + grid_net_ability;

sum_net_ability_direct_connections = zeros(size_bus,1);

for i = 1:size_bus
    for j = 1:line_max
        if bus(i,1) == branch(j,1) 
            sum_net_ability_direct_connections(i) = sum_net_ability_direct_connections(i) + na_sum(branch(j,2),2);
        elseif bus(i,1) == branch(j,2)
            sum_net_ability_direct_connections(i) = sum_net_ability_direct_connections(i) + na_sum(branch(j,1),2);
        end
    end
end

display(sum_net_ability_direct_connections)

size(na_sum)
size(sum_net_ability_direct_connections)
% scatter(na_sum(:,2),sum_net_ability_direct_connections)


net_ability_order = na_sum(:,2);

net_ability = sortrows(na_sum, -2);
% net_ability(:,1) = fix(net_ability(:,1));

for i = 1:size_bus
    fprintf('%6.0f %12.1f\n',net_ability(i,1),net_ability(i,2))
end



% for i = 1:10
%     fprintf('%0.0f\n',net_ability(i,1))
% end
% 
% fprintf('\n')
% 
% for i = 1:10
%     fprintf('%0.1f\n',net_ability(i,2))
% end
% 
% fprintf('\n')
% 
% for i = node_max-9:node_max
%     fprintf('%0.0f\n',net_ability(i,1))
% end
% 
% fprintf('\n')
% 
% for i = node_max-9:node_max
%     fprintf('%0.1f\n',net_ability(i,2))
% end

sum_connections_primary = zeros(size_bus,1);

for i = 1:size_bus
    for j = 1:line_max
        if bus(i,1) == branch(j,1) 
            sum_connections_primary(i) = sum_connections_primary(i) + 1;
        elseif bus(i,1) == branch(j,2)
            sum_connections_primary(i) = sum_connections_primary(i) + 1;
        end
    end
end

sum_connections_secondary = zeros(size_bus,1);

for i = 1:size_bus
    for j = 1:line_max
        if bus(i,1) == branch(j,1)
            k = branch(j,2);
            sum_connections_secondary(i) = sum_connections_primary(k) + sum_connections_secondary(i);
        elseif bus(i,1) == branch(j,2)
            k = branch(j,1);
            sum_connections_secondary(i) = sum_connections_primary(k) + sum_connections_secondary(i);
        end
    end
end

n_degree = 1;

x1 = sum_net_ability_direct_connections;
y = na_sum(:,2);

mdl1 = fitlm(x1,y);
R_squared_1 = mdl1.Rsquared.Ordinary;
p1 = polyfit(x1,y,n_degree);

data_points1 = min(x1):(max(x1)-min(x1))/100:max(x1);

if n_degree == 1
    lbf1 = p1(1).*data_points1 + p1(2);
elseif n_degree == 2
    lbf1 = p1(1).*data_points1.^2 + p1(2).*data_points1 + p1(3);
elseif n_degree == 3
    lbf1 = p1(1).*data_points1.^3 + p1(2).*data_points1.^2 + p1(3).*data_points1 + p1(4);
end

q = 1e2;
R_squared_1 = round(R_squared_1*q)/q;

% figure('units','normalized','outerposition',[0.1 1 0.4 0.6])
% plot(x1,y,'o',data_points1,lbf1,'r-','Linewidth',1.5)
legend({'Buses',['Regression Line, R^2 =  ' num2str(R_squared_1)]}, 'Location', 'Northwest','FontSize', 12)
title('IEEE 300 - ES Net Ability vs. Sum of ES Net Ability for Directly Connected Buses','FontSize', 15)
set(gca,'FontSize',12)
xlabel('Sum of ES Net Ability for Directly Connected Buses ( MW/\Omega )','FontSize', 15)
ylabel('ES Net Ability ( MW/\Omega )','FontSize', 15)
% xlim([0 max(x1)+1])
% ylim([0 280])



n_degree = 1;
x1 = sum_connections_primary;
% x2 = sum_connections_secondary;
y = net_ability_order;

% x1 = sum_net_ability_direct_connections;
% y = na_sum(:,2);
% 
mdl1 = fitlm(x1,y);
R_squared_1 = mdl1.Rsquared.Ordinary;
% mdl2 = fitlm(x2,y);
% R_squared_2 = mdl2.Rsquared.Ordinary;
p1 = polyfit(x1,y,n_degree);
% p2 = polyfit(x2,y,n_degree);



data_points1 = min(x1):(max(x1)-min(x1))/100:max(x1);
% data_points2 = min(x2):(max(x2)-min(x2))/100:max(x2);

if n_degree == 1
    lbf1 = p1(1).*data_points1 + p1(2);
%     lbf2 = p2(1).*data_points2 + p2(2);
elseif n_degree == 2
    lbf1 = p1(1).*data_points1.^2 + p1(2).*data_points1 + p1(3);
%     lbf2 = p2(1).*data_points2.^2 + p2(2).*data_points2 + p2(3);
elseif n_degree == 3
    lbf1 = p1(1).*data_points1.^3 + p1(2).*data_points1.^2 + p1(3).*data_points1 + p1(4);
%     lbf2 = p2(1).*data_points2.^3 + p2(2).*data_points2.^2 + p2(3).*data_points2 + p2(4);
end

q = 1e2;
R_squared_1 = round(R_squared_1*q)/q;
% R_squared_2 = round(R_squared_2*q)/q;

figure('units','normalized','outerposition',[0.1 1 0.4 0.6])
% subplot(1,2,1)
% hold on
plot(x1,y,'o',data_points1,lbf1,'r-','Linewidth',1.5)
legend({'Buses',['Regression Line, R^2 =  ' num2str(R_squared_1)]}, 'Location', 'northwest','FontSize', 12)
title({'IEEE 300 - Identical Line Impedances','ES Net Ability vs. Connections'},'FontSize', 15)
set(gca,'FontSize',12)
xlabel('Number of Connections','FontSize', 15)
ylabel('ES Net Ability ( MW/\Omega )','FontSize', 15)
xlim([0 max(x1)+1])
% ylim([0 350])
xticks([0:1:(max(x1)+1)])

% sum_connections_primary = [1 2 3 4 5 6 1 2 3 2 5 6 7]';

[x,~] = size(sum_connections_primary);

number = zeros(max(sum_connections_primary),1);

for i = 1:x
    test = sum_connections_primary(i);
    number(test) = number(test) + 1;
end

list(:,2) = number./x;
list(:,1) = 1:max(sum_connections_primary);

sum(list(:,2))

list;

% subplot(1,2,2)
% plot(x2,y,'o',data_points2,lbf2)
% legend('Buses',['Regression Line, R^2 =  ' num2str(R_squared_2)]','Location', 'Northwest')
% title('Primary and Secondary Connections')
% xlabel('Number of Connections')
% ylabel('Net Ability ( MW/\Omega )')

average_net_ability = sum(net_ability(:,2))/node_max