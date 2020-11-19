function [bus, branch] = c5

% bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin

bus = [1 3 1 1; 
        2 2 1 1;  
        3 1 1 1; 
        4 1 1 1; 
        5 1 1 1];
    
branch =   [1 2 0.02 0.06;
            1 3 0.08 0.24;
            2 3 0.06 0.25;
            2 4 0.06 0.18;
            2 5 0.04 0.12;
            3 4 0.01 0.03;
            4 5 0.08 0.24];
            
return 