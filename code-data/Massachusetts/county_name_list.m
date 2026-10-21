function [state, m, name_origin_list, time_start, time_end, top_five_pop_flow_nodes, time_step_size , zeta_vec, zeta,  Nv, gamma, beta_1_vec, phi_0] = county_name_list()

%% variable parameters 

state = 25;

m = 14;

name_origin_list = {'Barnstable County', 'Berkshire County', 'Bristol County', 'Dukes County', 'Essex County', 'Franklin County', 'Hampden County', 'Hampshire County', 'Middlesex County', 'Nantucket County', 'Norfolk County', 'Plymouth County', 'Suffolk County', 'Worcester County'};

time_start = 1;

time_end = 84;

top_five_pop_flow_nodes = [9,13,14,5,11];

%% constant parameters which we may change 

time_step_size = 0.10;                                                     % can compute N_time_steps

zeta_vec = [5, 6, 7, 8, 9, 10];                                            % put the diferent zeta values 
             
zeta = zeta_vec(1) ; % Dimensionless multiplicative factor that modifies the infectivity at the heterogeneous node.

Nv = 3;

gamma = 0.1429;                                                            % 1/7 wekly recovered

beta_1_vec = linspace(0,0.1428,1);                                         % simulate the infection rates

phi_0 = 0.004;                                                             % mean flux value in the first 81 days i.e. 0.006795

end