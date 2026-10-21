function td_matrix_phi = TDFluxMatrix(N_o,N_d,time_step_size,time_start,time_end,name_origin_list);

m = N_o;
tq = time_start:time_step_size:time_end;                                   % query points                                                                                                                                                 % temperature query points

time = time_start:1:time_end;

N_time_steps = ((time_end-time_start)/time_step_size)+1;

% Define the flux functions for each origin-destination pair

td_matrix_phi = zeros(N_o, N_d, N_time_steps);                             % preallocate for time dependent flux

% o is the origin
% d is the orifgin


for o = 1:m
    
    origin_to_destination = strcat(name_origin_list(o), '_Mobility_Flux');
    
    origin_to_destination_csv = strcat(name_origin_list(o), '_Mobility_Flux','.csv');
    
    new_origin_to_destination = csvread(string(origin_to_destination_csv) ,1,0);
    
    origin_to_destination_new = categorical(origin_to_destination);
    
    origin_to_destination = origin_to_destination_new;
    
    origin_to_destination = new_origin_to_destination;
    
    for d = 1:m                   % since d=1 is the week numbers so we don't need to read this data
        % Load the flux data for origin o and destination d for all 52 weeks
        %flux_data = load_flux_data(o, d); % assume there is a function to load data for each pair
        %time=Blaine_77(:,1)
        o_d = origin_to_destination(time_start:time_end,d);
        flux_data = interp1(time,o_d, tq);
        %Store the flux data in the 3D array phi
        td_matrix_phi(o,d, :) = flux_data;
    end
end


%% Save the time dependent flux matrix

%writematrix(td_matrix_phi,'td_matrix_phi-od.csv');  % export a matrix as a CSV file


end