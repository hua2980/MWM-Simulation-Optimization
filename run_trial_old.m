% takes weights, PC coordinates & parameters and returns new weights and performance measures (single trial)
function [weights, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants, latency] = ...
    run_trial_old (weights0, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, Vdecay, ac_const, beta, etdecay, ...
    lrate, discf, noise, platform_x, platform_y, starting_x, starting_y, speed, wall_pun)
% speed in meters per second (!) taken as input - latency computed accordingly

% FIXED PARAMETRES OF THE EXPERIMENT

pool_diameter = 1.4; %Maze diameter (\phi) in metres (m)
platform_radius = 0.06; %Platform radius

N_pc = 211; %Population of place cells
N_ac = 36; %Population of action cells
which = 0;

dist = 0;
wall_zone = 0;
quadrants = zeros(4,1); %Percentage spent on each quadrant 

weights = weights0; %Initialize modifiable weights 

el_tr = zeros(N_pc, N_ac); %Initialize eligibility traces matrix

%Initialize trajectories
track_x = starting_x; %Current position of trajectory is equal to the 
                      %starting location of the animal
track_y = starting_y;
vel_x = 0;
vel_y = 0;

% NAVIGATION LOOP
while ((track_x(end) - platform_x)^2 + (track_y(end) - platform_y)^2 > ... 
    platform_radius^2)

    weights = weights*(1-noise) + rand(N_pc, N_ac)*Wmult*noise;
    
    %Calculate PC activation 
    PC_activation = zeros(1,N_pc); 
    for i = 1:N_pc
        PC_activation(i) = exp(-((track_x(end) - PC_x(i))^2 + (track_y(end)...
        - PC_y(i))^2)/(2*sigma_pc^2));
    end
    
    %Calculate AC activation (i.e. value of the action, Q)
    if (length(track_x) > 1)
        prevQ = AC_activation(which); %Displays the Q value before movement
    end
    
    AC_activation = PC_activation * weights; 
    
    %Make an action
    ACsel = AC_activation.^beta;
    ACsel = ACsel / sum(ACsel);
    ASrand = rand;
    which = 1; ASsum = ACsel(1);
    while (which < N_ac && ASsum < ASrand)
        which = which + 1;
        ASsum = ASsum + ACsel(which);
    end
    
    %Eligibility traces
    el_tr = el_tr * etdecay;
    for j = 1:N_ac
        itmp = min(abs(j-which),N_ac-abs(j-which));
        actgaus = exp(-(itmp*itmp)/(2*sigma_ac*sigma_ac));
        el_tr(:,j) = el_tr(:,j) + actgaus*AC_activation(j)*(PC_activation');
    end
    
    vel_x = [vel_x (vel_x(end)+ac_const*cos(which/N_ac*2*pi))*Vdecay];
    vel_y = [vel_y (vel_y(end)+ac_const*sin(which/N_ac*2*pi))*Vdecay];
    %velocity per time step (not second)
    track_x = [track_x track_x(end)+vel_x(end)];
    track_y = [track_y track_y(end)+vel_y(end)];
    
    %Check if not out of bounds, reset location & speed if so
    if (track_x(end)^2 + track_y(end)^2 > (pool_diameter/2)^2)
        ratio = (track_x(end)^2 + track_y(end)^2)/((pool_diameter/2)^2);
        track_x(end) = track_x(end)/sqrt(ratio);
        track_y(end) = track_y(end)/sqrt(ratio);
        vel_x(end) = track_x(end) - track_x(end-1);
        vel_y(end) = track_y(end) - track_y(end-1);
    end
    
    if (length(track_x) > 2)
       if ((track_x(end) - platform_x)^2 + (track_y(end) - platform_y)^2 < ... 
    platform_radius^2)
           rew = 10; %found platform - reward
       elseif (track_x(end)^2+track_y(end)^2 > (0.99*pool_diameter/2)^2)
           rew = -wall_pun; %hit wall - punishment
       else
           rew = 0; %didn't find - no reward
       end
       currQ = AC_activation(which);
       %disp(['Qs: ',num2str(currQ),' - ',num2str(prevQ)]);
       tderr = rew + discf*currQ - prevQ; %temporal difference error
       weights = max(0,weights + lrate*tderr*el_tr);
%         for i = 1:N_ac
%             for j = 1:N_pc
%                 weights(j,i) = max(weights(j,i) + lrate*tderr*el_tr(j,i),0);
%             end
%         end
    end
    %disp([num2str(track_x(end)),',',num2str(track_y(end))]); %display the 
    %coordinates of each action cell within the loop
    laststep = sqrt((track_x(end)-track_x(end-1))^2 + (track_y(end)-track_y...
        (end-1))^2);
    dist = dist + laststep;
    if (track_x(end)^2 + track_y(end)^2 > 0.8*(pool_diameter/2)^2)
        wall_zone = wall_zone + 1;
    elseif (track_x(end) > 0 && track_y(end) > 0)
        quadrants(1) = quadrants(1) + 1;
    elseif (track_x(end) < 0 && track_y(end) > 0)
        quadrants(2) = quadrants(2) + 1;
    elseif (track_x(end) < 0 && track_y(end) < 0)
        quadrants(3) = quadrants(3) + 1;
    else
        quadrants(4) = quadrants(4) + 1;
    end
    if (length(track_x) > 100) % evaluate latency after 100+ steps
        speed_ts = mean((vel_x(2:end).^2+vel_y(2:end).^2).^0.5); % speed in meters/time step
        latency = (length(track_x)-1) * speed_ts / speed; % convert to seconds
        if (latency > 60) % if more than a minute, stop
            break;
        end
    end
end

latency = length(track_x)-1;  % latency in time steps
wall_zone = wall_zone/latency;
quadrants = quadrants/latency;
speed_ts = mean((vel_x(2:end).^2+vel_y(2:end).^2).^0.5); % speed in meters/time step

latency = latency * speed_ts / speed; % latency in seconds