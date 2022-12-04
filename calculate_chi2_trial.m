% calculate chi new version
function [chi2] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
    platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
    Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y)

PMs2 = zeros(5,Nruns); 
for rep = 1:Nruns
    %Generate initial weights for each run
    %weights = weights1;

    %Generate place cells for each run
%     PC_x = zeros(1,N_pc); %1xN_pc matrix containing the x = 0 coordinate for each place cell
%     PC_y = zeros(1,N_pc); %1xN_pc matrix containing the y = 0 coordinate for each place cell
%     for i = 1:N_pc %For each place cell:
%         PC_x(i) = (rand - 0.5)*pool_diameter; %Random positions of place cells
%         PC_y(i) = (rand - 0.5)*pool_diameter;
%         while (PC_x(i)^2 + PC_y(i)^2 > (pool_diameter/2)^2) %Checks for out of bounds
%             PC_x(i) = (rand - 0.5)*pool_diameter;
%             PC_y(i) = (rand - 0.5)*pool_diameter;
%         end
%     end



            idx = randi(4); %randomly choose one of 4 starting locations
            starting_x = starting_xs(idx);
            starting_y = starting_ys(idx);

             [wres, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants, ...
            latency, speed_std,speed_ps, mean_angle,time_step] = ...
            run_trial (weights, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, ...
            Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, ...
            platform_x, platform_y, starting_x, starting_y, speed, hitwall);
            %run trial
            weights = wres;

            PMs2(1,set,rep) = latency;
            PMs2(2,set,rep) = dist*100;
            PMs2(3,set,rep) = quadrants(4);
            PMs2(4,set,rep) = quadrants(2);
            PMs2(5,set,rep) = wall_zone;
            % new parameters
            PMs2(6,day,trial) = speed_std*100;
            PMs2(7,day,trial) = mean_angle;
            PMs2(8,day,trial) = time_step;
           
            %record performance measures
         
 end
    PMmod= mean(PMs2,2);
    % a= zeros(5,1);
    chi2=0;
    for i= 1:5
        a = (PMexp(i,6)-PMmod(i)).^2 / PMexp(i,7).^2;
        chi2= chi2+a;
    end
    
    chi2=chi2;
    
% Calculate Chi_Old version


% function [chi2] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
%     platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
%     Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y)
% 
% PMs2 = zeros(5,Nruns); 
% for rep = 1:Nruns
%     %Generate initial weights for each run
%     %weights = weights1;
% 
%     %Generate place cells for each run
% %     PC_x = zeros(1,N_pc); %1xN_pc matrix containing the x = 0 coordinate for each place cell
% %     PC_y = zeros(1,N_pc); %1xN_pc matrix containing the y = 0 coordinate for each place cell
% %     for i = 1:N_pc %For each place cell:
% %         PC_x(i) = (rand - 0.5)*pool_diameter; %Random positions of place cells
% %         PC_y(i) = (rand - 0.5)*pool_diameter;
% %         while (PC_x(i)^2 + PC_y(i)^2 > (pool_diameter/2)^2) %Checks for out of bounds
% %             PC_x(i) = (rand - 0.5)*pool_diameter;
% %             PC_y(i) = (rand - 0.5)*pool_diameter;
% %         end
% %     end
% 
% 
% 
%             idx = randi(4); %randomly choose one of 4 starting locations
%             starting_x = starting_xs(idx);
%             starting_y = starting_ys(idx);
% 
%             [wres, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants, latency] = ...
%             run_trial (weights1, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, ...
%                 Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, ...
%                 platform_x, platform_y, starting_x, starting_y, speed, hitwall);
%             %run trial
%             
% 
%             PMs2(1,rep) = latency;
%             PMs2(2,rep) = dist*100;
%             PMs2(3,rep) = quadrants(4);
%             PMs2(4,rep) = quadrants(2);
%             PMs2(5,rep) = wall_zone;
%             %record performance measures
%         
%     end
%     PMmod= mean(PMs2,2);
%     % a= zeros(5,1);
%     chi2=0;
%     for i= 1:5
%         a = (PMexp(i,6)-PMmod(i)).^2 / PMexp(i,7).^2;
%         chi2= chi2+a;
%     end
%     
%     chi2=chi2;
%     
