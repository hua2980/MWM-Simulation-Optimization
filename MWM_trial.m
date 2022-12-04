%expmean = readmatrix('expmean.csv');
%expmean = readmatrix('PMexp_trial.csv');
expmean = csvread('PMexp_trial.csv');
Strain=1;
Temperature=18;

N_pc = 211; %Population of place cells [100..300]
N_ac = 36; %Population of action cells [25..50]

plot_trajectories = 0; %yes - 1, no - 0
plot_cognitive_maps = 0; %yes - 1, no - 0
pln = plot_trajectories + plot_cognitive_maps;
Nruns = 40; %how many runs to run if not plotting anything
Nsets = 50; % take 50 times random alpha, beta,...
Nmice = 1; % estimate the number of mice

pool_diameter = 1.4; %Maze diameter (\phi) in metres (m)
platform_radius = 0.06; %Platform radius (m)
sigma_pc = 0.1;  %place cell sigma (standard deviation), in meters [0.05..0.2]
sigma_ac = 2;  %action cell sigma (standard deviation), in action cells [1..3]

% etdecay = 0.83; %Eligibility trace decay (lambda) [0.75..0.95] LESS THAN GAMMA!
% beta = 6;  %Exploration-exploitation factor (\beta) [0.5..12]
% alpha = 0.01;  %Learning rate (\alpha) [0.005..0.02]
% gamma = 0.85;  %Discount factor (\gamma) [0.75..0.95]

Vdecay = 0.82;  %velocity decay [0.75..0.95]
ac_const = 0.02;  %acceleration const [0.01..0.03]
Wnoise = 0.0004;  %Weight noise [0.0001..0.0007]
Wmult = 0.1;   %Weight multiplier [0.05..0.15]
hitwall = 5;  %punishment for hitting the wall [0..1]
speed = 0.175;  %mouse speed (m/s) [0.1..0.25]

Ntrials = 4;  %number of trials per day
Ndays = 5;  %number of days

PMs = zeros(5,Ndays,Ntrials,Nmice);  %multiple runs

%Platform coordinates:
platform_x = cos(-pi/4)*pool_diameter/4; %x coordinate
platform_y = sin(-pi/4)*pool_diameter/4; %y coordinate

%Starting locations of the modeled animal (4 different ones):
strad = pool_diameter/2*0.90; %15% of maze radius to the wall
starting_xs = strad * [cos(pi/6) cos(pi/3) cos(7*pi/6) cos(4*pi/3)]; %x coordinates
starting_ys = strad * [sin(pi/6) sin(pi/3) sin(7*pi/6) sin(4*pi/3)]; %y coordinates

th = 0:pi/50:2*pi; %for plotting circles :)



% run multiple times without plotting!

allpars=[];
% allweight = zeros(N_pc, N_ac, Nmice);
% allPC_x = zeros(Nmice,N_pc);

%%

for mice = 1:Nmice
    %Generate initial weights for each run
    weights1 = rand(N_pc,N_ac)*Wmult;

    %Generate place cells for each run
    PC_x = zeros(1,N_pc); %1xN_pc matrix containing the x = 0 coordinate for each place cell
    PC_y = zeros(1,N_pc); %1xN_pc matrix containing the y = 0 coordinate for each place cell
    for i = 1:N_pc %For each place cell:
        PC_x(i) = (rand - 0.5)*pool_diameter; %Random positions of place cells
        PC_y(i) = (rand - 0.5)*pool_diameter;
        while (PC_x(i)^2 + PC_y(i)^2 > (pool_diameter/2)^2) %Checks for out of bounds
            PC_x(i) = (rand - 0.5)*pool_diameter;
            PC_y(i) = (rand - 0.5)*pool_diameter;
        end
    end
   

    
    for day = 1:Ndays 
        day
%    for rep=1:Nmice
%         %Generate initial weights for each run
%         if day==1
%             weights = rand(N_pc,N_ac)*Wmult;
%             allweight(N_pc,N_ac,rep)= weights;
% 
%             %Generate place cells for each run
%             PC_x = zeros(1,N_pc); %1xN_pc matrix containing the x = 0 coordinate for each place cell
%             PC_y = zeros(1,N_pc); %1xN_pc matrix containing the y = 0 coordinate for each place cell
%             for i = 1:N_pc %For each place cell:
%                 PC_x(i) = (rand - 0.5)*pool_diameter; %Random positions of place cells
%                 PC_y(i) = (rand - 0.5)*pool_diameter;
%                 while (PC_x(i)^2 + PC_y(i)^2 > (pool_diameter/2)^2) %Checks for out of bounds
%                     PC_x(i) = (rand - 0.5)*pool_diameter;
%                     PC_y(i) = (rand - 0.5)*pool_diameter;
%                 end
%             end
%             allPC_x(rep,:) = PC_x;
%             allPC_y(rep,:) = PC_y;
%         else
%             weights= allweight(:,:,rep);
%             PC_x = allPC_x(rep,:);
%             PC_y = allPC_y(rep,:);
            
        for trial = 1:Ntrials
           trial
           
            [alpha, beta, gamma, etdecay, chi2, output]= optimize_chi_trial (...
            weights1,th,starting_ys,starting_xs,strad, platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,...
            pool_diameter,platform_radius,sigma_pc,sigma_ac,...
            Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PC_x, PC_y,...
            Strain, Temperature, day,trial, expmean);
        % may consider weight
%%
            allpars = [allpars;output];
            csvwrite('allpars.csv',allpars);
            
            PMs0=zeros(5,Nruns);
            allweights= zeros(N_pc,N_ac,Nruns);
         for rep=1:Nruns
            idx = randi(4); %randomly choose one of 4 starting locations
            starting_x = starting_xs(idx);
            starting_y = starting_ys(idx);

             [wres, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants, ...
            latency, speed_std,speed_ps, mean_angle,time_step] = ...
            run_trial (weights, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, ...
            Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, ...
            platform_x, platform_y, starting_x, starting_y, speed, hitwall);
            %run trial
            allweights(:,:,rep) = wres;

            PMs0(1,rep) = latency;
            PMs0(2,rep) = dist*100;
            PMs0(3,rep) = quadrants(4);
            PMs0(4,rep) = quadrants(2);
            PMs0(5,rep) = wall_zone;
            % new PMs
            PMs0(6,day,trial) = speed_std*100;
            PMs0(7,day,trial) = mean_angle;
            PMs0(8,day,trial) = time_step;
            
            %record performance measures
            
         end
            
         PMs(:,day,trial,mice) = mean(PMs0,2);
         weights1 = mean(allweights,3);

%             PMs(1,day,trial,mice) = latency;
%             PMs(2,day,trial,mice) = dist*100;
%             PMs(3,day,trial,mice) = quadrants(4);
%             PMs(4,day,trial,mice) = quadrants(2);
%             PMs(5,day,trial,mice) = wall_zone;
            %record performance measures
            
            dims = size(PMs);
            outf = [];
           if (length(dims) == 3)
            for i = 1:dims(1)
                for j = 1:dims(2)
                    for k = 1:dims(3)
                        outf = [outf; i j k PMs(i,j,k)];
                    end
                end
            end
        elseif (length(dims) == 4)     
            for i = 1:dims(1)
                for j = 1:dims(2)
                    for k = 1:dims(3)
                        for l = 1:dims(4)
                            outf = [outf; i j k l PMs(i,j,k,l)];
                        end
                    end
                end
            end
           end
        csvwrite('output.csv',outf);
            
        end
end
    
end