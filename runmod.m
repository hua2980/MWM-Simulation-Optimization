N_pc = 211; %Population of place cells [100..300]
N_ac = 36; %Population of action cells [25..50]

plot_trajectories = 0; %yes - 1, no - 0
plot_cognitive_maps = 0; %yes - 1, no - 0
pln = plot_trajectories + plot_cognitive_maps;
Nruns = 10; %how many runs to run if not plotting anything

pool_diameter = 1.4; %Maze diameter (\phi) in metres (m)
platform_radius = 0.06; %Platform radius (m)
sigma_pc = 0.1;  %place cell sigma (standard deviation), in meters [0.05..0.2]
sigma_ac = 2;  %action cell sigma (standard deviation), in action cells [1..3]

Vdecay = 0.82;  %velocity decay [0.75..0.95]
ac_const = 0.02;  %acceleration const [0.01..0.03]
Wnoise = 0.0004;  %Weight noise [0.0001..0.0007]
Wmult = 0.1;   %Weight multiplier [0.05..0.15]
hitwall = 0.5;  %punishment for hitting the wall [0..1]
speed = 0.175;  %mouse speed (m/s) [0.1..0.25]

Ntrials = 4;  %number of trials per day
Ndays = 25;  %number of days


for alpha = [0.005,0.01,0.015,0.02]
    etdecay = 0.75; %Eligibility trace decay (lambda) [0.75..0.95] LESS THAN GAMMA!
    beta = 6;  %Exploration-exploitation factor (\beta) [0.5..12]
    %alpha = 0.01;  %Learning rate (\alpha) [0.005..0.02]
    gamma = 0.95;  %Discount factor (\gamma) [0.75..0.95]

    alpha
    
if (pln > 0.5) %if any plots
    clf
    PMs = zeros(5,Ndays,Ntrials);  %performance measures to compute: latency, distance
    %time in target quadrant, opposite quadrant, and wall zone
else
    PMs = zeros(5,Ndays,Ntrials,Nruns);  %multiple runs
end


%Platform coordinates:
platform_x = cos(-pi/4)*pool_diameter/4; %x coordinate
platform_y = sin(-pi/4)*pool_diameter/4; %y coordinate

%Starting locations of the modeled animal (4 different ones):
strad = pool_diameter/2*0.85; %15% of maze radius to the wall
starting_xs = strad * [cos(pi/6) cos(pi/3) cos(7*pi/6) cos(4*pi/3)]; %x coordinates
starting_ys = strad * [sin(pi/6) sin(pi/3) sin(7*pi/6) sin(4*pi/3)]; %y coordinates

th = 0:pi/50:2*pi; %for plotting circles :)

if (pln > 0.5)

%Generate initial weights
weights = rand(N_pc,N_ac)*Wmult;

%Generate place cells
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
    for trial = 1:Ntrials
        
        idx = randi(4); %randomly choose one of 4 starting locations
        starting_x = starting_xs(idx);
        starting_y = starting_ys(idx);

        [wres, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants, latency] = ...
        run_trial (weights, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, ...
            Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, ...
            platform_x, platform_y, starting_x, starting_y, speed, hitwall);
        %run trial
        weights = wres;
               
        PMs(1,day,trial) = latency;
        PMs(2,day,trial) = dist;
        PMs(3,day,trial) = quadrants(4)*100;
        PMs(4,day,trial) = quadrants(2)*100;
        PMs(5,day,trial) = wall_zone*100;
        %record performance measures
        
        if (plot_trajectories)
        subplot(Ndays, Ntrials*pln,(day-1)*Ntrials*pln+(trial-1)*pln+1);
        hold on
        %plot the trajectory
        for i = 1:(length(track_x))-1
            line(track_x(i:i+1),track_y(i:i+1),'Color',[i/length(track_x),0,1-i/length(track_x)]);
        end
        %plot the maze and platform
        plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),'k');
        plot(platform_x+platform_radius*cos(th),platform_y+platform_radius*sin(th),'k');
        end
        
        if (plot_cognitive_maps)
        subplot(Ndays, Ntrials*pln,(day-1)*Ntrials*pln+trial*pln);
        hold on
        %plot the cognitive map
        for x = -(pool_diameter/2):(pool_diameter/6):(pool_diameter/2)
            for y = -(pool_diameter/2):(pool_diameter/6):(pool_diameter/2)
                if (x^2 + y^2 <= (pool_diameter/2)^2)
                    for k = 1:N_ac
                        PC_activation = zeros(1,N_pc);
                        for i = 1:N_pc
                            PC_activation(i) = exp(-((x - PC_x(i))^2 + (y...
                            - PC_y(i))^2)/(2*sigma_pc^2));
                        end
                %Calculate AC activation (i.e. value of the action)
                        AC_activation = zeros(1,N_ac);
                        for i = 1:N_ac
                            for j = 1:N_pc
                                AC_activation(i) = AC_activation(i) + ...
                                    PC_activation(j)*weights(j,i);
                            end
                        end
                        x2 = x + (AC_activation(k)/10)*cos(k/N_ac*2*pi);
                        y2 = y + (AC_activation(k)/10)*sin(k/N_ac*2*pi);
                        hold on;
                        line([x x2],[y y2],'Color',[k/N_ac 0 1-k/N_ac]);
                    end
                end
            end
        end
        %plot the maze and platform
        plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),'k');
        plot(platform_x+platform_radius*cos(th),platform_y+platform_radius*sin(th),'k'); 
        end

    end
end
else
% run multiple times without plotting!

for rep = 1:Nruns
%Generate initial weights for each run
weights = rand(N_pc,N_ac)*Wmult;

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
    
    
        for trial = 1:Ntrials
        
        idx = randi(4); %randomly choose one of 4 starting locations
        starting_x = starting_xs(idx);
        starting_y = starting_ys(idx);

        [wres, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants, latency] = ...
        run_trial (weights, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, ...
            Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, ...
            platform_x, platform_y, starting_x, starting_y, speed, hitwall);
        %run trial
        weights = wres;
               
        PMs(1,day,trial,rep) = latency;
        PMs(2,day,trial,rep) = dist*100;
        PMs(3,day,trial,rep) = quadrants(4);
        PMs(4,day,trial,rep) = quadrants(2);
        PMs(5,day,trial,rep) = wall_zone;
        %record performance measures
    end
end
end
end



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
                    outf = [outf; i j k l PMs(i,j,k,l)];% Pars(j,1:4)];
                end
            end
        end
    end
end

name= ['aplha',mat2str(alpha),'.csv'];
csvwrite(name,outf);

end
