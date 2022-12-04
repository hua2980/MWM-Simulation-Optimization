function[alpha, beta, gamma, etdecay, chi2, output]= optimize_chi_trial (...
    weights1, th,starting_ys,starting_xs,strad, platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,...
    pool_diameter,platform_radius,sigma_pc,sigma_ac,...
    Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays, PC_x, PC_y,...
    Strain, Temperature, day,trial, expmean)
% function to return the best metaparameter every day

% expmean = readmatrix('expmean.csv');
% expmean = csvread('expmean2.csv');
% Strain=1;
% Temperature=18;
% day=1;
PMexp = expmean(find((expmean(:,2)==Strain) & (expmean(:,3)== Temperature) & (expmean(:,4)==day)& (expmean(:,5)==trial)),:);
idx=[1 4 2 5 3];
PMexp = PMexp(idx,:);

PMs = zeros(5,Nsets,Nruns);  %multiple runs
Pars = zeros(Nsets,5);

% get 50 random metaparameters
for set = 1:Nsets
%% 
% PMs = zeros(5,Nsets,Nruns);  %multiple runs
% Pars = zeros(Nsets,5);

    alpha= 0.005+(0.02-0.005).*rand;   %Learning rate (\alpha) [0.005..0.02]
    beta = 0.5+ (12-0.5).*rand;  %Exploration-exploitation factor (\beta) [0.5..12]
    gamma = 0.75+(0.95-0.75).*rand;  %Discount factor (\gamma) [0.75..0.95]
    etdecay = 0.75+(0.95-0.75).*rand; %Eligibility trace decay (lambda) [0.75..0.95] LESS THAN GAMMA!
    
    while etdecay >= gamma
        etdecay = 0.75+(0.95-0.75).*rand;
        gamma = 0.75+(0.95-0.75).*rand;
    end
    
    Pars(set,1) = alpha;
    Pars(set,2) = beta;
    Pars(set,3) = gamma;
    Pars(set,4) = etdecay;
%%
      % test=[];
    
    for rep = 1:Nruns
    %Generate initial weights for each run
%     if day==1
%         weights = rand(N_pc,N_ac)*Wmult;
%     else
%         weights = weights1;
%     end
    
         weights = weights1;
    

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

            PMs(1,set,rep) = latency;
            PMs(2,set,rep) = dist*100;
            PMs(3,set,rep) = quadrants(4);
            PMs(4,set,rep) = quadrants(2);
            PMs(5,set,rep) = wall_zone;
            
            PMs(6,day,trial) = speed_std*100;
            PMs(7,day,trial) = mean_angle;
            PMs(8,day,trial) = time_step;
            
            %record performance measures
          
            
%         PMmod= sum(PMs,3)./rep;
%         chi2=0;
%         for i= 1:5
%             a = (PMexp(i,6)-PMmod(i,set)).^2 / PMexp(i,7).^2;
%             chi2= chi2+a;
%         end
%         test = [test;chi2];
    
    end
    
    PMmod= mean(PMs,3);
    % a= zeros(5,1);
    chi2=0;
    for i= 1:5
        a = (PMexp(i,6)-PMmod(i,set)).^2 / PMexp(i,7).^2;
        chi2= chi2+a;
    end
    Pars(set,5) = chi2;
    
    % plot(test)
end




%%

% best estimation
finalpar= zeros(5,15);
chi33= [];

chi2min = mink(Pars(:,5),15);
chi2find= zeros(15,5);
for i= 1:15
    idx= find(Pars(:,5)== chi2min(i));
    chi2find(i,:)= Pars(idx,:);
end

for set=1:15
    Pars2= [];
    
    alpha= chi2find(set,1);
    beta= chi2find(set,2);
    gamma= chi2find(set,3);
    etdecay= chi2find(set,4);
    chi2 = chi2find(set,5);
    
    %%
    % calculate better alpha
    step1 = 0.04.*(0.02-0.005);
    step2 = 0.2.*(0.02-0.005);
    
    while 1
        alpha0=alpha;
        alpha= alpha+step1;
        if alpha > 0.02
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
        
        chi33= [chi33;chi2t];
    
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha0 beta gamma etdecay chi2];
    Pars2 = [Pars2;aa];
    
    %alpha-step1
    alpha= chi2find(set,1);
    chi2= chi2find(set,5);
    while 1
        alpha0=alpha;
        alpha= alpha-step1;
        if alpha < 0.005
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
        
        chi33= [chi33;chi2t];
        
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha0 beta gamma etdecay chi2];
    Pars2 = [Pars2;aa];
    
    % alpha+step2
    alpha= chi2find(set,1);
    chi2= chi2find(set,5);
    while 1
        alpha0=alpha;
        alpha= alpha+step2;
        if alpha > 0.02
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
        
        chi33= [chi33;chi2t];
        
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha0 beta gamma etdecay chi2];
    Pars2 = [Pars2;aa];
    
    %alpha-step2
    alpha= chi2find(set,1);
    chi2= chi2find(set,5);
    while 1
        alpha0=alpha;
        alpha= alpha-step2;
        if alpha < 0.005
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
        
        chi33= [chi33;chi2t];
        
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha0 beta gamma etdecay chi2];
    Pars2 = [Pars2;aa];
    
    %get best alpha
    idx= find(Pars2(:,5)==min(Pars2(:,5)));
    alpha= unique(Pars2(idx,1));
    beta= unique(Pars2(idx,2));
    gamma= unique(Pars2(idx,3));
    etdecay= unique(Pars2(idx,4));
    chi2 = unique(Pars2(idx,5));
    
    finalpar(1,set) = alpha;
    %%
    
    
    
    idx= idx(1);
    
    % calculate the best beta
    step1 = 0.04.*(12-0.5);
    step2 = 0.2.*(12-0.5);
    
    while 1
        beta0=beta;
        beta= beta+step1;
        if beta > 12
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
        
        chi33 = [chi33;chi2t];
    
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha beta0 gamma etdecay chi2];
    Pars2 = [Pars2;aa];
    
    
    %beta-step1
    beta= Pars2(idx,2);
    chi2 = Pars2(idx,5);
    while 1
        beta0=beta;
        beta= beta-step1;
        if beta < 0.5
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
        
        chi33 = [chi33;chi2t];
    
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha beta0 gamma etdecay chi2];
    Pars2 = [Pars2;aa];
    
    % beta-step2
    beta= Pars2(idx,2);
    chi2 = Pars2(idx,5);
    while 1
        beta0=beta;
        beta= beta+step2;
        if beta > 12
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
          
        chi33= [chi33;chi2t];
    
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha beta0 gamma etdecay chi2];
    Pars2 = [Pars2;aa];
    
    %beta-step2
    beta= Pars2(idx,2);
    chi2 = Pars2(idx,5);
    while 1
        beta0=beta;
        beta= beta-step2;
        if beta < 0.5
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
          
        chi33= [chi33;chi2t];
    
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha beta0 gamma etdecay chi2];
    Pars2 = [Pars2;aa];
    
    %get best beta
    idx= find(Pars2(:,5)==min(Pars2(:,5)));
    alpha= unique(Pars2(idx,1));
    beta= unique(Pars2(idx,2));
    gamma= unique(Pars2(idx,3));
    etdecay= unique(Pars2(idx,4));
    chi2 = unique(Pars2(idx,5));
    
    finalpar(2,set) = beta;
    
    
    
    %%
    
    
    
    idx= idx(1);
    
    % calculate the best gamma 0.75~0.95
    step1 = 0.04.*(0.95-0.75);
    step2 = 0.2.*(0.95-0.75);
    
    while 1
        gamma0=gamma;
        gamma= gamma+step1;
        if gamma > 0.95
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
        
        chi33= [chi33;chi2t];
    
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha beta gamma0 etdecay chi2];
    Pars2 = [Pars2;aa];
    
    
    %gamma-step1
    gamma= Pars2(idx,3);
    chi2 = Pars2(idx,5);
    while 1
        gamma0=gamma;
        gamma= gamma-step1;
        if gamma < 0.75
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
        
        chi33= [chi33;chi2t];
    
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha beta gamma0 etdecay chi2];
    Pars2 = [Pars2;aa];
    
    % gamma+step2
    gamma= Pars2(idx,3);
    chi2 = Pars2(idx,5);
    while 1
        gamma0=gamma;
        gamma= gamma+step2;
        if gamma > 0.95
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
          
        chi33= [chi33;chi2t];
    
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha beta gamma0 etdecay chi2];
    Pars2 = [Pars2;aa];
    
    %gamma-step2
    gamma= Pars2(idx,3);
    chi2 = Pars2(idx,5);
    while 1
        gamma0=gamma;
        gamma= gamma-step2;
        if gamma < 0.75
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
          
        chi33= [chi33;chi2t];
    
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha beta gamma0 etdecay chi2];
    Pars2 = [Pars2;aa];
    
    %get best gamma
    idx= find(Pars2(:,5)==min(Pars2(:,5)));
    alpha= unique(Pars2(idx,1));
    beta= unique(Pars2(idx,2));
    gamma= unique(Pars2(idx,3));
    etdecay= unique(Pars2(idx,4));
    chi2 = unique(Pars2(idx,5));
    
    finalpar(3,set) = gamma;
    
    %%
    
    
    
    idx= idx(1);
    
    % calculate the best etdecay [0.75..0.95]
    step1 = 0.04.*(0.95-0.75);
    step2 = 0.2.*(0.95-0.75);
    
    while 1
        etdecay0=etdecay;
        etdecay= etdecay+step1;
        if etdecay > 0.95
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
        
        chi33= [chi33;chi2t];
    
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha beta gamma0 etdecay0 chi2];
    Pars2 = [Pars2;aa];
    
    
    %etdecay-step1
    etdecay= Pars2(idx,4);
    chi2 = Pars2(idx,5);
    while 1
        etdecay0=etdecay;
        etdecay= etdecay-step1;
        if etdecay < 0.75
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
        
        chi33= [chi33;chi2t];
    
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha beta gamma etdecay0 chi2];
    Pars2 = [Pars2;aa];
    
    % etdecay-step2
    etdecay= Pars2(idx,4);
    chi2 = Pars2(idx,5);
    while 1
        etdecay0=etdecay;
        etdecay= etdecay+step2;
        if etdecay > 0.95
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
          
        chi33= [chi33;chi2t];
    
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha beta gamma etdecay0 chi2];
    Pars2 = [Pars2;aa];
    
    %etdecay-step2
    etdecay= Pars2(idx,4);
    chi2 = Pars2(idx,5);
    while 1
        etdecay0=etdecay;
        etdecay= etdecay-step2;
        if etdecay < 0.75
            break
        end
        
        [chi2t] = calculate_chi2_trial (alpha, beta, gamma, etdecay,th,starting_ys,starting_xs,strad,...
        platform_x,platform_y,N_pc,N_ac,Nruns,Nsets,pool_diameter,platform_radius,sigma_pc,sigma_ac,...
        Vdecay,ac_const,Wnoise,Wmult,hitwall,speed,Ntrials,Ndays,PMexp,weights1,PC_x, PC_y);
          
        chi33= [chi33;chi2t];
    
        if chi2t <= chi2
            chi2 = chi2t;
        else
            break
        end  
    end
    
    aa = [alpha beta gamma etdecay0 chi2];
    Pars2 = [Pars2;aa];
    
    %get best etdecay
    idx= find(Pars2(:,5)==min(Pars2(:,5)));
    alpha= unique(Pars2(idx,1));
    beta= unique(Pars2(idx,2));
    gamma= unique(Pars2(idx,3));
    etdecay= unique(Pars2(idx,4));
    chi2 = unique(Pars2(idx,5));
    
    finalpar(4,set) = etdecay;
    finalpar(5,set) = chi2;
end

%%
% the best estimation
set= find(finalpar(5,:)==min(finalpar(5,:)));
alpha= finalpar(1,set);
beta= finalpar(2,set);
gamma= finalpar(3,set);
etdecay= finalpar(4,set);
chi2 = finalpar(5,set);
output= [alpha beta gamma etdecay chi2]
