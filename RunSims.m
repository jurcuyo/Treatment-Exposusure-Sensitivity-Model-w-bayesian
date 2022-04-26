%% RunsIndivFits.m

%  This code determines parameter values that fit our BLI model to 
%  experimental data, and is based on Working_Mar2019.m from the old 
%  directory prior to codebase initial commit. 
clear variables

%% User Input to Determine Which Data to Analyze

PDXline=input('specify PDX line (in format ''GBM39'' with no spaces): ');
Site=input('specify ''IC'' for intracranial tumor or ''Flank'': ');

CodeFolder   = '~/ProjectFolder/Codes';
InputFolder  = '~/ProjectFolder/Data';
ResultsFolder= '~/ProjectFolder/Results';
PlotsFolder  = '~/ProjectFolder/Plots';

%% Load Data (variables listed below)
% tpoints: all time points in total
% group1: sham control data
% group2: treated data
% DrugStart: when does drug pulsing initate
% labelstr: for plot titles and such

cd(InputFolder)
data_untx = xlsread('data.xlsx',[PDXline,' ',Site],'C2:Z7');
data_tx = xlsread('data.xlsx',[PDXline,' ',Site],'C8:Z14');
if data_untx(1,:) == data_tx(1,:)
    tpoints = data_untx(1,:);
else
    error("Check excel- tpoints do not match.")
end
group1 = data_untx(2:end,:);
group2 = data_tx(2:end,:);

if strcmp(PDXline,'GBM6')==1 
    if strcmp(Site,'Flank')==1
        labelstr='G6 Flank 11.03.17';
        DrugStart=21;
    elseif strcmp(Site,'IC')==1
        labelstr='G6 IC 01.23.18';
        DrugStart=7;
    end
elseif strcmp(PDXline,'GBM12')==1
    if strcmp(Site,'Flank')==1
        labelstr='G12 Flank 03.13.18';
        DrugStart=14;
    elseif strcmp(Site,'IC')==1
        labelstr='G12 IC 03.13.18';
        DrugStart=7;
    end
elseif strcmp(PDXline,'GBM39')==1 
    alteredDose = 0.01; %Drug stability
    if strcmp(Site,'Flank')==1
        labelstr='G39 Flank 09.07.17';
        DrugStart=14;
    elseif strcmp(Site,'IC')==1
        labelstr='G39 IC 11.21.17';
        DrugStart=7;
    end
end
cd(CodeFolder)

%% Data & Analysis
errorBars = input('Enter ''Y'' for error bars: ');
data_vis = input('Enter ''Y'' for data visualization: ');
for treatType=[1 2] %runs group1 for untreated, then group2 for treated
cd(CodeFolder)

%% Defining Data for Groups
FitResults=struct; %Create Structure to Hold fit results
group(group<=0) = NaN; %Replaces all 0's & negatives with NaN's

% Check for groups that terminated before the end of the experiment; 
% adjust time to match the existence of the group:
if size(tpoints,2)>size(group,2)
    t_thisExpt = tpoints(:,1:size(group,2));
else
    t_thisExpt = tpoints;
end

% Removes unnecessary NaN's for the time points and flux measurements
t_thisExpt(:,find(all(isnan(group),1)))=[];
group(:,all(isnan(group),1))=[];

% For running model simulations
stepsize=.1;
simT=0:stepsize:max(t_thisExpt);
simT_short = 0:ceil(max(t_thisExpt)/300/0.25)*0.25:max(t_thisExpt);

%% Plot BLI flux for Experiment
figure(1);clf
FigureSize(25,20,'centimeters')
for mouseNum = 1:size(group,1)
    h=semilogy(t_thisExpt(~isnan(group(mouseNum,:))),group(mouseNum,~isnan(group(mouseNum,:))),'x--','LineWidth',2, 'MarkerSize',8);
    hold on
end
xlabel('Time (days)','FontSize',20)
ylabel('Log of Total Flux (p/s)','FontSize',20)
if strcmp(data_vis,'Y')==0
    title([labelstr ': Group ' num2str(treatType)])
end
ax=gca;
ax.FontSize = 20;
grid on

%% Run Fits and Simulations
if treatType == 1 && strcmp(data_vis,'Y')==0 % Sham-Control // ADC
    %% Fit Rho (and C0) to Group 1 data
    
    % Initial guesses for parameters C0 and rho:
    C0_preFit = 1*10^5;
    rho_preFit = log(2)/(60/24); % ~60 divisions/day

    shamControl_fits = NaN(size(group,1),2);
    unTx_total = NaN(size(group,1),length(simT));
    
    % Function call by mouse - removes NaNs, runs linear regression, & obtains parameter estimates
    % Estimates used for model simulations, saved for plots later
    for mouseNum=1:size(group,1)
        shamControl_fits(mouseNum,:) = runfit_noTx([C0_preFit rho_preFit],t_thisExpt,group(mouseNum,:));
        t_indvMouse=t_thisExpt(~isnan(group(mouseNum,:)));
        simT_indvMouse=0:stepsize:(max(t_indvMouse));
        unTx_total(mouseNum,1:length(simT_indvMouse)) = runModelSims(simT_indvMouse,[shamControl_fits(mouseNum,2) 10^0 0 0 0 shamControl_fits(mouseNum,1)],zeros());
    end

    rho_dist_noTrunc = fitdist(shamControl_fits(:,2),'Normal');
    rho_distribution = truncate(rho_dist_noTrunc,0,inf); % Truncate since biologically must be >=0
    
elseif treatType==grp414 && strcmp(data_vis,'Y')==0
    %% Load Flank fits
    cd(ResultsFolder)
    if strcmp(Site,'IC')
        if exist('alteredDose','var')
            flankResults = load([PDXline,'_Dose',num2str(alteredDose),'_Flank_BayesianFit_Output_',char(month(datetime('today'),'shortname')),num2str(year(datetime('today'))),'.mat'],'FitResults');
        else
            flankResults = load([PDXline,'_Flank_BayesianFit_Output_',char(month(datetime('today'),'shortname')),num2str(year(datetime('today'))),'.mat'],'FitResults');
        end
        
        muS_fromFlank = NaN(length(flankResults.FitResults),1);
        for mouseNum=1:length(flankResults.FitResults)
            muS_fromFlank(mouseNum) = flankResults.FitResults(mouseNum).Parameters.muS;
        end
        
        flankPost = zeros(size(flankResults.FitResults(mouseNum).Process.muS_grid));
        for m=1:5
            marg = squeeze(sum(sum(sum(flankResults.FitResults(m).Process.posterior_grid,3),2),1));
            flankPost = flankPost + marg/sum(marg(:));
        end
        
        clear flankResults
    end

    %% Drug Dynamics
    doseInterval=7; % Time between pulses in days
    totalpulses=floor((t_thisExpt(end))/doseInterval); % Total number of pulses
    decayRate = log(2)/7;
    dose=0.1; %mg - this should be a concentration but not using actual doses
    doses=ones(1,totalpulses)*dose; % vector of pulses 
    
    % If dosing doesn't start at day 7, add zeros in the doses vector:
    if DrugStart==14
        doses(1)=0;
    elseif DrugStart==21
        doses(1:2)=[0 0];
    elseif DrugStart==28
        doses(1:3)=[0 0 0];
    end

    % In the case of drug stability issue
    if exist('alteredDose','var') && strcmp(PDXline,'GBM39')==1 && strcmp(Site,'Flank')
        doses(2:3)=alteredDose; % First 2 drug doses
    end

    doses = repmat(doses,[size(group,1) 1]);
    
    % Initialize solution vectors of drug level and integral of drug at each time point:
    A=zeros(size(simT)); % initializing A
    int_A=zeros(size(t_thisExpt)); % initializing the integral of D(t) values vector
    int_Asim=zeros(size(simT));
    int_Asim_short = zeros(size(simT_short));

    for p=1:length(doses(1,:)) % loop through pulses
        % Compute drug vs time due to that pulse:
        simTHeaviside = zeros(size(simT));
        simTshortHeaviside = zeros(size(simT_short));
        thisExpHeaviside = zeros(size(t_thisExpt));
        
        simTHeaviside((simT-7*p>0)) = 1;
        simTshortHeaviside((simT_short-7*p>0)) = 1;
        thisExpHeaviside((t_thisExpt-7*p>0))=1;
        
        A=A+doses(:,p).*2^p*exp(-decayRate.*simT).*simTHeaviside;%heaviside(simT-7*p); % adding the effects of sequential pulses on the later time points 

        % Compute integral of drug vs time:
        int_A=int_A+ doses(:,p).*(2^p).*((exp(-decayRate*7*p)-exp(-decayRate*t_thisExpt))/decayRate).*thisExpHeaviside;%heaviside(t_thisExpt-7*p);
        int_Asim=int_Asim+ doses(:,p).*(2^p).*((exp(-decayRate*7*p)-exp(-decayRate*simT))/decayRate).*simTHeaviside;%heaviside(simT-7*p);
        int_Asim_short = int_Asim_short + doses(:,p).*(2^p).*((exp(-decayRate*7*p)-exp(-decayRate*simT_short))/decayRate).*simTshortHeaviside;%heaviside(simT_short-7*p);
    end
    
    %% Fit muS/gamma (and remaining parameters) to Group 3/4 (treated) data
    % Bayesian approach
    for mouseNum = 1:length(group(:,1))
        t_indvMouse=t_thisExpt(~isnan(group(mouseNum,:)));
        simT_indvMouse=0:stepsize:(max(t_indvMouse));
        simTshort_indvMouse=0:ceil(max(t_thisExpt)/300/0.25)*0.25:(max(t_indvMouse));
        
        if length(simT_indvMouse) < length(simTshort_indvMouse)
            simT_useBayesian = simT_indvMouse;
            intA_useBayesian = int_Asim;
        else
            simT_useBayesian = simTshort_indvMouse;
            intA_useBayesian = int_Asim_short;
        end
        q_grid = linspace(-10,-2,20);
        z_grid = linspace(0,1,20);
        
        % Use distribution for rho obtained from sham control (group 1)
        rho_grid = linspace(max(rho_distribution.mean - 4*rho_distribution.std,0),rho_distribution.mean + 4*rho_distribution.std,20);
        p_rho = pdf(rho_distribution,rho_grid) * 1/(1-cdf(rho_dist_noTrunc,0));
        
        if strcmp(Site,'Flank')
            muS_grid = linspace(0,10,40);
            p_muS = ones(length(muS_grid),1); %Since the prior for muS is a uniform in Flank
            gamma_sim = 1;
            M = single(NaN(length(t_indvMouse),length(q_grid),length(z_grid),length(rho_grid),length(muS_grid))); % to hold all simulations at data timepoints
            N = single(NaN(length(simT_useBayesian),length(q_grid),length(z_grid),length(rho_grid),length(muS_grid))); % to hold all simulations at all timepoints
            sse = NaN(length(q_grid),length(z_grid),length(rho_grid),length(muS_grid)); % to hold SSE
        elseif strcmp(Site,'IC')
            muS_grid = linspace(0,10,40);
            gamma_grid = linspace(0,1,40); 

            M = single(NaN(length(t_indvMouse),length(q_grid),length(z_grid),length(rho_grid),length(muS_grid),length(gamma_grid))); % to hold all simulations at data timepoints
            N = single(NaN(length(simT_useBayesian),length(q_grid),length(z_grid),length(rho_grid),length(muS_grid),length(gamma_grid))); % to hold all simulations at all timepoints
            sse = single(NaN(length(q_grid),length(z_grid),length(rho_grid),length(muS_grid),length(gamma_grid))); % to hold SSE
        end
        
        if strcmp(Site,'Flank')
            squaredError_eachT = NaN(size(M));
            pLikelihood_eachT = NaN(size(M));
            posterior_eachT = NaN(size(M));
            p_likelihood = NaN(size(sse)); % to hold probability likelihood
            posterior = NaN(size(sse)); % to hold posterior
        elseif strcmp(Site,'IC')
            squaredError_eachT = single(NaN(size(M)));
            pLikelihood_eachT = single(NaN(size(M)));
            posterior_eachT = single(NaN(size(M)));
            p_likelihood = single(NaN(size(sse))); % to hold probability likelihood
            posterior = single(NaN(size(sse))); % to hold posterior
        end
        
        var_data = 1;  % How much variance do we expect to see in the BL imaging itself??
        alpha = 1;      % Some normalizing constant, currently set to 1
        
        tic
        %posterior(q,z,rho,muS)
        %posterior(q,z,rho,muS,gamma)
        for a = 1:length(q_grid)
            q_sim = 10^(q_grid(a));
            for b = 1:length(z_grid)
                z_sim = z_grid(b);
                for c = 1:length(rho_grid)
                    rho_sim = rho_grid(c);
                    % Calculate C_0 based on the value for rho & first data point
                    C0_sim = group(mouseNum,1)*exp(-rho_sim*t_indvMouse(1));
                    
                    for d = 1:length(muS_grid)
                        muS_sim = muS_grid(d);
                        if strcmp(Site,'Flank')
                            
                            M(:,a,b,c,d) = runModelSims(t_indvMouse,[rho_sim,q_sim,muS_sim,z_sim,gamma_sim,C0_sim],int_A(mouseNum,1:length(t_indvMouse)));
                            N(:,a,b,c,d) = runModelSims(simT_useBayesian,[rho_sim,q_sim,muS_sim,z_sim,gamma_sim,C0_sim],intA_useBayesian(mouseNum,1:length(simT_useBayesian)));

                            sse(a,b,c,d) = sum((log10(group(mouseNum,1:length(t_indvMouse))) - log10(M(:,a,b,c,d)')).^2);
                            p_likelihood(a,b,c,d) = alpha*exp(-1/(2*var_data) * sse(a,b,c,d));
                            posterior(a,b,c,d) = p_rho(c) * p_muS(d) * p_likelihood(a,b,c,d);
                            
                            % Note: p_z=1, p_q=1 b/c all uniform
                        elseif strcmp(Site,'IC')
                            for e = 1:length(gamma_grid)
                                
                                gamma_sim = gamma_grid(e);
                                
                                M(:,a,b,c,d,e) = runModelSims(t_indvMouse,[rho_sim,q_sim,muS_sim,z_sim,gamma_sim,C0_sim],int_A(mouseNum,1:length(t_indvMouse)));
                                N(:,a,b,c,d,e) = runModelSims(simT_useBayesian,[rho_sim,q_sim,muS_sim,z_sim,gamma_sim,C0_sim],intA_useBayesian(mouseNum,1:length(simT_useBayesian)));
                                
                                sse(a,b,c,d,e) = sum((log10(group(mouseNum,1:length(t_indvMouse))') - log10(M(:,a,b,c,d,e))).^2,'omitnan');
                                p_likelihood(a,b,c,d,e) = alpha*exp(-1/(2*var_data) * sse(a,b,c,d,e));
                                posterior(a,b,c,d,e) = p_rho(c) * flankPost(d) * p_likelihood(a,b,c,d,e);                                
                                % Note: p_gamma=1, p_z=1, p_q=1 b/c all uniform
                            end
                        end
                    end
                end
            end
        end
        toc
        
        [maxVal,maxIndx] = max(posterior(:));
        
        if strcmp(Site,'Flank')
            FitResults(mouseNum).Process = cell2struct({q_grid,z_grid,rho_grid,muS_grid,sse,p_likelihood,posterior}',{'q_grid','z_grid','rho_grid','muS_grid','sse_grid','likelihood_grid','posterior_grid'},1);
            [fit_a,fit_b,fit_c,fit_d] = ind2sub(size(posterior),maxIndx);
            best_sse = sse(fit_a,fit_b,fit_c,fit_d);
            best_likelihood = p_likelihood(fit_a,fit_b,fit_c,fit_d);
            best_posterior = posterior(fit_a,fit_b,fit_c,fit_d);
        elseif strcmp(Site,'IC')
            FitResults(mouseNum).Process = cell2struct({q_grid,z_grid,rho_grid,muS_grid,gamma_grid,sse,p_likelihood,posterior}',{'q_grid','z_grid','rho_grid','muS_grid','gamma_grid','sse_grid','likelihood_grid','posterior_grid'},1);
            [fit_a,fit_b,fit_c,fit_d,fit_e] = ind2sub(size(posterior),maxIndx);
            best_sse = sse(fit_a,fit_b,fit_c,fit_d,fit_e);
            best_likelihood = p_likelihood(fit_a,fit_b,fit_c,fit_d,fit_e);
            best_posterior = posterior(fit_a,fit_b,fit_c,fit_d,fit_e);
        end
        FitResults(mouseNum).Process.All = M;
        
        q_fit = q_grid(fit_a);
        z_fit = z_grid(fit_b);
        rho_fit = rho_grid(fit_c);
        C0_fit = group(mouseNum,1)*exp(-rho_fit*t_indvMouse(1));
        muS_fit = muS_grid(fit_d);
        if strcmp(Site,'IC')
            gamma_fit = gamma_grid(fit_e);
        else
            gamma_fit = 1;
        end
        
        if strcmp(errorBars,'Y')
            [dN1,dN2]=size(N);
            longN = NaN(dN2,dN1);
            boundsN = NaN(length(simT_useBayesian),2);
            
            hist_minExp = floor(log10(min(min(group))))-1;
            hist_maxExp = ceil(log10(max(max(group))))+1;
            
            histBinNum = (hist_maxExp - hist_minExp + 1) * 10;
            
            histPost = NaN(length(simT_useBayesian),histBinNum,4);
            
            histPost(1:length(simT_useBayesian),:,1) = repmat(logspace(hist_minExp,hist_maxExp,histBinNum),length(simT_useBayesian),1); %Bins
            
            tic
            for tt = 1:length(simT_useBayesian)
                longN(:,tt) = N(tt,:);
                
                for iii = 1:histBinNum
                    if iii == 1 %First Case
                        histPost(tt,iii,2) = sum(posterior(longN(:,tt) < histPost(tt,iii,1)));
                    else
                        histPost(tt,iii,2) = sum(posterior(longN(:,tt) >= histPost(tt,iii-1,1) & longN(:,tt) < histPost(tt,iii,1)));
                    end
                    histPost(tt,iii,3) = histPost(tt,iii,2) / sum(posterior(:)); %pdf
                    histPost(tt,iii,4) = sum(histPost(tt,1:iii,3)); %cdf
                end
                
                [~,closest_lb] = min(abs(histPost(tt,:,4)-.025));
                if histPost(tt,closest_lb,4) == .025 || closest_lb==1
                    boundsN(tt,1) = histPost(tt,closest_lb,1); %lower bounds
                else
                    if histPost(tt,closest_lb,4) > .025
                        leftpoint_lb_Post = squeeze(histPost(tt,closest_lb - 1,[1 4]));
                        rightpoint_lb_Post = squeeze(histPost(tt,closest_lb,[1 4]));
                    elseif histPost(tt,closest_lb,4) < .025
                        leftpoint_lb_Post = squeeze(histPost(tt,closest_lb,[1 4]));
                        rightpoint_lb_Post = squeeze(histPost(tt,closest_lb + 1,[1 4]));
                    end
                    lb_slope_Post = (rightpoint_lb_Post(2)-leftpoint_lb_Post(2))/(rightpoint_lb_Post(1)-leftpoint_lb_Post(1));
                    lb_est = (0.025 - histPost(tt,closest_lb,4))/lb_slope_Post + histPost(tt,closest_lb,1);
                    boundsN(tt,1) = lb_est; %lower bounds
                end
                
                [~,closest_ub] = min(abs(histPost(tt,:,4)-.975));
                if histPost(tt,closest_ub,4) == .975 || closest_ub==histBinNum
                    boundsN(tt,2) = histPost(tt,closest_ub,1); %uppere bounds
                else
                    if histPost(tt,closest_ub,4) > .975
                        leftpoint_ub_Post = squeeze(histPost(tt,closest_ub - 1,[1 4]));
                        rightpoint_ub_Post = squeeze(histPost(tt,closest_ub,[1 4]));
                    elseif histPost(tt,closest_ub,4) < .975
                        leftpoint_ub_Post = squeeze(histPost(tt,closest_ub,[1 4]));
                        rightpoint_ub_Post = squeeze(histPost(tt,closest_ub + 1,[1 4]));
                    end
                    ub_slope_Post = (rightpoint_ub_Post(2)-leftpoint_ub_Post(2))/(rightpoint_ub_Post(1)-leftpoint_ub_Post(1));
                    ub_est = (0.975 - histPost(tt,closest_ub,4))/ub_slope_Post + histPost(tt,closest_ub,1);
                    boundsN(tt,2) = ub_est; %upper bounds
                end
            end
            toc
            
            %Intervals on posterior -> parameter
            norm_post = posterior(:)/sum(posterior(:));
            [sort_npost, sort_nidx] = sort(norm_post,'descend');
            cum_post = zeros(size(sort_npost));
            cum_post(1) = sort_npost(1);
            for i=2:length(sort_npost)
                cum_post(i) = cum_post(i-1)+sort_npost(i);
            end
            
            hpd95_post = sort_npost(1:find(cum_post>0.95,1,'first'));
            hpd95_idx = sort_nidx(1:find(cum_post>0.95,1,'first'));
            if strcmp(Site,'Flank')
                [hpda, hpdb, hpdc, hpdd] = ind2sub(size(posterior),hpd95_idx);
            elseif strcmp(Site,'IC')
                [hpda, hpdb, hpdc, hpdd, hpde] = ind2sub(size(posterior),hpd95_idx);
            end %[q,z,rho,muS,(gamma)]
        end
        
        FitResults(mouseNum).Parameters = cell2struct({rho_fit,q_fit,muS_fit,z_fit,gamma_fit,C0_fit}',{'rho','q','muS','z','gamma','C0'},1);
        
        mse = best_sse ./ length(t_thisExpt);
        rmse = sqrt(mse);
        FitResults(mouseNum).Statistics = cell2struct({best_sse,mse,rmse, best_likelihood, best_posterior}',{'sse','mse','rmse','p_likelihood','posterior'},1);
        
        FitResults(mouseNum).Process.rhoDist = rho_distribution;
        params=[FitResults(mouseNum).Parameters.rho,10^FitResults(mouseNum).Parameters.q,FitResults(mouseNum).Parameters.muS,FitResults(mouseNum).Parameters.z,FitResults(mouseNum).Parameters.gamma,FitResults(mouseNum).Parameters.C0];
        
        FitResults(mouseNum).SimBLI = runModelSims(simT_indvMouse,params,int_Asim(mouseNum,1:length(simT_indvMouse)));
        
        if strcmp(errorBars,'Y')
            FitResults(mouseNum).lb_BLI = boundsN(:,1);
            FitResults(mouseNum).ub_BLI = boundsN(:,2);

            FitResults(mouseNum).Parameters.q_hist = hpda;
            FitResults(mouseNum).Parameters.z_hist = hpdb;
            FitResults(mouseNum).Parameters.rho_hist = hpdc;
            FitResults(mouseNum).Parameters.muS_hist = hpdd;
            if strcmp(Site,'IC')
                FitResults(mouseNum).Parameters.gamma_hist = hpde;
            end
        end
    end
end % end loop of fits for sham control and treated groups

%% Save results
if strcmp(data_vis,'N')
    if exist('alteredDose','var') && strcmp(PDXline,'GBM39')
        paramsavestr=[PDXline,'_Dose',num2str(alteredDose),'_',Site,'_','BayesianFit_Output_late',char(month(datetime('today'),'shortname')),num2str(year(datetime('today'))),'.mat'];
    else
        paramsavestr=[PDXline,'_',Site,'_','BayesianFit_Output_late',char(month(datetime('today'),'shortname')),num2str(year(datetime('today'))),'.mat'];
    end
    cd(ResultsFolder)
    save(paramsavestr,'FitResults');
    cd(CodeFolder) % return to parent directory containing this code file
end

%% Generate Plots
cd(CodeFolder)
axLimits=[0 7*(floor(max(t_thisExpt/7)+1)) ax.YLim(1) ax.YLim(2)];

% Plot Antibody Dynamics:
if treatType~=1 && strcmp(data_vis,'Y')==0
    figure(2);clf
    FigureSize(25,20,'centimeters')
    plot(simT,A)
    xlabel('Time (days)')
    ylabel('Amount of Drug (mg)')
    title(['Drug Decay ' labelstr])
    
    % Create text box to show which parameters generated the drug curve:
    dim2=[0.15 0.15 0.13 0.12];
    str_f2={'to get this plot:',['dose = ',num2str(dose)],['\lambda = ',num2str(decayRate)]};
    h2_text=annotation(figure(2),'textbox',dim2,'String',str_f2,'FitBoxToText','on','BackgroundColor',[1 1 1]);
    h2_text.FontSize=12;
    
    xlim(axLimits(1:2))
end

% Plot BLI Dynamics:
figure(1);
hold on
ax.ColorOrderIndex = 1;
if strcmp(data_vis,'Y')==0
    for mouseNum=1:size(group,1)
        if treatType==1
            semilogy(simT,unTx_total(mouseNum,:),'LineWidth',2)
            
        elseif treatType==grp414
            t_indvMouse=t_thisExpt(~isnan(group(mouseNum,:)));
            simT_indvMouse=0:stepsize:(max(t_indvMouse));
            simTshort_indvMouse=0:ceil(max(t_thisExpt)/300/0.25)*0.25:(max(t_indvMouse));
            semilogy(simT_indvMouse,FitResults(mouseNum).SimBLI,'LineWidth',2)
            
            if strcmp(errorBars,'Y')
                if length(simT_indvMouse) < length(simTshort_indvMouse)
                    semilogy(simT_indvMouse,FitResults(mouseNum).lb_BLI,'-')
                    semilogy(simT_indvMouse,FitResults(mouseNum).ub_BLI,'-')
                else
%                     ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
                    semilogy(simTshort_indvMouse,FitResults(mouseNum).lb_BLI,'-')
%                     ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
                    semilogy(simTshort_indvMouse,FitResults(mouseNum).ub_BLI,'-')
                end
            end
        end
    end
    
    if treatType==grp414
        % Shade region of plot before treatment was implemented:
        xxx=linspace(0,DrugStart,20);yyy=logspace(12,12,20); 
        gcf; hold on; area(xxx,yyy,'FaceColor','k','FaceAlpha',0.1)
    end
end
hold on
axis(axLimits)
xticks(0:7:axLimits(2))

%% Save Plots
cd(PlotsFolder) % change file path for saving plots
if strcmp(data_vis,'Y')==1
    saveas(figure(1),['data_',PDXline,'_',Site,'_Group',num2str(treatType),'_',char(month(datetime('today'),'shortname')),num2str(year(datetime('today'))),'.png'])
    saveas(figure(1),['data_',PDXline,'_',Site,'_Group',num2str(treatType),'_',char(month(datetime('today'),'shortname')),num2str(year(datetime('today'))),'.fig'])
elseif exist('alteredDose','var') && strcmp(PDXline,'GBM39')
    saveas(figure(1),['bayesianfit_',PDXline,'_Dose',num2str(alteredDose),'_',Site,'_Group',num2str(treatType),'_',char(month(datetime('today'),'shortname')),num2str(year(datetime('today'))),'.png'])
    saveas(figure(1),['bayesianfit_',PDXline,'_Dose',num2str(alteredDose),'_',Site,'_Group',num2str(treatType),'_',char(month(datetime('today'),'shortname')),num2str(year(datetime('today'))),'.fig'])
else
    saveas(figure(1),['bayesianfit_',PDXline,'_',Site,'_Group',num2str(treatType),'_',char(month(datetime('today'),'shortname')),num2str(year(datetime('today'))),'.png'])
    saveas(figure(1),['bayesianfit_',PDXline,'_',Site,'_Group',num2str(treatType),'_',char(month(datetime('today'),'shortname')),num2str(year(datetime('today'))),'.fig'])
end
cd(CodeFolder) % return to parent directory containing this code file
end % ends the for loop (groups 1 & 3)



%%% ------------ END OF PROGRAM ------------- %%%
%%% ----------------------------------------- %%%
%%% ----------- FUNCTIONS BELOW ------------- %%%
%%% ----------------------------------------- %%%

%% Functions- No Treatment: Run Fitting, embeded logged model function
function [fitted_noTx] = runfit_noTx(initparams,tdata,ydata)
    % Define lower and upper bounds for lsqcurvefit
    %    [C0    rho] (as passed in by initparams)
    lb = [10^2  0];
    ub = [10^10 inf];
    
    % lsqcurvefit requires column vectors with no NaNs
    toEst_tdata = tdata(:);
    toEst_tdata = toEst_tdata(~isnan(ydata));
    
    toEst_ydata = ydata(:);
    toEst_ydata = toEst_ydata(~isnan(ydata));
    
    % Run lsqcurvefit for the untreated model
    fitted_noTx = lsqcurvefit(@log_untreated_model,initparams,toEst_tdata,log(toEst_ydata),lb,ub,optimset('Display','iter'));
    function Cu_logged = log_untreated_model(params,t)
        C0 = params(1);
        rho = params(2);
        Cu_logged = log(C0) + rho*t; % log of model C0*exp(rho*t)
    end
end

%% Functions- ADC Treatment: Run Simulation
function output = runModelSims(t,params,int_A)
    rho = params(1);
    q = params(2);
    muS = params(3);
    z = params(4);
    gamma = params(5);
    C0 = params(6);
    
    output = C0 * exp(rho.*t) .* (q*exp(-gamma * z * muS .* int_A) + (1 - q)*exp(-gamma * muS .* int_A));
    output(output<0)=0;
end
