
function outputFitPowerLaw = mind_fitPowerLaw_Function(fnameStruct, taskType, toggleTrialNormalization, expPlot)

% If toggleTrialNormalization is 1, it will also make the plot with the
% trials normalized to the mean trial length

%% Set up input variables, RNG, and config
outputFitPowerLaw.fnameStruct = fnameStruct;
outputFitPowerLaw.taskType    = taskType;

inputRNG = 42;
rng(inputRNG);
outputFitPowerLaw.RNG         = inputRNG;

fitRange = [0:0.001:1];
outputFitPowerLaw.fitRange    = fitRange;

outputFitPowerLaw.toggleTrialNormalization = toggleTrialNormalization;

% For the flag of the linear and polynomial fitting
flagRange = [3 6];
config.flagRange = flagRange;

% For the trial length normalization
xInputs = logspace(log10(0.0001),log10(0.4), 100);
config.xInputs = xInputs;

% expPlot = 5;
config.expPlot = expPlot;

outputFitPowerLaw.config = config;

%% Get outLeaves data

for i=1:length(fnameStruct)
    load(fnameStruct(i).fname_mani);
    tic
    outLeaves{i} = mind_findMinLeaves(outMind.dat, fitRange);
    toc
    clear outMind
    i
end

outputFitPowerLaw.outLeaves = outLeaves;

%% Linear and Polynomial Fit

options = optimoptions('fsolve'); 
options.MaxIterations = 10000;
options.MaxFunctionEvaluations = 50000;
outputFitPowerLaw.fitOptions = options;

exponents_l = [];
exponents_c = [];


for i=1:length(fnameStruct)
    
    ydat = outLeaves{i}.Cs;
    
    % Linear Portion
    xdata = log10(fitRange);
    ydata = log10(ydat');
    
    flag = (ydata>flagRange(1)) & (ydata<flagRange(2));
    fun = @(x,xdata) x(1)*xdata + x(2);
    x0 = [1,1];
    xfit = lsqcurvefit(fun,x0,xdata(flag),ydata(flag), [], [], options);
    xx = xdata(flag);
    
    ydataAll(:,i) = ydata;
    xdataAll(:,i) = xdata;
    funFitAll{i}  = fun(xfit,xx);
    xxAll{i}      = xx;
    xFitAll(:,i)   = xfit;
    exponents_l = [exponents_l, xfit(1)];
    
    
    % Polynomial Portion
    xdata2 = fitRange;
    ydata2 = ydat';
    
    flag2 = (ydata2>10^4) & (ydata2<10^7);
    fun2 = @(x,xdata2) x(1)* (xdata2.^x(2));
    x02 = [11,1];
    xfit2 = lsqcurvefit(fun2,x02,xdata2(flag2),ydata2(flag2), [], [], options);
    xx2 = xdata2(flag2);
    
    ydataAll2(:,i) = ydata2;
    xdataAll2(:,i) = xdata2;
    funFitAll2{i}  = fun2(xfit2,xx2);
    xxAll2{i}      = xx2;
    xFitAll2(:,i)  = xfit2;
    flag2All{i}    = flag2;
    exponents_c = [exponents_c, xfit2(2)];
    i
end

mean_l = mean(exponents_l);
sem_l  = nieh_sem(exponents_l);
ci_l   = nieh_ci(exponents_l',1000);
mean_c = mean(exponents_c);
sem_c  = nieh_sem(exponents_c);
ci_c   = nieh_ci(exponents_c',1000);

figure;
subplot(1,2,1)
hold on;
for i=1:length(fnameStruct)
    plot(xxAll{i}, funFitAll{i},'k-','linewidth',1)
    plot(xdataAll(:,i),ydataAll(:,i))
end
xlabel('log10(d)')
ylabel('log10(N_neighbors')
title(['mean: ', num2str(mean_l), ' SEM: ', num2str(sem_l), ' CI: [' , num2str(ci_l'), ']'])

subplot(1,2,2)
hold on;
for i=1:length(fnameStruct)
    plot(xxAll2{i}, funFitAll2{i},'k--')
    plot(xxAll2{i}, ydataAll2(flag2All{i},i))
end
xlabel('d')
ylabel('N_neighbors')
title(['mean: ', num2str(mean_c), ' SEM: ', num2str(sem_c), ' CI: [' , num2str(ci_c'), ']'])
set(gcf, 'Position',  [1000, 900, 1300, 400])

% disp(['mean: ', num2str(mean(exponents_l)), 'SEM: ', num2str(nieh_sem(exponents_l))])
% disp(['mean: ', num2str(mean(exponents_c)), 'SEM: ', num2str(nieh_sem(exponents_c))])

outputFitPowerLaw.exponents_l = exponents_l;
outputFitPowerLaw.exponents_c = exponents_c;
outputFitPowerLaw.mean_l  = mean_l;
outputFitPowerLaw.sem_l   = sem_l;
outputFitPowerLaw.ci_l    = ci_l;
outputFitPowerLaw.mean_c  = mean_c;
outputFitPowerLaw.sem_c   = sem_c;
outputFitPowerLaw.ci_c    = ci_c;

%% Calculate N as function of distance for all landmarks

if toggleTrialNormalization==1
    
    figure
    hold on;
    for i=1:length(fnameStruct)
        
        if i==1
            c = [0, 0.4470, 0.7410];
        elseif i==2
            c = [0.8500, 0.3250, 0.0980];
        elseif i==3
            c = [0.9290, 0.6940, 0.1250];
        elseif i==4
            c = [0.4940, 0.1840, 0.5560];
        elseif i==5
            c = [0.4660, 0.6740, 0.1880];
        elseif i==6
            c = [0.3010, 0.7450, 0.9330];
        elseif i==7
            c = [0.6350, 0.0780, 0.1840];
        end
        
        load(fnameStruct(i).fname_mani);
        data = outMind.dat.forestdat.rwd.Dg;
        [triallength, ~] = TrialLengthFromD(fnameStruct(i).fname, taskType, outMind);
        
        
        Cs = zeros(length(xInputs),length(data));
        for p_idx = 1:length(data)
            distances = data(p_idx,:);
            distances(p_idx) = inf; %ignore distance to itself
            for idx = 1:length(xInputs)
                Cs(idx,p_idx) = sum( distances<xInputs(idx) );
            end
            disp(p_idx/length(data))
        end
        
        for idx = 1:length(xInputs)
            Nneighbors = Cs(idx,:);
            m = bootstrp(1000,@mean,Nneighbors);
            m = sort(m);
            best = mean(Cs(idx,:));
            if m(25) > 0
                le = best-m(25);
            else
                le = best-1e-6;
            end
            ue = m(975)-best;
            errorbar(xInputs(idx)./(triallength/2.), best, le, ue, '.', 'MarkerFaceColor', c, 'MarkerEdgeColor', c , 'Color',c, 'LineWidth', 1, 'CapSize', 0)
            disp(xInputs(idx))
        end
        
        tempCs{i} = mean(Cs,2);
        temptriallength{i} = triallength;
        plot(xInputs./(triallength/2.), mean(Cs,2), 'Color', c, 'LineWidth', 2)
    end
        set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    
%     figure
%     hold on;
%     for i=1:length(tempCs)
%         
%         plot(xInputs./(temptriallength{i}/2.), tempCs{i}, 'LineWidth', 2)
%         set(gca, 'YScale', 'log')
%         set(gca, 'XScale', 'log')
%     end

    for s = 10.^[-4,-2,0,2,4,6,8,10,12,14]
        loglog([1e-3,1e-2,1,10,100],s.*[1e-3,1e-2,1,10,100].^expPlot,'--','Color',[0.5,0.5,0.5])
    end
    xlim([0.05,4])
    ylim([5e-5,1e5])
    axis square
    xlabel('D/Do');
    % Distance between points in units of length of trial
    ylabel('Average number of neighboring landmarks')
    title(['Dotted Lines are exponent = ' num2str(expPlot)])
    
    outputFitPowerLaw.tempCs = tempCs;
    outputFitPowerLaw.temptriallength = temptriallength;
    
end

