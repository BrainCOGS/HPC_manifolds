
%% Get the data

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 5, [11 4], fnameStruct(7).fname, 'none', 'towers', 1, 1);
behavioralVariables = nic_output.behavioralVariables;
trials = behavioralVariables.Trial;
evi_raw = behavioralVariables.Evidence;
ch = behavioralVariables.Choice;
pos_raw = behavioralVariables.Position;
trialtype = behavioralVariables.TrialType;


%% Plot trajectories

rng(1);

figure
hold on;
xx = 0:1:300;
sourceData_S8a_x = xx';
sourceCount=1;
for tt = unique(trials)'
    pos = pos_raw(tt == trials);
    evi = evi_raw(tt == trials);
    if evi(end)<0
        cc = [0.3010, 0.7450, 0.9330];
    else
        cc = [0.8500, 0.3250, 0.0980];
    end
    f = fit(pos, evi, 'smoothingspline', 'SmoothingParam', 0.0001);
    plot(xx,f(xx), 'LineWidth', 1, 'Color', cc)
    sourceData_S8a_y(:,sourceCount) = f(xx);
    sourceCount=sourceCount+1;
end

flagL = trialtype==0;
f = fit(pos_raw(flagL), evi_raw(flagL), 'smoothingspline', 'SmoothingParam', 0.0001);
plot(xx,f(xx),'LineWidth',2, 'Color', [0, 0.4470, 0.7410])
sourceData_S8a_y_allL = f(xx);

flagR = trialtype==1;
f = fit(pos_raw(flagR), evi_raw(flagR), 'smoothingspline', 'SmoothingParam', 0.0001);
plot(xx,f(xx),'LineWidth', 2, 'Color', [0.6350, 0.0780, 0.1840])
sourceData_S8a_y_allR = f(xx);

xlim([0,300])
ylim([-16,16])
set(gca, 'box', 'off')
axis square


%% Plot flow field

rng(1);

figure;

flag = abs(evi_raw)>=0;

a = round(pos_raw(flag)/10) + 1;
b = evi_raw(flag)  + 16;
Mx = zeros(max(a),max(b));
My = zeros(max(a),max(b));
M = zeros(max(a),max(b));
for idx = 1:(length(a)-1)
    M(a(idx),b(idx))  = M(a(idx),b(idx))+1;
    My(a(idx),b(idx)) = My(a(idx),b(idx)) + (b(idx+1) - b(idx));
    Mx(a(idx),b(idx)) = Mx(a(idx),b(idx)) + (a(idx+1) - a(idx));
end
Mx(end,:) = 0;
My(end,:) = 0;
Xx = 6*Mx./M;
Yy = 1*My./M;
Xx(isnan(Xx)) = 0;
Yy(isnan(Yy)) = 0;
Xx_smo = imgaussfilt(Xx',1);
Yy_smo = imgaussfilt(Yy',1);
[X,Y] = meshgrid(0:10:300, -15:1:13);
Xx_smo(Xx' == 0) =0;
Yy_smo(Yy' == 0) =0;
s = 2;
q = quiver(X(1:s:end,1:s:end), Y(1:s:end,1:s:end), Xx_smo(1:s:end,1:s:end), Yy_smo(1:s:end,1:s:end));
sourceData_S8b_x_coordinate = X(1:s:end,1:s:end);
sourceData_S8b_y_coordinate = Y(1:s:end,1:s:end);
sourceData_S8b_x_direction  = Xx_smo(1:s:end,1:s:end);
sourceData_S8b_y_direction  = Yy_smo(1:s:end,1:s:end);
q.MaxHeadSize = .02;

xlim([0,300])
ylim([-16,16])
set(gca, 'box', 'off')
axis square

