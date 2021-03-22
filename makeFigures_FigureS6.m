
%% The view angle manifold plot is taken from the code for Figure 3

%% Randomize one dimension and decode

% First run this on matlab on spock
% mind_runEncodingJobs_testVar('towers',"{'ViewAngle', 'Evidence'}", '/jukebox/tank/enieh/mind/FINAL/Towers/Encoding/');
% Do the same for Position + Y Velocity
% Do the same for Position + Time

filelocation = 'D:\Encoding\';
% filelocation = 'M:\enieh\mind\FINAL\Towers\Encoding\';
animals = {'E22', 'E39', 'E43', 'E44', 'E47', 'E48', 'E65'};

% View angle and evidence
for i=1:length(animals)
    load([filelocation animals{i} '_ViewAngleEvidence_testVar_standardizeF.mat']);
    outputRegressOutViewAngle2_All(i) = outputRegressOutViewAngle2;
    all_coefs(i,:) = outputRegressOutViewAngle2.all_coefs;
end

disp(signrank(all_coefs(:,2), all_coefs(:,1)));  % significance across animals

disp(mean( (all_coefs(:,2) - all_coefs(:,1))./all_coefs(:,1) ));

disp(std( (all_coefs(:,2) - all_coefs(:,1))./all_coefs(:,1) ) / sqrt(length(all_coefs)));

figure;
hold on;
nieh_barSEMpaired(all_coefs(:,1), all_coefs(:,2))
ylabel("Decoding Index")
xticks([1 2]);
xticklabels({'VA+R_E', 'VA+E'});
title(['Sign rank p-value is: ' num2str(signrank(all_coefs(:,2), all_coefs(:,1)))]);


% Position and Y velocity
for i=1:length(animals)
    load([filelocation animals{i} '_PositionYvelocity_testVar_standardizeF.mat']);
    outputRegressOutViewAngle2_All(i) = outputRegressOutViewAngle2;
    all_coefs(i,:) = outputRegressOutViewAngle2.all_coefs;
end

disp(signrank(all_coefs(:,2), all_coefs(:,1)));  % significance across animals

disp(mean( (all_coefs(:,2) - all_coefs(:,1))./all_coefs(:,1) ));

disp(std( (all_coefs(:,2) - all_coefs(:,1))./all_coefs(:,1) ) / sqrt(length(all_coefs)));

figure;
hold on;
nieh_barSEMpaired(all_coefs(:,1), all_coefs(:,2))
ylabel("Decoding Index")
xticks([1 2]);
xticklabels({'Y+R_yV', 'Y+yV'});
title(['Sign rank p-value is: ' num2str(signrank(all_coefs(:,2), all_coefs(:,1)))]);

% Position and Time
for i=1:length(animals)
    load([filelocation animals{i} '_PositionTime_testVar_standardizeF.mat']);
    outputRegressOutViewAngle2_All(i) = outputRegressOutViewAngle2;
    all_coefs(i,:) = outputRegressOutViewAngle2.all_coefs;
end

disp(signrank(all_coefs(:,2), all_coefs(:,1)));  % significance across animals

disp(mean( (all_coefs(:,2) - all_coefs(:,1))./all_coefs(:,1) ));

disp(std( (all_coefs(:,2) - all_coefs(:,1))./all_coefs(:,1) ) / sqrt(length(all_coefs)));

figure;
hold on;
nieh_barSEMpaired(all_coefs(:,1), all_coefs(:,2))
ylabel("Decoding Index")
xticks([1 2]);
xticklabels({'Y+R_T', 'Y+T'});
title(['Sign rank p-value is: ' num2str(signrank(all_coefs(:,2), all_coefs(:,1)))]);


%% Decode the PC1 and PC2 of E and VA from 5-dim manifold

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

outputRegressOutViewAngle1_EandVA = mind_decodePCAVariable(fnameStruct,'towers', {'Evidence', 'ViewAngle'});


%% Comparing decoding of View Angle in Alternation Task vs Towers task

mind_script_compareTowersAlternateViewAngleSLIM;


%% Distributions of View Angle values (histogram version)

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');
for i=1:length(fnameStruct)
    
    fname = fnameStruct(i).fname;
    nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 5, [11 4], fname,'none','towers',1,1);
    behavioralVariables = nic_output.behavioralVariables;
    VA{i} = rad2deg(behavioralVariables.ViewAngle);
    
end

fnameStruct_Alternation = mind_makeFnameStruct('Edward','Alternation','none');
for i=1:length(fnameStruct_Alternation)
    
    fname = fnameStruct_Alternation(i).fname;
    nic_output = extractVariables('all', 2, 'goodTrials', 2, 0, 0, 5, [11 4], fname,'none','alternation',1,1);
    behavioralVariables = nic_output.behavioralVariables;
    VA_ALT{i} = rad2deg(behavioralVariables.ViewAngle);
    
end


for i=1:length(VA)
h = histogram(VA{i}, [-82.5:5:82.5], 'Normalization', 'probability');
prob1(i,:) = h.BinCounts/sum(h.BinCounts);
end

for i=1:size(prob1,2)
   sem1(i) = nieh_sem(prob1(:,i)); 
end

centers = [-82.5:5:77.5]+2.5;
mean1 = mean(prob1);


for i=1:length(VA_ALT)
h = histogram(VA_ALT{i}, [-82.5:5:82.5], 'Normalization', 'probability');
prob1A(i,:) = h.BinCounts/sum(h.BinCounts);
end

for i=1:size(prob1A,2)
   sem1A(i) = nieh_sem(prob1A(:,i)); 
end

mean1A = mean(prob1A);

tow_CI = nieh_ci(prob1,1000);
alt_CI = nieh_ci(prob1A,1000);

close all;

figure; 
hold on;
shadedErrorBar_CI(centers,mean1,tow_CI,'lineProps', 'b'); 
shadedErrorBar_CI(centers,mean1A,alt_CI,'lineProps', 'r'); 
plot(centers,prob1A','g')
plot(centers,prob1','y')
legend('Towers','Alternation');
xlabel('View Angle');
ylabel('Probability');
title('Blue is Towers, Red is Alt, Green is Alt ind, Yellow is Tow Ind');


%% View Angle Distributions, per position bin, individual lines

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');
fnameStruct_Alternation = mind_makeFnameStruct('Edward','Alternation','none');

outputMakeVAPlot = mind_makeVAplot(fnameStruct,'towers');
outputMakeVAPlot_Alternation = mind_makeVAplot(fnameStruct_Alternation,'Alternation');

figure;
for i=1:length(outputMakeVAPlot.fnameStruct)
   hold on;
   plot([1:31],outputMakeVAPlot.meanL(i,:),'Color','g','LineWidth', .1)
   plot([1:31],outputMakeVAPlot.meanR(i,:),'Color','g','LineWidth', .1)
end
for i=1:length(outputMakeVAPlot_Alternation.fnameStruct)
   hold on;
   plot([1:31],outputMakeVAPlot_Alternation.meanL(i,:),'Color','r','LineWidth', .1)
   plot([1:31],outputMakeVAPlot_Alternation.meanR(i,:),'Color','r','LineWidth', .1)
end
plot([1:31],mean(outputMakeVAPlot.meanL),'Color', 'g','LineWidth', 3);
plot([1:31],mean(outputMakeVAPlot.meanR),'Color', 'g','LineWidth', 3);
plot([1:31],mean(outputMakeVAPlot_Alternation.meanL),'Color', 'r','LineWidth', 3);
plot([1:31],mean(outputMakeVAPlot_Alternation.meanR),'Color', 'r','LineWidth', 3);

xlim([1 31])
ylim([-80 80])
xticks([1 11 21 31])
xticklabels({'0' , '100', '200', '300'});
xlabel('Position (cm)');
ylabel('View Angle (degrees)');
title('Shaded area is 95% CI, Green=Towers, Red=Alternation');


%% Binary Decoding


fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

% Can run this on the laptop OR use the spock code below
outputNonlinearDecoding_VarGroup_Binary_dim2 = mind_nonlinearDecoding_VarGroup(fnameStruct, 5, 'GP', {'Choice', 'PriorChoice','PriorCorrect'}, 'towers', 0, 1, 2);
outputNonlinearDecoding_VarGroup_Binary_dim3 = mind_nonlinearDecoding_VarGroup(fnameStruct, 5, 'GP', {'Choice', 'PriorChoice','PriorCorrect'}, 'towers', 0, 1, 3);
save('outputNonlinearDecoding_BinaryVariables_withPosVal.mat', '-v7.3');
outputNonlinearDecoding_VarGroup_Binary_dim4 = mind_nonlinearDecoding_VarGroup(fnameStruct, 5, 'GP', {'Choice', 'PriorChoice','PriorCorrect'}, 'towers', 0, 1, 4);
outputNonlinearDecoding_VarGroup_Binary_dim5 = mind_nonlinearDecoding_VarGroup(fnameStruct, 5, 'GP', {'Choice', 'PriorChoice','PriorCorrect'}, 'towers', 0, 1, 5);
save('outputNonlinearDecoding_BinaryVariables_withPosVal.mat', '-v7.3');
outputNonlinearDecoding_VarGroup_Binary_dim6 = mind_nonlinearDecoding_VarGroup(fnameStruct, 5, 'GP', {'Choice', 'PriorChoice','PriorCorrect'}, 'towers', 0, 1, 6);
outputNonlinearDecoding_VarGroup_Binary_dim7 = mind_nonlinearDecoding_VarGroup(fnameStruct, 5, 'GP', {'Choice', 'PriorChoice','PriorCorrect'}, 'towers', 0, 1, 7);
save('outputNonlinearDecoding_BinaryVariables_withPosVal.mat', '-v7.3');

% For the individual points figures, used the .sh file on spock
% sbatch --time=300 --mem=60000 -c 4  -o "/jukebox/tank/enieh/mind/logs/mind_nonlinearDecoding_VarGroup_dim7_%j.log" -p all /jukebox/tank/enieh/mind/OLD_Unused/mind_nonlinearDecoding_VarGroup_dim7.sh
load('M:\enieh\mind\outputNonlinearDecoding_VarGroup_Binary_dim4.mat')
load('M:\enieh\mind\outputNonlinearDecoding_VarGroup_Binary_dim5.mat')
load('M:\enieh\mind\outputNonlinearDecoding_VarGroup_Binary_dim2.mat')
load('M:\enieh\mind\outputNonlinearDecoding_VarGroup_Binary_dim3.mat')
load('M:\enieh\mind\outputNonlinearDecoding_VarGroup_Binary_dim7.mat')
load('M:\enieh\mind\outputNonlinearDecoding_VarGroup_Binary_dim6.mat')

% Make the plot
figure;
subplot(1,3,1)
nieh_barSEM(outputNonlinearDecoding_VarGroup_Binary_dim2(1).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim3(1).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim4(1).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim5(1).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim6(1).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim7(1).meancorPerList);
        hold on;
        scatter([ones(length(fnameStruct),1); ...
                 ones(length(fnameStruct),1)*2; ...
                 ones(length(fnameStruct),1)*3; ...
                 ones(length(fnameStruct),1)*4; ...
                 ones(length(fnameStruct),1)*5; ...
                 ones(length(fnameStruct),1)*6; ...
                 ],[outputNonlinearDecoding_VarGroup_Binary_dim2(1).meancorPerList ...
                    outputNonlinearDecoding_VarGroup_Binary_dim3(1).meancorPerList ...
                    outputNonlinearDecoding_VarGroup_Binary_dim4(1).meancorPerList ...
                    outputNonlinearDecoding_VarGroup_Binary_dim5(1).meancorPerList ...
                    outputNonlinearDecoding_VarGroup_Binary_dim6(1).meancorPerList ...
                    outputNonlinearDecoding_VarGroup_Binary_dim7(1).meancorPerList ...
                    ], '.');
        
xticklabels({'2', '3', '4', '5', '6', '7'});
xlabel('Dims Embedded');
ylabel('Accuracy')
ylim([.4 1]);
set(gca, 'box', 'off')
title('Choice');

subplot(1,3,2)
nieh_barSEM(outputNonlinearDecoding_VarGroup_Binary_dim2(2).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim3(2).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim4(2).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim5(2).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim6(2).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim7(2).meancorPerList);
        hold on;
        scatter([ones(length(fnameStruct),1); ...
                 ones(length(fnameStruct),1)*2; ...
                 ones(length(fnameStruct),1)*3; ...
                 ones(length(fnameStruct),1)*4; ...
                 ones(length(fnameStruct),1)*5; ...
                 ones(length(fnameStruct),1)*6; ...
                 ],[outputNonlinearDecoding_VarGroup_Binary_dim2(2).meancorPerList ...
                    outputNonlinearDecoding_VarGroup_Binary_dim3(2).meancorPerList ...
                    outputNonlinearDecoding_VarGroup_Binary_dim4(2).meancorPerList ...
                    outputNonlinearDecoding_VarGroup_Binary_dim5(2).meancorPerList ...
                    outputNonlinearDecoding_VarGroup_Binary_dim6(2).meancorPerList ...
                    outputNonlinearDecoding_VarGroup_Binary_dim7(2).meancorPerList ...
                    ], '.');
xticklabels({'2', '3', '4', '5', '6', '7'});
xlabel('Dims Embedded');
ylabel('Accuracy')
ylim([.4 1]);
set(gca, 'box', 'off')
title('Prior Choice');

subplot(1,3,3)
nieh_barSEM(outputNonlinearDecoding_VarGroup_Binary_dim2(3).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim3(3).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim4(3).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim5(3).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim6(3).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim7(3).meancorPerList);
        hold on;
        scatter([ones(length(fnameStruct),1); ...
                 ones(length(fnameStruct),1)*2; ...
                 ones(length(fnameStruct),1)*3; ...
                 ones(length(fnameStruct),1)*4; ...
                 ones(length(fnameStruct),1)*5; ...
                 ones(length(fnameStruct),1)*6; ...
                 ],[outputNonlinearDecoding_VarGroup_Binary_dim2(3).meancorPerList ...
                    outputNonlinearDecoding_VarGroup_Binary_dim3(3).meancorPerList ...
                    outputNonlinearDecoding_VarGroup_Binary_dim4(3).meancorPerList ...
                    outputNonlinearDecoding_VarGroup_Binary_dim5(3).meancorPerList ...
                    outputNonlinearDecoding_VarGroup_Binary_dim6(3).meancorPerList ...
                    outputNonlinearDecoding_VarGroup_Binary_dim7(3).meancorPerList ...
                    ], '.');
xticklabels({'2', '3', '4', '5', '6', '7'});
xlabel('Dims Embedded');
ylabel('Accuracy')
ylim([.4 1]);
set(gca, 'box', 'off')
title('Prior Correct');
