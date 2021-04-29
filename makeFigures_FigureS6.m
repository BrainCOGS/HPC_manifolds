
%% The view angle manifold plot is taken from the code for Figure 3

%% Randomize one dimension and decode

% First run this on matlab on spock
% mind_runEncodingJobs_testVar('towers',"{'ViewAngle', 'Evidence'}", '/jukebox/tank/enieh/mind/FINAL/Towers/Encoding/');
% Do the same for Position + Y Velocity
% Do the same for Position + Time

%filelocation = 'D:\Encoding\';
%filelocation = 'M:\enieh\mind\FINAL\Towers\Encoding\';
filelocation = 'C:\Neuroscience\imaging\FINAL\encoding_Data\';
animals = {'E22', 'E39', 'E43', 'E44', 'E47', 'E48', 'E65'};

% View angle and evidence
for i=1:length(animals)
    load([filelocation animals{i} '_ViewAngleEvidence_testVar_standardizeF.mat']);
    outputRegressOutViewAngle2_All(i) = outputRegressOutViewAngle2;
    all_coefs(i,:) = outputRegressOutViewAngle2.all_coefs;
end

disp(signrank(all_coefs(:,2), all_coefs(:,1)));  

figure;
hold on;
nieh_barSEMpaired(all_coefs(:,1), all_coefs(:,2))
ylabel("Decoding Index")
xticks([1 2]);
xticklabels({'VA+R_E', 'VA+E'});
title(['Sign rank p-value is: ' num2str(signrank(all_coefs(:,2), all_coefs(:,1)))]);

sourceData_S6b = all_coefs;

% Position and Y velocity
for i=1:length(animals)
    load([filelocation animals{i} '_PositionYvelocity_testVar_standardizeF.mat']);
    outputRegressOutViewAngle2_All(i) = outputRegressOutViewAngle2;
    all_coefs(i,:) = outputRegressOutViewAngle2.all_coefs;
end

disp(signrank(all_coefs(:,2), all_coefs(:,1)));  % significance across animals

figure;
hold on;
nieh_barSEMpaired(all_coefs(:,1), all_coefs(:,2))
ylabel("Decoding Index")
xticks([1 2]);
xticklabels({'Y+R_yV', 'Y+yV'});
title(['Sign rank p-value is: ' num2str(signrank(all_coefs(:,2), all_coefs(:,1)))]);

sourceData_S6c = all_coefs;

% Position and Time
for i=1:length(animals)
    load([filelocation animals{i} '_PositionTime_testVar_standardizeF.mat']);
    outputRegressOutViewAngle2_All(i) = outputRegressOutViewAngle2;
    all_coefs(i,:) = outputRegressOutViewAngle2.all_coefs;
end

disp(signrank(all_coefs(:,2), all_coefs(:,1)));  % significance across animals

figure;
hold on;
nieh_barSEMpaired(all_coefs(:,1), all_coefs(:,2))
ylabel("Decoding Index")
xticks([1 2]);
xticklabels({'Y+R_T', 'Y+T'});
title(['Sign rank p-value is: ' num2str(signrank(all_coefs(:,2), all_coefs(:,1)))]);

sourceData_S6d = all_coefs;


%% Decode the PC1 and PC2 of E and VA from 5-dim manifold

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');
varTypeList = {'Evidence', 'ViewAngle'};

load("C:\Neuroscience\imaging\FINAL\decoding_Data\decodePC1and2_EandVA.mat")
% Or run this code:
% outputRegressOutViewAngle1_EandVA = mind_decodePCAVariable(fnameStruct,'towers', varTypeList);

figure; 
nieh_barSEMpaired(outputRegressOutViewAngle1_EandVA.all_coefs(:,1), outputRegressOutViewAngle1_EandVA.all_coefs(:,2))
sourceData_S6e = outputRegressOutViewAngle1_EandVA.all_coefs;
ylabel('Decoding index (r)');
xlabel('PC');
set(gca, 'box', 'off')
title(['Manifold Decoding PC1 and PC2 of PCA on ', varTypeList{1}, ' and ', varTypeList{2} ])

%% Comparing decoding of View Angle in Alternation Task vs Towers task

fnameStruct_Alternation = mind_makeFnameStruct('Edward','Alternation','NONE');
fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

load("C:\Neuroscience\imaging\FINAL\decoding_Data\decodeVATowersAlternation.mat");
% Or run this code:
% outputCompareTowersAlternationVA = mind_compareTowersAlternationVA(fnameStruct, fnameStruct_Alternation);

figure;
nieh_barSEM(outputCompareTowersAlternationVA.meancorrAll_VA, outputCompareTowersAlternationVA.meancorrAll_VA_ALT);
sourceData_S6f = [outputCompareTowersAlternationVA.meancorrAll_VA' outputCompareTowersAlternationVA.meancorrAll_VA_ALT'];
hold on;
scatter([ones(length(outputCompareTowersAlternationVA.meancorrAll_VA),1); ones(length(outputCompareTowersAlternationVA.meancorrAll_VA_ALT),1)*2],[outputCompareTowersAlternationVA.meancorrAll_VA outputCompareTowersAlternationVA.meancorrAll_VA_ALT], '.');
xticklabels({'Towers - View Angle', 'Alternation - View Angle'});
xtickangle(45)
ylabel('Decoding index (r)');
set(gca, 'box', 'off')
title(['Rank sum p-value is: ' num2str(ranksum(outputCompareTowersAlternationVA.meancorrAll_VA, outputCompareTowersAlternationVA.meancorrAll_VA_ALT))]);


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

close all;

rng(1)
tow_CI = nieh_ci(prob1,1000);
alt_CI = nieh_ci(prob1A,1000);

% Flip it to make it work in shadedErrorBar
tow_CI_fixed = [tow_CI(2,:); tow_CI(1,:)];
alt_CI_fixed = [alt_CI(2,:); alt_CI(1,:)];

tow_CI_fixed(1,:) = tow_CI_fixed(1,:)-mean1;
tow_CI_fixed(2,:) = mean1-tow_CI_fixed(2,:);

alt_CI_fixed(1,:) = alt_CI_fixed(1,:)-mean1A;
alt_CI_fixed(2,:) = mean1A-alt_CI_fixed(2,:);

figure; 
hold on;
shadedErrorBar(centers,mean1,tow_CI_fixed,'lineProps', 'b'); 
shadedErrorBar(centers,mean1A,alt_CI_fixed,'lineProps', 'r'); 
plot(centers,prob1A','g')
plot(centers,prob1','y')
sourceData_S6h_Tow_mean = [centers' mean1'  tow_CI_fixed(2,:)' tow_CI_fixed(1,:)'];
sourceData_S6h_Tow_ind  = prob1';
sourceData_S6h_Alt_mean = [centers' mean1A' alt_CI_fixed(2,:)' alt_CI_fixed(1,:)'];
sourceData_S6h_Alt_ind  = prob1A';
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

sourceData_S6g_Tow_meanL =  mean(outputMakeVAPlot.meanL)';
sourceData_S6g_Tow_meanR =  mean(outputMakeVAPlot.meanR)';
sourceData_S6g_Alt_meanL =  mean(outputMakeVAPlot_Alternation.meanL)';
sourceData_S6g_Alt_meanR =  mean(outputMakeVAPlot_Alternation.meanR)';
sourceData_S6g_Position  = [0:10:300]';
sourceData_S6g_Tow_indL  = outputMakeVAPlot.meanL';
sourceData_S6g_Tow_indR  = outputMakeVAPlot.meanR';
sourceData_S6g_Alt_indL  = outputMakeVAPlot_Alternation.meanL';
sourceData_S6g_Alt_indR  = outputMakeVAPlot_Alternation.meanR';

xlim([1 31])
ylim([-80 80])
xticks([1 11 21 31])
xticklabels({'0' , '100', '200', '300'});
xlabel('Position (cm)');
ylabel('View Angle (degrees)');
title('Green=Towers, Red=Alternation');


%% Binary Decoding

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

load('C:\Neuroscience\imaging\FINAL\decoding_Data\decodeBinary_all.mat')

% Or run this code
% outputNonlinearDecoding_VarGroup_Binary_dim2 = mind_nonlinearDecoding_VarGroup(fnameStruct, 5, 'GP', {'Choice', 'PriorChoice','PriorCorrect'}, 'towers', 0, 1, 2);
% outputNonlinearDecoding_VarGroup_Binary_dim3 = mind_nonlinearDecoding_VarGroup(fnameStruct, 5, 'GP', {'Choice', 'PriorChoice','PriorCorrect'}, 'towers', 0, 1, 3);
% save('outputNonlinearDecoding_BinaryVariables_withPosVal.mat', '-v7.3');
% outputNonlinearDecoding_VarGroup_Binary_dim4 = mind_nonlinearDecoding_VarGroup(fnameStruct, 5, 'GP', {'Choice', 'PriorChoice','PriorCorrect'}, 'towers', 0, 1, 4);
% outputNonlinearDecoding_VarGroup_Binary_dim5 = mind_nonlinearDecoding_VarGroup(fnameStruct, 5, 'GP', {'Choice', 'PriorChoice','PriorCorrect'}, 'towers', 0, 1, 5);
% save('outputNonlinearDecoding_BinaryVariables_withPosVal.mat', '-v7.3');
% outputNonlinearDecoding_VarGroup_Binary_dim6 = mind_nonlinearDecoding_VarGroup(fnameStruct, 5, 'GP', {'Choice', 'PriorChoice','PriorCorrect'}, 'towers', 0, 1, 6);
% outputNonlinearDecoding_VarGroup_Binary_dim7 = mind_nonlinearDecoding_VarGroup(fnameStruct, 5, 'GP', {'Choice', 'PriorChoice','PriorCorrect'}, 'towers', 0, 1, 7);
% save('outputNonlinearDecoding_BinaryVariables_withPosVal.mat', '-v7.3');

% *************** OLD ***************************
% For the individual points figures, used the .sh file on spock
% sbatch --time=300 --mem=60000 -c 4  -o "/jukebox/tank/enieh/mind/logs/mind_nonlinearDecoding_VarGroup_dim7_%j.log" -p all /jukebox/tank/enieh/mind/OLD_Unused/mind_nonlinearDecoding_VarGroup_dim7.sh
% load('M:\enieh\mind\outputNonlinearDecoding_VarGroup_Binary_dim4.mat')
% load('M:\enieh\mind\outputNonlinearDecoding_VarGroup_Binary_dim5.mat')
% load('M:\enieh\mind\outputNonlinearDecoding_VarGroup_Binary_dim2.mat')
% load('M:\enieh\mind\outputNonlinearDecoding_VarGroup_Binary_dim3.mat')
% load('M:\enieh\mind\outputNonlinearDecoding_VarGroup_Binary_dim7.mat')
% load('M:\enieh\mind\outputNonlinearDecoding_VarGroup_Binary_dim6.mat')
% ****************** OLD END ************************

% Make the plot
figure;
subplot(1,3,1)
nieh_barSEM(outputNonlinearDecoding_VarGroup_Binary_dim2(1).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim3(1).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim4(1).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim5(1).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim6(1).meancorPerList, ...
            outputNonlinearDecoding_VarGroup_Binary_dim7(1).meancorPerList);
sourceData_S6i_Choice = [outputNonlinearDecoding_VarGroup_Binary_dim2(1).meancorPerList*100; ...
                         outputNonlinearDecoding_VarGroup_Binary_dim3(1).meancorPerList*100; ...
                         outputNonlinearDecoding_VarGroup_Binary_dim4(1).meancorPerList*100; ...
                         outputNonlinearDecoding_VarGroup_Binary_dim5(1).meancorPerList*100; ...
                         outputNonlinearDecoding_VarGroup_Binary_dim6(1).meancorPerList*100; ...
                         outputNonlinearDecoding_VarGroup_Binary_dim7(1).meancorPerList*100];
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
sourceData_S6i_PrevChoice = [outputNonlinearDecoding_VarGroup_Binary_dim2(2).meancorPerList*100; ...
                             outputNonlinearDecoding_VarGroup_Binary_dim3(2).meancorPerList*100; ...
                             outputNonlinearDecoding_VarGroup_Binary_dim4(2).meancorPerList*100; ...
                             outputNonlinearDecoding_VarGroup_Binary_dim5(2).meancorPerList*100; ...
                             outputNonlinearDecoding_VarGroup_Binary_dim6(2).meancorPerList*100; ...
                             outputNonlinearDecoding_VarGroup_Binary_dim7(2).meancorPerList*100];        
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
sourceData_S6i_PrevCorrect = [outputNonlinearDecoding_VarGroup_Binary_dim2(3).meancorPerList*100; ...
                              outputNonlinearDecoding_VarGroup_Binary_dim3(3).meancorPerList*100; ...
                              outputNonlinearDecoding_VarGroup_Binary_dim4(3).meancorPerList*100; ...
                              outputNonlinearDecoding_VarGroup_Binary_dim5(3).meancorPerList*100; ...
                              outputNonlinearDecoding_VarGroup_Binary_dim6(3).meancorPerList*100; ...
                              outputNonlinearDecoding_VarGroup_Binary_dim7(3).meancorPerList*100];                
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

%% Correlation between view angle and evidence

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

for i=1:length(fnameStruct)
    nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 5, [11 4], fnameStruct(i).fname,'none','towers',1,1);
    evidence = nic_output.behavioralVariables.Evidence;
    viewangle = nic_output.behavioralVariables.ViewAngle;
    corr1 = corrcoef(evidence, -viewangle);
    VA_E_corr(i) = corr1(1,2);
end
mean(VA_E_corr)
nieh_sem(VA_E_corr)