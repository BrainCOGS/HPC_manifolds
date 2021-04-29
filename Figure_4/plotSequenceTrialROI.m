function [sourceData_4a_Time, sourceData_4a_Position] = plotSequenceTrialROI(doubletNum, outputDoublets, fname, timestep)

%roi1 is in normal values, not global

if timestep>1
    warndlg('Are you sure that is timestep and not fs?');
end

aTrials = outputDoublets.saveAll_basics.data_trialnrs;

%% Extract rois and trials

doubletNum_1 = outputDoublets.saveAll_doublets(doubletNum).first_cell;
doubletNum_2 = outputDoublets.saveAll_doublets(doubletNum).second_cell;

roi1 = [doubletNum_1 doubletNum_2];
numtrial = outputDoublets.saveAll_doublets(doubletNum).trials_appear;

%% Get Data
nic_output = extractVariables(roi1, 2, aTrials(numtrial)', 2, 100, 0, 5, [0 0], fname,'none','towers',1,1);
ROI_trial = reshape(nic_output.ROIactivities,[101 length(numtrial) length(roi1)]);

nic_output2 = extractVariables(roi1, 2, aTrials(numtrial)',2,0,0,5,[0 0],fname,'none','towers',1,1);
ROI_time = nic_output2.ROIactivities;
behavioralVariables = nic_output2.behavioralVariables;
trialn = behavioralVariables.Trial;

%% Plot
figure;
subplot(1,2,1)
hold on;
sCount = 1;
for j=1:length(numtrial)
    hold on
    for i=1:length(roi1)
        plotdata = ROI_time(trialn==aTrials(numtrial(j)),i);%score.dataDFF(begin1:end1,roi1(i));
        plotdata_smooth = mind_preprocess(plotdata,5,0);
        plotdata_norm   = plotdata_smooth/(max(plotdata_smooth));
        xdata = [0:1:length(plotdata_norm)-1].*timestep;
        if i==1
            plot(xdata,plotdata_norm+(1*j),'Color', [.9 .5 .1])
            %             sourceData_4a_Time_cell1_y(1:length(plotdata_norm),j) = plotdata_norm+(1*j);
            %             sourceData_4a_Time_cell1_x(1:length(plotdata_norm),j) = xdata;
        else
            plot(xdata,plotdata_norm+(1*j),'Color', [0 .5 0])
            %             sourceData_4a_Time_cell2_y(1:length(plotdata_norm),j) = plotdata_norm+(1*j);
            %             sourceData_4a_Time_cell2_x(1:length(plotdata_norm),j) = xdata;
        end
        sourceData_4a_Time(1:length(plotdata_norm),sCount) = xdata;
        sCount = sCount + 1;
        sourceData_4a_Time(1:length(plotdata_norm),sCount) = plotdata_norm+(1*j);
        sCount = sCount + 1;
        
        ylim([0 length(numtrial)+1]);
        xlim([0 7]);
    end
end

subplot(1,2,2)
hold on
sCount = 1;
for j=1:length(numtrial)
    hold on
    for i=1:length(roi1)
        plotdataY(:,j,i) = squeeze(ROI_trial(:,j,i));
        plotdata_smoothY(:,j,i) = mind_preprocess(plotdataY(:,j,i),5,0);
        plotdata_normY = plotdata_smoothY(:,j,i)/(max(plotdata_smoothY(:,j,i)));
        if i==1
            plot([0:100],plotdata_normY+(1*j),'Color', [.9 .5 .1])
            %             sourceData_4a_Pos_cell1_y(1:length(plotdata_normY),j) = plotdata_normY+(1*j);
            %             sourceData_4a_Pos_cell1_x(1:length(plotdata_normY),j) = [0:100].*3;
        else
            plot([0:100],plotdata_normY+(1*j),'Color', [0 .5 0])
            %             sourceData_4a_Pos_cell2_y(1:length(plotdata_normY),j) = plotdata_normY+(1*j);
            %             sourceData_4a_Pos_cell2_x(1:length(plotdata_normY),j) = [0:100].*3;
        end
        sourceData_4a_Position(1:length(plotdata_normY),sCount) = [0:100].*3;
        sCount = sCount + 1;
        sourceData_4a_Position(1:length(plotdata_normY),sCount) = plotdata_normY+(1*j);
        sCount = sCount + 1;
        
        ylim([0 length(numtrial)+1]);
        xticks([0 33.33333 66.66667 100]);
        xticklabels({'0', '100' , '200', '300'})
    end
end
suptitle(num2str(roi1));
