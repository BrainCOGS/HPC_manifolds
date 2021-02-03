function [mouseID, date, isModeling] = parseModelingName(fname)
% isModeling is 1 if the filename is a modeling file

C = strsplit(fname,'\');

keepC = C{end};

D = strsplit(keepC,'_');

isModeling = ~isempty(strfind(fname,'modeling'));

if isModeling==1
    mouseID = D{1};
    mouseID = mouseID(2:end);
    mouseID = str2num(mouseID);
    
    date = D{2};
    date = str2num(date);
else
    mouseID = D{4};
    mouseID = mouseID(2:end);
    mouseID = str2num(mouseID);
    
    date = D{6};
    date = strsplit(date,'.');
    date = date{1};
    date = str2num(date);
end