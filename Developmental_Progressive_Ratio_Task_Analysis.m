% Last modification: 4/3/2023
% Last modifier: Arturo Torres-Herraez
% General description: Main program that runs through the raw data,
% extracts it and performs the data analysis

% -- % General starting settings % -- %
close all; clear variables; clc; % Clear the workspace

% -- % Define path to operant box data % -- %
corepth = 'Z:\Emily C\PN FLX -  Reward Processing\Development Fall\Adult\MatLab Data\Group 2\Day 10 PR 3\'; % Fill the quotations with the path to your data. E.g., 'Z:\MedPC_Data\Experiment 1\Progressive Ratio\' 

% -- % Look for different data folders in your path % -- %
% Within your data folder you may have more than one day/experiment folder containing all your mice data 
files2analyze = dir(corepth);
for i = length(files2analyze):-1:1
    if strcmp(files2analyze(i).name(1),'.') == 1 || files2analyze(i).isdir == 0 ...
            ||  contains(files2analyze(i).name,'Template') == 1
        files2analyze(i) = [];
    end
end

% -- % Loop through the folders found in your path % -- %
for i = 1:length(files2analyze)
    localpth= [files2analyze(i).name,'\'];
    pth = strcat(corepth,localpth);

    % -- % Function that extracts MedPc data to matlab format % -- %
    rawlist = Raw_Data(pth);
    Data = rawlist.data;

    % -- % Create a variable containing the Animal IDs for all individuals % -- %
    tmp = dir(pth);
    for ii = length(tmp):-1:1
        if strcmp(tmp(ii).name(1),'.') == 1 || tmp(ii).isdir == 0
            tmp(ii) = [];
        end
    end
    
    Animal_Id = cell(length(tmp),1);
    for ii = 1:length(tmp)
        Animal_Id{ii} = tmp(ii).name;
    end

    % -- % Define the folder where you want to save your analysis % -- %
    Path2savingfolder = strcat(corepth,localpth);
    
    % -- % Function that analyzes the extracted data % -- %
    [OPERANT_DATA] = Analysis_Functions(Data);
    OPERANT_DATA.MiceId = Animal_Id;
    
    % -- % Save the analized data  % -- %
    save([pth,'Developmental_Progressive_Ratio_Task_Analysis'],'OPERANT_DATA')
end