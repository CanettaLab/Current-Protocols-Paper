function [Computed_data] = Analysis_Functions(Data)

% Last modification: 4/3/2023
% Last modifier: Arturo Torres-Herraez
% General description: Analysis_Functions is a list of processes used to compute raw data into final outputs.

% Inputs:
    % Data: Cells with matrices containing timestamps in the first
    %       column and code representing behavioral events from the operant
    %       box system for each animal

%  Outputs:
    % Computed_data: structure containing the processed data, e.g., total
    % number of lever presses/nosepokes, total number of reinforcers,
    % latencies to first lever press/nosepoke, etc

%% Define variables
SessionDur = ones(length(Data),1)*nan;
Reinforcers = ones(length(Data),1)*nan;
missed_dips = ones(length(Data),1)*nan;
responses = ones(length(Data),1)*nan;
n_responses_per_trial = cell(length(Data),1);
responserate = ones(length(Data),1)*nan;
responserate_per_trial = cell(length(Data),1);
meanLatency2reward = ones(length(Data),1)*nan;
meanLatency2firstresponse = ones(length(Data),1)*nan;
AllLatency2reward = cell(length(Data),1);
AllLatency2firstresponse =  cell(length(Data),1);

%% Compute data
% Loop through the data of all individuals included in Data
for i = 1:length(Data)
    %% Defining time of structural and behavioral events
    %--% Time Lever Extension %--%
    if length(Data{i}(Data{i}(:,2) == 28,1)) == length(Data{i}(Data{i}(:,2) == 30,1))
        t_lever_right = [Data{i}(Data{i}(:,2) == 28,1) Data{i}(Data{i}(:,2) == 30,1)];
    elseif length(Data{i}(Data{i}(:,2) == 28,1)) > length(Data{i}(Data{i}(:,2) == 30,1))
        if length(Data{i}(Data{i}(:,2) == 28,1)) == 2 
            idx = find(Data{i}(:,2) == 28);
            Data{i}(idx(1),:) = [];
            t_lever_right = [Data{i}(Data{i}(:,2) == 28,1) Data{i}(Data{i}(:,2) == 30,1)];
        else
            t_lever_right = ones(length(Data{i}(Data{i}(:,2) == 28,1)),2)*nan;
            t_lever_right(:,1) = Data{i}(Data{i}(:,2) == 28,1);
            t_lever_right(1:length(Data{i}(Data{i}(:,2) == 30,1)),2) = Data{i}(Data{i}(:,2) == 30,1);
        end
    end
    if ~isempty(t_lever_right)
        if isnan(t_lever_right(end,2))
            t_lever_right(end,2) = Data{i}(Data{i}(:,1) == max(Data{i}(:,1)),1);
        end
        if t_lever_right(end,2) == 0
            t_lever_right(end,2) = Data{i}(end,1);
        end
    end
    
    %--% Time of the response (either lever presses or nose pokes) %--%
    if length(Data{i}(Data{i}(:,2) == 1016,1)) == length(Data{i}(Data{i}(:,2) == 1018,1))
        t_response = [Data{i}(Data{i}(:,2) == 1016,1) Data{i}(Data{i}(:,2) == 1018,1)];
    elseif length(Data{i}(Data{i}(:,2) == 1016,1)) > length(Data{i}(Data{i}(:,2) == 1018,1))
        t_response = ones(length(Data{i}(Data{i}(:,2) == 1016,1)),2)*nan;
        t_response(:,1) = Data{i}(Data{i}(:,2) == 1016,1);
        t_response(1:length(Data{i}(Data{i}(:,2) == 1018,1)),2) = Data{i}(Data{i}(:,2) == 1018,1);
    end
    if ~isempty(t_response)
        if isnan(t_response(end,2))
            t_response(end,:) = [];
        end
    end
    
    %--% Time Dipper Up %--%
    if length(Data{i}(Data{i}(:,2) == 25,1)) == length(Data{i}(Data{i}(:,2) == 26,1))
        t_Dipper = [Data{i}(Data{i}(:,2) == 25,1) Data{i}(Data{i}(:,2) == 26,1)];
    elseif length(Data{i}(Data{i}(:,2) == 25,1)) > length(Data{i}(Data{i}(:,2) == 26,1))
        t_Dipper = ones(length(Data{i}(Data{i}(:,2) == 25,1)),2)*nan;
        t_Dipper(:,1) = Data{i}(Data{i}(:,2) == 25,1);
        t_Dipper(1:length(Data{i}(Data{i}(:,2) == 26,1)),2) = Data{i}(Data{i}(:,2) == 26,1);
    end
    if ~isempty(t_Dipper)
        if isnan(t_Dipper(end,2))
            t_Dipper(end,:) = [];
        end
    end
    
    %--% Time Head Entry %--%
    t_head_entry = Data{i}(Data{i}(:,2) == 1011,1);
    
    %% Computing behavioral data
    SessionDur(i) = (Data{1,i}(end,1) - Data{1,i}(1,1));
    
    %--% Create a trial list %--%
    if ~isempty(t_Dipper)
        if Data{1,i}(end,1) > t_Dipper(end,2)
            A = [Data{1,i}(1,1);t_Dipper(:,2)];
            B = [t_Dipper(:,2);SessionDur(i)];
            t_trials = [A B];
        else
            t_trials = t_Dipper;
            t_trials(1,1) = Data{1,i}(1,1);
        end
    else
        t_trials = [Data{1,i}(1,1) Data{1,i}(end,1)];
    end
    
    %--% number of responses before reinforcement %--% 
    n_responses_per_trial{i} = ones(size(t_trials,1),1)*nan;
    if ~isempty(t_Dipper)
        for dp = 1:size(t_trials,1)
            if dp == 1
                tmp = t_response(t_response(:,1) <= t_Dipper(dp,1));
            elseif dp == size(t_trials,1)
                tmp = t_response(t_response(:,1) >= t_Dipper(dp-1,2) & ...
                    t_response(:,1) <= t_trials(dp,2));
            else
                tmp = t_response(t_response(:,1) >= t_Dipper(dp-1,2) & ...
                    t_response(:,1) <= t_Dipper(dp,1));
            end
            n_responses_per_trial{i}(dp) = length(tmp);
        end
    else
        tmp = t_response(t_response(:,1) >= t_trials(1,1) & ...
                    t_response(:,1) <= t_trials(1,2));
        n_responses_per_trial{i}(1) = length(tmp);
    end

    %--% Latency to first response per trial %--%
    latency2firstresponse = ones(size(t_trials,1),1)*nan;
    if ~isempty(t_Dipper)
        for dp = 1:size(t_trials,1)
            if dp == 1
                tmp = t_response(t_response(:,1) <= t_Dipper(dp,1));
                latency2firstresponse(dp) = tmp(1,1) - Data{i}(1,1);
            elseif dp == size(t_trials,1)
                tmp = t_response(t_response(:,1) >= t_Dipper(dp-1,2) & ...
                    t_response(:,1) <= t_trials(dp,2));
                if ~isempty(tmp)
                    latency2firstresponse(dp) = tmp(1,1) - t_Dipper(dp-1,2);
                else
                    latency2firstresponse(dp) = nan;
                end
            else
                tmp = t_response(t_response(:,1) >= t_Dipper(dp-1,2) & ...
                    t_response(:,1) <= t_Dipper(dp,1));
                latency2firstresponse(dp) = tmp(1,1) - t_Dipper(dp-1,2);
            end
        end
    else
        tmp = t_response(t_response(:,1) >= t_trials(1,1) & ...
                    t_response(:,1) <= t_trials(1,2));
        if ~isempty(tmp)
            latency2firstresponse(dp) = tmp(1,1) - t_trials(1,1);
        end
    end

    %--% response rate %--% 
    DipperUpDur = t_Dipper(:,2)-t_Dipper(:,1);
    responserate(i) = length(t_response)./(SessionDur(i) - sum(DipperUpDur));
    responserate_per_trial{i} = ones(size(t_trials,1),1)*nan;
    if ~isempty(t_Dipper)
        for dp = 1:size(t_trials,1)
            if dp == 1
                tmp = t_response(t_response(:,1) <= t_Dipper(dp,1));
                responserate_per_trial{i}(dp) = length(tmp)./(t_Dipper(dp,1) - Data{1,i}(1,1));
            elseif dp == size(t_trials,1)
                tmp = t_response(t_response(:,1) >= t_Dipper(dp-1,2) & ...
                    t_response(:,1) <= t_trials(dp,2));
                if ~isempty(tmp)
                    responserate_per_trial{i}(dp) = length(tmp)./(t_trials(dp,2) - t_Dipper(dp-1,2));
                else
                    responserate_per_trial{i}(dp) = 0;
                end
            else
                tmp = t_response(t_response(:,1) >= t_Dipper(dp-1,2) & ...
                    t_response(:,1) <= t_Dipper(dp,1));
                responserate_per_trial{i}(dp) = length(tmp)./(t_Dipper(dp,1) - t_Dipper(dp-1,2));
            end
        end
    else
        tmp = t_response(t_response(:,1) >= t_trials(1,1) & ...
            t_response(:,1) <= t_trials(1,2));
        if ~isempty(tmp)
            responserate_per_trial{i}(dp) = length(tmp)./(t_trials(1,2) - t_trials(1,1));
        else
            responserate_per_trial{i}(1) = 0;
        end
    end
    
    %--% reinforcements %--%
    Reinforcers(i) = size(t_Dipper,1);
    
    %--% Rewarded trials and latency to reward %--%
    Rewarded = 0;
    Latency2reward = ones(size(t_Dipper,1),1)*nan;
    for dp = 1:size(t_Dipper,1)
        tmp = t_head_entry(t_head_entry >= t_Dipper(dp,1) & t_head_entry <= t_Dipper(dp,2));
        if ~isempty(tmp)
            Rewarded = Rewarded + 1;
            Latency2reward(dp) = tmp(1) - t_Dipper(dp,1);
        end
    end
    
    %--% Non reinforced trials %--%
    missed_dips(i) = Reinforcers(i) - Rewarded;
    
    AllLatency2firstresponse{i} = latency2firstresponse;
    meanLatency2firstresponse(i) = mean(latency2firstresponse,'omitnan');
    AllLatency2reward{i} = Latency2reward;
    meanLatency2reward(i) = mean(Latency2reward,'omitnan');
    responses(i) = length(t_response);

end

%% ======== Generating the Data Structure ==========================%
Computed_data.SessionDur = SessionDur; % Session duration (in seconds)
Computed_data.Reinforcers = Reinforcers; % Total number of dipper presentations
Computed_data.MissedDips = missed_dips; % Number of trials which mouse did not make a head entry during the time the dipper was up
Computed_data.Responses = responses; % Total number of responses (lever presses or nose pokes)
Computed_data.number_responses_per_trial = n_responses_per_trial; % Number of responses per trial (lever presses or nose pokes)
Computed_data.MeanResponseRate = responserate; % Total number of responses divided by the duration of the session (excluding time of dipper up) 
Computed_data.ResponseRatePerTrial = responserate_per_trial; % Number of responses per trial divided by the duration of the trial (excluding time of dipper up) 
Computed_data.Latency2rewardPerTrial = AllLatency2reward; % Time from dipper up to the first head entry into the dipper port for each trial
Computed_data.meanLatency2reward = meanLatency2reward; % Average of the latency to reward per trial
Computed_data.Latency2firstresponsePerTrial = AllLatency2firstresponse; % Time from dipper down to the first response in the next trial
Computed_data.meanLatency2firstresponse = meanLatency2firstresponse; % Average of the latency to first response per trial

end