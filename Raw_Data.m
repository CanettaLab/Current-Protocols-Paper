function rawlist = Raw_Data(pth)

% Last modification: 4/3/2023
% Last modifier: Arturo Torres-Herraez
% General description: Raw_Data extracts the data from MedPC boxes and
% translates it into a matlab meaningful structure.

% Inputs:
    % pth: string indicating the path to the location of the MedPC raw data.
    % E.g., 'Z:\MedPC_Data\Experiment 1\Progressive Ratio\'

%  Outputs:
    % Rawlist: structure containing: 
        % name: Cells with the animal IDs
        % data: Cells with matrices containing timestamps in the first
        %       column and code representing behavioral events from the operant
        %       box system for each of the specified animals (matched with name)
        % program: Cells with the name of the program run in the operant
        %          box


rawlist = struct('name',{{}},'data',{{}},'program',{{}});

files = dir(pth);
elim = [];
foldr = [];
for i = 1:length(files)
	if files(i).name(1) == '.' || files(i).name(1) == '~'
		elim = [elim i];
	elseif isdir(strcat(pth,files(i).name)) == 1
		elim = [elim i];
		foldr = [foldr i];
	elseif files(i).bytes == 0
		elim = [elim i];
	elseif strcmp(files(i).name,'outfile.BbB')
		elim = [elim i];
    elseif contains(files(i).name,'.mat')
        elim = [elim i];
	end
end
foldr = files(foldr);
files(elim) = [];

ndex = 0;
for i = 1:length(files)
	rawlist.name(ndex+1) = {files(i).name};
    rawlist.program(ndex+1) = {local_MedPCprog(strcat(pth,files(i).name))};
    
    rawlist.data(ndex+1) = {local_randycode(strcat(pth,files(i).name))};
	
	ndex = ndex+1;
end

for i = 1:length(foldr)
    out = Raw_Data(strcat(pth,foldr(i).name,'/'));
    odex = length(out.name);
    rawlist.name(ndex+1:ndex+odex) = out.name;
	rawlist.data(ndex+1:ndex+odex) = out.data;
	rawlist.program(ndex+1:ndex+odex) = out.program;
	ndex = ndex+odex;
end


end

function prog = local_MedPCprog(pth)
%LOCAL_MEDPCPROG Designed to handle generic randycode data

	fid = fopen(pth);
	q = textscan(fid, '%s');
    fclose(fid);
    ddex = find(ismember(q{1,1}, 'MSN:')==1) + 1;
    prog = q{1,1}(ddex);
end

function data = local_randycode(pth)
%LOCAL_RANDYCODE Designed to handle generic randycode data

	fid = fopen(pth);
    q = textscan(fid, '%s');
    fclose(fid);
    ddex = find(ismember(q{1,1}, 'W:')==1) + 1;
    d = q{1,1}(ddex:length(q{1,1}));
    
    for i = 1:length(d)
        if d{i}(length(d{i})) == ':'
            ph(i) = (1==0);
        else
            ph(i) = (1==1);
        end
    end

    ph = find(ph);
    for i = 1:length(ph)
        data(i,1) = str2num(d{ph(i)});
        data(i,2) = mod(data(i,1),10000); % A way to get the last four digits
        data(i,1) = (data(i,1)-data(i,2))/10000000; % Subtraction is done to remove the event code
        % The denominator has 10000 * 10 * time Used to be 10000000
    end
end
