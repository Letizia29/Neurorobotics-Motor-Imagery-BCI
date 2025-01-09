function subj = loadSubjData(folder_name)

% This function loads all the data of a subject
file_list = dir(strcat(folder_name, '/*.gdf'));

%file_list = dir('*.gdf'); % Gets the list of all subjects folders

subj.offline = struct('filename', {});
subj.online = struct('filename', {});

% Sort files into offline and online
for i = 1:length(file_list)
    if contains(file_list(i).name, 'offline', 'IgnoreCase', true)
        subj.offline(end+1).filename = file_list(i).name;
    elseif contains(file_list(i).name, 'online', 'IgnoreCase', true)
        subj.online(end+1).filename = file_list(i).name;
    end
end

subj.offline_data.s = struct('data', {});
subj.offline_data.h = struct('data', {});

for i = 1:length(subj.offline)
        filename = subj.offline(i).filename;

        file_to_load = strcat(folder_name, '/', filename);

        % Load the file
        [s, h] = sload(file_to_load);

        subj.offline_data.s(end+1).data = s; % Append data to s
        subj.offline_data.h(end+1).data = h; % Append data to h
end

subj.online_data.s = struct('data', {});
subj.online_data.h = struct('data', {});

for i = 1:length(subj.online)
        filename = subj.online(i).filename;

        file_to_load = strcat(folder_name, '/', filename);

        % Load the file
        [s, h] = sload(file_to_load);

        subj.online_data.s(end+1).data = s; % Append data to s
        subj.online_data.h(end+1).data = h; % Append data to h
end

end