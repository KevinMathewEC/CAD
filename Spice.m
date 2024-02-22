filepath = 'C:\Users\Kevin\OneDrive - Indian Institute of Science\IISc sem 2\CAD\HW1\Spice_q1.txt';

% Open the file
fileID = fopen(filepath, 'r');

% Initialize cell arrays to store the extracted fields
component = {};%component 'R' or 'C' ...
row1 = [];
col1 = [];
value = {};% value of component
node1=0;
node2=0;
node3=0;
node4=0;
tokens='';
% Initialize the hash table as a cell array
global hashTable;
hashTable=zeros(8,8);

global RHS_table;
RHS_table=zeros(8,1);

RHS_stampValue=zeros(8,1);
positions={};%positions in the array where the stamp needs to be added
stampValue=[];

max_node=0;

Ib_node =8;% node for modified nodal analysis

Ib_node_used=0;%keep track if modified nodal analysis used
h=100e-6;% time step for simulation

soln_matrix=zeros(8,1);
global time1;
global time2;% declare as global to keep track of time period of simulation

global pulse_signal;
global sine_wave;
global pwl_signal;

global CapacitorStruct ;
CapacitorStruct = struct('node1', {}, 'node2', {}, 'value', {});
% Read lines from the file
while ~feof(fileID)
    addToRHS=0;
    addToMatrix=0;
    valid_text=0;
    % Read a line from the file
    line = fgetl(fileID);

    % Parse the line using regular expression
    tokens = regexp(line, '\s+', 'split');

    % Check if the line contains '.TRAN'
    if contains(line, '.TRAN', 'IgnoreCase', true)

        % Extract the values using regular expressions
        % tokens = regexp(line, '\.TRAN\s+(\w+)\s+(\w+)', 'tokens');
        if ~isempty(tokens)
            % Extract the time values
            time1 = tokens{2};
            time2 = tokens{3};

            time1=convertUnit(time1(1:end-1));
            time2=convertUnit(time2(1:end-1));
            % Run Transient response with the extracted values
            Transient_run(max_node,Ib_node_used,h,Ib_node,"Trans");  %
        end
    end

    %    if regexp(line, '^\s*\.alter')
    if contains(line, 'PULSE', 'IgnoreCase', true)
        % Define regular expression pattern to match PULSE parameters
        pattern = 'PULSE\(\s*(\S+)V\s*(\S+)V\s*(\S+)s\s*(\S+)s\s*(\S+)s\s*(\S+)s\s*(\S+)s\)';

        % Match pattern in the line
        match = regexp(line, pattern, 'tokens', 'once');

        % Check if a match was found
        if isempty(match)
            error('Failed to parse PULSE signal parameters.');
        end

        % Convert matched values to numerical format
        pulse_params = cellfun(@convertUnit, match);

        % Assign the extracted values to variables
        V1 = pulse_params(1);
        V2 = pulse_params(2);
        Tdelay = pulse_params(3);
        Trise = pulse_params(4);
        Tfall = pulse_params(5);
        Ton = pulse_params(6);
        Period = pulse_params(7);

        % Time vector
        t = 0:h:time2;

        % Initialize array to store the pulse signal values
        pulse_signal = zeros(size(t));
        disp(Ton);
        disp(Period);

        % Generate the pulse signal
        for i = 1:numel(t)
            % Check if current time is within the pulse window
            if mod(t(i), Period) <= Ton
                pulse_signal(i) = V2;
            else
                pulse_signal(i) = V1;
            end
        end

        %  Transient_run(max_node,Ib_node_used,h,Ib_node,"Pulse");
        % % Plot the pulse signal
        % plot(t, pulse_signal);
        % xlabel('Time');
        % ylabel('Amplitude');
        % title('Pulse Signal');
        %
        % % Example: Get the value of the pulse signal at time t = 1.3
        % index_t = find(t == 1.3);
        % value_at_t = pulse_signal(index_t);
        % disp(['Value of pulse signal at time t = 1.3: ', num2str(value_at_t)]);
    end

    %%If input voltage is sine wave
    if contains(line, 'SIN', 'IgnoreCase', true)

        % Define regular expression pattern to match the parameters
        pattern = '(sin|SIN)\(\s*(\S+)v\s*(\S+)v\s*(\S+)Hz\s*(\S+)s\s*(\S+)\s*(\S+)\s*\)';

        % Extract parameters using regular expressions
        matches = regexp(line, pattern, 'tokens', 'once', 'ignorecase');

        if ~isempty(matches)
            % Convert extracted tokens to numeric values
            Voffset = convertUnit(matches{2}); % Offset voltage
            Vamp = convertUnit(matches{3}); % Amplitude
            freq = convertUnit(matches{4}); % Frequency
            td = convertUnit(matches{5}); % Time delay
            theta = convertUnit(matches{6}); % Phase angle

            % Generate time vector
            t = 0:h:time2; % Adjust time range as needed

            % Generate sine wave
            sine_wave = Voffset + Vamp * sin(2*pi*freq*t + deg2rad(theta));

            %   Transient_run(max_node,Ib_node_used,h,Ib_node,"Sine");


            % Plot the sine wave
            % plot(t, sine_wave);
            % xlabel('Time (s)');
            % ylabel('Voltage (V)');
            % title('Generated Sine Wave');
            % grid on;
        end
    end

    %%If input voltage is sine wave
    if contains(line, 'PWL', 'IgnoreCase', true)
        pattern = '(?i)PWL\(([^)]+)\)';

        % Extract parameters using regular expressions
        matches = regexp(line, pattern, 'tokens', 'once');

        if ~isempty(matches)
            % Extract time and voltage pairs
            pairs = strsplit(matches{1}, {' ', ','});

            % Initialize arrays to store time and voltage values
            numPairs = numel(pairs) / 2;
            time = zeros(1, ceil(numPairs));
            voltage = zeros(1, ceil(numPairs));

            % Extract time and voltage values from pairs
            for i = 1:numPairs
                time(i)=convertUnit(pairs{(i-1)*2+1}(1:end-1));
                voltage(i) =convertUnit(pairs{i*2}(1:end-1));

            end
            % Remove duplicate time points
            [time, uniqueIdx] = unique(time);
            voltage = voltage(uniqueIdx);

            % Interpolate data points to create PWL signal
            interpTime = 0:h:time(end); % Time points at h intervals
            pwl_signal = interp1(time, voltage, interpTime, 'linear'); % Linear interpolation

            % Plot the piecewise linear graph
            plot(interpTime, pwl_signal);
            xlabel('Time (ms)');
            ylabel('Voltage (V)');
            title('Piecewise Linear Graph');
            grid on;
            time2=time(i)
            Transient_run(max_node,Ib_node_used,h,Ib_node,"Pwl")
            % Plot the piecewise linear graph
        else
            disp('No PWL pattern found in the input line.');
        end
    end


    % Check if the line matches the required format
    if numel(tokens) == 4

        % value = [value, tokens{4}];
        component =  tokens{1};
        node1 =  str2double(tokens{2});
        node2 =  str2double(tokens{3});
        value =  tokens{4};
        valid_text=1;

    elseif numel(tokens) == 6
        component =  tokens{1};
        node1 =  str2double(tokens{2});
        node2 =  str2double(tokens{3});
        node3 =  str2double(tokens{4});
        node4 =  str2double(tokens{5});
        value =  tokens{6};
        valid_text=1;
    end
    max_node=max(max_node,max(node1,node2));

    if(valid_text)

        firstChar = lower(component(1));
        % Check if the first character is 'r'
        switch firstChar
            case 'r'
                value=convertUnit(value);
                value=double(1/value);
                if ((node1==0)|| (node2==0))% if any of the nodes are ground
                    if(node1)
                        positions={[node1,node1]};
                        stampValue={value};
                    else
                        positions={[node2,node2]};
                        stampValue={-value};
                    end
                else
                    positions = {[node1, node1], [node1, node2], [node2, node1], [node2, node2]};
                    stampValue={value,-value,-value,value};
                end

                addToMatrix=1;
                addToRHS=0;
            case 'i'
                value=convertUnit(value);
                addToMatrix=0;
                addToRHS=1;
                if ((node1==0)|| (node2==0))% if any of the nodes are ground
                    if(node1)
                        positions={[node1,1]};
                        RHS_stampValue={-value};
                    else
                        positions={[node2,1]};
                        RHS_stampValue={value};
                    end
                else
                    positions = {node1, node2};
                    RHS_stampValue={-value,value};
                end

            case 'c'
                value=convertUnit(value);
                addToMatrix=1;
                addToRHS=1;
                if ((node1==0)|| (node2==0))% if any of the nodes are ground
                    if(node1)
                        positions={[node1,node1]};
                        RHS_stampValue(node1)=(value/h).*soln_matrix(node1);
                        stampValue={value/h};
                    else
                        positions={[node2,node2]};
                        RHS_stampValue(node2)=(-value/h).*soln_matrix(node2);
                        stampValue={-value/h};
                    end

                else
                    positions = {[node1, node1], [node1, node2], [node2, node1], [node2, node2]};
                    RHS_stampValue(node1)=(value/h)*soln_matrix(node2);
                    RHS_stampValue(node2)= (-value/h)*soln_matrix(node1);
                    stampValue={value/h,-value/h,-value/h,value/h};
                end

                % Store node1, node2, and value in the structure array
                CapacitorStruct(end+1).node1 = node1;
                CapacitorStruct(end).node2 = node2;
                CapacitorStruct(end).value = value;




            case 'v'
                Ib_node_used=1;
                value=value(1:end-1);
                value=convertUnit(value);
                addToMatrix=1;
                addToRHS=1;
                if ((node1==0)|| (node2==0))% if any of the nodes are ground
                    if(node1)
                        positions={[node1,Ib_node],[Ib_node,node1]};
                    else
                        positions={[node2,Ib_node],[Ib_node,node2]};
                    end

                    stampValue={1,1};
                else
                    positions = {[node1,node1],[node1,Ib_node],[node2,node2],[node2,Ib_node],[Ib_node,node1],[Ib_node,node2]};
                    stampValue={0,1,0,-1,+1,-1};
                end
                RHS_stampValue(Ib_node)=value;

            case 'g'
                value=convertUnit(value);
                addToMatrix=1;
                addToRHS=0;
                if ((node3==0)|| (node4==0))% if any of the nodes are ground
                    if(node3)
                        positions={[node3,node3]};
                        stampValue={value};
                    else
                        positions={[node4,node4]};
                        stampValue={-value};
                    end
                else
                    positions = {[node3, node3], [node3, node4], [node4, node3], [node4, node4]};
                    stampValue={value,-value,-value,value};
                end

        end


        % Assign the value to the corresponding positions in the hash table


        for i = 1:size(positions,2)

            row = positions{i}(1);
            col = positions{i}(2);
            if(addToMatrix)
                hashTable(row, col) = hashTable(row, col)+stampValue{i};
            end
        end
        if(addToRHS)
            current_row=0;

            for i = 1:size(positions,2)
                row = positions{i}(1);
                if (current_row ~= row)
                    RHS_table(row)= RHS_table(row)+RHS_stampValue(row);
                    current_row=row; %avoid duplicate entry on same row, allow reuse of same posittion array

                end
            end
        end
    end

    %%Invert Matrix


end

% Close the file
fclose(fileID);
%disp(RHS_table)

% Print the hash table
% disp('Hash Table:');
for row = 1:size(hashTable,1)
    for col = 1:size(hashTable, 2)
        if ~isempty(hashTable(row, col))
            %         fprintf('(%d, %d): %d\t', row, col, hashTable(row, col));
        end
    end
    % fprintf('\n');
end





function convertedValue = convertUnit(inputStr)
% Define unit multipliers
multipliers = containers.Map({'k', 'M', 'u', 'n','m'}, {1e3, 1e6, 1e-6,1e-9,1e-3});

% Check if the unit suffix is empty and extract it differently
numStr = inputStr(1:end-1);  % Extract numerical value
unitStr = inputStr(end);     % Extract unit suffix

if(unitStr~='M')
    unitStr=lower(unitStr);
end
% Get the multiplier corresponding to the unit suffix
if multipliers.isKey(unitStr)
    % Convert numerical value to a number
    numericalValue = str2double(numStr);

    % Get the multiplier corresponding to the unit suffix
    multiplier = multipliers(unitStr);

    % Calculate the converted value
    convertedValue = numericalValue * multiplier;
elseif isnumeric(str2double(unitStr))
    % Convert numerical value to a number
    convertedValue = str2double(inputStr);

end
end

% Define your function to handle the time values
function Transient_run(max_node,Ib_node_used,h,Ib_node,str)
global time1;
global time2;
global hashTable;
global RHS_table;
global CapacitorStruct;
global pulse_signal;
global sine_wave;
global pwl_signal;

soln_matrix=[];
% Extract the submatrix of size max_node*max_node
submatrix = hashTable(1:max_node, 1:max_node);
stampTable = submatrix;
if Ib_node_used
    % Append the last column and last row to the submatrix
    last_column = hashTable(1:max_node, end);
    last_row = hashTable(end, 1:max_node);

    % Construct the new square matrix
    stampTable = [submatrix, last_column; last_row, hashTable(end, end)];
end

% Modify RHS table to include only rows 1 through max_node

if Ib_node_used
    last_row = RHS_table(end, 1);
    RHS_table =[RHS_table(1:max_node);last_row];

else
    RHS_table = RHS_table(1:max_node);
end


% Find inverse of stamp table
invStampTable = inv(stampTable);
%  disp('Inverse of stamp table:');
%  disp(invStampTable);
execution_cycle=1;%keep track of execution cycle
time_array =[];

for t = 0:h:time2

    soln_column = zeros(size(invStampTable,2),1);
    row_sum=0;
    for i = 1:size(invStampTable,1)
        for j = 1:size(invStampTable,2)
            row_sum=row_sum+invStampTable(i,j).*RHS_table(j,1);
        end
        soln_column(i) = row_sum;
        row_sum=0;
    end

    soln_matrix=[soln_matrix,soln_column];

    %time_array = [time_array, t];
    % Convert the input string to lowercase to make the comparison case-insensitive
    strLower = lower(str);

    % Check if the string equals "pulse"
    isPulse = strcmp(strLower, 'pulse');
    if(isPulse)
        RHS_table(end)= pulse_signal(execution_cycle);

    end

    isSine = strcmp(strLower, 'sine');
    if(isSine)

        RHS_table(end)= sine_wave(execution_cycle);

    end

    isPwl=strcmp(strLower,'pwl');
    if(isPwl)
        RHS_table(end)=pwl_signal(execution_cycle);
    end



    for i=1:size(CapacitorStruct,2)
        if(CapacitorStruct(i).node1 && CapacitorStruct(i).node2)
            RHS_table(CapacitorStruct(i).node1) = (CapacitorStruct(i).value/h)*((soln_matrix(CapacitorStruct(i).node1,execution_cycle))-(soln_matrix(CapacitorStruct(i).node2,execution_cycle)));
            RHS_table(CapacitorStruct(i).node2) = (CapacitorStruct(i).value/h)*((soln_matrix(CapacitorStruct(i).node2,execution_cycle))-(soln_matrix(CapacitorStruct(i).node1,execution_cycle)));
        elseif(CapacitorStruct(i).node1)
            RHS_table(CapacitorStruct(i).node1) = (CapacitorStruct(i).value/h)*(soln_matrix(CapacitorStruct(i).node1,execution_cycle));
        elseif (CapacitorStruct(i).node2)
            RHS_table(CapacitorStruct(i).node2) = (CapacitorStruct(i).value/h)*(soln_matrix(CapacitorStruct(i).node2,execution_cycle));
        end
    end

    execution_cycle=execution_cycle+1;


end
figure;
for i = 1:max_node
    subplot(max_node, 1, i); % Create subplot for each row
    plot(soln_matrix(i, :)); % Plot the ith row of soln_matrix
    title(['Voltage node ', num2str(i)]); % Set title for the subplot
    xlabel('time'); % Set x-axis label
    ylabel('Volts(V)'); % Set y-axis label
end
end
