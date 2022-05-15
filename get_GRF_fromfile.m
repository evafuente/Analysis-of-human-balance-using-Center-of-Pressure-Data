function GRF = get_GRF_fromfile(filename, n)
% GRF = get_GRF_fromfile(filename, n)
% nuskaito duomenis iðsaugotus .txt failiukuose
% inputs: 
%       filename: failo, kurá reikia nuskaityti, pavadinimas
%       n: duomenø ilgis
% Outputs:
%       GRF: A 2D array of doubles containing the row column data in filename

% %% Setup the Import Options and import the data
% opts = delimitedTextImportOptions("NumVariables", 19);
% 
% % Specify range and delimiter
% opts.DataLines = [1, Inf];
% opts.Delimiter = "\t";
% 
% % Specify column names and types
% opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19"];
% opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
% 
% % Specify file level properties
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% 
% % Import the data
% GRF_data = readtable(filename, opts);
% 
% %% Convert to output type
% GRF = table2array(GRF_data(:,18:19));
% 
% %% Clear temporary variables
% clear opts

fid=fopen(filename); % atidarom .txt bylà
next_line = fgetl(fid); % nuskaitom pirmà eilutæ
GRF = zeros(n, 2); % sukuriame masyvà atmintyje
i=1; 

while next_line ~= -1
    row_data = sscanf(next_line, '%f');
    GRF(i,:) = row_data(8:9);
    next_line = fgetl(fid);
    i=i+1;
    if i>n
        break
    end
end
fclose(fid);