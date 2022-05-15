%% Clear workspace
clearvars
close all
clc

%%
addpath SubjectsMeasurements\
[~, ~, SubjData] = xlsread('BDSinfo.xlsx','Sheet1','A1:BL1931');
SubjData(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),SubjData)) = {''}; 
filename = SubjData{2,1};
TempBDS_file = importBalanceTXTfile(filename);
ds = dataset('File',fullfile('BDSinfo.txt'));

%% read existing .txt files of measuremnts and save file names
BalanceFilesID = dir('SubjectsMeasurements\*.txt');
fid = fopen('Bylu_sarasas.txt', 'w');
for i=1:length(BalanceFilesID)
    FailoVardas = BalanceFilesID(i,1).name;
    fprintf(fid, '%s \r\n', FailoVardas);
end
fclose(fid);
clearvars

%% read excel 
BDSinfo = readtable('BDSinfo.xlsx');

%% read Balance txt files
addpath SubjectsMeasurements\
matavimo_bylos = importdata('Bylu_sarasas.txt'); %file with all name data files
NumFiles = length(matavimo_bylos);

PATIENTS=160;
DATA_PACIENT=12;

for i=1:(DATA_PACIENT*PATIENTS) 
    s1 = 'SubjectsMeasurements\';
    s2 = matavimo_bylos{i}; %BDS with number.txt
    s = strcat(s1,s2);
    filename = s;
    TempBDS_file = importBalanceTXTfile(filename); %import data
%     COPxcm=importBalanceTXTfile(filename, [1,3]) %me saca las filas
    ID = strsplit(matavimo_bylos{i, 1},'.'); %BDS with number
    BalanceDB(i).ID =  ID{1,1}; %BDS with number
    BalanceDB(i).BDSfile = TempBDS_file; %BDS with number.txt

    if strcmp(BDSinfo.Surface(i), 'Firm')
           COPxcm(i,:)=TempBDS_file.COPxcm;
           COPycm(i,:)=TempBDS_file.COPycm;
           
        if strcmp(BDSinfo.Vision(i), 'Open')
            COPxcm_open(i,:)=TempBDS_file.COPxcm; 
            COPycm_open(i,:)=TempBDS_file.COPycm;
        else
            COPxcm_closed(i,:)=TempBDS_file.COPxcm; 
            COPycm_closed(i,:)=TempBDS_file.COPycm;
        end
    end
end

%% Remove zeros from matrix. Those zeros represent the rows that are not firm surface.

COPxcm = remove_zeros(COPxcm); %COPx patients by rows
COPycm = remove_zeros(COPycm);

% Matrices AUX
COPxcm_open = remove_zeros(COPxcm_open);
COPxcm_closed = remove_zeros(COPxcm_closed);
COPycm_open = remove_zeros(COPycm_open);
COPycm_closed = remove_zeros(COPycm_closed);

% remove first 10 sec
for i=1:1000
   COPxcm(:,i)=[];
   COPycm(:,i)=[];
   COPxcm_open(:,1)=[];
   COPxcm_closed(:,1)=[];
   COPycm_open(:,1)=[];
   COPycm_closed(:,1)=[];
end

%% Mean and std of COP for each patient
COPxcm_mean_pac = mean(COPxcm,2);
COPycm_mean_pac = mean(COPycm,2);

COPx_mean_open = mean(COPxcm_open, 'all');
COPy_mean_open = mean(COPycm_open, 'all');
COPx_mean_closed = mean(COPxcm_closed, 'all');
COPy_mean_closed = mean(COPycm_closed, 'all');

COPx_std_open = std2(COPxcm_open);
COPx_std_closed = std2(COPxcm_closed);
COPy_std_open = std2(COPycm_open);
COPy_std_closed = std2(COPycm_closed);

%% RANGE
[COPxcm_range_pac] = get_rango(COPxcm);
[COPycm_range_pac] = get_rango(COPycm);

COPxO_range_mean = mean(get_rango(COPxcm_open), 'all');
COPyO_range_mean = mean(get_rango(COPycm_open), 'all');
COPxC_range_mean = mean(get_rango(COPxcm_closed), 'all');
COPyC_range_mean = mean(get_rango(COPycm_closed), 'all');

COPxO_range_std = std2(get_rango(COPxcm_open));
COPyO_range_std = std2(get_rango(COPycm_open));
COPxC_range_std = std2(get_rango(COPxcm_closed));
COPyC_range_std = std2(get_rango(COPycm_closed));

%% RANGE RATIO = RANGE ML/RANGE AL
ratio_open_mean = COPyO_range_mean/COPxO_range_mean;
ratio_closed_mean = COPyC_range_mean/COPxC_range_mean;

%% TOTAL LENGHT PATH 
COPxcm_traj = get_length_traj(COPxcm);
COPycm_traj = get_length_traj(COPycm);

COPxO_traj_mean = mean(get_length_traj(COPxcm_open), 'all');
COPyO_traj_mean = mean(get_length_traj(COPycm_open), 'all');
COPxC_traj_mean = mean(get_length_traj(COPxcm_closed), 'all');
COPyC_traj_mean = mean(get_length_traj(COPycm_closed), 'all');

COPxO_traj_std = std2(get_length_traj(COPxcm_open));
COPyO_traj_std = std2(get_length_traj(COPycm_open));
COPxC_traj_std = std2(get_length_traj(COPxcm_closed));
COPyC_traj_std = std2(get_length_traj(COPycm_closed));

total_traj_open = sqrt(COPxO_traj_mean*COPxO_traj_mean+COPyO_traj_mean*COPyO_traj_mean);
total_traj_closed = sqrt(COPxC_traj_mean*COPxC_traj_mean+COPyC_traj_mean*COPyC_traj_mean);

%% VELOCITY AND ACEL STUDY
[COPxcm_vel, COPxcm_acel] = get_vel(COPxcm);
[COPycm_vel, COPycm_acel] = get_vel(COPycm);

[COPxO_vel, COPxO_acel] = get_vel(COPxcm_open);
[COPyO_vel, COPyO_acel] = get_vel(COPycm_open);
[COPxC_vel, COPxC_acel] = get_vel(COPxcm_closed);
[COPyC_vel, COPyC_acel] = get_vel(COPycm_closed);

COPxO_vel_mean = mean(COPxO_vel, 'all');
COPyO_vel_mean = mean(COPyO_vel, 'all');
COPxC_vel_mean= mean(COPxC_vel, 'all');
COPyC_vel_mean= mean(COPyC_vel, 'all');

COPxO_acel_mean = mean(COPxO_acel, 'all');
COPyO_acel_mean = mean(COPyO_acel, 'all');
COPxC_acel_mean = mean(COPxC_acel, 'all');
COPyC_acel_mean = mean(COPyC_acel, 'all');

COPxO_vel_std = std2(COPxO_vel);
COPyO_vel_std = std2(COPyO_vel);
COPxC_vel_std= std2(COPxC_vel);
COPyC_vel_std= std2(COPyC_vel);

COPxO_acel_std = std2(COPxO_acel);
COPyO_acel_std = std2(COPyO_acel);
COPxC_acel_std = std2(COPxC_acel);
COPyC_acel_std = std2(COPyC_acel);

total_vel_open = sqrt(COPxO_vel_mean*COPxO_vel_mean+COPyO_vel_mean*COPyO_vel_mean);
total_vel_closed = sqrt(COPxC_vel_mean*COPxC_vel_mean+COPyC_vel_mean*COPyC_vel_mean);

total_accel_open = sqrt(COPxO_acel_mean*COPxO_acel_mean+COPyO_acel_mean*COPyO_acel_mean);
total_accel_closed = sqrt(COPxC_acel_mean*COPxC_acel_mean+COPyC_acel_mean*COPyC_acel_mean);

%% SWAY AREA
COPcm2_SA = get_sway_area(COPxcm, COPycm);

openSA_mean = mean(get_sway_area(COPxcm_open, COPycm_open), 'all');
closedSA_mean = mean(get_sway_area(COPxcm_closed, COPycm_closed), 'all');

openSA_std = std2(get_sway_area(COPxcm_open, COPycm_open));
closedSA_std = std2(get_sway_area(COPxcm_closed, COPycm_closed));

%% CREATE NEW TABLE FOR MODEL
load('BDSinfo_LEARN2.mat', 'BDSinfo22')
BDSinfo_model = [head(BDSinfo22, size(COPxcm,1)) array2table(COPxcm_range_pac) array2table(COPycm_range_pac) array2table(COPxcm_mean_pac) array2table(COPycm_mean_pac) array2table(COPxcm_traj) array2table(COPycm_traj) array2table(COPxcm_vel) array2table(COPycm_vel) array2table(COPxcm_acel) array2table(COPycm_acel) array2table(COPcm2_SA)];

%% PLOT
subplot(2,2,1)
plot_axis = reordercats(categorical({'x axis OE vs CE', 'y axis OE vs CE'}),{'x axis OE vs CE', 'y axis OE vs CE'});
plot_range_x = [COPxO_range_mean COPxC_range_mean];
plot_range_y = [COPyO_range_mean COPyC_range_mean];
bar(plot_axis,[plot_range_x;plot_range_y])
title('Range mean');

subplot(2,2,2)
plot_traj_x = [COPxO_traj_mean COPxC_traj_mean];
plot_traj_y = [COPyO_traj_mean COPyC_traj_mean];
bar(plot_axis,[plot_traj_x;plot_traj_y])
title('Trajectory lenght mean');

subplot(2,2,3)
plot_vel_x = [COPxO_vel_mean COPxC_vel_mean];
plot_vel_y = [COPyO_vel_mean COPyC_vel_mean];
bar(plot_axis,[plot_vel_x;plot_vel_y])
title('Velocity mean');

subplot(2,2,4)
plot_acel_x = [COPxO_acel_mean COPxC_acel_mean];
plot_acel_y = [COPyO_acel_mean COPyC_acel_mean];
bar(plot_axis,[plot_acel_x;plot_acel_y])
title('Acceleration mean');

figure()
subplot(1,2,1)
plot_axis_range = categorical({'Range ratio OE vs CE'});
plot_RR = [ratio_open_mean ratio_closed_mean];
bar(plot_axis_range,plot_RR)
title('Range ratio mean X/Y');

subplot(1,2,2)
plot_axis_SA = categorical({'Sway area OE vs CE'});
plot_SA = [openSA_mean closedSA_mean];
bar(plot_axis_SA,plot_SA)
title('Sway area mean');

% % show trajectory
% figure()
% plot(COPxcm(5,:),COPycm(5,:));
% title('COP Trajectory');
% xlabel('COPxcm')
% ylabel('COPycm')
% %amplitude
% t=0.01:0.01:50;
% figure()
% subplot(121)
% plot(t,COPxcm(5,:));
% title('AP movement');
% xlabel('time')
% ylabel('COPxcm')
% subplot(122)
% plot(t,COPycm(5,:));
% title('ML movement');
% ylabel('COPycm')
% xlabel('time')

%% FUNTIONS
function [COP] = remove_zeros(COP)
    BC=sum(abs(COP),1)==0;
    COP(sum(abs(COP),2)==0,:)=[ ];
    COP(:,BC)=[ ];
end

function [COPcm_range_pac] = get_rango(COP)
    COPcm_range_aux=sort(transpose(COP));
    for i=1:size(COP,1)
        COPcm_range_pac(i,1) = abs(COPcm_range_aux(1,i)-COPcm_range_aux(size(COP,2),i));
    end
end

function [traj_pac] = get_length_traj(COP)
    for m=1:size(COP,1)
        traj_pto_pto = 0;
        for i=1:size(COP,2) 
            if i==1
                pos1 = 0;
            else
                pos1 = pos2;
            end
            pos2 = COP(m,i);
            traj_pto_pto = abs(pos2-pos1)+ traj_pto_pto;
        end
        traj_pac(m,1) = traj_pto_pto;
    end
end

function [vel_avg_pac, acel_avg_pac] = get_vel(COP)
    for m=1:size(COP,1)
        for i=1:size(COP,2) 
            if i==1
                pos1 = 0;
            else
                pos1 = pos2;
            end
            pos2 = COP(m,i);
            vel_temp_pac(i,1) = abs(pos2-pos1);
            
            if i==1
                vel1 = 0;
            else
                vel1 = vel_temp_pac(i,1);
            end
            
            acel_temp_pac(i,1) = abs(vel_temp_pac(i,1)-vel1);
        end
        vel_avg_pac(m,1) = mean(vel_temp_pac, 'all');
        acel_avg_pac(m,1) = mean(acel_temp_pac, 'all');
    end
end

function SA = get_sway_area(COPx, COPy)
    for m=1:size(COPx,1)
        sum = 0;
        for i=1:size(COPx,2)-1
            x1 = COPx(m,i);
            x2 = COPx(m,i+1);
            y1 = COPy(m,i);
            y2 = COPy(m,i+1);
            sum = 0.5*(abs(x2*y1-x1*y2)+sum);
        end
        SA(m,1) = sum;
    end
end

