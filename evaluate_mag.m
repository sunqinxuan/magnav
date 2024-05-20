clear;clc;
close all;

addpath('.\data')
addpath('.\src')

%% input
data_original_filename = 'Flt1003_train.h5';
line_number = 1003.08; % 1003.02 1003.08 1003.04

% load('.\data\save_nn.mat')
% add_channel_num= [30 62 37 90 27 36 32 26];

%% organize to data_x, data_y
% no embedding yet
TL_filename_name = './build/data_TL.h5';

mag_3_c = h5read(TL_filename_name,'/mag_3_c');
mag_4_c = h5read(TL_filename_name,'/mag_4_c');
mag_5_c = h5read(TL_filename_name,'/mag_5_c');

%%
data_info = h5info(data_original_filename);
data_line = h5read(data_original_filename,'/line');
i1 = find(data_line==line_number, 1 );
i2 = find(data_line==line_number, 1, 'last' );

tt = h5read(data_original_filename,'/tt');
tt = tt(i1:i2);

%%

mag_3_c = mag_3_c(i1:i2,:);
mag_4_c = mag_4_c(i1:i2,:);
mag_5_c = mag_5_c(i1:i2,:);
% data_x = [mag_3_c,mag_4_c,mag_5_c];

mag_1_uc = h5read(data_original_filename,'/mag_1_uc');
mag_1_c = h5read(data_original_filename,'/mag_1_c');
mag_1_dc = h5read(data_original_filename,'/mag_1_dc');
mag_1_igrf = h5read(data_original_filename,'/mag_1_igrf');

mag_1_uc = mag_1_uc(i1:i2,:);
mag_1_c = mag_1_c(i1:i2,:);
mag_1_dc = mag_1_dc(i1:i2,:);
mag_1_igrf = mag_1_igrf(i1:i2,:);

%%
% half_input_wid=2;
% y_real = slg';
% y_real = y_real(half_input_wid+1:end-half_input_wid);
y_real = mag_1_uc';
y_val_3=mag_3_c';
y_val_4=mag_4_c';
y_val_5=mag_5_c';

y_real=detrend(y_real);
y_val_3=detrend(y_val_3);
y_val_4=detrend(y_val_4);
y_val_5=detrend(y_val_5);

fprintf('mag 3 - ensemble rmse on %f = %f\n',line_number,sqrt(mean((y_val_3 - y_real).^2)));
fprintf('mag 4 - ensemble rmse on %f = %f\n',line_number,sqrt(mean((y_val_4 - y_real).^2)));
fprintf('mag 5 - ensemble rmse on %f = %f\n',line_number,sqrt(mean((y_val_5 - y_real).^2)));

%% plot

figure()
plot(tt, y_real,'k')
hold on
plot(tt, y_val_3,'r')
hold on
plot(tt, y_val_4,'g')
hold on
plot(tt, y_val_5,'b')
legend('SGL Mag 1','predicted 3','predicted 4','predicted 5')
xlabel('Time [s]')
ylabel('Magnetic Field [nT]')
set(gcf,'color','white')
hold off

figure()
plot(tt, zeros(length(y_real),1),'k--')
hold on
plot(tt, y_val_3 - y_real,'r')
hold on
plot(tt, y_val_4 - y_real,'g')
hold on
plot(tt, y_val_5 - y_real,'b')
legend('benchmark','predicted 3','predicted 4','predicted 5')
xlabel('Time [s]')
ylabel('Absolute Error [nT]')
set(gcf,'color','white')
hold off