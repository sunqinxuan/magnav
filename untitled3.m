clear;clc;
% close all;

addpath('.\data')
addpath('.\src')

%% input
data_original_filename = 'Flt1002_train.h5';
line_number = 1002.02; % 1003.02 1003.04 1003.08

data_info = h5info(data_original_filename);
data_line = h5read(data_original_filename,'/line');
i1 = find(data_line==line_number, 1 );
i2 = find(data_line==line_number, 1, 'last' );

tt = h5read(data_original_filename,'/tt');
tt = tt(i1:i2);

flux_b_t = h5read(data_original_filename,'/flux_b_t');
flux_c_t = h5read(data_original_filename,'/flux_c_t');
flux_d_t = h5read(data_original_filename,'/flux_d_t');
mag_1_uc = h5read(data_original_filename,'/mag_1_uc');
mag_3_uc = h5read(data_original_filename,'/mag_3_uc');
mag_4_uc = h5read(data_original_filename,'/mag_4_uc');
mag_5_uc = h5read(data_original_filename,'/mag_5_uc');

flux_b_t = flux_b_t(i1:i2);
flux_c_t = flux_c_t(i1:i2);
flux_d_t = flux_d_t(i1:i2);
mag_1_uc = mag_1_uc(i1:i2);
mag_3_uc = mag_3_uc(i1:i2);
mag_4_uc = mag_4_uc(i1:i2);
mag_5_uc = mag_5_uc(i1:i2);

figure;
plot(tt,flux_b_t,'--');
hold on;
plot(tt,flux_c_t,'--');
hold on;
plot(tt,flux_d_t,'--');
hold on;
plot(tt,mag_1_uc);
hold on;
plot(tt,mag_3_uc);
hold on;
plot(tt,mag_4_uc);
hold on;
plot(tt,mag_5_uc);
hold on;

legend('flux b','flux c','flux d','mag 1','mag 3','mag 4','mag 5')


