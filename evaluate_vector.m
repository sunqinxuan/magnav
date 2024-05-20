clear;clc;
close all;

addpath('.\data')
addpath('.\src')

%% input
data_original_filename = 'Flt1003_train.h5';
line_number = 1003.02; % 1003.02 1003.04 1003.08

data_info = h5info(data_original_filename);
data_line = h5read(data_original_filename,'/line');
i1 = find(data_line==line_number, 1 );
i2 = find(data_line==line_number, 1, 'last' );

tt = h5read(data_original_filename,'/tt');
tt = tt(i1:i2);

% load('.\data\save_nn.mat')
% add_channel_num= [30 62 37 90 27 36 32 26];

%% organize to data_x, data_y
% no embedding yet
TL_filename_name = './build/data_TL.h5';


flux_c_t=h5read(TL_filename_name,'/flux_c_t');
flux_c_x=h5read(TL_filename_name,'/flux_c_x');
flux_c_y=h5read(TL_filename_name,'/flux_c_y');
flux_c_z=h5read(TL_filename_name,'/flux_c_z');

flux_d_t=h5read(TL_filename_name,'/flux_d_t');
flux_d_x=h5read(TL_filename_name,'/flux_d_x');
flux_d_y=h5read(TL_filename_name,'/flux_d_y');
flux_d_z=h5read(TL_filename_name,'/flux_d_z');

% flux_c_t=h5read(data_original_filename,'/flux_c_t');
% flux_c_x=h5read(data_original_filename,'/flux_c_x');
% flux_c_y=h5read(data_original_filename,'/flux_c_y');
% flux_c_z=h5read(data_original_filename,'/flux_c_z');
% 
% flux_d_t=h5read(data_original_filename,'/flux_d_t');
% flux_d_x=h5read(data_original_filename,'/flux_d_x');
% flux_d_y=h5read(data_original_filename,'/flux_d_y');
% flux_d_z=h5read(data_original_filename,'/flux_d_z');

% slg = slg(i1:i2,:);
% mag_3_c = mag_3_c(i1:i2,:);
% mag_4_c = mag_4_c(i1:i2,:);
% mag_5_c = mag_5_c(i1:i2,:);
flux_c_t=flux_c_t(i1:i2,:);
flux_c_x=flux_c_x(i1:i2,:);
flux_c_y=flux_c_y(i1:i2,:);
flux_c_z=flux_c_z(i1:i2,:);

flux_d_t=flux_d_t(i1:i2,:);
flux_d_x=flux_d_x(i1:i2,:);
flux_d_y=flux_d_y(i1:i2,:);
flux_d_z=flux_d_z(i1:i2,:);

%%

% mag_1_uc = h5read(data_original_filename,'/mag_1_uc');
% mag_1_c = h5read(data_original_filename,'/mag_1_c');
% mag_1_dc = h5read(data_original_filename,'/mag_1_dc');
% mag_1_igrf = h5read(data_original_filename,'/mag_1_igrf');
flux_b_x = h5read(data_original_filename,'/flux_b_x');
flux_b_y = h5read(data_original_filename,'/flux_b_y');
flux_b_z = h5read(data_original_filename,'/flux_b_z');
flux_b_t = h5read(data_original_filename,'/flux_b_t');

% flux_b_mod=flux_b_t;
% for i=1:size(flux_d_t,1)
%     flux_b_mod(i)=sqrt(flux_b_x(i)*flux_b_x(i)+flux_b_y(i)*flux_b_y(i)+flux_b_z(i)*flux_b_z(i));
% end

% mag_1_uc = mag_1_uc(i1:i2,:);
% mag_1_c = mag_1_c(i1:i2,:);
% mag_1_dc = mag_1_dc(i1:i2,:);
% mag_1_igrf = mag_1_igrf(i1:i2,:);
flux_b_x=flux_b_x(i1:i2,:);
flux_b_y=flux_b_y(i1:i2,:);
flux_b_z=flux_b_z(i1:i2,:);
flux_b_t=flux_b_t(i1:i2,:);

%%
% half_input_wid=2;
% y_real = slg';
% y_real = y_real(half_input_wid+1:end-half_input_wid);
y_real = flux_b_t';
y_val_c=flux_c_t';
y_val_d=flux_d_t';

y_real=detrend(y_real);
y_val_c=detrend(y_val_c);
y_val_d=detrend(y_val_d);

fprintf('flux c - ensemble rmse on %f = %f\n',line_number,sqrt(mean((y_val_c - y_real).^2)));
fprintf('flux d - ensemble rmse on %f = %f\n\n',line_number,sqrt(mean((y_val_d - y_real).^2)));

% plot_error(tt, y_real, y_val_c, y_val_d)

%% plot

% figure()
% plot(tt, y_real,'k')
% hold on
% plot(tt, y_val_c,'r')
% hold on
% plot(tt, y_val_d,'g')
% hold on
% legend('flux b','predicted flux c','predicted flux d')
% xlabel('Time [s]')
% ylabel('Magnetic Field [nT]')
% set(gcf,'color','white')
% hold off
% 
% figure()
% plot(tt, zeros(length(y_real),1),'k--')
% hold on
% plot(tt, y_val_c - y_real,'r')
% hold on
% plot(tt, y_val_d - y_real,'g')
% hold on
% legend('benchmark','predicted flux c','predicted flux d')
% xlabel('Time [s]')
% ylabel('Absolute Error [nT]')
% set(gcf,'color','white')
% hold off

%%
y_real_x = flux_b_x';
y_val_c_x=flux_c_x';
y_val_d_x=flux_d_x';

y_real_x=detrend(y_real_x);
y_val_c_x=detrend(y_val_c_x);
y_val_d_x=detrend(y_val_d_x);

fprintf('flux c x - ensemble rmse on %f = %f\n',line_number,sqrt(mean((y_val_c_x - y_real_x).^2)));
fprintf('flux d x - ensemble rmse on %f = %f\n\n',line_number,sqrt(mean((y_val_d_x - y_real_x).^2)));

% plot_error(tt, y_real_x, y_val_c_x, y_val_d_x)

%%
y_real_y = flux_b_y';
y_val_c_y=flux_c_y';
y_val_d_y=flux_d_y';

y_real_y=detrend(y_real_y);
y_val_c_y=detrend(y_val_c_y);
y_val_d_y=detrend(y_val_d_y);

fprintf('flux c y - ensemble rmse on %f = %f\n',line_number,sqrt(mean((y_val_c_y - y_real_y).^2)));
fprintf('flux d y - ensemble rmse on %f = %f\n\n',line_number,sqrt(mean((y_val_d_y - y_real_y).^2)));

% plot_error(tt, y_real_y, y_val_c_y, y_val_d_y)

%%
y_real_z = flux_b_z';
y_val_c_z=flux_c_z';
y_val_d_z=flux_d_z';

y_real_z=detrend(y_real_z);
y_val_c_z=detrend(y_val_c_z);
y_val_d_z=detrend(y_val_d_z);

fprintf('flux c z - ensemble rmse on %f = %f\n',line_number,sqrt(mean((y_val_c_z - y_real_z).^2)));
fprintf('flux d z - ensemble rmse on %f = %f\n\n',line_number,sqrt(mean((y_val_d_z - y_real_z).^2)));

% plot_error(tt, y_real_z, y_val_c_z, y_val_d_z)

%%
plot_figure(tt, y_real, y_val_c, y_val_d, y_real_x, y_val_c_x, y_val_d_x, y_real_y, y_val_c_y, y_val_d_y, y_real_z, y_val_c_z, y_val_d_z)

