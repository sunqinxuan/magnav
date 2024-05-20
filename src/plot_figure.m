function plot_figure(tt, y_real, y_val_c, y_val_d, y_real_x, y_val_c_x, y_val_d_x, y_real_y, y_val_c_y, y_val_d_y, y_real_z, y_val_c_z, y_val_d_z)

figure('Position', [400, 100, 1000, 1200])

plot_error(tt, y_real, y_val_c, y_val_d,4,2,1);
plot_error(tt, y_real_x, y_val_c_x, y_val_d_x,4,2,3);
plot_error(tt, y_real_y, y_val_c_y, y_val_d_y,4,2,5)
plot_error(tt, y_real_z, y_val_c_z, y_val_d_z,4,2,7)