function plot_error(tt, y_real, y_val_c, y_val_d, rows, cols, idx)

subplot(rows, cols, idx)
plot(tt, y_real,'k')
hold on
plot(tt, y_val_c,'r')
hold on
plot(tt, y_val_d,'g')
hold on
legend('flux b','predicted flux c','predicted flux d')
xlabel('Time [s]')
ylabel('Magnetic Field [nT]')
set(gcf,'color','white')
hold off

subplot(rows, cols, idx+1)
plot(tt, zeros(length(y_real),1),'k--')
hold on
plot(tt, y_val_c - y_real,'r')
hold on
plot(tt, y_val_d - y_real,'g')
hold on
legend('benchmark','predicted flux c','predicted flux d')
xlabel('Time [s]')
ylabel('Absolute Error [nT]')
set(gcf,'color','white')
hold off