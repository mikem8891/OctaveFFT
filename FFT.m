# clear commandline and figures
clc;
clear all;
close all;

# open file
offset = 1; % first row of data
nop = 2^8;  % number of points taken from the file
end_row = offset + nop;
data = dlmread('accel.csv', '', [offset 0 end_row 1]);

# parse data
time_raw  = data(:, 1); % 1st column is time
accel_raw = data(:, 2); % 2nd column is accelerations

f1 = figure('Name', 'Time domain accelerations');
plot(time_raw, accel_raw, 'o');
hold on

# create new data set with evenly spaced data point
L = length(time_raw);     % Length of signal
T = time_raw(L)/(L-1);    % Sample period
Fs = 1/T;                % Sampling frequency (hz)
time = (0:L-1)*T;        % Time
accel = interp1(time_raw, accel_raw, time, 'cubic');

plot(time, accel, 'x');
hold off
xlabel("t (seconds)");
ylabel("acceleration (G)");

# Fourier Transform

C = fft(accel);

# Fourier Coefficents
A =  2*real(C(1:L/2+1))/L; % cos coeffients
B = -2*imag(C(1:L/2+1))/L; % sin coeffients
A(1) = C(1)/L;             % offset

C_mag = sqrt(A.^2 + B.^2); % Magnitudes
C_ang = atan2(B, A);       % Phase angles
f = Fs*(1:(L/2))/L;        % frequency domain (Fs/L to Fs/2)

f2 = figure('Name', 'Frequency domain accelerations');
semilogx(f, C_mag(2:end));

# Fourier Series (for debugging)
j=1;
t_plot_max = min(100*T, time_raw(L));
time_i = (0:T/5:t_plot_max);

for t = time_i
  FS(j) = 0;
  k = (L-1) * t / time_raw(L);
  dtheta = 2*pi*k/L;
  theta = 0;
  for i = 1:L/2+1
    theta = (i-1)*dtheta;
    FS(j) = FS(j) + C_mag(i)*cos(theta - C_ang(i));
  end
  j = j+1;
end

figure(f1);
hold on
plot(time_i, FS);
hold off
xlim([0 t_plot_max]);



