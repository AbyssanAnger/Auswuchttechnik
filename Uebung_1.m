load("awt_UE1_Aufgabe.mat");

accel_1 = awt_UE1(:,1);
accel_2 = awt_UE1(:,2);
trigger = awt_UE1(:,3);
time = awt_UE1(:,4); % in µs

t_end = sum(time)*10e-6; % in s

delta_t = mean(awt_UE1(:,4));

f_s = 1/delta_t; % in Hz

t = linspace(0, t_end, 3275);

%% Plotting raw data

tiledlayout(3,1)

nexttile
fig1 = figure(1);

plot(t, accel_1);
title("Beschleunigung 1")

nexttile
plot(t, accel_2);
title("Beschleunigung 2")

nexttile
plot(t, trigger);
title("Triggerfunktion")

%% Zuschnitt Triggerfunktion
diff_trigg = diff(trigger); % diff function for peak finding
first_index = find(diff_trigg == 1, 1, "first") + 1;
last_index = find(diff_trigg == 1, 1, "last"); % index of last non zero element

trigger_cut = trigger(first_index:last_index);
t_cut = t(1: length(trigger_cut)); %

%% Zuschnitt Beschleunigung
accel_1_cut = accel_1(first_index:last_index);
accel_2_cut = accel_2(first_index:last_index);

%% Plot der zugeschnittenen Signale
tiledlayout(3, 1)

fig2 = figure(2);

nexttile
plot(t_cut, accel_1_cut);
title("geschnittene Beschleunigung 1")

nexttile
plot(t_cut, accel_2_cut);
title("geschnittene Beschleunigung 2")

nexttile
plot(t_cut, trigger_cut);
title("geschnittenes Triggersignal")

%% FFT der Signale
L = length(t_cut); % Length of signal

fft_trigger = 2*fft(trigger_cut)/L;

fig_7 = figure(7);
plot(f_s/L*(0:L-1), abs(fft_trigger))

trigger_amplitude = max(abs(fft_trigger(2:end/2)));
trigger_index = find(abs(fft_trigger) == trigger_amplitude);
trigger_index = trigger_index(1);

fft_accel_1 = 2*fft(accel_1_cut)/L;
accel_1_frequency = fft_accel_1(trigger_index+1);
accel_1_amplitude = max(abs(fft_accel_1(2:end/2)));
accel_1_angle = angle(fft_accel_1(trigger_index))*180/pi;

fft_accel_2 = 2*fft(accel_2_cut)/L;
accel_2_frequency = fft_accel_2(trigger_index+1);
accel_2_amplitude = max(abs(fft_accel_2(2:end/2)));
accel_2_angle = angle(fft_accel_2(trigger_index))*180/pi;

%% 3D Plot der Fourier Signale



%% Plot Polarcoordinates

fig_8 = figure(8);
polarplot([0 accel_1_angle/180*pi], [0 accel_1_amplitude], "black-o", "DisplayName", "Beschleunigungssensor 1");
hold on
polarplot([0 accel_2_angle/180*pi], [0 accel_2_amplitude], "magenta-o", "DisplayName", "Beschleunigungssensor 2");
hold off
legend show

%% instationäre Drehzahlen
stft(trigger_cut, f_s, Window=kaiser(256,5),OverlapLength=220,FFTLength=512)