load("awt_UE1_Aufgabe.mat");
%% Daten aufnehmen
Null_Lauf=awt_messen("COM5", 5);
%% Daten aufnehmen
Ein_Lauf=awt_messen("COM5", 5); % Testgewicht angebracht, position in Winkel und Masse dokumentiert

%% Daten Laden

null_lauf = load("Null_Lauf.mat", "data");
null_lauf = null_lauf.data;

ein_lauf = load("Testlauf_Eins.mat", "data");
ein_lauf = ein_lauf.data_test_1;

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
function [freq_index, amp, angle_deg] = signal_fft(signal, f_s)
    L = length(signal);
    fft_signal = 2*fft(signal)/L;

    % Max. Amplitude (ignoriere DC-Anteil)
    amp = max(abs(fft_signal(2:end/2)));

    % Index des Maximums
    freq_index = find(abs(fft_signal) == amp, 1);

    % Phasenwinkel in Grad
    angle_deg = angle(fft_signal(freq_index)) * 180/pi;
end

%% 3D Plot der Fourier Signale
plot()

%% Abfrage nach positivem Massenausgleich
% answer = questdlg('Möchten Sie einen positiven Massenausgleich?', ...
% 	'Dessert Menu', ...
% 	'Ice cream','Cake','No thank you','No thank you');
% % Handle response
% switch answer
%     case 'Ja'
%         disp([answer ' Es wird ein positiver Massenausgleich verwendet.'])
%         dessert = 1;
%     case 'Cake'
%         disp([answer ' coming right up.'])
%         dessert = 2;
%     case 'No thank you'
%         disp('I''ll bring you your check.')
%         dessert = 0;
% end

%% Plot Polarcoordinates
% Option für Einstellen von Winkeln an denen tatsächlich ein Setzungsgewicht anbringbar ist, z.B. bei vorgegebenen Löchern.
% Option für Wiedergabe von negativem Massenausgleich => 180 Grad versetzt

fig_8 = figure(8);
polarplot([0 accel_1_angle/180*pi], [0 accel_1_amplitude], "black-o", "DisplayName", "Beschleunigungssensor 1");
hold on
polarplot([0 accel_2_angle/180*pi], [0 accel_2_amplitude], "magenta-o", "DisplayName", "Beschleunigungssensor 2");
hold off
legend show

%% instationäre Drehzahlen
stft(trigger_cut, f_s, Window=kaiser(256,5),OverlapLength=220,FFTLength=512)

%% U_test in komplexe Zahl umrechnen
gewicht = 3,6; % in gramm
winkel = 90; % in Grad

betrag = gewicht*cosd(winkel);
komplexer_winkel = gewicht*sind(winkel);

u_test = betrag + 1i*komplexer_winkel;
u_test;

% b_hat = (accel_1_frequency - (-0,1128+1i*0,036))/u_test;
% u_wuchtsetzung = - (-0,1128+1i*0,036/b_hat);


%% Daten Laden
tmp = load("Null_Lauf.mat", "data");
null_lauf = tmp.data;

tmp = load("Testlauf_Eins.mat", "data_test_1");
ein_lauf = tmp.data_test_1;

%% Zeitwerte kalkulieren
[t_0, f_s_0] = calculate_time(null_lauf);
[t_1, f_s_1] = calculate_time(ein_lauf);

%% Zuschnitt
[accel_1_cut, accel_2_cut, trigger_cut, t_cut] = ...
    cut_signals(null_lauf(:,1), null_lauf(:,2), null_lauf(:,3), t_0);

%% Plotting raw data
plot_raw_data(null_lauf, t_0);
plot_raw_data(ein_lauf, t_1);

%% ====== FUNKTIONEN ======
function [t, f_s] = calculate_time(data)
    time = data(:,4);                 % Zeitspalte (µs)
    delta_t = mean(time) * 1e-6;      % Abtastschritt [s]
    f_s = 1/delta_t;                  % Samplingfrequenz [Hz]
    t_end = length(time) * delta_t;   % Gesamtdauer [s]
    t = linspace(0, t_end, length(time));
end

function [accel_1_cut, accel_2_cut, trigger_cut, t_cut] = cut_signals(accel_1, accel_2, trigger, t)
    diff_trigg = diff(trigger);
    first_index = find(diff_trigg > 0.5, 1, "first") + 1;
    last_index  = find(diff_trigg > 0.5, 1, "last");

    trigger_cut = trigger(first_index:last_index);
    accel_1_cut = accel_1(first_index:last_index);
    accel_2_cut = accel_2(first_index:last_index);
    t_cut = t(1:length(trigger_cut));
end

function plot_raw_data(Lauf, t)
    figure;
    tiledlayout(3,1);

    nexttile;
    plot(t, Lauf(:,1));
    title("Beschleunigung 1");

    nexttile;
    plot(t, Lauf(:,2));
    title("Beschleunigung 2");

    nexttile;
    plot(t, Lauf(:,3));
    title("Triggerfunktion");
end

function Testsetzung(gewicht, winkel)
    betrag = gewicht*cos(winkel)