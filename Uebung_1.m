load("awt_UE1_Aufgabe.mat");

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



%% instationäre Drehzahlen
stft(trigger_cut, f_s, Window=kaiser(256,5),OverlapLength=220,FFTLength=512)

% b_hat = (accel_1_frequency - (-0,1128+1i*0,036))/u_test;
% u_wuchtsetzung = - (-0,1128+1i*0,036/b_hat);

%% Daten aufnehmen
null_lauf = awt_messen("COM5", 5);

%% Daten Laden
% tmp = load("null_lauf.mat", "data");
% null_lauf = tmp.data;

%% Zeitwerte kalkulieren
[t_0, f_s_0] = calculate_time(null_lauf);

%% Zuschnitt
[accel_1_cut_0, accel_2_cut_0, trigger_cut_0, t_cut_0] = ...
    cut_signals(null_lauf(:,1), null_lauf(:,2), null_lauf(:,3), t_0); % Null-Lauf geschnitten

%% Plotting raw data
plot_raw_data(null_lauf, t_0, "(0-Lauf)");

%% Plotting cut data
plot_cut_data(accel_1_cut_0, accel_2_cut_0, trigger_cut_0, t_cut_0, "(0-Lauf) geschnitten");

%% FFT
% FFT für alle Signale erstellen
% Amplitude von Triggerfunktion aus dem Frequenzbereich bestimmen
% Index der Amplitude der Triggerfunktion finden
% Den Index auf die Beschleunigungen im Frequenzbereich anwenden und dort
% die Werte auslesen
[trigger_0_fft, trigger_cut_0_index] = trigger_data(trigger_cut_0);
[accel_1_0_fft, accel_1_0_amp, accel_1_0_angle] = create_fft(accel_1_cut_0, trigger_cut_0_index);
[accel_2_0_fft, accel_2_0_amp, accel_2_0_angle] = create_fft(accel_2_cut_0, trigger_cut_0_index);

%% FFT plotting
[accel_1_0_fft_plot] = plot_fft(accel_1_0_fft, f_s_0, "Beschleunigung 1 (0-Lauf)");
[accel_2_0_fft_plot] = plot_fft(accel_2_0_fft, f_s_0, "Beschleunigung 2 (0-Lauf)");
[trigger_0_fft_plot] = plot_fft(trigger_0_fft, f_s_0, "Triggerfunktion (0-Lauf)");

%% Komplexe Zahl bilden

accel_1_0 = accel_1_0_amp + 1i*accel_1_0_angle;
accel_2_0 = accel_2_0_amp + 1i*accel_1_0_angle;
%% 3D Plot der Fourier Signale
% tiled(trigger_fft_plot, accel_1_0_fft_plot, accel_2_0_fft_plot, "0-Lauf");
%% Komplexe Zahlen bilden

accel_1_0 = accel_1_0_amp + 1i*accel_1_0_angle;
accel_2_0 = accel_2_0_amp + 1i*accel_2_0_angle;
%% Polar coordinate plot
plot_polarcoordinates(accel_1_0_amp, accel_1_0_angle, accel_2_0_amp, accel_2_0_angle);

%% Umrechnung in komplexe Werte
testsetzung_gewicht = 3.38;
testsetzung_winkel = 0;

u_test = testsetzung(testsetzung_gewicht, testsetzung_winkel);

%% Testsetzung Daten aufnehmen
ein_lauf=awt_messen("COM5", 5); % Testgewicht angebracht, position in Winkel und Masse dokumentiert

%% Zeitwerte kalkulieren
[t_1, f_s_1] = calculate_time(ein_lauf);

%% Zuschnitt
[accel_1_cut_1, accel_2_cut_1, trigger_cut_1_1, t_cut_1] = ...
    cut_signals(ein_lauf(:,1), ein_lauf(:,2), ein_lauf(:,3), t_1); % Ein-Lauf geschnitten

%% Plotting raw data
plot_raw_data(ein_lauf, t_1,  "(1-Lauf)");

%% Plotting cut data
plot_cut_data(accel_1_cut_1, accel_2_cut_1, trigger_cut_1_1, t_cut_1, "(1-Lauf) geschnitten");

%% FFT
% FFT für alle Signale erstellen
% Amplitude von Triggerfunktion aus dem Frequenzbereich bestimmen
% Index der Amplitude der Triggerfunktion finden
% Den Index auf die Beschleunigungen im Frequenzbereich anwenden und dort
% die Werte auslesen
[trigger_1_1_fft, trigger_cut_1_1_index] = trigger_data(trigger_cut_1_1);
[accel_1_1_fft, accel_1_1_amp, accel_1_1_angle] = create_fft(accel_1_cut_1, trigger_cut_1_1_index);
[accel_2_1_fft, accel_2_1_amp, accel_2_1_angle] = create_fft(accel_2_cut_1, trigger_cut_1_1_index);

%% FFT plotting
[accel_1_1_fft_plot] = plot_fft(accel_1_1_fft, f_s_1, "Beschleunigung 1 (1-Lauf)");
[accel_2_1_fft_plot] = plot_fft(accel_2_1_fft, f_s_1, "Beschleunigung 2 (1-Lauf)");
[trigger_1_1_fft_plot] = plot_fft(trigger_1_1_fft, f_s_1, "Triggerfunktion (1-Lauf)");

%% Komplexe Zahl bilden

accel_1_1 = accel_1_1_amp + 1i*accel_1_1_angle;
accel_2_1 = accel_2_1_amp + 1i*accel_1_1_angle;
%% Plot FFT 0-Lauf zu 1-Lauf
% Zuschnitt 1 Lauf auf gleiche Größe wie 0 Lauf!

%% Berechnung b
b_hat = (accel_1_1 - accel_1_0)/u_test;

%% Berechnung Wuchsetzung u_hat
u_hat = - (accel_1_0/b_hat);

%% Von Komplex zu real
[gewicht_1, winkel_1] = complex_to_real(u_hat);

%% plot Polar 
fig8 = figure(8);

polarplot([0 accel_1_0_angle/180*pi], [0 accel_1_0_amp], "black-o", "DisplayName", "a_0");
hold on
polarplot([0 accel_1_1_angle/180*pi], [0 accel_1_1_amp], "magenta-o", "DisplayName", "a_test");
hold on
polarplot([accel_1_0_angle/180*pi accel_1_1_angle/180*pi], [accel_1_0_amp accel_1_1_amp], "magenta-o", "DisplayName", "differenz");
hold off
legend show
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
    last_index  = find(diff_trigg == 1, 1, "last");

    trigger_cut = trigger(first_index:last_index);
    accel_1_cut = accel_1(first_index:last_index);
    accel_2_cut = accel_2(first_index:last_index);
    t_cut = t(1:length(trigger_cut));
end

function plot_raw_data(Lauf, t, overtitle)
    figure;
    tlc = tiledlayout(3,1);

    nexttile;
    plot(t, Lauf(:,1));
    title("Beschleunigung 1");

    nexttile;
    plot(t, Lauf(:,2));
    title("Beschleunigung 2");

    nexttile;
    plot(t, Lauf(:,3));
    title("Triggerfunktion");

    title(tlc, overtitle);
end

function plot_cut_data(accel_1, accel_2, trigger, t, overtitle)
    figure;
    tlc = tiledlayout(3,1);

    nexttile;
    plot(t, accel_1);
    title("Beschleunigung 1");

    nexttile;
    plot(t, accel_2);
    title("Beschleunigung 2");

    nexttile;
    plot(t, trigger);
    title("Triggerfunktion");

    title(tlc, overtitle)
end

function [fft_signal, index, amp] = trigger_data(signal)
    L = length(signal);
    fft_signal = 2*fft(signal)/L;

    amp = max(abs(fft_signal(2:end/2)));

    index = find(abs(fft_signal) == amp);
    index = index(1);
end

function [fft_signal, amp, angle_deg] = create_fft(signal, trigger_index)
    L = length(signal);
    fft_signal = 2*fft(signal)/L;

    freq = fft_signal(trigger_index + 1);
    amp = max(abs(fft_signal(2:end/2)));
    angle_deg = angle(fft_signal(trigger_index))*180/pi;
end

function [fig] = plot_fft(fft_signal, f_s,  fft_signal_name)
    fig = figure;
    L = length(fft_signal(1:end/2));
    f = f_s/L*(1:L);
    signal = abs(fft_signal(1:end/2));
    plot(f, signal);
    title(fft_signal_name)
end

function tiled(fig1, fig2, fig3, titel)
    figure;
    tlc = tiledlayout(3,1);

    nexttile;
    fig1;
    title("Beschleunigung 1");

    nexttile;
    fig2;
    title("Beschleunigung 2");

    nexttile;
    fig3;
    title("Triggerfunktion");

    title(tlc, titel)

end

function [complex] = testsetzung(gewicht, winkel)
    real = gewicht*cosd(winkel);
    imaginary = gewicht*sind(winkel);
    complex = real + 1i*imaginary;
end

function [gewicht, winkel] = complex_to_real(complex)
    gewicht = abs(complex);
    winkel = rad2deg(angle(complex));
end

function plot_polarcoordinates(amp1, angle1, amp2, angle2)
% Option für Einstellen von Winkeln an denen tatsächlich ein Setzungsgewicht anbringbar ist, z.B. bei vorgegebenen Löchern.
% Option für Wiedergabe von negativem Massenausgleich => 180 Grad versetzt
    figure();
    polarplot([0 angle1/180*pi], [0 amp1], "black-o", "DisplayName", "Beschleunigungssensor 1");
    hold on
    polarplot([0 angle2/180*pi], [0 amp2], "magenta-o", "DisplayName", "Beschleunigungssensor 2");
    hold off
    legend show
end