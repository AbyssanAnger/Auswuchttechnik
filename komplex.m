function [komplex_1, komplex_2] = komplex(X)

beschleunigung_1 = X(:,1);

beschleunigung_2 = X(:,2);

trigger = X(:,3);

t = linspace(0,sum(X(:,4))*1e-6,length(trigger));

delta_t = mean(X(:,4))*1e-6;



differenz = diff(trigger);

indices = find(differenz==1);


t_cut = t(1:indices(end));

beschleunigung_1_cut = beschleunigung_1(indices(1)+1:indices(end),1);

beschleunigung_2_cut = beschleunigung_2(indices(1)+1:indices(end),1);

trigger_cut = trigger(indices(1)+1:indices(end),1);


L = length(trigger_cut);
Fs = 1/delta_t;


trigger_fft = fft(trigger_cut)/L*2;
f = Fs/L*(0:L-1);
peak = max(abs(trigger_fft(2:end)));
index_max = find(abs(trigger_fft)==peak);

%plot(f(1:end/2),abs(trigger_fft(1:end/2)))

beschleunigung_1_fft = fft(beschleunigung_1_cut)/L*2;

beschleunigung_2_fft = fft(beschleunigung_2_cut)/L*2;


%plot(f(1:end/2),abs(beschleunigung_1_fft(1:end/2)))

%plot(f(1:end/2),abs(beschleunigung_2_fft(1:end/2)))


komplex_1 = beschleunigung_1_fft(index_max(1));

komplex_2 = beschleunigung_2_fft(index_max(1));


end