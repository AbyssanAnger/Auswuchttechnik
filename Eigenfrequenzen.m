load("awt_UE1_Zusatz1.mat");

delta_t = mean(awt_UE1_Zusatz1(:,4))*1e-6;
Fs = 1/delta_t;
beschleunigung_1 = awt_UE1_Zusatz1(:,1);

[s,f,t] = stft(beschleunigung_1,Fs,Window=kaiser(1500,5),OverlapLength=1400,FFTLength=1500);

mesh(t,f,abs(s)./f);
view(2)
ylim([0 50])
clim([0 10])
colorbar
colormap turbo