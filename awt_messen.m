function awt_messen = awt_messen(port,t)
    device = serialport(port,500000);

%%
    flush(device);
%%
    pause(t)
%%
    data = read(device,device.NumBytesAvailable,'char');
    
%    data

    x=find(data==newline);

    data_cut = str2num(data(x(1):x(end)));

    data_V = data_cut;
    data_V(:,1:2)=(data_cut(:,1:2)-325)./1024.*5;

    data_nobias=data_V;
    data_nobias(:,1:2)=data_V(:,1:2);


    data_g=data_nobias;
    data_g(:,1:2)=data_nobias(:,1:2)/0.3;

    awt_messen=data_g;
end
