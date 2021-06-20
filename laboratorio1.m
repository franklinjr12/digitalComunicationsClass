close all;

Fsampling = 50;
Tsymbol = 1;
Vamp = 0.1;
Nbits = 100;

t = 0:1/Fsampling:Nbits-1/Fsampling;
bit1 = Vamp*ones(1,Fsampling);
bit0 = -Vamp*ones(1,Fsampling);
% y = [bit1 bit0 bit1 bit0 bit1];
% Nbits = length(t)/length(bit1);
y = zeros(1,length(t));
for i=1:length(t)/Nbits:length(t)
    if (randi(100,1) <= 50)
        y(i:i+length(bit1)-1) = bit0;
    else
        y(i:i+length(bit1)-1) = bit1;
    end
end

 
figure;
plot(t,y);
xlabel('tempo (s)');
title('generator');
axis([0 Nbits -1.2*Vamp 1.2*Vamp]);

noise = randn(1,length(y));
yr = y+noise;

figure;
plot(t,yr);
xlabel('tempo (s)');
title('signal');

yfiltered = lowpass(yr,Fsampling/2.2, Fsampling);
figure;
plot(yfiltered);
title('yfiltered');

ytruncated = yfiltered;

step = length(y)/Nbits;
for i = 1:step:length(y)
    received_bit =  mean(ytruncated(i:i+step-1));
    if received_bit > 0 
        ytruncated(i:i+step-1) = Vamp;
    else
        ytruncated(i:i+step-1) = -Vamp;
    end 
end
figure;
plot(t, ytruncated);
axis([0 Nbits -1.2*Vamp 1.2*Vamp]);
xlabel('tempo (s)');
title('receptor');

n_wrong_bits = 0;
for i = 1:step:length(y)
    generated_bit = mean(y(i:i+step-1));
    if generated_bit > 0 
        generated_bit = Vamp;
    else
        generated_bit = -Vamp;
    end
    received_bit =  mean(ytruncated(i:i+step-1));
    if received_bit > 0 
        received_bit = Vamp;
    else
        received_bit = -Vamp;
    end    
    if generated_bit ~= received_bit
        n_wrong_bits = n_wrong_bits+1;
    end
end
%=================================================
%=================================================

%initial values and arrays
bit_error_rate = [];
wave_energy = [];
k=4;
max_bit_err = 10^-k; 
Nbits=100*10^k;
% Nbits=20;

Fsampling = 50;
Tsymbol = 1;
Vamp = 0.1;

t = 0:1/Fsampling:Nbits-1/Fsampling;
y = zeros(1,length(t));

profile on

%iterate until 10dB
while(10*log10((Vamp^2)*Fsampling*Tsymbol)<10)

    %generate y signal with random bits
    bit1 = Vamp*ones(1,Fsampling);
    bit0 = -Vamp*ones(1,Fsampling);
    for i=1:length(t)/Nbits:length(t)
        if (randi(100,1) <= 50)
            y(i:i+length(bit1)-1) = bit0;
        else
            y(i:i+length(bit1)-1) = bit1;
        end
    end 
    
    %add transmission noise
    noise = randn(1,length(y));
    yr = y+noise;

    %filter received signal
    yfiltered = lowpass(yr,Fsampling/2.2, Fsampling);

    %digitalize
    ytruncated = yfiltered;
    for i = 1:step:length(y)
        received_bit =  mean(ytruncated(i:i+step-1));
        if received_bit > 0 
            ytruncated(i:i+step-1) = Vamp;
        else
            ytruncated(i:i+step-1) = -Vamp;
        end 
    end

    %calculate number of wrong bits
    n_wrong_bits = 0;
    step = length(y)/Nbits;
    for i = 1:step:length(y)
        generated_bit = mean(y(i:i+step-1));
        if generated_bit > 0 
            generated_bit = Vamp;
        else
            generated_bit = -Vamp;
        end
        received_bit =  mean(ytruncated(i:i+step-1));
        if received_bit > 0 
            received_bit = Vamp;
        else
            received_bit = -Vamp;
        end    
        if generated_bit ~= received_bit
            n_wrong_bits = n_wrong_bits+1;
        end
    end

    %calculate statistics for making final graphs
    bit_error_rate = [bit_error_rate, (n_wrong_bits/Nbits)];
    wave_energy = [wave_energy , 10*log10(Vamp^2*Fsampling*Tsymbol)];
    10*log10(Vamp^2*Fsampling*Tsymbol)
    
    %update Amplitude for next iteraction    
    Vamp = Vamp+0.1;
end

figure();
semilogy(wave_energy, bit_error_rate);
axis([0 max(wave_energy) 0 1.2*max(bit_error_rate)]);
xlabel('wave energy dB');
ylabel('bit error rate log');
disp('Final amplitude: ');
disp(Vamp);