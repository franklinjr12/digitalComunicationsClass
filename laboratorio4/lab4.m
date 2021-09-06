close all;
Fs = 24000;
T = 1/8000;
r = Fs*T; % Fator de oversampling
t = -5*T:1/Fs:5*T;
t = t + 1e-10;
% alfa = 0.5;
alfa = 0.00001;

% Filtro raiz de cosseno levantado
h = rcosine(1/T,Fs,'sqrt',alfa,5);
figure;plot(t,h)
xlabel('Tempo (intervalo de símbolo)')
title('Raiz de cosseno levantado')
legend('\alpha=0.5');
grid           


%generate y signal with random bits
Nbits = 10;
Vamp = 1;
bits = zeros(1,Nbits);
for i=1:Nbits
    bits(i) = rand > 0.5;
end
s=2*bits-1;
s_up = zeros(1,length(s)*r);
s_up(1:r:r*length(s)) = s*Vamp; % sequência com oversampling
x=conv(s_up,h);
%add transmission noise
noise = (rand(1,length(x))-0.5)*2;

s_trans = x+noise;
figure
eixo = (0:length(x)-1)*1/Fs;
points = zeros(1,length(x));
t_start_symbol = 0.000625;
for i=1:length(eixo)
    if (eixo(i) >= t_start_symbol && eixo(i) < T*Nbits + t_start_symbol)
        if (mod(eixo(i),T)==0)
            if (x(i) > 0)
                points(i)=Vamp;
            else
                points(i)=-Vamp;
            end
        end
    end
end    
plot(eixo,x,'r');
hold on
stem(eixo, points);
grid
xlabel('Tempo (intervalo de símbolo)')
title('Sinal transmitido')
%0.00075-0.000625

figure
eixo = (0:length(s_trans)-1)*1/Fs;
plot(eixo,s_trans,'r');
grid
xlabel('Tempo (intervalo de símbolo)')
title('Sinal transmitido com ruido')

receiver = conv(s_trans,h);
figure;
eixo = (0:length(receiver)-1)*1/Fs;
points = zeros(1,length(receiver));
plot(eixo,receiver);
hold on
grid
xlabel('Tempo (intervalo de símbolo)')
title('Sinal recebido')

bits_received = zeros(1,Nbits);
j=1;
t_start_symbol = 0.00125;
for i=1:length(receiver)    
    if (eixo(i) >= t_start_symbol && eixo(i) < T*Nbits + t_start_symbol ...
            && j <= Nbits)
        eixo(i)
        if (mod(eixo(i),T)==0)
            if (receiver(i) > 0)
                points(i)=Vamp;
                bits_received(j)=1;
                j=j+1;
            else
                points(i)=-Vamp;
                bits_received(j)=0;
                j=j+1;                
            end
        end
    end
end
stem(eixo,points);
wrong_bits = sum( bits ~= bits_received)