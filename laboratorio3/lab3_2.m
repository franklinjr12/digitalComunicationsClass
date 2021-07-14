close all;
clear all;

rand('state',0);
randn('state',0);

bits=1e6; %Número de bits a serem simulados.
% bits=10; %Número de bits a serem simulados.

M=4;

 b=(rand(1,bits)); %Geração de bits 0 e 1 equiprováveis.

% temp=b;
% for i=M:-1:1
%     temp
%     temp(b>((i-1)/M))=2*i-1;
% end
% b=temp;
% disp(b);

for i=1:bits
    if (b(i) < 1/4) 
        b(i) = 1;
    elseif (b(i) < 2/4) 
        b(i) = 3;
    elseif (b(i) < 3/4) 
        b(i) = 5;
    else
        b(i) = 7;
    end
end

s=b;
N0=1; %N0 será fixa em 1.
%EMF foi suposto como 1, portanto a variância do ruído após
%o filtro casado é N0/2 apenas.
%ruido proporcional ao máximo do meu sinal
n=(M*2-1)*randn(1,bits)*sqrt(N0/2); 

ber_array = [];
EbN0_array = [];

EbN0dB=0; %Valor de Eb/N0 a ser considerado na simulação.
while (EbN0dB < 10)
    
    EbN0=10^(EbN0dB/10); %Eb/N0 em escala linear.
    Eb=EbN0*N0; %Eb requerido para atingir a razão Eb/N0 de interesse.
    Es=Eb*log2(M); %Es calculado a partir de Eb. Como a modulação é binária Es=Eb.

    y=(sqrt(Es)*s+n);

    b_est=zeros(1,bits); %Decisor.
    for i=1:bits
        if (y(i) < 3*sqrt(Es)) 
            b_est(i) = 1;
        elseif (y(i) < 5*sqrt(Es)) 
            b_est(i) = 3;
        elseif (y(i) < 7*sqrt(Es)) 
            b_est(i) = 5;
        else
            b_est(i) = 7;
        end
    end
%     n
%     y
%     b_est
    erros=sum(b~=b_est); %Contagem de erros.
    ber=erros/bits; %Cálculo da BER.
    Pb=qfunc(sqrt(2*EbN0)); %BER teórica.   
    
    ber_array = [ber_array ber];
    EbN0_array = [EbN0_array EbN0];
    
    EbN0dB = EbN0dB + 1;
end

plot(EbN0_array, ber_array);
title('EbN0 vs BER');
xlabel('EbN0');
ylabel('BER');