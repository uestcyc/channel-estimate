clear all;

close all;

fprintf('OFDM信道估计仿真\n\n');

carrier_count=64;%-----------载波数目

num_symbol=50;%--------------OFDM符号个数

Guard=8;%--------------------循环前缀

pilot_Inter=8;%--------------导频间隔

modulation_mode=4;%---------调制方式

SNR=[1:1:12];%-------------信噪比取值

NumLoop=10;%-----------------循环次数

num_bit_err=zeros(length(SNR),NumLoop);

num_bit_err_dft=zeros(length(SNR),NumLoop);

num_bit_err_ls=zeros(length(SNR),NumLoop);

num_bit_err_mmse=zeros(length(SNR),NumLoop);

MSE=zeros(length(SNR),NumLoop);

MSE1=zeros(length(SNR),NumLoop);

MSE2=zeros(length(SNR),NumLoop);

%%%%%%%%%%%%%%%%%%%%%%%主程序循环%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for c1=1:length(SNR)

fprintf('\n\n\n仿真信噪比为%f\n\n',SNR(c1));

for num1=1:NumLoop

%---------------产生发送的随机序列——————————————-

BitsLen=carrier_count*num_symbol;

BitsTx=randi([0 3],1,BitsLen);

%---------------符号调制---------------------------------------

Modulated_Sequence=qammod(BitsTx,modulation_mode);

%---------------导频格式---------------------------------------

% pilot_len=carrier_count;
% 
% pilot_symbols=round(rand(1,pilot_len));
% 
% for i=1:pilot_len
% 
% if pilot_symbols(1,i)==0
% 
% pilot_symbols(1,i)=pilot_symbols(1,i)-1;
% 
% else
% 
% pilot_symbols(1,i)=pilot_symbols(1,i);
% 
% end
% 
% end
% 
% pilot_symbols=pilot_symbols';
 pilot_symbols=qammod(randi([0 3],64,1),4);
%----------------计算导频和数据数目----------------------------

num_pilot=ceil(num_symbol/pilot_Inter);

if rem(num_symbol,pilot_Inter)==0

num_pilot=num_pilot+1;

end

num_data=num_symbol+num_pilot;

%----------------导频位置计算----------------------------------

pilot_Indx=zeros(1,num_pilot);

Data_Indx=zeros(1,num_pilot*(pilot_Inter+1));

for i=1:num_pilot-1

pilot_Indx(1,i)=(i-1)*(pilot_Inter+1)+1;

end

pilot_Indx(1,num_pilot)=num_data;

for j=0:num_pilot

Data_Indx(1,(1+j*pilot_Inter):(j+1)*pilot_Inter)=(2+j*(pilot_Inter+1)):((j+1)*(pilot_Inter+1));

end

Data_Indx=Data_Indx(1,1:num_symbol);

%----------------导频插入-------------------------------------

piloted_ofdm_syms=zeros(carrier_count,num_data);

piloted_ofdm_syms(:,Data_Indx)=reshape(Modulated_Sequence,carrier_count,num_symbol);

piloted_ofdm_syms(:,pilot_Indx)=repmat(pilot_symbols,1,num_pilot);

%----------------IFFT变换———————————————————

time_signal=sqrt(carrier_count)*ifft(piloted_ofdm_syms);

%----------------加循环前缀------------------------------------

add_cyclic_signal=[time_signal((carrier_count-Guard+1:carrier_count),:);time_signal];

Tx_data_trans=reshape(add_cyclic_signal,1,(carrier_count+Guard)*num_data);

%----------------信道处理--------------------------------------

Tx_signal_power=sum(abs(Tx_data_trans(:)).^2)/length(Tx_data_trans(:));

noise_var=Tx_signal_power/(10^(SNR(c1)/10));

Rx_data=awgn(Tx_data_trans,SNR(c1),'measured');

%----------------信号接收、去循环前缀、FFT变换-----------------

Rx_signal=reshape(Rx_data,(carrier_count+Guard),num_data);

Rx_signal_matrix=zeros(carrier_count,num_data);

Rx_signal_matrix=Rx_signal(Guard+1:end,:);

Rx_carriers=fft(Rx_signal_matrix)/sqrt(carrier_count);

%----------------导频和数据提取--------------------------------

Rx_pilot=Rx_carriers(:,pilot_Indx);

Rx_fre_data=Rx_carriers(:,Data_Indx);

%----------------导频位置信道响应LS估计------------------------

pilot_patt=repmat(pilot_symbols,1,num_pilot);

pilot_esti=Rx_pilot./pilot_patt;

%----------------LS估计的线性插值———————————————

int_len=pilot_Indx;

len=1:num_data;

for ii=1:carrier_count

channel_H_ls(ii,:)=interp1(int_len,pilot_esti(ii,1:(num_pilot)),len,'spline');

end

channel_H_data_ls=channel_H_ls(:,Data_Indx);

%----------------LS估计中发送数据的估计值----------------------

Tx_data_estimate_ls=Rx_fre_data.*conj(channel_H_data_ls)./(abs(channel_H_data_ls).^2);

%----------------DFT估计--------------------------------------

Tx_pilot_estimate_ifft=ifft(pilot_esti);

padding_zero=zeros(64,num_pilot);

Tx_pilot_estimate_ifft_padding_zero=[Tx_pilot_estimate_ifft;padding_zero];

Tx_pilot_estimate_dft=fft(Tx_pilot_estimate_ifft_padding_zero);

%----------------DFT估计的线性插值———————————————

int_len=pilot_Indx;

len=1:num_data;

for ii=1:carrier_count

channel_H_dft(ii,:)=interp1(int_len,Tx_pilot_estimate_dft(ii,1:(num_pilot)),len,'linear');

end

channel_H_data_dft=channel_H_dft(:,Data_Indx);

%----------------DFT估计中发送数据的估计值----------------------

Tx_data_estimate_dft=Rx_fre_data.*conj(channel_H_data_dft)./(abs(channel_H_data_dft).^2);
%----------------MMSE估计--------------------------------------
beta=1;
trms=1e-6;
t_max=8;
N=carrier_count;
Rhh=zeros(N,N);
for k=1:N
    for l=1:N
        Rhh(k,l)=(1-exp((-1)*t_max*((1/trms)+1i*2*pi*(k-l)/N)))./(trms*(1-exp((-1)*t_max/trms))*((1/trms)+1i*2*pi*(k-l)/N));
    end
end

Q=Rhh*(eye(N)/(Rhh+(beta/SNR(c1))*eye(N)));
Hmmse=Q*pilot_esti;

%----------------MMSE估计中发送数据的估计值----------------------
int_len=pilot_Indx;

len=1:num_data;

for ii=1:carrier_count

channel_H_mmse(ii,:)=interp1(int_len,Hmmse(ii,1:(num_pilot)),len,'linear');

end

channel_H_data_mmse=channel_H_mmse(:,Data_Indx);

Tx_data_estimate_mmse=Rx_fre_data.*conj(channel_H_data_mmse)./(abs(channel_H_data_mmse).^2);
%----------------DFT符号解调------------------------------------

demod_in_dft=Tx_data_estimate_dft(:).';

demod_out_dft=qamdemod(demod_in_dft,modulation_mode);

%----------------LS符号解调------------------------------------

demod_in_ls=Tx_data_estimate_ls(:).';

demod_out_ls=qamdemod(demod_in_ls,modulation_mode);


%----------------LS符号解调------------------------------------

demod_in_mmse=Tx_data_estimate_mmse(:).';

demod_out_mmse=qamdemod(demod_in_mmse,modulation_mode);
%----------------误码率的计算----------------------------------

for i=1:length(BitsTx)

if demod_out_dft(i)~=BitsTx(i)

num_bit_err_dft(c1,num1)=num_bit_err_dft(c1,num1)+1;

end

if demod_out_ls(i)~=BitsTx(i)

num_bit_err_ls(c1,num1)=num_bit_err_ls(c1,num1)+1;

end

if demod_out_mmse(i)~=BitsTx(i)

num_bit_err_mmse(c1,num1)=num_bit_err_mmse(c1,num1)+1;

end

end

end

end

BER_dft=mean(num_bit_err_dft.')/length(BitsTx);

BER_ls=mean(num_bit_err_ls.')/length(BitsTx);

BER_mmse=mean(num_bit_err_mmse.')/length(BitsTx);
%%%%%%%%%%%%%%%%%%%a%%%%%%%%主程序循环换结束%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

semilogy(SNR,BER_dft,'-go',SNR,BER_ls,'-k+',SNR,BER_mmse,'-b*');
legend('dft','ls','mmse');
title('OFDM系统的LS和DFT和MMSE信道估计');
xlabel('SNR'),ylabel('BER')