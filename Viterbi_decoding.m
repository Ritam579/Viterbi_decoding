% Convolutional Code hard decoding

N = 10^6;
Eb_N0_dB = [0:1:10];
Ec_N0_dB = Eb_N0_dB - 10*log10(2);
ref_bits_h = [0 0 ; 0 1; 1 0; 1 1];
ref_tab = [0 0 0 0 ;0 0 0 0 ;1 1 0 0 ;0 0 1 1];

for yy = 1:length(Eb_N0_dB)
    x = rand(1,N) > 0.5;
    conv1 = mod(conv(x,[111]),2);
    conv2 = mod(conv(x,[101]),2);
    conv_ap = [conv1;conv2];
    conv_ap = conv_ap(:).'; % Convert into a row variable

    s = 2*conv_ap - 1; %BPSK Modulator

    n = 1/sqrt(2)*[randn(size(conv_ap))+1j*randn(size(conv_ap))];%Noise

    y = s + 10^(-Ec_N0_dB(yy)/20)*n; %AWGN channel

    %Receiver Demodulation
    dem_h = real(y) > 0; %thresholding

    %Hard Decision Decoder

    %Path Matrix and Selected Path
    pm_h = zeros(4,1); %hard path metric
    sp_h = zeros(4,length(y)/2); %Selected Path for hard decoding

    %Decision by taking two bits together
    for b_x = 1:length(y)/2;
        dem_h1 = dem_h(2*b_x-1:2*b_x);

        %Hamming distance
        comb_h = repmat(dem_h1,4,1); %repeat h1 4 times
        d_hamming = sum(xor(comb_h,ref_bits_h),2);

        if (b_x == 1) || (b_x == 2)

            %Path metric and selection path for state 0
            pm_hard1 = pm_h(1,1) + d_hamming(1);
            pm_h_n(1,1) = pm_hard1;
            sp_h1(1,1) = 1;

            %Path metric and selection path for state 1
            pm_hard1 = pm_h(3,1) + d_hamming(3);
            pm_h_n(2,1) = pm_hard1;
            sp_h1(2,1) = 1;

            %Path metric and selection path for state 2
            pm_hard1 = pm_h(1,1) + d_hamming(4);
            pm_h_n(3,1) = pm_hard1;
            sp_h1(3,1) = 1;

            %Path metric and selection path for state 3
            pm_hard1 = pm_h(3,1) + d_hamming(2);
            pm_h_n(4,1) = pm_hard1;
            sp_h1(4,1) = 3;

        else
            %path metric and selection path for state 0
            pm_hard1 = pm_h(1,1) + d_hamming(1);
            pm_hard2 = pm_h(2,1) + d_hamming(4);
            [pm_h_n(1,1) idx] = min([pm_hard1,pm_hard2]);
            sp_h1(1,1) = idx;

            %path metric and selection path for state 1
            pm_hard1 = pm_h(3,1) + d_hamming(3);
            pm_hard2 = pm_h(4,1) + d_hamming(2);
            [pm_h_n(2,1) idx] = min([pm_hard1,pm_hard2]);
            sp_h1(2,1) = idx+2;


            %path metric and selection path for state 2
            pm_hard1 = pm_h(1,1) + d_hamming(4);
            pm_hard2 = pm_h(2,1) + d_hamming(1);
            [pm_h_n(3,1) idx] = min([pm_hard1,pm_hard2]);
            sp_h1(3,1) = idx;


            %path metric and selection path for state 3
            pm_hard1 = pm_h(3,1) + d_hamming(2);
            pm_hard2 = pm_h(4,1) + d_hamming(3);
            [pm_h_n(4,1) idx] = min([pm_hard1,pm_hard2]);
            sp_h1(4,1) = idx+2;

        end

        pm_h = pm_h_n;
        sp_h(:,b_x) = sp_h1;
    end

    %Trace back unit
    state_cur_h = 1;
    X_hat_hard = zeros(1,length(y)/2);
    for jj = length(y)/2:-1:1
        state_prev_h = sp_h(state_cur_h,jj);
        X_hat_hard(jj) = ref_tab(state_cur_h,state_prev_h);
        state_cur_h = state_prev_h;
    end

    %Counting the errors
    err_h(yy) = size(find([x - X_hat_hard(1:N)]),2);
end

ber_h = err_h/N;
theoryBer = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10)));
semilogy(Eb_N0_dB,theoryBer,'bd-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,ber_h,'mp-','LineWidth',2);
axis([0 10 10^(-5) 0.5]);
grid on
legend('Theory-Uncoded','Simulation-Hard-Viterbi');





