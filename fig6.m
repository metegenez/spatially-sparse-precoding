clear all
snrValueDB = 0; 
SNR = 10.^(snrValueDB./10);
tryNumber = 300; 
%% Fig6 Parameters
parameters_Fig6 = containers.Map('KeyType','char','ValueType','any');
parameters_Fig6("numberTransmitAntennas") = 64; % Number of transmit antennas
parameters_Fig6("numberRecieveAntennas") = 16; % Number of receive antennas
parameters_Fig6("numberDataStreams") = 1; % Number of data streams
parameters_Fig6("numberRFChains") = 4; % Number of RF chains for precoding and combining
parameters_Fig6("numberCluster") = 8; % Number of clusters
parameters_Fig6("numberRayPerCluster") = 10; % Number of rays per cluster

%% Angle Spread vs Spectral Eff
angleSpreadValues = 0:1:15;
spectralEffOptimal = zeros(tryNumber,length(angleSpreadValues));
spectralEffHybrid = zeros(tryNumber,length(angleSpreadValues));
for a = 1:length(angleSpreadValues)
    parameters_Fig6("angularSpread") = angleSpreadValues(a);
    for i = 1:tryNumber
        channel = ChannelGeneration(parameters_Fig6);
        tempObj = OptimalUnconstraint(SNR,channel);
        spectralEffOptimal(i,a) = tempObj.spectralEfficiency;
        tempObj = HybridSparsePrecoding(SNR,channel);
        spectralEffHybrid(i,a) = tempObj.spectralEfficiency;
    end
end
% Averaging Tries
spectralEffOptimalSNR = mean(spectralEffOptimal,1); 
spectralEffHybridSNR = mean(spectralEffHybrid,1);
% figure(); 
hold on
color = rand(1,3);
l1 = plot(angleSpreadValues,spectralEffOptimalSNR,'-s','Color',color,'LineWidth',2.0,'MarkerSize',8.0, 'DisplayName', sprintf("Optimal Uns. %dx%d, N_s=%d, N_{RF}=%d",parameters_Fig6("numberTransmitAntennas"),parameters_Fig6("numberRecieveAntennas"),parameters_Fig6("numberDataStreams"),parameters_Fig6("numberRFChains") ));
l2 = plot(angleSpreadValues,spectralEffHybridSNR,'-o','Color',color,'LineWidth',2.0,'MarkerSize',8.0,'DisplayName', sprintf("Hybrid Comb. %dx%d, N_s=%d, N_{RF}=%d",parameters_Fig6("numberTransmitAntennas"),parameters_Fig6("numberRecieveAntennas"),parameters_Fig6("numberDataStreams"),parameters_Fig6("numberRFChains") ));
%% v2
parameters_Fig6("numberTransmitAntennas") = 256; % Number of transmit antennas
parameters_Fig6("numberRecieveAntennas") = 64; % Number of receive antennas
spectralEffOptimal = zeros(tryNumber,length(angleSpreadValues));
spectralEffHybrid = zeros(tryNumber,length(angleSpreadValues));
for a = 1:length(angleSpreadValues)
    parameters_Fig6("angularSpread") = angleSpreadValues(a);
    for i = 1:tryNumber
        channel = ChannelGeneration(parameters_Fig6);
        tempObj = OptimalUnconstraint(SNR,channel);
        spectralEffOptimal(i,a) = tempObj.spectralEfficiency;
        tempObj = HybridSparsePrecoding(SNR,channel);
        spectralEffHybrid(i,a) = tempObj.spectralEfficiency;
    end
end
% Averaging Tries
spectralEffOptimalSNR = mean(spectralEffOptimal,1); 
spectralEffHybridSNR = mean(spectralEffHybrid,1);
hold on
color = rand(1,3);
l3 = plot(angleSpreadValues,spectralEffOptimalSNR,'-s','Color',color,'LineWidth',2.0,'MarkerSize',8.0, 'DisplayName', sprintf("Optimal Uns. %dx%d, N_s=%d, N_{RF}=%d",parameters_Fig6("numberTransmitAntennas"),parameters_Fig6("numberRecieveAntennas"),parameters_Fig6("numberDataStreams"),parameters_Fig6("numberRFChains") ));
l4 = plot(angleSpreadValues,spectralEffHybridSNR,'-o','Color',color,'LineWidth',2.0,'MarkerSize',8.0,'DisplayName', sprintf("Hybrid Comb. %dx%d, N_s=%d, N_{RF}=%d",parameters_Fig6("numberTransmitAntennas"),parameters_Fig6("numberRecieveAntennas"),parameters_Fig6("numberDataStreams"),parameters_Fig6("numberRFChains") ));
%% v3
parameters_Fig6("numberTransmitAntennas") = 64; % Number of transmit antennas
parameters_Fig6("numberRecieveAntennas") = 16; % Number of receive antennas
parameters_Fig6("numberDataStreams") = 2; % Number of data streams
spectralEffOptimal = zeros(tryNumber,length(angleSpreadValues));
spectralEffHybrid = zeros(tryNumber,length(angleSpreadValues));
for a = 1:length(angleSpreadValues)
    parameters_Fig6("angularSpread") = angleSpreadValues(a);
    for i = 1:tryNumber
        channel = ChannelGeneration(parameters_Fig6);
        tempObj = OptimalUnconstraint(SNR,channel);
        spectralEffOptimal(i,a) = tempObj.spectralEfficiency;
        tempObj = HybridSparsePrecoding(SNR,channel);
        spectralEffHybrid(i,a) = tempObj.spectralEfficiency;
    end
end
% Averaging Tries
spectralEffOptimalSNR = mean(spectralEffOptimal,1); 
spectralEffHybridSNR = mean(spectralEffHybrid,1);
hold on
color = rand(1,3);
l5 = plot(angleSpreadValues,spectralEffOptimalSNR,'-s','Color',color,'LineWidth',2.0,'MarkerSize',8.0, 'DisplayName', sprintf("Optimal Uns. %dx%d, N_s=%d, N_{RF}=%d",parameters_Fig6("numberTransmitAntennas"),parameters_Fig6("numberRecieveAntennas"),parameters_Fig6("numberDataStreams"),parameters_Fig6("numberRFChains") ));
l6 = plot(angleSpreadValues,spectralEffHybridSNR,'-o','Color',color,'LineWidth',2.0,'MarkerSize',8.0, 'DisplayName', sprintf("Hybrid Comb. %dx%d, N_s=%d, N_{RF}=%d",parameters_Fig6("numberTransmitAntennas"),parameters_Fig6("numberRecieveAntennas"),parameters_Fig6("numberDataStreams"),parameters_Fig6("numberRFChains") ));
%% v4
parameters_Fig6("numberTransmitAntennas") = 64; % Number of transmit antennas
parameters_Fig6("numberRecieveAntennas") = 16; % Number of receive antennas
parameters_Fig6("numberDataStreams") = 2; % Number of data streams
parameters_Fig6("numberRFChains") = 6;
spectralEffOptimal = zeros(tryNumber,length(angleSpreadValues));
spectralEffHybrid = zeros(tryNumber,length(angleSpreadValues));
for a = 1:length(angleSpreadValues)
    parameters_Fig6("angularSpread") = angleSpreadValues(a);
    for i = 1:tryNumber
        channel = ChannelGeneration(parameters_Fig6);
        tempObj = OptimalUnconstraint(SNR,channel);
        spectralEffOptimal(i,a) = tempObj.spectralEfficiency;
        tempObj = HybridSparsePrecoding(SNR,channel);
        spectralEffHybrid(i,a) = tempObj.spectralEfficiency;
    end
end
% Averaging Tries
spectralEffOptimalSNR = mean(spectralEffOptimal,1); 
spectralEffHybridSNR = mean(spectralEffHybrid,1);
hold on
color = rand(1,3);
l7 = plot(angleSpreadValues,spectralEffOptimalSNR,'-s','Color',color,'LineWidth',2.0,'MarkerSize',8.0, 'DisplayName', sprintf("Optimal Uns. %dx%d, N_s=%d, N_{RF}=%d",parameters_Fig6("numberTransmitAntennas"),parameters_Fig6("numberRecieveAntennas"),parameters_Fig6("numberDataStreams"),parameters_Fig6("numberRFChains") ));
l8 = plot(angleSpreadValues,spectralEffHybridSNR,'-o','Color',color,'LineWidth',2.0,'MarkerSize',8.0,'DisplayName', sprintf("Hybrid Comb. %dx%d, N_s=%d, N_{RF}=%d",parameters_Fig6("numberTransmitAntennas"),parameters_Fig6("numberRecieveAntennas"),parameters_Fig6("numberDataStreams"),parameters_Fig6("numberRFChains") ));
legend
xlabel("Angle Spread (degrees)")
ylabel("Spectral Efficiency(bits/s/Hz)")