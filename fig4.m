clear all
snrValuesDB = -40:5:0; %in Figures
snrValues = 10.^(snrValuesDB./10);
tryNumber = 100; 
%% Fig4 Parameters
parameters_Fig4 = containers.Map('KeyType','char','ValueType','any');
parameters_Fig4("numberTransmitAntennas") = 256; % Number of transmit antennas
parameters_Fig4("numberRecieveAntennas") = 64; % Number of receive antennas
parameters_Fig4("numberDataStreams") = 1; % Number of data streams
parameters_Fig4("numberRFChains") = 6; % Number of RF chains for precoding and combining
parameters_Fig4("numberCluster") = 8; % Number of clusters
parameters_Fig4("numberRayPerCluster") = 10; % Number of rays per cluster
parameters_Fig4("angularSpread") = 7.5; % Angular spread of 7.5 degree
spectralEffOptimal = zeros(tryNumber,length(snrValues));
spectralEffHybrid = zeros(tryNumber,length(snrValues));
spectralEffBeam = zeros(tryNumber,length(snrValues));
for s = 1:length(snrValues)
    SNR = snrValues(s);
    for i = 1:tryNumber
        channel = ChannelGeneration(parameters_Fig4);
        tempObj = OptimalUnconstraint(SNR,channel);
        spectralEffOptimal(i,s) = tempObj.spectralEfficiency;
        tempObj = HybridSparsePrecoding(SNR,channel);
        spectralEffHybrid(i,s) = tempObj.spectralEfficiency;
        tempObj = BeamSteering(SNR,channel);
        spectralEffBeam(i,s) = tempObj.spectralEfficiency;
        
    end
end
% Averaging Tries
spectralEffOptimalSNR = mean(spectralEffOptimal,1); 
spectralEffHybridSNR = mean(spectralEffHybrid,1);
spectralEffBeamSNR = mean(spectralEffBeam,1); 
% figure(); 
hold on
l1 = plot(snrValuesDB,spectralEffOptimalSNR,'-s','Color',[0 0.5 0],'LineWidth',2.0,'MarkerSize',8.0);
l2 = plot(snrValuesDB,spectralEffHybridSNR,'-o','Color',[0 0.45 0.74],'LineWidth',2.0,'MarkerSize',8.0);hold on;
l3 = plot(snrValuesDB,spectralEffBeamSNR,'-d','Color',[0.85 0.33 0.1],'LineWidth',2.0,'MarkerSize',8.0);hold on;
%% Data stream
parameters_Fig4("numberDataStreams") = 2; % Number of data streams
spectralEffOptimal = zeros(tryNumber,length(snrValues));
spectralEffHybrid = zeros(tryNumber,length(snrValues));
for s = 1:length(snrValues)
    SNR = snrValues(s);
    for i = 1:tryNumber
        channel = ChannelGeneration(parameters_Fig4);
        tempObj = OptimalUnconstraint(SNR,channel);
        spectralEffOptimal(i,s) = tempObj.spectralEfficiency;
        tempObj = HybridSparsePrecoding(SNR,channel);
        spectralEffHybrid(i,s) = tempObj.spectralEfficiency;

        
    end
end
% Averaging Tries
spectralEffOptimalSNR = mean(spectralEffOptimal,1); 
spectralEffHybridSNR = mean(spectralEffHybrid,1);
spectralEffBeamSNR = mean(spectralEffBeam,1); 
l4 = plot(snrValuesDB,spectralEffOptimalSNR,'-s','LineWidth',2.0,'MarkerSize',8.0);
l5 = plot(snrValuesDB,spectralEffHybridSNR,'-o','LineWidth',2.0,'MarkerSize',8.0);hold on;
legend([l1 l2 l3 l4 l5],'Optimal unconstrained precoding N_s = 1','Hybrid precoding and combining N_s = 1','Beam steering N_s = 1','Optimal unconstrained precoding N_s = 2','Hybrid precoding and combining N_s = 2','Location','northwest','FontSize', 15);
xlabel("SNR (dB)",'FontSize', 20)
ylabel("Spectral Efficiency(bits/s/Hz)",'FontSize', 20)
