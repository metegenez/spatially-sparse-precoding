classdef BeamSteering < handle
    
    properties
        channelSimulation
        SNR
        spectralEfficiency
    end
    
    methods
        function obj = BeamSteering(SNR,Channel)
            obj.channelSimulation = Channel;
            obj.SNR = SNR; %linear
            obj.launch();
        end
        
        function obj = launch(obj)
            %%%% A. Optimal unconstrained precoding and combining
            [indexArrResponseTx,indexArrResponseRx] = obj.findSteeringVector();
            H = obj.channelSimulation.channelMatrix;
            F_BS = [];
            W_BS = [];
            
            for n = 1:obj.channelSimulation.numberDataStreams
                F_BS = [F_BS obj.channelSimulation.arrayResponseTx(:,indexArrResponseTx)];
                W_BS = [W_BS obj.channelSimulation.arrayResponseRx(:,indexArrResponseRx)];
            end
            RnBS = W_BS'*W_BS;
            obj.spectralEfficiency = abs(log2(det(eye(obj.channelSimulation.numberDataStreams)...
                + (obj.SNR/obj.channelSimulation.numberDataStreams)...
                * (RnBS\(W_BS'*H*F_BS*F_BS'*H'*W_BS)))));
               
        end
        function [indexArrResponseTx,indexArrResponseRx] = findSteeringVector(obj)
            % This is the function to obtain the beam steering vectors at the transmitter and receiver
            % by finding the array response vectors corresponding to the
            % largest effective channel gain. This algorithm is brute-force search. 
            
            numberPath = size(obj.channelSimulation.arrayResponseTx,2);
            
            effectiveGain = zeros(numberPath);
            H = obj.channelSimulation.channelMatrix;
            if obj.channelSimulation.numberDataStreams == 1
                for t = 1:numberPath
                    for r = 1:numberPath
                        % Compute the absolute effective channel gain.
                        effectiveGain(t,r) = abs(obj.channelSimulation.arrayResponseRx(:,r)'*H*obj.channelSimulation.arrayResponseTx(:,t));
                    end              
                end
            else               
                error('Check Transmit Data Streams Number');   
            end
            
            [indexArrResponseTx,indexArrResponseRx] = find(effectiveGain == max(max(effectiveGain)));

            
        end
    end
end

