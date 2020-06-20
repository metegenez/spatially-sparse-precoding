classdef OptimalUnconstraint < handle

    
    properties
        channelSimulation
        SNR
        spectralEfficiency
    end
    
    methods
        function obj = OptimalUnconstraint(SNR,Channel)
            obj.channelSimulation = Channel;
            obj.SNR = SNR; %linear
            obj.launchObj();
        end
        
        function obj = launchObj(obj)         
            %%%% A. Optimal unconstrained precoding and combining
            H = obj.channelSimulation.channelMatrix;
            [~,~,V] = svd(H);
            fOpt = V(:,1:obj.channelSimulation.numberDataStreams); % Nt x Ns, optimal precoder
            wOpt = ((1/sqrt(obj.SNR))*(fOpt'* H'* H* fOpt+...
                obj.channelSimulation.numberDataStreams/obj.SNR * ...
                eye(obj.channelSimulation.numberDataStreams)) \ (fOpt' * H'))'; % Nr x Ns, optimal combiner
            rOpt = wOpt'*wOpt;
            obj.spectralEfficiency = abs(log2(det(eye(obj.channelSimulation.numberDataStreams)+...
                (obj.SNR/obj.channelSimulation.numberDataStreams)*...
                (rOpt\(wOpt'*H*fOpt*fOpt'*H'*wOpt)))));   
        
        end
        
        
    end
end

