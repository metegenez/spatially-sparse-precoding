classdef ChannelGeneration < handle
    % Simulation for environment and optimiation types.
    
    properties
        numberTransmitAntennas  % Number of transmit antennas
        numberRecieveAntennas % Number of receive antennas
        numberDataStreams % Number of data streams
        numberRFChains% Number of RF chains for precoding and combining 
        numberCluster % Number of clusters
        numberRayPerCluster % Number of rays per cluster
        angularSpread % Angular spread of 7.5 degree
    end
    properties (Access=public)
        channelMatrix
        arrayResponseTx
        arrayResponseRx
        alpha
    end
    
    methods
        function obj = ChannelGeneration(parameters)
            obj.numberTransmitAntennas = parameters("numberTransmitAntennas");
            obj.numberRecieveAntennas = parameters("numberRecieveAntennas");
            obj.numberDataStreams = parameters("numberDataStreams");
            obj.numberRFChains = parameters("numberRFChains");
            obj.numberCluster = parameters("numberCluster");
            obj.numberRayPerCluster = parameters("numberRayPerCluster");
            obj.angularSpread = parameters("angularSpread");
            obj.launch();
        end
        
        function obj = launch(obj)
            minAOD = -30;
            maxAOD = 30;
            minEOD = 80;
            maxEOD = 100; 
            clusterAOD = rand(obj.numberCluster,1)*(maxAOD-minAOD)+minAOD; % [-30,30] degree
            clusterEOD = rand(obj.numberCluster,1)*(maxEOD-minEOD)+minEOD; % [80,100] degree
            % Receiver
            clusterAOA = (rand(obj.numberCluster,1)-0.5)*360; % [-179,180] degree
            clusterEOA = rand(obj.numberCluster,1)*180; % [0,179] degree
            b = obj.angularSpread/sqrt(2); % Scaling parameter, degree          
            rndm = rand(obj.numberRayPerCluster*obj.numberCluster,1)-0.5;          
            % Dimension of AOD, EOD, AOA, EOA : (n_ray x n_cluster) x 1
            rayAOD = repelem(clusterAOD,obj.numberRayPerCluster,1)-b*sign(rndm).*log(1-2.*abs(rndm));
            rayEOD = repelem(clusterEOD,obj.numberRayPerCluster,1)-b*sign(rndm).*log(1-2.*abs(rndm));
            rayAOA = repelem(clusterAOA,obj.numberRayPerCluster,1)-b*sign(rndm).*log(1-2.*abs(rndm));
            rayEOA = repelem(clusterEOA,obj.numberRayPerCluster,1)-b*sign(rndm).*log(1-2.*abs(rndm));
            
            [txPosition, rxPosition] = obj.obtainPositionVectors();
            
            SphericalUnitVecTx = obj.getSphericalUnitVector(rayEOD,rayAOD); % 3 x n_ray*n_cluster
            SphericalUnitVecRx = obj.getSphericalUnitVector(rayEOA,rayAOA); % 3 x n_ray*n_cluster
            
            obj.arrayResponseTx = (1/sqrt(obj.numberTransmitAntennas))*exp(1i*pi*txPosition.'*SphericalUnitVecTx); % n_t x n_ray*n_cluster
            obj.arrayResponseRx = (1/sqrt(obj.numberRecieveAntennas))*exp(1i*pi*rxPosition.'*SphericalUnitVecRx); % n_r x n_ray*n_cluster
            
            obj.generateComplexGainPath();
            obj.generateChannelMatrix();
        end
        
        function [txPosition, rxPosition] = obtainPositionVectors(obj)
            %% Obtain antenna element position vectors (normalized by half of the wavelength)
           % Transmitter
           tempNtH = sqrt(obj.numberTransmitAntennas);
           tempNtV = sqrt(obj.numberTransmitAntennas);
           
           X_Tx = zeros(1,obj.numberTransmitAntennas);
           [Y_Tx,Z_Tx] = meshgrid(0:tempNtH-1,0:tempNtV-1);
           txPosition = [X_Tx;Y_Tx(:).';Z_Tx(:).']; % 3 x n_t
           
           % Receiver
           tempNrH = sqrt(obj.numberRecieveAntennas);
           tempNrV = sqrt(obj.numberRecieveAntennas);
           
           xRx = zeros(1,obj.numberRecieveAntennas);
           [yRx,zRx] = meshgrid(0:tempNrH-1,0:tempNrV-1);
           rxPosition = [xRx;yRx(:).';zRx(:).']; % 3 x n_r
        end
        
        function obj = generateComplexGainPath(obj)
            obj.alpha = sqrt(1/2)*(randn(1,obj.numberRayPerCluster*obj.numberCluster)+1i*randn(1,obj.numberRayPerCluster*obj.numberCluster));
        end
        function sphericalUnitVector = getSphericalUnitVector(obj,theta,phi)
            sphericalUnitVector = [(sind(theta).*cosd(phi)).';...
                (sind(theta).*sind(phi)).';...
                (cosd(theta)).'];
        end
        function obj = generateChannelMatrix(obj)
        obj.channelMatrix = sqrt((obj.numberTransmitAntennas*obj.numberRecieveAntennas)...
            /(obj.numberRayPerCluster*obj.numberCluster))*obj.arrayResponseRx*diag(obj.alpha)*obj.arrayResponseTx';

       
        end
    end
end

