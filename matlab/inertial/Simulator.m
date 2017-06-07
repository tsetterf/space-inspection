classdef Simulator < handle
    %SIMULATOR Performs a SPHERES simulation and then outputs all of the files (except the images)
    % necessary to run TP_1511 in Linux into the GSdata folder
    
    properties
        type;               % simulation type, either 'TA' for triaxial, or 'AS1' for axis-symm 1
        simDir;             % path to SPHERES simulation directory
        cfg;                % SPHERES simulation configuration
        res;                % SPHERES simulation results
        blenderSimData;     % the blender simulation data for rendering
        imageTimeTags;      % the image numbers and times
        navStateData;       % the navigation state data with star tracker noise
        globMetData;        % global metrology data
        thrusterFiringTimesData;    % the on and off times of the inspector thrusters
        baseDir;            % the base directory where the TS00AS1 and TS00TA folders will be written
        resultsDir;         % the path the TS00AS1 or TS00TA folders
        gsDataDir;          % path to GSdata folder
        imageDir;           % path to images folder
        tStart = 0;         % SPHERES time at which to start the data [ms]
        tStartBlend = 329000; % start blender data as required
        tEnd;               % SPHERES time at which the data ends [ms]
        t;                  % SPHERES times simulated region [ms]
        nT;                 % number of timesteps
        inspState;          % true state for inspector
        tgtState;           % true state for target
    end
    
    methods
        %% Constructor
        function this = Simulator(type)
            
            this.type = type;
            if strcmp(this.type,'TA')
                error('TA case not setup yet');
            elseif strcmp(this.type,'AS1')
                this.simDir = 'D:\PhD\svns\TestProjects\VERTIGO_Wifi\MIT_P3420\Simulation';
                this.baseDir = 'D:\PhD\Thesis\Matlab\v0.8\';
                this.resultsDir = [this.baseDir 'TS00AS1\'];
                this.gsDataDir = [this.resultsDir 'GSdata\'];
                this.imageDir = [this.resultsDir 'OpticsMount_0_Images\'];
            end
                                                       
        end
        %% Simulate the scenario
        function simulate(this)
                                               
            % Go to the simulation directory, then configure and build the SPHERES sim
            currDir = pwd;
            cd(this.simDir);
            this.cfg = configureSim();
            BuildSimulation(this.cfg);
            
            % Create a diary to log the output imu data, then run the simulator
            delete([this.gsDataDir 'imu_data.txt']);
            diary([this.gsDataDir 'imu_data.txt']);
            this.res = RunSimulation(this.cfg);
            diary off;
            
            cd(currDir);
            
            % Get the end time from the true state times
            this.tEnd = min(this.res.trueStateT{1}(end),this.res.trueStateT{2}(end));
            
            % Get the vector of true times            
            this.t = this.res.trueStateT{1}(this.res.trueStateT{1} > this.tStart ...
                                                & this.res.trueStateT{1} < this.tEnd);
            this.nT = length(this.t);
            
            % Get the true states of the inspector and their times
            this.inspState  = this.res.trueStates{1}(:, this.res.trueStateT{1} > this.tStart ...
                                                    & this.res.trueStateT{1} < this.tEnd);                                                
            this.tgtState   = this.res.trueStates{2}(:, this.res.trueStateT{2} > this.tStart ...
                                                    & this.res.trueStateT{2} < this.tEnd);
            
            
            % Create the required data
            this.createData();
            
            disp('Done!');
                        
        end
        %% Create the required data files
        function createData(this)
                        
            % Delete old data            
            delete([this.imageDir 'imageTimetags.csv']);
            delete([this.gsDataDir 'navState_data.txt']);
            delete([this.gsDataDir 'globMet_data.txt']);
            delete([this.gsDataDir 'thrusterFiringTimes_data.txt']);
            delete([this.gsDataDir 'blenderSim_data.txt']);
            delete([this.resultsDir 'spheresData.mat']);
            delete([this.resultsDir 'simulatorObject.mat']);
                        
            this.getImageTimeTags();
            this.getNavStateData();
            this.getglobMetData();
            this.getThrusterFiringTimesData();
            this.getBlenderSimData();
            
            % Store a data variable to spheresData.mat
            cfg = this.res;
            save([this.resultsDir 'spheresData.mat'],'cfg');
            
            % Store the entire SPHERES simulator object, including the ground truth in this.res
            s = this;
            save([this.resultsDir 'simulatorObject.mat'],'s');
            
        end
        %% Get the timetags for all of the images
        function getImageTimeTags(this)
            
            % Clear image time tags
            this.imageTimeTags = [];
            
            % Get image time tags
            this.imageTimeTags(:,1) = 0:this.nT-1;
            this.imageTimeTags(:,2) = this.t;
            
            % Write data to file
            dlmwrite([this.imageDir 'imageTimetags.csv'],this.imageTimeTags,'precision',7);
            
        end
        %% Get the navigation state data with added star tracker noise
        function getNavStateData(this)
            
            % Clear nav state data
            this.navStateData = [];
            
            % Gyro and acceleration bias
            bg = [0 0 0]'; ba = [0 0 0]';
            
            % Star tracker noise. Assume 10x worse that 1-sigma 20 arcsecond resolution
            sigmaStarTracker = 10 * (20*pi/648000);
            covRBtoW = diag([sigmaStarTracker^2 sigmaStarTracker^2 sigmaStarTracker^2]);
            
            % Loop through states and get nav state data
            for i = 1:this.nT
                
                % Get inspector rotation, translation, and velocity
                RBtoW = quat2rot(this.inspState(7:10,i))';
                tWtoB_W = this.inspState(1:3,i);
                vB_W = this.inspState(4:6,i);
                
                % Add noise to the rotation
                dTheta = mvnrnd(zeros(3,1),covRBtoW);
                RBtoW = RBtoW*Exp(dTheta);
                
                % Combine results into a pose matrix
                PWtoB = [RBtoW tWtoB_W; zeros(1,3) 1];
                
                % Store the nav state data as an array
                this.navStateData(i,:) = [ this.t(i) reshape(PWtoB',1,16) vB_W' bg' ba' ];
                
            end
                        
            % Write data to file
            dlmwrite([this.gsDataDir 'navState_data.txt'],this.navStateData,'precision',7);
            
        end
        %% Get the global metrology data closest to the true times
        function getglobMetData(this)
            
            % Clear the global metrology data
            this.globMetData = zeros(this.nT,15);
            
            % Looop through the timesteps getting the closest global metrology
            for i = 1:this.nT
                
                % Get index of closest global metrology and the time and state
                [~,iGm] = min(abs(double(this.res.data{1}.BackTel.StdTime) - this.t(i)));
                tGm     = double(this.res.data{1}.BackTel.StdTime(iGm));
                stateGm = double(this.res.data{1}.BackTel.StdState(:,iGm));
                                                   
                % Store the global metrology data
                this.globMetData(i,:) = [ this.t(i) tGm stateGm' ];
                
            end
            
            % Write data to file
            dlmwrite([this.gsDataDir 'globMet_data.txt'],this.globMetData,'precision',7);            
            
        end
        %% Get the thruster firing times data from the SPHERES debug packets
        function getThrusterFiringTimesData(this) 
            

            % Get debug data
            debugData = this.res.data{1}.DebugUnsignedVal;
            
            % Clear the thruster firing times data
            this.thrusterFiringTimesData = [];
            
            % Pre-allocate space
            thrTimes  = zeros(size(debugData,2),26);
            
            % No time difference between SPHERES time and VERTIGO in simulation
            timeDiff = 0;
                        
            % Get thruster data
            indThrTimes = 1;
            for i = 1:size(debugData,2) 

                % Exclude non-thruster related packets
                if sum(debugData(3:end,i)) ~= 0

                    % Time of the control period
                    ctrlTime = debugData(1,i)*100;
                    thrTimes(indThrTimes,1) = ctrlTime - timeDiff;
                    thrTimes(indThrTimes,2) = ctrlTime;

                    % Get activated thrusters (==0 for thrusters 0-5, ==1 for thrusters 6-11)
                    actThrust = dec2bin(debugData(14,i),6);

                    % Populate the thruster times
                    for j = 1:6
                        if str2double(actThrust(7-j)) == 0 ...
                           && (debugData(2*j,i) ~= debugData(2*j+1,i))
                                thrTimes(indThrTimes,   (j+2)) = debugData(2*j  ,i);
                                thrTimes(indThrTimes,12+(j+2)) = debugData(2*j+1,i);
                        elseif (debugData(2*j,i) ~= debugData(2*j+1,i))        
                            thrTimes(indThrTimes,   (j+2)+6) = debugData(2*j  ,i);
                            thrTimes(indThrTimes,12+(j+2)+6) = debugData(2*j+1,i);
                        end
                    end
                    indThrTimes = indThrTimes + 1;        
                end
            end

            % Remove rows containing all zeros
            thrTimes = thrTimes(any(thrTimes~=0,2),:);

            % Store thruster times
            this.thrusterFiringTimesData = thrTimes;
            
            % Write to file
            dlmwrite([this.gsDataDir 'thrusterFiringTimes_data.txt'],thrTimes,'precision',7);

        end
        %% Get the blender simulation data, where the scale of translations is increased by 100x
        function getBlenderSimData(this)
            
            % Clear blender sim data
            this.blenderSimData = [];
            
            % Get data in the format:
            %   x_lc,y_lc,z_lc,x_rc,y_rc,z_rc,qv1_c,qv2_c,qv3_c,qw_c,qv1_e,qv2_e,qv3_e,qw_e
            % 		where l=left, r=right, c=camera, and e=endurance spacecraft which is 
            %		always positioned at the origin.
                             
            % Rotation from camera to body (vertigo)
            RCvtoB = [0 0 1; -1 0 0; 0 -1 0];
            
            % Rotation from blender camera to vertigo camera
            RCbtoCv = [1 0 0; 0 -1 0; 0 0 -1];
            
            % Rotation from blender camera to body
            RCbtoB = RCvtoB * RCbtoCv;
            
            % Translation from body to left and right camera frames
            tBtoCl_B = [23.4-2.5  4.5 3.0]';
            tBtoCr_B = [23.4-2.5 -4.5 3.0]';
            
            % Major axis of endurance is the y axis, so need a transform from spheres to endurance
            REtoBt = [0 1 0; 1 0 0; 0 0 -1];
                             
            % Loop through the inspector data to properly position the left and right cameras
            this.blenderSimData = zeros(this.nT,14);
            for i = 1:this.nT
               
                if this.t(i) < this.tStartBlend
                    continue;
                end
                
                % Get inspector rotation and camera rotation
                RBtoW  = quat2rot(this.inspState(7:10,i))';
                RCbtoW = RBtoW * RCbtoB;
                qWtoCb = rot2quat(RCbtoW');
                                
                % Get left and right camera position
                tWtoB_W  = this.inspState(1:3,i)*100;
                tWtoCl_W = tWtoB_W + RBtoW * tBtoCl_B;
                tWtoCr_W = tWtoB_W + RBtoW * tBtoCr_B;
                
                % Get body rotation of target and then of endurance
                RBttoW = quat2rot(this.tgtState(7:10,i))';
                REtoW  = RBttoW * REtoBt;
                qWtoE  = rot2quat(REtoW');
                                            
                % Record the data
                this.blenderSimData(i,:) = [tWtoCl_W' tWtoCr_W' qWtoCb' qWtoE'];  
                  
                % For testing
%                 disp(['bpy.data.objects[''Camera''].rotation_quaternion = (' num2str(qWtoCb(4)) ',' ...
%                         num2str(qWtoCb(1)) ',' num2str(qWtoCb(2)) ',' num2str(qWtoCb(3)) ')']);                  
%                 disp(['bpy.data.objects[''Camera''].location = (' ...
%                         num2str(tWtoCl_W(1)) ',' num2str(tWtoCl_W(2)) ',' num2str(tWtoCl_W(3)) ')']);            
%                 disp(['bpy.data.objects[''0parent''].rotation_quaternion = (' num2str(qWtoE(4)) ',' ...
%                         num2str(qWtoE(1)) ',' num2str(qWtoE(2)) ',' num2str(qWtoE(3)) ')']);
                                                
            end
                        
            % Write data to file
            dlmwrite([this.gsDataDir 'blenderSim_data.txt'],this.blenderSimData,'precision',7);
                        
        end
        
    end
    
end

