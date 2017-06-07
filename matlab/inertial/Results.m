classdef Results < handle
    %RESULTSPrepares the results of TP_1511_DynamicPoseSLAM for use (uses vertigo time)
    
    properties
        % Poses vio = visual inertial odom, gm = global metrology, isam = iSAM2, bsam = batch SAM
        PWtoBvio    = [];       % world to inspector
        PWtoBgm     = [];       % world to inspector
        PWtoBtgm    = [];       % world to target
        PWtoGgm     = [];       % world to geometric
        PGtoBgm     = [];       % geometric to inspector
        PWtoBbsam   = [];       % world to inspector (GTSAM)
        PWtoGbsam   = [];       % world to geometric (GTSAM)
        PGtoBbsam   = [];       % geometric to inspector (GTSAM)
        PWtoBisam   = [];       % world to inspector (GTSAM)
        PWtoGisam   = [];       % world to geometric (GTSAM)
        PGtoBisam   = [];       % geometric to inspector (GTSAM)
        PWtoBtrue   = [];       % world to inspector (simulated)
        PWtoBttrue  = [];       % world to target (from truth if simulated, PA corrected if not)
        % Translation from geometric frame and world frame to CM
        tGtoBt_Ggm    = [];
        tGtoBt_Gbsam  = [];
        tGtoBt_Gisam  = [];
        tGtoBt_Gtrue  = [];
        tWtoBt_Wbsam  = [];
        tWtoBt_Wisam  = [];
        tWtoBt_Wtrue  = [];
        % Rotations and their covariances
        RBtoGisam     = [];
        RBtoWisam     = [];
        % Velocities
        vB_Wvio     = [];       % velocity of inspector
        vB_Wgm      = [];       % velocity of inspector
        vB_Wbsam    = [];       % velocity of inspector
        vB_Wisam    = [];       % velocity of inspector
        vB_Wtrue    = [];       % velocity of inspector
        % Covariances
        covPWtoBisam   = [];
        covPGtoBisam   = [];
        covtGtoBt_Gisam = [];
        covRBtoGisam   = [];
        covRBtoWisam   = [];
        covVB_Wisam    = [];
        % Angular velocities, their times, and state numbers
        omegaBt_Btgm    = [];
        omegaBt_Gisam   = [];
        omegaBt_Bttrue  = [];
        tOmegaGisam;
        covOmegaGisam   = [];
        omegaStateNums  = [];
        % Rotation from target body to geometric (for estimated vals, one for each state)
        RBttoGisam   = []; 
        RBttoGbsam   = [];
        RBttoGgm     = [];
        RBttoGtrue   = [];
        % Check for simulated data
        isSimulated;            % true if is simulated data        
        % Directories
        dataDir     = '';       % top level test session data directory
        resultsDir  = '';       % results directory with gtsam results etc.
        % Numbers
        imageNums  = [];
        stateNums  = [];
        nStates;
        % Times
        tInspState   = [];      % inspector state times (vertigo times) [s]
        tInspGlobMet = [];      % inspector glob met times (vertigo times) [s]
        tInspTgtTrue = [];      % inspector & target simulated times (spheres & vertigo times) [s]
        tInspImu     = [];      % inspector imu times (vertigo times)  [s]
        tTgtGlobMet  = [];      % target glob met times (vertigo times) [s]
        t0InspSph;              % initial spheres time on the inspector [s]
        t0TgtSph;               % initial spheres time on the target [s]
        % Feature maps (in frame Gc) and stereo frames (in frame B)
        featureMap      = [];   % inliers in centered target geometric frame Gc
        featureMapAll   = [];   % all pts in centered target geometric frame Gc
        firstInliers    = [];   % first inliers in inspector body frame B
        stereoFrames    = [];   % inliers in inspector body frame B
        stereoFramesAll = [];   % all pts in inspector body frame B
        % Spheres data
        spheresData = [];
        % Run info
        frameRate;              % programmed camera frame rate [fps]
        dtFrame;                % programmed time between frames [s]
        spinRate;               % approx magnitude of programmed spin rate [rad/s]
        startTime;              % start time from parameter file (vertigo time) [s]
        endTime;                % end time from parameter file (vertigo time) [s] 
        % Optional array of inertial properties optimizers
        ipo;                    % cell array of optimized inertial properties objects
        % For real data, use the gyro-determined consensus value for Jtrue and RBttoGstrue
        % For simulated data, overwrite with true simulated inertia ratios in GSdata (data dir)
        Jtrue = [1.220680149970564 1.167197687274214]'; % gyro-determined spheres inertia ratios TS53T7R3
        RBttoGstrue  = [0.997678278401242 0.002931380968306 0.068040133832561;   % gyro-determined 
                       -0.002964735103303 0.999995529408951 0.000389240129122;   % TS53T7R3
                       -0.068038688641841 -0.000590057395114 0.997682508957719]; 
        % Execution times
        equalizeTime;           % time to equalize images [s]
        blobTrackTime;          % time to do blob tracking [s]
        stereoFrameTime;        % time to get stereo frame [s]
        visOdomTime;            % time for vis odom [s]
        loopClosureTime;        % time for loop closure [s]
        isamTime;               % time for isam [s]
    end
    
    methods
        %% Constructor with frame rate and aproximate spin rate in rpm
        function this = Results(dataDir, resultsDir, frameRate, spinRpm, startTime, endTime)
            
            % Store frame and spin rates
            this.frameRate      = frameRate;
            this.dtFrame        = 1/frameRate;
            this.spinRate       = spinRpm/60*2*pi;
            this.startTime      = startTime/1000;
            this.endTime        = endTime/1000;
            
            % Load data
            this.dataDir        = dataDir;
            this.resultsDir     = resultsDir;
            spheresData_        = load([dataDir '/spheresData.mat']);
            if isfield(spheresData_,'data')     % real results
                this.isSimulated = 0;
                this.spheresData = spheresData_.data;
            else                                % simulated results
                this.isSimulated = 1;
                this.spheresData = spheresData_.cfg.data;
                trueStateInsp   = spheresData_.cfg.trueStates{1};
                trueStateTgt    = spheresData_.cfg.trueStates{2};
                tTrue           = spheresData_.cfg.trueStateT{1}/1000;
                this.Jtrue      = csvread([dataDir '/GSdata/inertiaRatios_data.txt'])';
                this.Jtrue      = this.Jtrue(1:2);
            end
            globMetInsp         = double(this.spheresData{1}.BackTel.StdState);
            globMetTgt          = double(this.spheresData{2}.BackTel.StdState);
            imuData             = csvread([dataDir '/GSdata/imu_data.txt']);
            vioData             = csvread([dataDir '/GSdata/navState_data.txt']);
            dtSph               = mean(imuData(:,1)-imuData(:,2))/1000;
            tInsp               = double(this.spheresData{1}.BackTel.StdTime)/1000 ... % conv time                                       
                                    - double(this.spheresData{1}.TestsTime)/1000 + dtSph; 
            tTgt                = double(this.spheresData{2}.BackTel.StdTime)/1000 ... % conv time  
                                    - double(this.spheresData{2}.TestsTime)/1000 + dtSph;                   
            stateNumData        = csvread([resultsDir '/stateNums_data.txt']);
           	bsamData            = csvread([resultsDir '/gtsamBatch_data.txt']);
            bsamGeomData        = csvread([resultsDir '/gtsamGeometricBatch_data.txt']);
            transGeomData       = csvread([resultsDir '/gtsamTransGeometric_data.txt']);
            this.featureMap     = csvread([resultsDir '/featureMapIncrementalForwardInliers.csv']);
            this.featureMapAll  = csvread([resultsDir '/featureMapIncremental.csv']);
            this.firstInliers   = csvread([resultsDir '/firstInliers.csv']);
            blobPosition        = csvread([resultsDir '/gtsamBlobPosition_data.txt']);        
        
            % Crop data to relevant time period
            this.tInspState     = stateNumData(1:end-1,1)/1000;     % 1 unused image in sim and iss
            this.stateNums      = stateNumData(1:end-1,2);          % 1 unused image in sim and iss
            this.imageNums      = stateNumData(1:end-1,3);          % 1 unused image in sim and iss
            t0                  = this.tInspState(1);       t0ms = t0*1000;
            tf                  = this.tInspState(end);     tfms = tf*1000;
            globMetInsp         = globMetInsp(:,tInsp>=t0 & tInsp<=tf);
            globMetTgt          = globMetTgt(:,tTgt>=t0 & tTgt<=tf);
            ind0InspSph         = find(tInsp>=t0,1);
            ind0TgtSph          = find(tTgt>=t0,1);
            tInsp               = tInsp(tInsp>=t0 & tInsp<=tf);
            tTgt                = tTgt(tTgt>=t0 & tTgt<=tf);
            vioData             = vioData(vioData(:,1)>=t0ms & vioData(:,1)<=tfms,:);
            vioData             = [(this.imageNums(1):this.imageNums(end))' vioData];
            vioData             = vioData(ismember(vioData(:,1),this.imageNums),2:end);
            this.tInspGlobMet   = tInsp;
            this.tTgtGlobMet    = tTgt;
            this.t0InspSph      = double(this.spheresData{1}.BackTel.StdTime(ind0InspSph))/1000 ... 
                                    - double(this.spheresData{1}.TestsTime)/1000;
            this.t0TgtSph       = double(this.spheresData{2}.BackTel.StdTime(ind0TgtSph))/1000 ... 
                                    - double(this.spheresData{2}.TestsTime)/1000;
            if this.isSimulated
                trueStateInsp       = trueStateInsp(:,tTrue>=t0 & tTrue<=tf);
                trueStateTgt        = trueStateTgt(:,tTrue>=t0 & tTrue<=tf);
                tTrue               = tTrue(tTrue>=t0 & tTrue<=tf);
                this.tInspTgtTrue   = tTrue;
            end
            
            % Make time t = 0 at state 0
            this.tInspGlobMet   = this.tInspGlobMet - this.tInspState(1);
            this.tTgtGlobMet    = this.tTgtGlobMet - this.tInspState(1);
            this.tInspTgtTrue   = this.tInspTgtTrue - this.tInspState(1);
            this.tInspState     = this.tInspState' - this.tInspState(1);
            
            % Get execution time data
            execTimeData        = csvread([resultsDir '/execTimes.txt']);
            this.equalizeTime        = execTimeData(:,2);
            this.blobTrackTime       = execTimeData(:,3);
            this.stereoFrameTime     = execTimeData(:,4);
            this.visOdomTime         = execTimeData(:,5);
            this.loopClosureTime     = execTimeData(:,6);
            this.isamTime            = execTimeData(:,7);
            
            % Get intial conditions and use them to create a reasonable position for the first 
            % geometric frame            
            rTgt        = 0.1;  % [m]
            RWtoB0      = quat2rot(globMetInsp(7:10,1));
            RWtoBt0     = quat2rot(globMetTgt(7:10,1));
            tWtoB0_W    = globMetInsp(1:3,1);
            tWtoG0_W    = tWtoB0_W/norm(tWtoB0_W) * rTgt;     % approx target as a sphere
            RWtoG0      = RWtoB0;                             % G frame rot same as B frame rot
            RB0toG0     = eye(3);                             % G frame rot same as B frame rot
            RBttoG      = RB0toG0 * RWtoB0 * RWtoBt0';        % rot from tgt B to tgt G
            this.tGtoBt_Ggm = -RWtoG0 * tWtoG0_W;             % geometric frame estimate
            this.RBttoGgm = RBttoG;                           % target body to geometric
            
            % Get the position of the target body in the world frame  
            this.tWtoBt_Wisam  = blobPosition(1,:)'; 
            this.tWtoBt_Wbsam  = blobPosition(2,:)';
            if this.isSimulated
                this.tWtoBt_Wtrue  = [0 0.2 0]';
            else
                this.tWtoBt_Wtrue  = [0 0 0]';                
            end
            % Correction terms for mis-estimated blob position
            tWtoWbsam_W = this.tWtoBt_Wtrue - this.tWtoBt_Wbsam; 
            tWtoWisam_W = this.tWtoBt_Wtrue - this.tWtoBt_Wisam;
            
            % Get the poses of the inspector, target, and geometric frame at all image times
            this.nStates    = length(this.stateNums);
            this.PWtoBvio   = zeros(4,4,this.nStates);        % world to inspector
            this.PWtoBgm    = zeros(4,4,this.nStates);        % world to inspector
            this.PWtoBtgm   = zeros(4,4,this.nStates);        % world to target
            this.PWtoGgm    = zeros(4,4,this.nStates);        % world to geometric
            this.PGtoBgm    = zeros(4,4,this.nStates);        % geometric to inspector
            this.PWtoBbsam  = zeros(4,4,this.nStates);        % world to inspector
            this.PWtoGbsam  = zeros(4,4,this.nStates);        % world to geometric
            this.PGtoBbsam  = zeros(4,4,this.nStates);        % geometric to inspector
            this.PWtoBisam  = zeros(4,4,this.nStates,this.nStates); % world to inspector
            this.PWtoGisam  = zeros(4,4,this.nStates,this.nStates); % world to geometric
            this.PGtoBisam  = zeros(4,4,this.nStates,this.nStates); % geometric to inspector
            this.PWtoBtrue  = zeros(4,4,this.nStates);        % world to inspector 
            this.PWtoBttrue = zeros(4,4,this.nStates);        % world to target
            this.covPWtoBisam = zeros(6,6,this.nStates,this.nStates); % world to inspector
            this.covPGtoBisam = zeros(6,6,this.nStates,this.nStates); % geometric to inspector
            % Rotations
            this.RBtoGisam      = zeros(3,3,this.nStates,this.nStates); % body to geometric
            this.RBtoWisam      = zeros(3,3,this.nStates,this.nStates); % body to world
            this.covRBtoGisam   = zeros(3,3,this.nStates,this.nStates); % cov body to geometric
            this.covRBtoWisam   = zeros(3,3,this.nStates,this.nStates); % cov body to world
            % Velocities
            this.vB_Wvio    = zeros(3,this.nStates);          % velocity of inspector
            this.vB_Wgm     = zeros(3,this.nStates);          % velocity of inspector
            this.vB_Wbsam   = zeros(3,this.nStates);          % velocity of inspector
            this.vB_Wisam   = zeros(3,this.nStates,this.nStates);  % velocity of inspector
            this.vB_Wtrue   = zeros(3,this.nStates);          % velocity of inspector
            this.covVB_Wisam = zeros(3,3,this.nStates);       % velocity of inspector
            % Angular velocities
            this.omegaBt_Btgm   = zeros(3,this.nStates);      % angular velocity of target
            this.omegaBt_Bttrue = zeros(3,this.nStates);      % angular velocity of target
            % Translations
            this.tGtoBt_Gisam    = zeros(3,this.nStates);     % geometric to target body
            this.covtGtoBt_Gisam  = zeros(3,3,this.nStates);  % cov geometric to target body
            % Stereo frames
            this.stereoFrames    = cell(1,this.nStates);      % inliers in inspector body frame
            this.stereoFramesAll = cell(1,this.nStates);      % all pts in inspector body frame
            
            for i = 1:this.nStates
                [~,iInsp]       = min(abs(this.tInspGlobMet - this.tInspState(i)));  % index in insp glob met
                [~,iTgt]        = min(abs(this.tTgtGlobMet  - this.tInspState(i)));  % index in tgt glob met
                [~,iInspTgtTrue]= min(abs(this.tInspTgtTrue - this.tInspState(i)));  % index in simulated

                % Poses
                this.PWtoBvio(:,:,i)    = reshape(vioData(i,2:17), 4, 4)';
                
                this.PWtoBgm(:,:,i)     = ...
                    [ quat2rot(globMetInsp(7:10,iInsp))' globMetInsp(1:3,iInsp); 0 0 0 1 ];
                this.PWtoBtgm(:,:,i)    = ...
                	[ quat2rot(globMetTgt(7:10,iTgt))' globMetTgt(1:3,iTgt); 0 0 0 1 ]; 
                
                this.PWtoBbsam(:,:,i)   = reshape(bsamData(i,2:17), 4, 4)' ...
                                            + [zeros(3) tWtoWbsam_W; 0 0 0 1]; 
                                
                RWtoGi                  = RBttoG * quat2rot(globMetTgt(7:10,iTgt));
                tWtoGi_W                = RWtoGi'*RWtoG0 * tWtoG0_W;     % logbook #3 pp 41
                this.PWtoGgm(:,:,i)     = [ RWtoGi' tWtoGi_W; 0 0 0 1 ];
                this.PGtoBgm(:,:,i)     = ...
                    [ RWtoGi  -RWtoGi*tWtoGi_W; 0 0 0 1 ] * this.PWtoBgm(:,:,i);

                this.PGtoBbsam(:,:,i)   = reshape(bsamGeomData(i,2:17), 4, 4)';
                this.PWtoGbsam(:,:,i)   = this.PWtoBbsam(:,:,i) * ... % PWtoG = PWtoB * (PGtoB}^{-1}
                    [ this.PGtoBbsam(1:3,1:3,i)' -this.PGtoBbsam(1:3,1:3,i)'*this.PGtoBbsam(1:3,4,i); 0 0 0 1 ];
                
                if this.isSimulated
                    this.PWtoBtrue(:,:,i)   = [ quat2rot(trueStateInsp(7:10,iInspTgtTrue))' ...
                                                 trueStateInsp(1:3,iInspTgtTrue); 0 0 0 1 ]; 
                    this.PWtoBttrue(:,:,i)  = [ quat2rot(trueStateTgt(7:10,iInspTgtTrue))' ...
                                                 trueStateTgt(1:3,iInspTgtTrue); 0 0 0 1 ];
                else
                    this.PWtoBtrue(:,:,i)   = this.PWtoBvio(:,:,i);
                    this.PWtoBttrue(:,:,i)  = [ ...
                        this.PWtoBtgm(1:3,1:3,i)*Exp([pi 0 0])*this.RBttoGstrue ...
                        this.PWtoBtgm(1:3,4,i); 0 0 0 1 ];
                end
                          
                % Velocities
                this.vB_Wvio(:,i)       = vioData(i,18:20)';
                this.vB_Wgm(:,i)        = globMetInsp(4:6,iInsp);
                this.vB_Wbsam(:,i)      = bsamData(i,18:20)';
                if this.isSimulated
                    this.vB_Wtrue(:,i)  = trueStateInsp(4:6,iInspTgtTrue);
                else
                    this.vB_Wtrue(:,i)  = this.vB_Wvio(:,i);
                end
                
                % Angular velocities; in case of experimental data, mitigate flipped axes effect          
                this.omegaBt_Btgm(:,i)       = globMetTgt(11:13,iTgt);
                if this.isSimulated
                    this.omegaBt_Bttrue(:,i) = trueStateTgt(11:13,iInspTgtTrue);
                else
                    this.omegaBt_Bttrue(:,i) = Exp([pi 0 0]') * this.RBttoGstrue' ...
                        * this.omegaBt_Btgm(:,i);
                end
                
                % Isam data will be 4 dimensional. Do not collect the data until 2nd timestep
                if i == 1
                    continue;
                end
                stateNum        = this.stateNums(i);
                isamData        = csvread([resultsDir '/gtsamIncremental' ...
                                                            num2str(stateNum) '.txt']);
                isamGeomData    = csvread([resultsDir '/gtsamGeometricIncremental'...
                                                            num2str(stateNum) '.txt']);
                isamCovData     = csvread([resultsDir '/gtsamCov' ...
                                                            num2str(stateNum) '.txt']);
                isamCovGeomData = csvread([resultsDir '/gtsamCovGeometric' ...
                                                            num2str(stateNum) '.txt']);
                for j = 1:i                   
                                                        
                    this.PWtoBisam(:,:,j,i) = reshape(isamData(j,2:17), 4, 4)' ....
                                                + [zeros(3) tWtoWisam_W; 0 0 0 0 ];                                          
                    this.PGtoBisam(:,:,j,i) = reshape(isamGeomData(j,2:17), 4, 4)';
                    this.PWtoGisam(:,:,j,i) = this.PWtoBisam(:,:,j,i) * ... 
                        [ this.PGtoBisam(1:3,1:3,j,i)' ...        % PWtoG = PWtoB * (PGtoB}^{-1}
                            -this.PGtoBisam(1:3,1:3,j,i)'*this.PGtoBisam(1:3,4,j,i); 0 0 0 1 ];
                    this.covPWtoBisam(:,:,j,i) = reshape(isamCovData(j,2:37), 6, 6)';
                    this.covPGtoBisam(:,:,j,i) = reshape(isamCovGeomData(j,2:37), 6, 6)';
                
                    this.vB_Wisam(:,j,i)        = isamData(j,18:20)';
                    this.covVB_Wisam(:,:,j,i)   = reshape(isamCovData(j,38:46), 3, 3)';
                        
                    % Rotations
                    this.RBtoWisam(:,:,j,i)    = this.PWtoBisam(1:3,1:3,j,i);
                    this.RBtoGisam(:,:,j,i)    = this.PGtoBisam(1:3,1:3,j,i);
                    this.covRBtoWisam(:,:,j,i) = this.covPWtoBisam(1:3,1:3,j,i);
                    this.covRBtoGisam(:,:,j,i) = this.covPGtoBisam(1:3,1:3,j,i);
                end
                
                if i >= 2
                    this.tGtoBt_Gisam(:,i)        = transGeomData(i-1,2:4);
                    this.covtGtoBt_Gisam(:,:,i)   = reshape(transGeomData(i-1,5:13), 3, 3)';
                end                                  
                
            end
            
            
            % Get stereo frames by reverse engineering from final isam result
            for i = 1:this.nStates
                                
                RBtoG       = this.PGtoBisam(1:3,1:3,i,end);    % insp body to tgt geom
                tGtoB_G     = this.PGtoBisam(1:3,4,i,end);      % tgt geom to insp body
                tGtoBt_G    = this.tGtoBt_Gisam(:,end);         % tgt geom to tgt body
                
                tBttoP_Gc    = this.featureMap(this.featureMap(:,1) == i-1,2:4)';       % inliers in Gc
                tBttoP_GcAll = this.featureMapAll(this.featureMapAll(:,1) == i-1,2:4)'; % all pts in Gc
                Nin          = size(tBttoP_Gc,2);               % number of inliers
                Nall         = size(tBttoP_GcAll,2);            % number of tot pts
                
                tBtoP_B     = RBtoG' * ( -repmat(tGtoB_G,1,Nin) + repmat(tGtoBt_G,1,Nin) ...
                                            + tBttoP_Gc );
                tBtoP_Ball  = RBtoG' * ( -repmat(tGtoB_G,1,Nall) + repmat(tGtoBt_G,1,Nall) ...
                                            + tBttoP_GcAll );
                
                % Add to stereo frames
                this.stereoFrames{i} = tBtoP_B;
                this.stereoFramesAll{i} = tBtoP_Ball;                
                
            end
            
            % Get the translation from the geometric frame to center of mass (logbook #3 pp 155)
            % tGtoB_G = RWtoG0 * (tWtoBt0_W - tWtoB0_W) - RB0toG0 * (1/N Sum_{i=1}^N tBtoPi_B)
%             this.tGtoBt_Gbsam          = transGeomData(2,:)';
            this.tGtoBt_Gtrue          = this.PWtoBtrue(1:3,1:3,1)' * ...
                ( this.PWtoBttrue(1:3,4,1) - this.PWtoBtrue(1:3,4,1) ) ...
                - eye(3) * mean(this.firstInliers(:,2:4))';
            
            % Get the rotation from target body to geometric (init geom frame G0 === insp body
            % frame B0); correct for flipped axes in experimental case
            if this.isSimulated  
                RB0toW          = this.PWtoBtrue(1:3,1:3,1);
                RBt0toW         = this.PWtoBttrue(1:3,1:3,1);                
                this.RBttoGtrue = RB0toW' * RBt0toW;
            else
                RB0toW          = this.PWtoBgm(1:3,1:3,1);
                RGs0toW         = this.PWtoBtgm(1:3,1:3,1);
                this.RBttoGtrue = RB0toW' * RGs0toW * this.RBttoGstrue * Exp([-pi 0 0]');
            end
            
            % Get the angular velocities of the target, their covariance, and their time
            this.omegaBt_Gisam  = zeros(3,this.nStates-1,this.nStates-1);   % ang vel of target
            this.covOmegaGisam  = zeros(3,3,this.nStates-1,this.nStates-1); % cov of ang vel
            this.tOmegaGisam    = zeros(this.nStates-1,this.nStates-1);     % times for ang vels
            this.omegaStateNums = this.stateNums(1:this.nStates-1);
            Nbuff               = 5;                            % buffer length for ang vel avg
            for i = 1:this.nStates-1
                
                omegaBt_Gbuff       = zeros(3,Nbuff);           % buffer of last few valid ang vels
                imageNumBuff        = zeros(1,Nbuff);           % image nums of last few valid vels
                
                for j = 1:i-1
                
                    RBitoGi = this.RBtoGisam(:,:,j,i);        RBitoW  = this.RBtoWisam(:,:,j,i);
                    covRGi  = this.covRBtoGisam(:,:,j,i);     covRWi  = this.covRBtoWisam(:,:,j,i);
                    RBjtoW  = this.RBtoWisam(:,:,j+1,i);      RBjtoGj = this.RBtoGisam(:,:,j+1,i);
                    covRGj  = this.covRBtoGisam(:,:,j+1,i);   covRWj  = this.covRBtoWisam(:,:,j+1,i);
                    
                    % Calculate omegaBt_G and its covariance
                    omegaBt_G   = 1/this.dtFrame * Log(RBitoGi*RBitoW'*RBjtoW*RBjtoGj');
                    covOmegaG   = 1/this.dtFrame * ( RBitoGi*(covRGi + covRWi)*RBitoGi' ...
                        + RBitoGi*RBitoW'*RBjtoW*(covRWj + covRGj)*RBjtoW'*RBitoW*RBitoGi' );

                    % Get average of nonzero buffer angular velocities 
                    omegaBt_Gavg = omegaBt_Gbuff;
                    omegaBt_Gavg(:,all(~any(omegaBt_Gavg),1)) = [];     % remove zero columns
                    omegaBt_Gavg = mean(omegaBt_Gavg,2);
                    
                    % Perform logical checks, checking that the buffer is full, the images are
                    % adjacent, and that unless there is a large skip since the image buffer, that
                    % the angular velocity does not have a large jump in magnitude or vector value.
                    isFirstNbuffImages   = j < Nbuff;
                    areImagesAdjacent    = this.imageNums(j+1)-this.imageNums(j) == 1;
                    isLargeSkipSinceBuff = this.imageNums(j) - imageNumBuff(end) > Nbuff;
                    isMagJumpInLimits    = norm(omegaBt_G)/norm(omegaBt_Gavg) > 0.8 ...
                                            && norm(omegaBt_G)/norm(omegaBt_Gavg) < 1.25;
                    isVecJumpInLimits    = norm(omegaBt_G - omegaBt_Gavg)/norm(omegaBt_Gavg) < 0.6;
                    isValid              = areImagesAdjacent && ( isFirstNbuffImages || ...
                        isLargeSkipSinceBuff || (isMagJumpInLimits && isVecJumpInLimits) );
%                         (isMagJumpInLimits && isVecJumpInLimits) ); % old logic was more strict

                    if isValid

                        this.omegaBt_Gisam(:,j,i)     = omegaBt_G;
                        this.covOmegaGisam(:,:,j,i)   = covOmegaG;

                        omegaBt_Gbuff(:,1)          = [];
                        omegaBt_Gbuff(:,end+1)      = omegaBt_G;
                        imageNumBuff(:,1)           = [];
                        imageNumBuff(:,end+1)       = this.imageNums(j);

                        % Find the closest time to when a frame should have been taken
                        this.tOmegaGisam(j,i) = this.tInspState(1) ...
                          + round((this.tInspState(j)-this.tInspState(1))/this.dtFrame) * this.dtFrame;

                    else
                                                
                        this.tOmegaGisam(j,i)      = Inf;          % invalid
%                         disp(['Rejected angular velocity at iSAM state ' num2str(this.stateNums(i)) ...
%                                 ', image ' num2str(this.imageNums(j))]);
%                         disp(['  Image adjacency: '....
%                             num2str(this.imageNums(j+1)-this.imageNums(j))]);
%                         disp(['  Magnitude jump:  '...
%                             num2str(norm(omegaBt_G)/norm(omegaBt_Gavg))]);
%                         disp(['  Vector jump:     '...
%                             num2str(norm(omegaBt_G-omegaBt_Gavg)/norm(omegaBt_Gavg))]);
                        
                    end
                    
                end
                     
            end
                       
            % Create empty array of inertial properties optimizations, one at each stateNum
            this.ipo            = cell(1,this.nStates);
                        
        end
        
        %% Add an inertial properties optimization object and use it to calculate the optimization
        %  results. Under type, include 'isam', or 'bsam'
        function addInertialPropertiesOpt(this, ipo, type, stateNum)
                        
            % Store inertial properties optimization
            this.ipo{stateNum+1} = ipo;
            
            % Get the estimated rotation from target body to geometric
            eval(['this.RBttoG' type '(:,:,' num2str(stateNum+1) ') = ipo.RBttoG;']);
                        
            
        end
    
    end

end

