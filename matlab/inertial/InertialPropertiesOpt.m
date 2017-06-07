classdef InertialPropertiesOpt < handle
    %INERTIALPROPERTIESOPT Estimates the RigidBodyRotation object given a set of results including
    % rotation matrices RBtoG, RBtoW, their covariances covRBtoG, covRBtoW, and an array of times t
    % See SetterfieldPhD, 'Estimating Inertial Properties' chapter
    
    properties
        RBtoG;              % 3x3xnT rotation matrices from inspector body to target geometric frame
        RBtoGtrue;          % 3x3 true rotation matrix from results
        RBtoW;              % 3x3xnT rotation matrices from inspector body to world frame
        covRBtoG;           % 3x3xnT covariances of rotation RBtoG [rad^2]
        covRBtoW;           % 3x3xnT covariances of rotation RBtoW [rad^2]
        dt;                 % average time between measurements [s]
        tOmega;             % 1xnT-1 vector of times corresponding to angular velocities [s]
        nT;                 % number of timesteps 
        omegaB_G;           % 3xnT-1 target angular velocity in geometric frame [rad/s]
        omegaB_B;           % 3xnT-1 target angular velocity in body frame [rad/s]
        omegaB_Btrue;       % 3xnT-1 target true angular velocity in body frame (from sim or aligned gyros) [rad/s]
        covOmegaG;          % 3x3xnT-1 covariance of omegaB_G [rad^2]
        covOmegaB;          % 3x3xnT-1 covariance of omegaB_B [rad^2]
        isSingleAxis;       % boolean indicating if the body is undergoing single axis rotation
        principalAxesOpt;   % the PrincipalAxesOpt object to use for axis alignment
        RBttoG;             % estimated rotation from the target body to geometric frame
        principalAxesCost;  % cost of principal axes optimization
        RBttoW;             % 3x3xnT-1 rotation matrices from target body to world frame
        inertiaRatiosOpt;   % the InertiaRatiosOpt object to use for inertia ratio estimation
        inertiaRatiosCost;  % cost of the inertia ratios optimization
        rigidBodyRotationOpt;  % the RigidBodyRotation object to store the optimal solution
        rigidBodyRotationInit; % the RigidBodyRotation object to store the ang mom initialization
        constraintIndOpt;   % the index of the optimal equality constraint for inertia ratio opt
        omegaStateNums;     % the states for which valid angular velocities exist
        stateNum;           % the state number being optimized
        execTime;           % the execution time in milliseconds
        paoExecTime;        % principal axis opt execution time
        iroExecTime;        % inertia ratio optimization times
        r;                  % results
    end
    
    methods
        %% Constructor
        function this = InertialPropertiesOpt(results,stateNum)
            
            % Record stateNum being optimized
            this.stateNum = stateNum;
            
            % Narrow the results down to be for the closest state number
            [~,ind]         = min(abs(results.omegaStateNums-stateNum));
            
            this.tOmega     = results.tOmegaGisam(:,ind)';
            this.RBtoG      = results.RBtoGisam(:,:,:,ind);
            this.RBtoW      = results.RBtoWisam(:,:,:,ind);
            this.covRBtoG   = results.covRBtoGisam(:,:,:,ind);
            this.covRBtoW   = results.covRBtoWisam(:,:,:,ind);
            this.omegaB_G   = results.omegaBt_Gisam(:,:,ind);
            this.covOmegaG  = results.covOmegaGisam(:,:,:,ind);
            this.dt         = results.dtFrame;
            
            % Keep only valid data where time is not Inf and ang vel is not zero
            validOmegaInds  = ~isinf(this.tOmega) & [1 this.tOmega(2:end) ~= 0] ...
                                & sum(this.omegaB_G.^2,1) ~= 0;
            this.tOmega     = this.tOmega(validOmegaInds);
            this.RBtoG      = this.RBtoG(:,:,validOmegaInds);
            this.RBtoW      = this.RBtoW(:,:,validOmegaInds);
            this.covRBtoG   = this.covRBtoG(:,:,validOmegaInds);
            this.covRBtoW   = this.covRBtoW(:,:,validOmegaInds);
            this.omegaB_G   = this.omegaB_G(:,validOmegaInds);
            this.covOmegaG  = this.covOmegaG(:,:,validOmegaInds);
            this.omegaStateNums = results.omegaStateNums(validOmegaInds);
            this.nT         = length(this.tOmega);
            
            this.omegaB_Btrue   = results.omegaBt_Bttrue(:,validOmegaInds);
            
            % Create principal axes optimizer
            this.principalAxesOpt = PrincipalAxesOpt(this.omegaB_G);
            
            % Store results
            this.r = results;
                        
        end
        %% Perform optimization, optionally only using data until a certain time step
        function optimize(this,timeStep)
           
            % Time the execution
            tic;
            
            % By default, optimize until the end of the data
            if nargin < 2 || timeStep > this.nT
                timeStep = this.nT;
            end
            
            % Check whether the last three uncertainty ellipsoids intersected to determine if 
            % single or multi-axis rotation
            this.isSingleAxis = 1;
            for i = 1:timeStep
                if  i >= 3 ...
                  && ~this.ellipsoidIntersection(this.omegaB_G(:,1), this.omegaB_G(:,i), ...
                        this.covOmegaG(:,:,1), this.covOmegaG(:,:,i)) ...
                  && ~this.ellipsoidIntersection(this.omegaB_G(:,1), this.omegaB_G(:,i-1), ...
                        this.covOmegaG(:,:,1), this.covOmegaG(:,:,i-1)) ...
                  && ~this.ellipsoidIntersection(this.omegaB_G(:,1), this.omegaB_G(:,i-2), ...
                        this.covOmegaG(:,:,1), this.covOmegaG(:,:,i-2))
                    
                    this.isSingleAxis = 0;
                    disp(['Axis is multi-axis based on timesteps: ' num2str(i-2) ',' num2str(i-1) ...
                        ',' num2str(i)]);
                    
                    % Plot ellipsoid evidence
%                     figure(5); clf;
%                     this.plotAxisDetermination(i-2:i); xlabel('{}^G\omega_1 [rad/s]');
%                     ylabel('{}^G\omega_2 [rad/s]'); zlabel('{}^G\omega_3 [rad/s]');
%                     title(['Evidence of multi-axis rotation using 1-\sigma covariance ellipses '...
%                         'of timesteps: ' num2str(i-2) ', ' num2str(i-1) ', and ' num2str(i)]);
%                     drawnow;

                    break;
                    
                end
            end
                        
            if this.isSingleAxis     % single-axis rotation -------------------------------------------
                
                disp('Single axis rotation detected');
                
                % Assume frame G == frame B
                this.RBttoG = eye(3);
                
                % Create fully-symmetric rotation with frame G == frame B
                RBt0toW   = this.RBtoW(:,:,1) * this.RBtoG(:,:,1)' * this.RBttoG;
                omegaB0_B = this.RBttoG' * this.omegaB_G(:,1); 
                this.rigidBodyRotationOpt = RigidBodyRotation([1 1]',RBt0toW,omegaB0_B,'omega0',0);
                this.rigidBodyRotationInit = this.rigidBodyRotationOpt;
                
                % Unable to perform optimization, assign infinite cost
                this.principalAxesCost = Inf;
                this.inertiaRatiosCost = Inf;
                
                this.paoExecTime = toc;
                this.iroExecTime = 0;
                
            else                % multi-axis rotation --------------------------------------------
                                
                % Perform principal axis optimization
                [this.RBttoG,this.principalAxesCost] = this.principalAxesOpt.optimize();
                
                this.paoExecTime = toc;
                tic;
                
%                 figure(1); clf;
%                 this.principalAxesOpt.plot();
%                 drawnow;
                                
                % Get estimated angular velocities in body frame
                this.omegaB_B = this.RBttoG' * this.omegaB_G;
                
                % Get the estimated covariance of angular velocities in the body frame (logbook #3 
                % pp 130) as well as the orientation of body frame
                this.covOmegaB = zeros(3,3,this.nT); this.RBttoW = zeros(3,3,this.nT);
                for i = 1:this.nT
                    this.covOmegaB(:,:,i) = this.RBttoG' * this.covOmegaG(:,:,i) * this.RBttoG;
                    RGitoW = this.RBtoW(:,:,i) * this.RBtoG(:,:,i)';
                    this.RBttoW(:,:,i) = RGitoW * this.RBttoG;
                end
                                                
                % Create inertia ratios optimization and then optimize
                this.inertiaRatiosOpt = InertiaRatiosOpt(this.omegaB_B,this.covOmegaB, ...
                    this.RBttoW,this.tOmega, this.principalAxesOpt.energyState, ...
                    this.principalAxesOpt.inertiaSymmetry);
                [~,~,this.constraintIndOpt] = this.inertiaRatiosOpt.optimize();
                this.inertiaRatiosCost = this.inertiaRatiosOpt.getCost();
                
                this.iroExecTime = toc;
                
                % Get optimal and initial rigid body rotation objects
                this.rigidBodyRotationOpt = this.inertiaRatiosOpt.getRigidBodyRotation('opt');
                this.rigidBodyRotationInit = this.inertiaRatiosOpt.getRigidBodyRotation('init');
            
            end
            
            this.execTime = toc;
                        
        end
        %% Plot the resuls of the optimization; if ground truth available, plot elsewhere
        function plot(this)
            

            % Create an extended timeframe into the future
            tExt = [this.tOmega]; % ...
%                         this.tOmega(end)+(this.dt:this.dt:(this.tOmega(end)-this.tOmega(1)))];
            
            % Get intial and optimal values for angular velocities
            omegaB_Binit = this.rigidBodyRotationInit.predictOmega(tExt);
            omegaB_Bopt = this.rigidBodyRotationOpt.predictOmega(tExt);
%             omegaB_Gopt = this.RBttoG * omegaB_Bopt;
                                    
            % Get truth
            if this.r.isSimulated
                thetaCorr = Log(this.r.RBttoGtrue'*this.RBttoG);
                omegaB_BtrueCorr = Exp([-thetaCorr(1) 0 0]) * this.omegaB_Btrue;
            else
                omegaB_BtrueCorr = this.omegaB_Btrue;
            end
            
            % Plot the measured body ang vels as well as the initial and optimized solution            
            subplot(3,2,2);
            plot(this.tOmega,this.omegaB_B(1,:),'Color',[0.5 0 0.5],'LineWidth',1); hold on; grid on;
            plot(tExt,omegaB_BtrueCorr(1,:),'-c','LineWidth',2);
            plot(tExt,omegaB_Binit(1,:),'--r');
            plot(tExt,omegaB_Bopt(1,:),'-.k','LineWidth',2);
            ylabel('{}^B\omega_1 [rad/s]'); xlim([this.tOmega(1) this.tOmega(end)]);
            if this.r.isSimulated
                legend('Aligned Meas.','True (In-Plane Corrected)','Initialization','Optimal Fit');
            end
            title(['Body Angular Velocities Fit. J^*=[' ...
                num2str(this.rigidBodyRotationOpt.J(1)) ',' ...
                num2str(this.rigidBodyRotationOpt.J(2)) ']']);
            subplot(3,2,4);
            plot(this.tOmega,this.omegaB_B(2,:),'Color',[0.5 0 0.5],'LineWidth',1); hold on; grid on;
            plot(tExt,omegaB_Binit(2,:),'--r');
            plot(tExt,omegaB_BtrueCorr(2,:),'-c','LineWidth',2);
            plot(tExt,omegaB_Bopt(2,:),'-.k','LineWidth',2);
            ylabel('{}^B\omega_2 [rad/s]'); xlim([this.tOmega(1) this.tOmega(end)]);
            if ~this.r.isSimulated
                legend('Aligned Meas.','True','Initialization','Optimal Fit');
            end
            subplot(3,2,6);
            plot(this.tOmega,this.omegaB_B(3,:),'Color',[0.5 0 0.5],'LineWidth',1); hold on; grid on;
            plot(tExt,omegaB_Binit(3,:),'--r');
            plot(tExt,omegaB_BtrueCorr(3,:),'-c','LineWidth',2);
            plot(tExt,omegaB_Bopt(3,:),'-.k','LineWidth',2);
            ylabel('{}^B\omega_3 [rad/s]');
            xlabel('Time [s]'); xlim([this.tOmega(1) this.tOmega(end)]);
            
            % Plot the principal axes alignment on the left hand side
            subplot(3,2,[1 3 5]);
            this.principalAxesOpt.plot(0,1);
            title('Principal Axes Alignment'); 
            xlabel('{}^B\omega_1 [rad/s]'); ylabel('{}^B\omega_2 [rad/s]'); 
            zlabel('{}^B\omega_3 [rad/s]');
            plot3(omegaB_Bopt(1,:),omegaB_Bopt(2,:),omegaB_Bopt(3,:),'-.k','LineWidth',2); hold on; grid on;
            plot3(this.omegaB_Btrue(1,:),this.omegaB_Btrue(2,:),this.omegaB_Btrue(3,:),'-c','LineWidth',2);
%             legend('Measured {}^G\omega_B','Estimated {}^Gx_E','Estimated {}^Gy_E', ...
%                 'Estimated {}^Gz_E','Estimated {}^Gx_B','Estimated {}^Gy_B', ...
%                 'Estimated {}^Gz_B','Final Fit Conic 1','Final Fit Conic 2','Final Fit Conic 3', ...
%                 'Estimated {}^G\omega_B','True {}^B\omega_B','Location','EastOutside');
%             view(140,30);
            
        end
        %% Plot logic used to determine single/multi-axis for specific set of timesteps
        function plotAxisDetermination(this,timeSteps)
                        
            col = jet(length(timeSteps));
            
            for i = timeSteps
                
                % Get color index
                j = find(timeSteps==i);
                                
                % Plot lines from origin to angular velocities
                plot3([0 this.omegaB_G(1,1)],[0 this.omegaB_G(2,1)],[0 this.omegaB_G(3,1)],'-g');
                hold on; grid on; axis equal;
                plot3([0 this.omegaB_G(1,i)],[0 this.omegaB_G(2,i)],[0 this.omegaB_G(3,i)], ...
                    'Color',col(j,:));

                % Plot covariance ellipsoids
                InertialPropertiesOpt.plotEllipsoids(this.omegaB_G(:,1),this.omegaB_G(:,i), ...
                    this.covOmegaG(:,:,1),this.covOmegaG(:,:,i),[0 1 0],col(j,:));
                
            end
            
            % Plot triad
            omegaMax = max(max(this.omegaB_G));
            plot3([0 omegaMax],[0 0],[0 0],'-r');
            plot3([0 0],[0 omegaMax],[0 0],'-g');
            plot3([0 0],[0 0],[0 omegaMax],'-b');
            
        end
    end
    
    methods(Static)
        %% Checks whether two 1-sigma uncertainty covariances intersect; returns 1 if they do
        function ellipsoidsIntersect = ellipsoidIntersection(mu1,mu2,cov1,cov2)

            % Set k to make x'*Cov^{-1}*x = k contain 0.6826 of the probability (1-sigma equiv)
%             k = 14.157; % 3-sigma [Maybeck1979 pp 366]
            k = 3.527;  % 1-sigma [Maybeck1979 pp 366]
            
            % Use technique from Alfono 2003 to see if they intersect [logbook #3 pp 115]
            S1 = [ inv(cov1)   zeros(3,1) ;
                   zeros(1,3) -k          ];
            S2 = [ inv(cov2)   zeros(3,1) ;
                   zeros(1,3) -k          ];
            T1 = [ eye(3)                   zeros(3,1) ;
                   -mu1(1) -mu1(2) -mu1(3)  1          ];
            T2 = [ eye(3)                   zeros(3,1) ;
                   -mu2(1) -mu2(2) -mu2(3)  1          ];
            A = T1*S1*T1';
            B = T2*S2*T2';
                      
            % Default to no intersection of ellipsoid
            ellipsoidsIntersect = 0;
            
            % Check to see if one mean lies inside the other ellipsoid [Alfano2003]
            if [mu2;1]'*A*[mu2;1] <= 0 || [mu1;1]'*B*[mu1;1] <= 0
                ellipsoidsIntersect = 1;        % mu2 inside ellipsoid 1 or mu1 inside ellipsoid 2
            else    

                % Attempt to find viable intersection points
                [vecs,~] = eig(A\B);

                % Remove all with scalar term of zero
                vecs = vecs(:,vecs(4,:)~=0); 

                if ~isempty(vecs)

                    % Normalize so that the scalar term is 1
                    vecs = vecs./repmat(vecs(4,:),4,1);

                    % Loop through scaled eigenvectors and see if any of them 
                    % satisfy both x'*A*x and x'*B*x
                    for i = 1:size(vecs,2)
                        if abs(vecs(:,i)'*A*vecs(:,i)) < 1e-6 && abs(vecs(:,i)'*B*vecs(:,i)) < 1e-6
                            ellipsoidsIntersect = 1;    % x''*A*x=0 and x''*B*x=0 both satisfied
                        end
                    end
                end
                
            end
            
        end
        %% Plots two 1-sigma covariance ellipsoids (x'-mu')*cov^{-1}*(x-mu) = k
        function plotEllipsoids(mu1,mu2,cov1,cov2,col1,col2)
            
            % Set k to make x'*Cov^{-1}*x = k contain 0.6826 of the probability (1-sigma equiv)
%             k = 14.157; % 3-sigma [Maybeck1979 pp 366]
            k = 3.527;  % 1-sigma [Maybeck1979 pp 366]
            
            % Get eigvecs and eigvals of the inverse covariance matrices
            [vec1,vals1] = eig(cov1);
            [vec2,vals2] = eig(cov2);

            % Get rotation matrices REtoW and radii for both cases
            RE1toW = [vec1(:,1) vec1(:,2) cross(vec1(:,1),vec1(:,2))]; 
            RE2toW = [vec2(:,1) vec2(:,2) cross(vec2(:,1),vec2(:,2))];
            radii1 = sqrt(k)*sqrt(diag(vals1)); 
            radii2 = sqrt(k)*sqrt(diag(vals2));
                                    
            % Get centered ellipsoid for both cases in ellipsoid frame
            [x1_E,y1_E,z1_E] = ellipsoid(0,0,0,radii1(1),radii1(2),radii1(3));
            [x2_E,y2_E,z2_E] = ellipsoid(0,0,0,radii2(1),radii2(2),radii2(3));

            % Perform point rotation to rotate the ellipsoid
            x1_W = zeros(size(x1_E)); y1_W = zeros(size(x1_E)); z1_W = zeros(size(x1_E));
            x2_W = zeros(size(x2_E)); y2_W = zeros(size(x2_E)); z2_W = zeros(size(x2_E));
            for i = 1:size(x1_E,1)
                for j = 1:size(x1_E,2)
                    p1 = RE1toW * [x1_E(i,j) y1_E(i,j) z1_E(i,j)]' + mu1;
                    p2 = RE2toW * [x2_E(i,j) y2_E(i,j) z2_E(i,j)]' + mu2;
                    x1_W(i,j) = p1(1); y1_W(i,j) = p1(2); z1_W(i,j) = p1(3);
                    x2_W(i,j) = p2(1); y2_W(i,j) = p2(2); z2_W(i,j) = p2(3);
                end
            end
            
            surf(x1_W,y1_W,z1_W,'FaceColor',col1,'FaceAlpha',0.2,'EdgeAlpha',0.1); hold on; grid on;
            surf(x2_W,y2_W,z2_W,'FaceColor',col2,'FaceAlpha',0.2,'EdgeAlpha',0.1);
            plot3(mu1(1),mu1(2),mu1(3),'o','Color',col1);
            plot3(mu2(1),mu2(2),mu2(3),'o','Color',col2);
            
        end
    end
    
end

