classdef InertiaRatiosOpt < handle
    %INERTIARATIOOPT InertiaRatio optimizer.
    %   Attempts to find the optimal inertia ratios to fit a given sequence of body angular
    %   velocities, with knowledge of whether it is low / high energy.
    
    properties
        omegaB_B;               % 3xnT matrix of rigid body angular velocity data in the body frame
        covOmegaB;              % 3x3xnT matrix of angular velocity covariances in the body frame
        RBtoW;                  % 3x3xnT matrix of rigid body orientation
        t;                      % a vector of times when omegaB_B was sampled
        energyState;            % string, 'LE' when h^2 > 2*Ek*J2, 'HE' when h^2 < 2*Ek*J2
        inertiaSymmetry;        % either tri-axial 'TA', or axis-symmetric 'AS1' or 'AS3'
        inertiaConstraints;     % an InertiaConstraints object containing inertia ratio constraints
        omegaMax;               % the fit maximum velocities in each direction
        p0;                     % 1x3 (TA) 1x1 (AS) vector of params from angular momentum init
        J0;                     % 2x3 (TA) 2x1 (AS) matrix of inertia ratios from ang momentum init
        cost0;                  % 1x3 (TA) 1x1 (AS) vector of initial costs from ang momentum init
        t00;                    % 1x3 (TA) 1x1 (AS) vector of initial times from ang momuntum init
        pOpt;                   % 1x3 (TA) 1x1 (AS) vector of optimal parameters from optimization
        Jopt;                   % 2x3 (TA) 2x1 (AS) matrix of optimal inertia ratios from opt
        costOpt;                % 1x3 (TA) 1x1 (AS) vector of optimal cost for optimizations
        t0opt;                  % 1x3 (TA) 1x1 (AS) vector of optimal initial times
        indOpt;                 % the index of the constraint along which the optimal result lies
    end
    
    methods
        %% Construct from body angular velocity data, energy state, and vector of times
        function this = InertiaRatiosOpt(omegaB_B,covOmegaB,RBtoW,t,energyState,inertiaSymmetry)
            
            % Store angular velocity, orientation, times, and energy state
            this.omegaB_B = omegaB_B; 
            this.covOmegaB = covOmegaB;
            this.RBtoW = RBtoW;
            this.t = t;
            this.energyState = energyState;
            this.inertiaSymmetry = inertiaSymmetry;
            
            % Create inertia ratio constraints
            this.inertiaConstraints = InertiaConstraints(omegaB_B,energyState,inertiaSymmetry);
            
            % Get the maximum angular velocities fit in each direction from the conic fits
            this.omegaMax = zeros(3,1);
            if strcmp(this.inertiaSymmetry,'TA')
                this.omegaMax(1) = this.inertiaConstraints.conicFits(1).aSemi;      % proj plane 1-2
                this.omegaMax(3) = this.inertiaConstraints.conicFits(2).bSemi;      % proj plane 2-3
                if strcmp(energyState,'LE')         % (LE) case
                    this.omegaMax(2) = this.inertiaConstraints.conicFits(2).aSemi;  % proj plane 2-3
                elseif strcmp(energyState,'HE');    % (HE) case
                    this.omegaMax(2) = this.inertiaConstraints.conicFits(1).bSemi;  % proj plane 1-2                
                end
            elseif strcmp(this.inertiaSymmetry,'AS1')
                this.omegaMax(1) = mean(omegaB_B(1,:));
                this.omegaMax(2) = this.inertiaConstraints.conicFits(2).aSemi;      % proj plane 2-3
                this.omegaMax(3) = this.inertiaConstraints.conicFits(2).bSemi;      % proj plane 2-3
            elseif strcmp(this.inertiaSymmetry,'AS3')
                this.omegaMax(1) = this.inertiaConstraints.conicFits(1).aSemi;      % proj plane 1-2
                this.omegaMax(2) = this.inertiaConstraints.conicFits(1).bSemi;      % proj plane 1-2
                this.omegaMax(3) = mean(omegaB_B(3,:));                
            end
            
        end
        %% Get the RigidBodyRotation object corresponding (either 'opt', or 'init')
        function rigidBodyRotation = getRigidBodyRotation(this,optOrInit)
            
            % Get inertia ratios and initial time
            if strcmp(optOrInit,'init')
                J  = this.inertiaConstraints.ratios(this.p0(this.indOpt),this.indOpt);
                t0 = this.t00(this.indOpt); 
            else
                J  = this.inertiaConstraints.ratios(this.pOpt(this.indOpt),this.indOpt);     
                t0 = this.t0opt(this.indOpt); 
            end
            J1 = J(1); J2 = J(2); J3 = 1;
            
            % Estimate kinetic energy and angular momentum
            omegaBarSq = 1/size(this.omegaB_B,2)*sum(this.omegaB_B.^2,2);
            Ek = 1/2*( J1  *omegaBarSq(1) + J2  *omegaBarSq(2) + J3  *omegaBarSq(3) );
            h  = sqrt( J1^2*omegaBarSq(1) + J2^2*omegaBarSq(2) + J3^2*omegaBarSq(3) );
            
            % Create the RigidBodyRotation object
            rigidBodyRotation = RigidBodyRotation(J,this.RBtoW(:,:,1),[Ek h t0],'Ekht0',1);
            
        end
        %% Get the nutation angular velocity omegaP, the modulus k, and the signs s
        function [omegaP,k,s,omegaA,omegaT] = getConstants(this,p,constraintNum)
            
            % Get the inertia ratios to use
            J = this.inertiaConstraints.ratios(p,constraintNum);
            J1 = J(1); J2 = J(2); J3 = 1;
            
            % Estimate kinetic energy and angular momentum
            omegaBarSq = 1/size(this.omegaB_B,2)*sum(this.omegaB_B.^2,2);
            Ek = 1/2*( J1  *omegaBarSq(1) + J2  *omegaBarSq(2) + J3  *omegaBarSq(3) );
            h  = sqrt( J1^2*omegaBarSq(1) + J2^2*omegaBarSq(2) + J3^2*omegaBarSq(3) );
            
            % Get body nutation rate, elliptical modulus, and signs of the angular velocities
            if strcmp(this.inertiaSymmetry,'TA')        % tri-axial
                omegaA = 0; omegaT = 0;                 % only calculate for axis-symm
                if strcmp(this.energyState,'LE')        % (LE) case
                    omegaP = sqrt( (h^2 - 2*Ek*J3)*(J1-J2) / (J1*J2*J3) );
                    k = sqrt( (2*Ek*J1 - h^2)/(J2*(J1 - J2)) * (J2*(J2 - J3))/(h^2 - 2*Ek*J3)  );
                    s = [1 -1 1]';
                elseif strcmp(this.energyState,'HE')    % (HE) case
                    omegaP = sqrt( (2*Ek*J1 - h^2)*(J2-J3) / (J1*J2*J3) );
                    k = sqrt( (h^2 - 2*Ek*J3)/(J2*(J2 - J3)) * (J2*(J1 - J2))/(2*Ek*J1 - h^2)  );
                    s = [-1 1 1]';
                end
            else                                        % axis-symmetric
                if strcmp(this.inertiaSymmetry,'AS1')  
                    Ja = J1; Jt = J3;                   % axial and transverse inertias                
                elseif strcmp(this.inertiaSymmetry,'AS3')
                    Ja = J3; Jt = J1;                   % axial and transverse inertias 
                end
                omegaA = sqrt( (h^2-2*Jt*Ek)/(Ja^2-Jt*Ja) ); % axial ang vel
                omegaT = sqrt( (h^2-2*Ja*Ek)/(Jt^2-Jt*Ja) ); % transverse ang vel
                     
                omegaP = (Jt-Ja)/Jt * omegaA;
                k = 0;                                  % not required for sin and cos
                s = [1 1 1]';
            end
            
        end
        %% Predict the angular velocities at times t given a parameter p along the selected
        %  equality constraint and a value for cycle start time t0 
        function [omegaB_B,omegaP,k,s] = predictOmega(this,p,t0,t,constraintNum)
            
            % Get the constants required in evaluation of omegaB_B
            [omegaP,k,s,omegaA,omegaT] = this.getConstants(p,constraintNum);
            
            % Get the angular velocities at times t
            if strcmp(this.inertiaSymmetry,'TA')        % tri-axial
                m = k^2;
                [sn,cn,dn] = ellipj(omegaP*(t-t0), m, 1e-10);
                if strcmp(this.energyState,'LE')        % (LE) case
                    omegaB_B = [ s(1) * this.omegaMax(1) * dn ;
                                 s(2) * this.omegaMax(2) * sn ;
                                 s(3) * this.omegaMax(3) * cn ];
                elseif strcmp(this.energyState,'HE')    % (HE) case
                    omegaB_B = [ s(1) * this.omegaMax(1) * cn ;
                                 s(2) * this.omegaMax(2) * sn ;
                                 s(3) * this.omegaMax(3) * dn ];
                end
            elseif strcmp(this.inertiaSymmetry,'AS1')   % axis-symmetric J2==J3                
                omegaB_B = [    repmat(s(1) * omegaA,1,length(t))  ;
                                s(2) * omegaT * sin(omegaP*(t-t0)) ;
                                s(3) * omegaT * cos(omegaP*(t-t0)) ];
            elseif strcmp(this.inertiaSymmetry,'AS3')   % axis-symmetric J1==J1
                omegaB_B = [    s(1) * omegaT * sin(omegaP*(t-t0)) ;
                                s(2) * omegaT * cos(omegaP*(t-t0)) ;
                                repmat(s(1) * omegaA,1,length(t))  ];
            end
            
        end
        %% Predict the angular velocities using the optimal fit (after optimize)
        function [omegaB_B,omegaP,k,s] = predictOmegaOpt(this,t)
            [omegaB_B,omegaP,k,s] = this.predictOmega(this.pOpt(this.indOpt), ...
                                        this.t0opt(this.indOpt),t,this.indOpt);
        end
        %% Predict the initial angular velocity along the optimal constraint (after optimize)
        function [omegaB_B,omegaP,k,s] = predictOmegaInit(this,t)
            [omegaB_B,omegaP,k,s] = this.predictOmega(this.p0(this.indOpt), ... 
                                        this.t00(this.indOpt),t,this.indOpt);
        end
        %% Calculate the time in the past at which the canonical versions of sn(t0) = 0, cn(t0) = 1,
        %  dn(t0) = 1, by averaging the result obtained from each data point (RigidBodyRotation.m)
        function t0 = calculateT0(this,p,constraintNum)
                        
            % Get the constants required in evaluation of omegaB_B
            [omegaP,k,s] = this.getConstants(p,constraintNum);
            
            % Get the canonical angular velocities
            omegaB_Bc = [   s(1)*this.omegaB_B(1,:)/this.omegaMax(1)    ;
                            s(2)*this.omegaB_B(2,:)/this.omegaMax(2)    ;
                            s(3)*this.omegaB_B(3,:)/this.omegaMax(3)    ];
            omegaB_Bc(omegaB_Bc > 1) = 1; omegaB_Bc(omegaB_Bc < -1) = -1;    % fix numerical issues
                  
            
            if strcmp(this.inertiaSymmetry,'TA')
            
                % Get canonical version of angular velocities in cn and sn-axes
                omegaSnC = omegaB_Bc(2,:);
                if strcmp(this.energyState,'LE')
                    omegaCnC = omegaB_Bc(3,:);
                elseif strcmp(this.energyState,'HE')
                    omegaCnC = omegaB_Bc(1,:);
                end

                % Get the quarter-period in terms of u (e.g. sn(u, m))
                m = k^2;
                
                % Get the values of asn, which will be between -Tu and Tu
                u = ellipticF(atan2(omegaSnC,omegaCnC),m);
                
                % Get the estimated value of u0 from each sample                
                u0 = u - omegaP*this.t;
                
                % Get the estimated value of phi0 from each sample
                [snU0,cnU0,~] = ellipj(u0,m,1e-10);
                phi0 = atan2(snU0,cnU0);
       
                % Get the mean angle and its std dev by averaging the complex unit vector
                % representation (see logbook #3 pp 122)
                rho0Bar   = mean(cos(phi0)) + mean(sin(phi0))*1i;
                phi0Bar   = angle(rho0Bar);
                sigmaPhi0 = sqrt(log(1/abs(rho0Bar)^2));
                
                % Get the indices of phi0 values within 2 sigma using dot product test
                inds2Sig = ([cos(phi0)' sin(phi0)']*...
                                [cos(phi0Bar); sin(phi0Bar)] >= cos(2*sigmaPhi0));
                
                % Retake mean without outliers                
                rho0Bar   = mean(cos(phi0(inds2Sig))) + mean(sin(phi0(inds2Sig)))*1i;
                phi0Bar   = wrapTo2Pi(angle(rho0Bar));
                
                % Get averaged u0, which will be in the proper quadrant (no asin)
                u0Bar = ellipticF(phi0Bar,m);
                
                % Get the estimated value of t0 (see logbook #3 pp 123)
                t0 = -u0Bar/omegaP;
                
                %{
                
                % Get the value of asin, which will be between -pi/2 and pi/2
                phi = atan2(omegaSnC,omegaCnC);
                
                % Get the mean angle and its std dev by averaging the complex unit vector
                % representation (see logbook #3 pp 122)
                rhoBar   = mean(cos(phi)) + mean(sin(phi))*1i;
                phiBar   = angle(rhoBar);
                sigmaPhi = sqrt(1 - abs(rhoBar));
                
                % Retake mean with angles outside of 2 sigma of mean value using dot prod test
                phi     = phi([cos(phi)' sin(phi)']*[cos(phiBar); sin(phiBar)] >= cos(sigmaPhi));
                phiBar  = mean(phi);
                
                % Get the values of asn, which will be between -T and T
                u = ellipticF(phiBar,m);
                
                % Determine which quadrant the resulting solution is in
                if omegaSnC(1) >= 0 && omegaCnC(1) >= 0         % Q1
                    u = u;
                elseif (omegaSnC(1) >= 0 && omegaCnC(1) < 0)    ...
                    || (omegaSnC(1) <  0 && omegaCnC(1) < 0)    % Q2 or Q3
                    u = 2*Tu - u;
                else                                            % Q4
                    u = 4*Tu + u;
                end
                %}
                
                %{
                % Get the values of asn, which will be between -T and T
                u = ellipticF(asin(omegaSnC), m);

                % Find out which angular velocties are in each quadrant
                indsQ1   = (omegaSnC >= 0 & omegaCnC >= 0);
                indsQ2Q3 = (omegaSnC >= 0 & omegaCnC <  0) | (omegaSnC <  0 & omegaCnC < 0);
                indsQ4   = (omegaSnC <  0 & omegaCnC >=  0);

                % Add the appropriate number of quarter periods depending on quadrant
                u(indsQ1)   =        u(indsQ1);
                u(indsQ2Q3) = 2*Tu - u(indsQ2Q3);
                u(indsQ4)   = 4*Tu + u(indsQ4);
                
                % The magnitude of t0 will then be u/omegaP, but make it negative, so that
                % our experiment starts at time=0, rather than at t0. Additionally, each data
                % point was taken at time t, so that time needs to be subtracted.
                t0   = -u/omegaP + this.t;
                %}
                
            elseif strcmp(this.inertiaSymmetry,'AS1') || strcmp(this.inertiaSymmetry,'AS3')
                
                % Determine which axis is driven by sin and which is driven by cos
                if strcmp(this.inertiaSymmetry,'AS1')   
                    omegaSinC = omegaB_Bc(2,:); omegaCosC = omegaB_Bc(3,:);
                elseif strcmp(this.inertiaSymmetry,'AS3')
                    omegaSinC = omegaB_Bc(1,:); omegaCosC = omegaB_Bc(2,:);                    
                end
                
                % Get the angle at the cananonical sin and cos curves
                phi = atan2(omegaSinC,omegaCosC);
                                
                % Get the estimated theta0 value from each sample
                phi0 = phi - omegaP*this.t;
                
                % Get the mean angle and its std dev by averaging the complex unit vector
                % representation (see logbook #3 pp 122)
                rho0Bar   = mean(cos(phi0)) + mean(sin(phi0))*1i;
                phi0Bar   = angle(rho0Bar);
                sigmaPhi0 = sqrt(log(1/abs(rho0Bar)^2));
                
                % Get the indices of phi0 values within 2 sigma using dot product test
                inds2Sig = ([cos(phi0)' sin(phi0)']* ...
                                [cos(phi0Bar); sin(phi0Bar)] >= cos(2*sigmaPhi0));
                
                % Retake mean without outliers 
                rho0Bar = mean(cos(phi0(inds2Sig))) + mean(sin(phi0(inds2Sig)))*1i;
                phi0Bar = wrapTo2Pi(angle(rho0Bar));
                
                % Get the estimated value of theta0 (see logbook #3 pp 124)
                t0 = -phi0Bar/omegaP;
                
                % The magnitude of t0 will then be theta/omegaP, but make it negative, so that
                % our experiment starts at time=0, rather than at t0. Additionally, each data
                % point was taken at time t, so that time needs to be subtracted.
%                 t0   = -theta/omegaP + this.t;
                
            end
            
        end
        %% Calculate the angular velocity component parallel to the angular momentum (i.e.
        %  perpendicular to the invariable plane)
        function omegaPar_W = calculateOmegaPar(this)
            
            % Get the number of timesteps
            nT = length(this.t);
            
            % Get the angular velocties in the world frame W
            omegaB_W = zeros(3,nT);
            for i = 1:nT
                omegaB_W(:,i) = this.RBtoW(:,:,i) * this.omegaB_B(:,i);
            end
            
            % Use technique from [Masutani1994] (logbook #2 pp. #207, #3 pp. 53) to determine 
            % omegaParallel=2*Ek/h and invariable plane normal vector. Solve the normal equations 
            % for the cost function below, where n = omegaPar/||omegaPar||^2
            %   argmin_{n} Sum_{i} || n^T*omega_i - 1 ||^2 
            % = argmin_{n} Sum_{i} || A*n - b||^2
            A = omegaB_W';  b = ones(nT,1);  
            n = (A'*A)\A'*b;
            omegaPar_W = n/norm(n)^2;
            
        end
        %% Get an initialization for the parameter p0 for the specified equality constraint
        function [p0,J0] = initialize(this,constraintNum)
           
            % Get the number of timesteps
            nT = length(this.t);
            
            % Calculate the angular velocity parallel to the angular momentum
            omegaPar_W = this.calculateOmegaPar();
            
            % Assign optimal inertia ratios that match the direction of omegaPar with that
            % of the angular momentum h (see logbook #3 pp. 59-60 and 77-78). Solve the normal
            % equations for the cost function below, where n = [I1 I2 I3]'*||omegaPar||^2/2*T,
            % I1, I2, I3, are conventional principal inertias, and T is kinetic energy.
            %   argmin_{n} Sum_{i} || RBitoW*diag(omegaBi_Bi) - omegaPar ||^2
            % = argmin_{n} Sum_{i} || A*n - b ||^2
            A = zeros(3*nT,3); 
            b = repmat(omegaPar_W,nT,1);
            for i = 1:nT
                A(3*(i-1)+1:3*(i-1)+3,:) = this.RBtoW(:,:,i) * diag(this.omegaB_B(:,i));
            end
            n = (A'*A)\(A'*b);
            J = n(1:2)/n(3);
            
            % The inertia ratios obtained above will not be optimal in the timeseries sense.
            % Find the parameter value at the closest point to J on the equality constraint.
            p0 = this.inertiaConstraints.findClosestPoint(J,constraintNum);
            J0 = this.inertiaConstraints.ratios(p0,constraintNum);
            
        end
        %% Get the cost of specific proposed parameter p and initial time t0, using the specified
        %  equality constraint
        function Cost = cost(this,p,constraintNum)
            
            % Get the best value of t0 for this parameter
            t0 = this.calculateT0(p,constraintNum);
                        
            % Get the predicted omegaB_B at the times of the sampled omegaB_B
            omegaB_Bp = this.predictOmega(p,t0,this.t,constraintNum);
            
            % TEMP, testing --------------------------- resulst are comparable but not better, even
            % with cost being the sum of both of them
%             J = this.inertiaConstraints.ratios(p,constraintNum);
%             J1 = J(1); J2 = J(2); J3 = 1;
%             omegaBarSq = 1/size(this.omegaB_B,2)*sum(this.omegaB_B.^2,2);
%             Ek = 1/2*( J1  *omegaBarSq(1) + J2  *omegaBarSq(2) + J3  *omegaBarSq(3) );
%             h  = sqrt( J1^2*omegaBarSq(1) + J2^2*omegaBarSq(2) + J3^2*omegaBarSq(3) );
%             t0 = this.calculateT0(p,constraintNum);
%             rbr = RigidBodyRotation([J1 J2]',this.RBtoW(:,:,1),[Ek h t0],'Ekht0',1);
%             RBtoWp = rbr.predictOrientation(this.t);
            % ----------------------------------------------------
            
            % Get the cost as sum_{i=1}^{N} || omegaB_Bp_{i} - omegaB_B_{i} ||^2_covOmegaB
            Cost = 0;
            for i = 1:size(this.omegaB_B,2)
                Cost = Cost + (omegaB_Bp(:,i)-this.omegaB_B(:,i))' * ( this.covOmegaB(:,:,i) \ ...
                                (omegaB_Bp(:,i)-this.omegaB_B(:,i)) );
                
                % TEMP, testing ---------------------------
%                 Cost = Cost + norm( Log(this.RBtoW(:,:,i)'*RBtoWp(:,:,i)) );
                % -----------------------------------------
                            
            end
            
        end
        %% Find the optimal p to minimize the cost function along the specified constraint.
        %  Return the optimal value and the residual of the fit
        function [pOpt,costOpt] = optimizeAlongConstraint(this,constraintNum)
                             
            % Find minimum and maximum viable parameters to use as lower and upper bounds
            [pMin,~] = this.inertiaConstraints.findClosestPoint([1.01, 1.01]',constraintNum);
            if strcmp(this.inertiaSymmetry,'AS1')
                pMax = 2;
            else
                [pMax,~] = this.inertiaConstraints.findClosestPoint([15 15]',constraintNum);
            end
                        
            % Find an initialization point using the angular velocity method
            [p0,J0] = this.initialize(constraintNum);
            this.p0(constraintNum) = p0;
            this.J0(:,constraintNum) = J0;
            
            % Perform optimization
            opts = optimset('TolFun',1e-13,'MaxFunEvals',5e3,'display','off');
            [pOpt,costOpt,~,~,~,~,~] = fmincon(@(x) this.cost(x,constraintNum), ...
                p0,[],[],[],[],pMin,pMax,[],opts);
            
            % Get object values for record of optimization
            this.cost0(constraintNum) = this.cost(p0,constraintNum);
            this.t00(constraintNum) = this.calculateT0(p0,constraintNum);
            this.pOpt(constraintNum) = pOpt;
            this.Jopt(:,constraintNum) = this.inertiaConstraints.ratios(pOpt,constraintNum);
            this.costOpt(constraintNum) = costOpt;
            this.t0opt(constraintNum) = this.calculateT0(pOpt,constraintNum);
            
            disp(['Constraint ' num2str(constraintNum) ' | costOpt: ' ...
                    num2str(this.costOpt(constraintNum)) ', Jopt: ' ...
                    num2str(this.Jopt(1,constraintNum)') ', ' ...
                    num2str(this.Jopt(2,constraintNum)) ', t0opt: ' ...
                    num2str(this.t0opt(constraintNum))]);
                        
        end
        %% Find the optimal p to minimize the cost function along all constraints.
        %  Return the optimal parameter value, the residual of the fit, and the index of the 
        %  optimal constraint.
        function [pOpt,costOpt,indOpt] = optimize(this)
            
            tic;
            
            % Optimize along each of the equality constraint(s)
            if strcmp(this.inertiaSymmetry,'TA')
                this.costOpt = [ Inf Inf Inf]';
                for i = 1:3
                    try
                        this.optimizeAlongConstraint(i);
                    catch err
                        disp(['Constraint ' num2str(i) ' | Error ' err.identifier ', skipping...']);
                    end
                end
            elseif strcmp(this.inertiaSymmetry,'AS1') || strcmp(this.inertiaSymmetry,'AS3')
                this.costOpt = [ Inf Inf Inf]';
                this.optimizeAlongConstraint(1);
            end
            
            % Get the best results
            [costOpt,indOpt] = min(this.costOpt);
            pOpt = this.pOpt(indOpt);
            this.indOpt = indOpt;
            
            disp('Inertia ratios optimization complete');
            disp(['Optimal Constraint: ' num2str(indOpt)]);
            
            toc
            
        end        
        %% Plot the cost function for different values of p for a given constraint
        function plotCostVsP(this,constraintNum,lineSpec)
            
            % Find minimum viable parameter
            [pMin,~] = this.inertiaConstraints.findClosestPoint([1.01, 1.01]',constraintNum);
            
            % Create a range of parameters to plot over
            if strcmp(this.inertiaSymmetry,'TA')
                pRange = pMin:0.1:pMin+pi/2;
            elseif strcmp(this.inertiaSymmetry,'AS1')
                pRange = pMin:0.01:2;
            elseif strcmp(this.inertiaSymmetry,'AS3')
                pRange = pMin:0.01:4;
            end
            
            % Get the cost at each parameter value p
            Cost = zeros(1,length(pRange));
            disp('Plotting Cost vs p ...');
            for i = 1:length(pRange)
                if strcmp(this.inertiaSymmetry,'TA')    % takes longer for TA, so show progress
                    disp(['Plotting Cost vs p Step ' num2str(i) ' of ' num2str(length(pRange))]);
                end
                Cost(i) = this.cost(pRange(i),constraintNum);
            end
            
            % Plot
            plot(pRange,Cost,lineSpec);
            
        end
        %% Plot the predicted magnitude of angular velocity parallel to angular momentum for 
        %  different values of parameter p and a given constraint.
        function plotOmegaParVsP(this,constraintNum,lineSpec)
            
            % Find the best estimate of actual parallel angular velocity
            omegaPar_W = this.calculateOmegaPar();
            omegaParMagMeas = norm(omegaPar_W);
            
            % Find minimum viable parameter
            [pMin,~] = this.inertiaConstraints.findClosestPoint([1.01, 1.01]',constraintNum);
            
            % Create a range of parameters to plot over
            if strcmp(this.inertiaSymmetry,'TA')
                pRange = pMin:0.01:pMin+pi/2;
            elseif strcmp(this.inertiaSymmetry,'AS1')
                pRange = pMin:0.01:2;
            elseif strcmp(this.inertiaSymmetry,'AS3')
                pRange = pMin:0.01:4;
            end
            
            % Get the omegaPar magnitude (= 2*Ek/h) at each parameter value p
            omegaParMag = zeros(1,length(pRange));
            omegaBarSq = 1/size(this.omegaB_B,2)*sum(this.omegaB_B.^2,2);
            for i = 1:length(pRange)
                disp(['Plotting omegaPar vs p Step ' num2str(i) ' of ' num2str(length(pRange))]);
                J = this.inertiaConstraints.ratios(pRange(i),constraintNum);
                J1 = J(1); J2 = J(2); J3 = 1;
                Ek = 1/2*( J1  *omegaBarSq(1) + J2  *omegaBarSq(2) + J3  *omegaBarSq(3) );
                h  = sqrt( J1^2*omegaBarSq(1) + J2^2*omegaBarSq(2) + J3^2*omegaBarSq(3) );
                omegaParMag(i) = 2*Ek/h;
            end
            
            % Plot
            plot(pRange,omegaParMag,lineSpec); hold on;
            plot([min(pRange) max(pRange)],[omegaParMagMeas omegaParMagMeas],'-k');
            
        end
        %% Get the optimal inertia ratios
        function J = getInertiaRatios(this)
            J = this.Jopt(:,this.indOpt);
        end
        %% Get the optimal cost
        function cost = getCost(this)
            cost = this.costOpt(this.indOpt);
        end
            
        
    end
    
end

