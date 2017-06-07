classdef RigidBodyRotation < handle
    %RIGIDBODYROTATION Analytic determination of the angular velocity and attitude of a passive 
    % rotating rigid body
    
    properties
        
        J;                  % inertia ratios [J1 J2 J3]'
        RB0toW;             % rotation matrix from initial body frame to world frame
        omegaB0_B;          % initial angular velcocity in the body frame [rad/s]
        Ek;                 % normalized kinetic energy [rad^2/s^2]
        h;                  % normalized angular momentum [rad/s]
        t0;                 % initial time when sn(t0) = 0, cn(t0) = 1, dn(t0) = 1 [s]. Will be -ve.
        createMode;         % either 'omega0' or 'Ekht0' depending on intialization type
        suppressPrint;      % boolean == 1 when it is desired to print nothing
        
        RHtoW;              % rotation matrix from angular momentum frame to world frame
        s;                  % signs of angular velocities [1 -1 1]' for low energy, [-1 1 1]' for high energy 
        T;                  % quarter-period of sn and cn, half-period of dn [s]
        omegaMax;           % the maximum angular velocities in each axis
        k;                  % the modulus of the elliptic functions
        omegaP;             % the body nutation rate
        energyState;        % string, 'LE' when h^2 > 2*Ek*J2, 'ME' when h^2 = 2*Ek*J2, 
                            % 'HE' when h^2 < 2*Ek*J2
        inertiaSymmetry;    % either 'TA' for tri-axial, 'AS1' or 'AS3' for axis-symmetric 1 and 2,
                            % 'FS' for full-symmetric (see logbook #3 pp 107)
        delA = 1e-4;        % tolerance in inertia ratios for symmetry and principal axes [-]
        isPrincipalRot;     % boolean == 1 when rotating about a principal axis
        delP = 0.0087;      % tolerance in angular velocity angle from a principal axis [rad]
        rotAxis;            % equals [0 0 0]' for multi-axis rotation, otherwise gives rotation axis
    end
    
    methods
        %% Constructor with createMode = 'omega0' or 'Ekht0'
        function this = RigidBodyRotation(ratios,RB0toW,omega0orEkht0,createMode,suppressPrint)
            this.J = [ratios; 1];
            this.RB0toW = RB0toW;
            if strcmp(createMode,'omega0')
                this.omegaB0_B = omega0orEkht0;
            else
                this.Ek = omega0orEkht0(1);
                this.h  = omega0orEkht0(2);
                this.t0 = omega0orEkht0(3);
            end
            this.createMode = createMode;
            this.suppressPrint = suppressPrint;
            this.updateInternalVariables();
        end
        %% Getters and setters
        function ratios = getInertiaRatios(this)
            ratios = this.J;
        end
        function RB0toW = getInitialOrientation(this)
            RB0toW = this.RB0toW;
        end
        function omegaB0_B = getInitialAngularVelocity(this)
            omegaB0_B = this.omegaB0_B;
        end
        function Ek = getKineticEnergy(this)
            Ek = this.Ek;
        end
        function h = getAngularMomentum(this)
            h = this.h;
        end
        function t0 = getInitialTime(this)
            t0 = this.t0;
        end
        function h_W = getAngularMomentumVector(this)
            h_W = this.RB0toW * diag(this.J) * this.omegaB0_B;
        end
        function setInertiaRatios(this, ratios)
            this.J = ratios;
            this.updateInternalVariables();
        end
        function setInitialOrientation(this,RB0toW)
            this.RB0toW = RB0toW;
            this.updateInternalVariables();            
        end
        function setInitialAngularVelocity(this,omegaB0_B)
            this.omegaB0_B = omegaB0_B;
            this.createMode = 'omega0';
            this.updateInternalVariables();
        end
        function setKineticEnergy(this,Ek)
            this.Ek = Ek;
            this.createMode = 'Ekht0';
            this.updateInternalVariables();
        end
        function setAngularMomentum(this,h)
            this.h = h;
            this.createMode = 'Ekht0';
            this.updateInternalVariables();
        end
        function setInitialTime(this,t0)
            this.t0 = t0;
            this.createMode = 'Ekht0';
            this.updateInternalVariables();
        end
        %% Predict the body angular velocities at time t
        %  Can accept row vectors t
        function omegaB_B = predictOmega(this, t)
            
            nT = length(t);
            
            if sum(this.rotAxis) ~= 0                   % single-axis
                omegaB_B = repmat(this.omegaB0_B,1,length(t));
            elseif strcmp(this.inertiaSymmetry,'AS1')   % multi-axis, axis-symmetric J2==J3
                omegaB_B = [    repmat(this.s(1) * this.omegaMax(1),1,nT)     ;
                                this.s(2) * this.omegaMax(2) * sin(this.omegaP*(t-this.t0)) ;
                                this.s(3) * this.omegaMax(3) * cos(this.omegaP*(t-this.t0)) ];
            elseif strcmp(this.inertiaSymmetry,'AS3')   % multi-axis, axis-symmetric J1==J1
                omegaB_B = [    this.s(1) * this.omegaMax(1) * sin(this.omegaP*(t-this.t0)) ;
                                this.s(2) * this.omegaMax(2) * cos(this.omegaP*(t-this.t0)) ;
                                repmat(this.s(3) * this.omegaMax(3),1,nT)     ];
            elseif strcmp(this.inertiaSymmetry,'TA')    % multi-axis, tri-axial
                % Get the elliptic integral terms
                m = this.k^2;
                [sn,cn,dn] = ellipj(this.omegaP*(t-this.t0), m);

                % Get the angular velocity at time t
                if strcmp(this.energyState,'LE') || strcmp(this.energyState,'ME')
                    omegaB_B = [ this.s(1) * this.omegaMax(1) * dn ;
                                 this.s(2) * this.omegaMax(2) * sn ;
                                 this.s(3) * this.omegaMax(3) * cn ];
                else        
                    omegaB_B = [ this.s(1) * this.omegaMax(1) * cn ;
                                 this.s(2) * this.omegaMax(2) * sn ;
                                 this.s(3) * this.omegaMax(3) * dn ];
                end
            end
            
        end
        %% Calculate the orientation in the world frame RBtoW at time t, where t is a 1xnT row vector
        function RBtoW = predictOrientation(this, t)
            
            % Get length of times and initialize rotation matrices
            nT = length(t);
            RBtoW = zeros(3,3,nT);
            
            if sum(this.rotAxis) ~= 0  % single-axis
                for i = 1:nT
                    RBtoW(:,:,i) = this.RB0toW * Exp(this.omegaB0_B * t(i));
                end
            else                       % multi-axis
                RHtoB = this.calculateRHtoB(t);
               
                for i = 1:nT
                    RBtoW(:,:,i) = this.RHtoW * RHtoB(:,:,i)';
                end
            end
            
        end        
        %% Calculate the rotation matrix RHtoB at time(s) t, where t is a 1xnT row vector
        function RHtoB = calculateRHtoB(this, t)
            
            % Number of timesteps
            nT = length(t);
                
            % Get the current and maximum angular velocity
            omegaB_B = this.predictOmega(t);
            omega1 = omegaB_B(1,:); omega2 = omegaB_B(2,:); omega3 = omegaB_B(3,:);
            omega1m = this.omegaMax(1); omega2m = this.omegaMax(2); omega3m = this.omegaMax(3);
            
            % Get inertia tensor, angular momentum, angular momentum components and their maxes
            J1 = this.J(1); J2 = this.J(2); J3 = this.J(3); h = this.h;
            h1 = J1*omega1; h2 = J2*omega2; h3 = J3*omega3;
            h1m = J1*omega1m; h2m = J2*omega2m; h3m = J3*omega3m;             
            
            if strcmp(this.inertiaSymmetry,'TA')
                % Determine number of cycles of phi for finding alpha
                numCycles = floor( (t-this.t0)./(4*this.T) );   

                % Get the modulus
                m = this.k^2;

                % Get elliptic functions for finding alpha
                [snUm,cnUm,~] = ellipj(this.omegaP*(t-this.t0), m);

                % Get values of phi at time=t for finding alpha, using quadrants to properly
                % determine phiM. Note that asin returns values [-pi/2, pi/2]
                phiM = zeros(1,length(t));
                for i = 1:length(t)
                    if snUm(i) >= 0 && cnUm(i) >= 0                                     % Q1
                        phiM(i) = numCycles(i)*2*pi + asin(snUm(i));
                    elseif snUm(i) >= 0 && cnUm(i) < 0 || snUm(i) < 0 && cnUm(i) < 0    % Q2 or Q3
                        phiM(i) = numCycles(i)*2*pi + pi - asin(snUm(i));
                    else                                                                % Q4
                        phiM(i) = numCycles(i)*2*pi + 2*pi + asin(snUm(i));
                    end
                end

                % Determine rotation matrix based on energy level
                if strcmp(this.energyState,'LE')                    % Low energy case --------------

                    % Get values for calculations using Hurtado 2011 equations
                    A  = sqrt( 1 - (h1./h).^2 );
                    a0 = (J2^2*J3 - J2*J3^2)/(J1*J2 - J1*J3);                
                    n  = (J2^2*omega2m^2 - J3^2*omega3m^2) / (J3^2*omega3m^2);
                    ellint3 = 1/this.omegaP*ellipticPi(-n,phiM,m);  % see testEllipticPi.m and formula sheets                
                    alpha   = h/(J2*J3) .* ( a0*(t-this.t0) + (J2-a0)*(ellint3) );

                    % Form the rotation matrix from angular momentum frame to body frame
                    % R1 = [  h1/h     A               0              ;
                    %         h2/h    -h1*h2/(A*h^2)   h3/(A*h)       ;
                    %         h3/h    -h1*h3/(A*h^2)  -h2/(A*h)       ];
                    % R2 = [  1        0              0               ;
                    %         0        cos(alpha)     sin(alpha)      ;
                    %         0       -sin(alpha)     cos(alpha)      ];
                    R1col1 = [  h1./h       ;  h2./h                ;  h3./h                ];
                    R1col2 = [  A           ; -h1.*h2./(A.*(h.^2)) 	; -h1.*h3./(A.*(h.^2))  ];
                    R1col3 = [  zeros(1,nT) ;  h3./(A.*h)           ; -h2./(A.*h)           ];
                    R2col1 = [  ones(1,nT)  ;  zeros(1,nT)          ;  zeros(1,nT)          ];
                    R2col2 = [  zeros(1,nT) ;  cos(alpha)           ; -sin(alpha)       	];
                    R2col3 = [  zeros(1,nT) ;  sin(alpha)           ;  cos(alpha)           ];

                elseif strcmp(this.energyState,'ME')                % Medium enenrgy case -----------

                    % Get values for calculations using Hurtado 2011 equations (logbook #3 pp 108)                             
                    alpha   = (J1*h3m^2 + J3*h1m^2) ./ (J1*J3*h) .* (t - this.t0);

                    % Form the rotation matrix from angular momentum frame to body frame
                    % R1 = [   h3m/h          h1/h           -h1m*h2/h^2  ;
                    %          0              h2/h            h1/h1m      ;
                    %          h1m/h          h3/h            h2*h3m/h^2  ];
                    % R2 = [   cos(alpha)     0              -sin(alpha)  ;
                    %          0              1               0           ;
                    %          sin(alpha)     0               cos(alpha)  ];   
                    R1col1 = [   repmat(h3m/h,1,nT)     ;  zeros(1,nT)          ;  repmat(h1m/h,1,nT) ];
                    R1col2 = [   h1./h                  ;  h2./h                ;  h3./h              ];
                    R1col3 = [  -h1m*h2/h^2             ;  h1/h1m               ;  h2*h3m/h^2         ];
                    R2col1 = [   cos(alpha)             ;  zeros(1,nT)          ;  sin(alpha)         ];
                    R2col2 = [   zeros(1,nT)            ;  ones(1,nT)           ;  zeros(1,nT)        ];
                    R2col3 = [  -sin(alpha)             ;  zeros(1,nT)          ;  cos(alpha)         ];

                elseif strcmp(this.energyState,'HE')                % High energy case --------------

                    % Get values for calculations using Hurtado 2011 equations
                    A  = sqrt( 1 - (h3./h).^2 ); 
                    a0 = (J1^2*J2 - J1*J2^2)/(J1*J3 - J2*J3);
                    n  = (J2^2*omega2m^2 - J1^2*omega1m^2) / (J1^2*omega1m^2);
                    ellint3 = 1/this.omegaP*ellipticPi(-n,phiM,m);  % see testEllipticPi.m and formula sheets    
                    alpha   = h/(J1*J2) .* ( a0*(t-this.t0) + (J2-a0)*(ellint3) );

                    % Form the rotation matrix from angular momentum frame to body frame
                    % R1 = [  -h1*h3/(A*h^2)  h2/(A*h)        h1/h    ;
                    %         -h2*h3/(A*h^2)  -h1/(A*h)       h2/h    ;
                    %          A              0               h3/h    ];
                    % R2 = [   cos(alpha)     sin(alpha)      0       ;
                    %         -sin(alpha)     cos(alpha)      0       ;
                    %         0               0               1       ];  
                    R1col1 = [  -h1.*h3./(A.*(h.^2))    ; -h2.*h3./(A.*(h.^2))  ;  A            ];
                    R1col2 = [   h2./(A.*h)             ; -h1./(A.*h)           ;  zeros(1,nT)  ];
                    R1col3 = [   h1./h                  ;  h2./h                ;  h3./h        ];
                    R2col1 = [   cos(alpha)             ; -sin(alpha)           ;  zeros(1,nT)  ];
                    R2col2 = [   sin(alpha)             ;  cos(alpha)           ;  zeros(1,nT)  ];
                    R2col3 = [   zeros(1,nT)            ; zeros(1,nT)           ;  ones(1,nT)   ]; 

                end

                % Calculate the desired rotation matrix                    
                R1     = reshape([R1col1; R1col2; R1col3], 3, 3, nT);  
                R2     = reshape([R2col1; R2col2; R2col3], 3, 3, nT);
                RHtoB = zeros(3,3,nT);
                for i = 1:nT
                    RHtoB(:,:,i) = R1(:,:,i)*R2(:,:,i);              
                end
                
            elseif strcmp(this.inertiaSymmetry,'AS1') || strcmp(this.inertiaSymmetry,'AS3')
                
                % Get the axial and transverse angular momenta, transverse inertia ratio, and 
                % quaternion index order
                if strcmp(this.inertiaSymmetry,'AS1')
                    ha0 = h1m; ht0 = h3m; Jt = J3; qInd = [2 3 1 4]';
                elseif strcmp(this.inertiaSymmetry,'AS3')
                    ha0 = h3m; ht0 = h1m; Jt = J1; qInd = [1 2 3 4]';
                end
                
                % Get angles [Hughes 2004 pp 96]
                mu = this.omegaP * (t-this.t0);
                gamma = atan(ht0/ha0);              % nutation angle
                
                % Calculate the quaternion solution [Hughes 2004 pp 104]. Define H frame to start 
                % with lambda0 = 0
                qHtoB = zeros(4,nT);
                RHtoB = zeros(3,3,nT);
                for i = 1:nT
                    lambda = 0.5 * (h/Jt - this.omegaP) * (t(i)-this.t0);                    
                    qHtoB(qInd(1),i) = sin(gamma/2) * cos(lambda);
                    qHtoB(qInd(2),i) = sin(gamma/2) * sin(lambda);
                    qHtoB(qInd(3),i) = cos(gamma/2) * sin(mu(i) + lambda);
                    qHtoB(qInd(4),i) = cos(gamma/2) * cos(mu(i) + lambda);                    
                	RHtoB(:,:,i) = quat2rot(qHtoB(:,i));
                end
                
            end
            
        end
        %% Calculate the rotation from the angular momentum frame H to the world frame W
        function RHtoW = calculateRHtoW(this)
            
            % Calculate as in [Hurtado2014 Sec V]
            RHtoW = this.RB0toW * this.calculateRHtoB(0);
            
        end        
        %% Calculate the time in the past at which sn(t0) = 0, cn(t0) = 1, dn(t0) = 1 [s]. Will be -ve.
        function t0 = calculateT0(this)
            
            % Put initial angular velocity into canonical form, with domain [-1, 1], and
            % signs appropriate to the canonical versions of sn, cn, and dn
            omega0_Bc = [ this.s(1)*this.omegaB0_B(1)/this.omegaMax(1); ...
                          this.s(2)*this.omegaB0_B(2)/this.omegaMax(2); ...
                          this.s(3)*this.omegaB0_B(3)/this.omegaMax(3)  ];
            omega0_Bc(omega0_Bc > 1) = 1; omega0_Bc(omega0_Bc < -1) = -1; % fix numerical issues
            
            if strcmp(this.inertiaSymmetry,'TA')
                
                % Invert the Jacobi elliptic functions. First determine which angular
                % velocities correspond to the sn and cn curves respectively
                omega0_Bc_sn = omega0_Bc(2);
                if strcmp(this.energyState,'LE')
                    omega0_Bc_cn = omega0_Bc(3);
                elseif strcmp(this.energyState,'ME')
                    omega0_Bc_cn = omega0_Bc(3); % could also use omega0_Bc(1) since both are sech
                elseif strcmp(this.energyState,'HE')
                    omega0_Bc_cn = omega0_Bc(1);
                end

                % Get the quarter-period of the elliptic functions in terms of u (i.e. sn(u, m))
                Tu = this.T * this.omegaP;

                % Get the value asn, which will be between -Tu and Tu
                m = this.k^2;
                u = ellipticF(asin(omega0_Bc_sn), m);
                
                % Determine which quadrant the resulting solution is in
                if omega0_Bc_sn >= 0 && omega0_Bc_cn >= 0         % Q1
                    u = u;
                elseif (omega0_Bc_sn >= 0 && omega0_Bc_cn < 0)    ...
                    || (omega0_Bc_sn <  0 && omega0_Bc_cn < 0)    % Q2 or Q3
                    u = 2*Tu - u;
                else                                              % Q4
                    u = 4*Tu + u;
                end
                
                % The magnitude of t0 will then be u/omegaP, but make it negative, so that
                % our experiment starts at time=0, rather than at t0
                t0 = -u/this.omegaP;                
                
            elseif strcmp(this.inertiaSymmetry,'AS1') || strcmp(this.inertiaSymmetry,'AS3')
                
                % Determine which axis is driven by sin and which is driven by cos
                if strcmp(this.inertiaSymmetry,'AS1')
                    omega0_Bc_sin = omega0_Bc(2); omega0_Bc_cos = omega0_Bc(3);
                elseif strcmp(this.inertiaSymmetry,'AS3')
                    omega0_Bc_sin = omega0_Bc(1); omega0_Bc_cos = omega0_Bc(2);                    
                end
                
                % Get the angle for the cananonical sin curve
                theta = asin(omega0_Bc_sin);
                
                % Correct for quadrant
                if omega0_Bc_sin >= 0 && omega0_Bc_cos >= 0         % Q1
                    theta = theta;
                elseif (omega0_Bc_sin >= 0 && omega0_Bc_cos < 0)    ...
                    || (omega0_Bc_sin <  0 && omega0_Bc_cos < 0)    % Q2 or Q3
                    theta = pi - theta;
                else                                                % Q4
                    theta = 2*pi + theta;
                end            
                
                t0 = -theta/this.omegaP;                
            end
            
        end
        %% Calculate the quarter-period of sn and cn, the half period of dn [s]. 
        %  I.e. Sn and cn are periodic every 4*T secs, and dn is periodic every 2*T secs
        function T = calculateT(this)
            
            if strcmp(this.inertiaSymmetry,'TA')
                
                % Get the modulus
                m = this.k^2;

                % Get the period in terms of u (e.g. sn(u, m))
                Tu = ellipticF(pi/2, m);

                % Get the period in terms of time t (e.g. sn(omegaP*t, m)) [s]
                T = Tu/this.omegaP;
                
            elseif strcmp(this.inertiaSymmetry,'AS1') || strcmp(this.inertiaSymmetry,'AS3')
                
                % Get the period as the inverse of frequency
                T = 2*pi*abs(1/this.omegaP);
                            
            end
            
        end
        %% Update all internal variables based on current inertia ratios, angular momentum, and 
        %  initial orientation
        function updateInternalVariables(this)
                        
            % Calculate kinetic energy and ang mom if initialized using initial angular velocity
            if strcmp(this.createMode,'omega0')
                this.Ek = 0.5 * this.omegaB0_B' * diag(this.J) * this.omegaB0_B;
                this.h  = norm( diag(this.J) * this.omegaB0_B );                
            end
            
            % Get the values of J1, J2, and J3
            J1 = this.J(1); J2 = this.J(2); J3 = this.J(3);
            
            % Determine whether the body is in the low, medium or high energy state
            if this.h^2 - 2*this.Ek*J2 > 1e-10
                this.energyState = 'LE';
            elseif abs(this.h^2 - 2*this.Ek*J2) <= 1e-10
                this.energyState = 'ME';                
            elseif this.h^2 - 2*this.Ek*J2 < -1e-10
                this.energyState = 'HE';
            end
            
            % Determine whether the inertias are tri-axial, axis-symmetric 1, axis-symmetric 3,
            % or fully-symmetric (see logbook #3 pp 107)
            if J1 >= J2 + this.delA && J2 >= J3 + this.delA                  
                this.inertiaSymmetry = 'TA';            % tri-axial
            elseif J1 >= J2 + this.delA && J1 >= J3 + this.delA && J2-J3 <= this.delA 
                this.inertiaSymmetry = 'AS1';           % axis-symmetric 1
                Ja = J1; Jt = J3;                       % axial and transverse inertias
                omegaA = sqrt( (this.h^2-2*Jt*this.Ek)/(Ja^2-Jt*Ja) ); % axial ang vel
                omegaT = sqrt( (this.h^2-2*Ja*this.Ek)/(Jt^2-Jt*Ja) ); % transverse ang vel
            elseif J1-J2 <= this.delA && J1 >= J3 + this.delA && J2 >= J3 + this.delA
                this.inertiaSymmetry = 'AS3';           % axis-symmetric 3
                Ja = J3; Jt = J1;                       % axial and transverse inertias 
                omegaA = sqrt( (this.h^2-2*Jt*this.Ek)/(Ja^2-Jt*Ja) ); % axial ang vel
                omegaT = sqrt( (this.h^2-2*Ja*this.Ek)/(Jt^2-Jt*Ja) ); % transverse ang vel
            elseif J1-J2 <= this.delA && J2-J3 <= this.delA && J1-J3 <= this.delA
                this.inertiaSymmetry = 'FS';            % fully-symmetric
            end
                        
            % Get the signs by convention
            if strcmp(this.inertiaSymmetry,'TA')
                if strcmp(this.energyState,'LE')
                    this.s = [1 -1 1]';                    
                elseif strcmp(this.energyState,'ME')
                    % Since sech > 0 and tanh > 0, sign must mirror omegaB0 when created with 
                    % omega0. If created with Ekht0, use ME convention from Hurtado 2011
                    if strcmp(this.createMode,'omega0')
                        this.s = sign(this.omegaB0_B);
                        this.s(this.s==0) = 1;    % default to +ve sign if omegaB0_B is 0 in an axis
                    elseif strcmp(this.createMode,'Ekht0')
                        this.s = [1 1 -1]';
                    end
                elseif strcmp(this.energyState,'HE')
                    this.s = [-1 1 1]';
                end
            else
                this.s = [1 1 1]';
            end
            
            % Get maximum angular velocities
            this.omegaMax = zeros(3,1);
            if strcmp(this.inertiaSymmetry,'TA')        % tri-axial
                this.omegaMax(1) = sqrt( (this.h^2-2*this.Ek*J3) / (J1*(J1-J3)) );
                if strcmp(this.energyState,'LE') || strcmp(this.energyState,'ME')
                    this.omegaMax(2) = sqrt( (2*this.Ek*J1-this.h^2) / (J2*(J1-J2)) );
                elseif strcmp(this.energyState,'HE');
                    this.omegaMax(2) = sqrt( (this.h^2-2*this.Ek*J3) / (J2*(J2-J3)) );
                end
                this.omegaMax(3) = sqrt( (2*this.Ek*J1-this.h^2) / (J3*(J1-J3)) );
            elseif strcmp(this.inertiaSymmetry,'AS1')   % J2==J3                    
                    this.omegaMax(1) = omegaA;
                    this.omegaMax(2) = omegaT;
                    this.omegaMax(3) = omegaT;
            elseif strcmp(this.inertiaSymmetry,'AS3')   % J1==J2                  
                    this.omegaMax(1) = omegaT;
                    this.omegaMax(2) = omegaT;
                    this.omegaMax(3) = omegaA;
            end

            % Get modulus for elliptic functions
            if strcmp(this.inertiaSymmetry,'TA')                    
                a = sqrt( (this.h^2-2*this.Ek*J3) / (J2*(J2-J3)) );
                b = sqrt( (2*this.Ek*J1-this.h^2) / (J2*(J1-J2)) );                
                if strcmp(this.energyState,'LE')
                    this.k = b/a;
                elseif strcmp(this.energyState,'ME')
                    this.k = 1;
                elseif strcmp(this.energyState,'HE')
                    this.k = a/b;
                end
                this.k = real(max(0,min(1,this.k)));
            end
            
             % Get body nutation rate
            if strcmp(this.inertiaSymmetry,'TA')
                if strcmp(this.energyState,'LE') || strcmp(this.energyState,'ME')
                    this.omegaP = sqrt( (this.h^2-2*this.Ek*J3)*(J1-J2) / (J1*J2*J3) );
                else
                    this.omegaP = sqrt( (2*this.Ek*J1-this.h^2)*(J2-J3) / (J1*J2*J3) );
                end
            elseif strcmp(this.inertiaSymmetry,'AS1') || strcmp(this.inertiaSymmetry,'AS3')
                this.omegaP = (Jt-Ja)/Jt * omegaA;
            end
            
            if ~strcmp(this.inertiaSymmetry,'FS')
                
                % Get quarter-period for sn and cn, half-period for dn, or one period of sin, cos (AS)
                this.T = this.calculateT();
                        
                % Calculate initial time if not previously specified
                if strcmp(this.createMode,'omega0')
                   this.t0 = this.calculateT0();
                end
                
            end
            
            % Detect principal axis rotations from when the Jt is close to the inertia ratios. 
            % In the case of intermediate axis, the time also must have progressed to the point 
            % where it is 'close' to the axis. Find rotation axis and initial angular velocity (if
            % required).
            if ~strcmp(this.inertiaSymmetry,'FS') && abs(J1-this.h^2/(2*this.Ek)) < this.delA ...
                    && ~strcmp(this.inertiaSymmetry,'AS3')
                this.isPrincipalRot = 1;    
                this.rotAxis = [1 0 0]';
                if strcmp(this.createMode,'Ekht0')
                    this.omegaB0_B = this.rotAxis * sqrt(this.Ek/J1);
                end
            elseif ~strcmp(this.inertiaSymmetry,'FS') ...
                    && abs(J2-this.h^2/(2*this.Ek)) < this.delA && abs(this.t0)/this.T > 0.1 
                this.isPrincipalRot = 1;    
                this.rotAxis = [0 1 0]';
                if strcmp(this.createMode,'Ekht0')
                    this.omegaB0_B = this.rotAxis * sqrt(this.Ek/J2);
                end
            elseif ~strcmp(this.inertiaSymmetry,'FS') && ~strcmp(this.inertiaSymmetry,'AS1') ...
                    && abs(J3-this.h^2/(2*this.Ek)) < this.delA                    
                this.isPrincipalRot = 1;    
                this.rotAxis = [0 0 1]';
                if strcmp(this.createMode,'Ekht0')
                    this.omegaB0_B = this.rotAxis * sqrt(this.Ek/J3);
                end
            else
                this.isPrincipalRot = 0;               
            end
            
            % When the rotation is not about a principal axis, determine rotation axis ([0 0 0]'
            % if multi-axis), and angular velocity (if required).
            if ~this.isPrincipalRot
                if strcmp(this.inertiaSymmetry,'TA')            % TA -----------------------------
                    this.rotAxis = [0 0 0]';
                    if strcmp(this.createMode,'Ekht0')
                        this.omegaB0_B = this.predictOmega(0);
                    end
                elseif strcmp(this.inertiaSymmetry,'AS1')       % AS1 J2=J3 ----------------------
                    if strcmp(this.energyState,'LE')
                        this.rotAxis = [0 0 0]';
                        if strcmp(this.createMode,'Ekht0')
                            this.omegaB0_B = this.predictOmega(0);
                        end
                    elseif strcmp(this.energyState,'ME')
                        if strcmp(this.createMode,'Ekht0')
                            disp(['Warning: Axial-Symmetry 1 ME rotation axis undefined,' ...
                                'using 1/sqrt(2)*[0 1 1]''.']);
                            this.rotAxis = 1/sqrt(2)*[0 1 1]';
                            this.omegaB0_B = this.rotAxis * sqrt(2*this.Ek/mean(this.J(2:3)));
                        elseif strcmp(this.createMode,'omega0')
                            this.rotAxis = this.omegaB0_B/norm(this.omegaB0_B);
                        end
                    elseif strcmp(this.energyState,'HE')
                        error('Error: High Energy Axial-Symmetry 1 is not a viable option.');
                    end
                elseif strcmp(this.inertiaSymmetry,'AS3')       % AS3 J1=J2 ----------------------
                    if strcmp(this.energyState,'LE')
                        error('Error: Low Energy Axial-Symmetry 2 is not a viable option.');
                    elseif strcmp(this.energyState,'ME')
                        if strcmp(this.createMode,'Ekht0')
                            disp(['Warning: Axial-Symmetry 2 ME rotation axis undefined,' ...
                                'using 1/sqrt(2)*[1 1 0]''.']);
                            this.rotAxis = 1/sqrt(2)*[1 1 0]';
                            this.omegaB0_B = this.rotAxis * sqrt(2*this.Ek/mean(this.J(1:2)));
                        elseif strcmp(this.createMode,'omega0')
                            this.rotAxis = this.omegaB0_B/norm(this.omegaB0_B);
                        end
                    elseif strcmp(this.energyState,'HE')
                        this.rotAxis = [0 0 0]';
                        if strcmp(this.createMode,'Ekht0')
                            this.omegaB0_B = this.predictOmega(0);
                        end
                    end
                elseif strcmp(this.inertiaSymmetry,'FS')        % FS -----------------------------
                    if strcmp(this.createMode,'Ekht0')
                        disp(['Warning: Fully-symmetric rotation axis undefined,' ...
                            'using 1/sqrt(3)*[1 1 1]''.']);
                        this.rotAxis = 1/sqrt(3)*[1 1 1]';
                        this.omegaB0_B = this.rotAxis * sqrt(2*this.Ek/mean(this.J));
                    elseif strcmp(this.createMode,'omega0')
                        this.rotAxis = this.omegaB0_B/norm(this.omegaB0_B);
                    end
                end
            end

            if ~strcmp(this.inertiaSymmetry,'FS')
                
                % Get rotation from angular momentum frame to world frame
                this.RHtoW = this.calculateRHtoW();
                
            end
            
            if ~this.suppressPrint
                disp('RigidBodyRotation -------------------------');
                disp(['energyState     = ' this.energyState]);
                disp(['inertiaSymmetry = ' this.inertiaSymmetry]);
                disp(['isPrincipalRot  = ' num2str(this.isPrincipalRot)]);
                disp(['rotAxis    = ' num2str(this.rotAxis')]);
                disp(['[J1 J2 J3] = ' num2str(this.J')]);
                disp(['omegaP     = ' num2str(this.omegaP) ' rad/s']);
                disp(['T          = ' num2str(this.T) ' s']);
                disp(['t0         = ' num2str(this.t0) ' s']);
                disp(['omegaMax   = ' num2str(this.omegaMax') ' rad/s']);
                disp(['omegaB0_B  = ' num2str(this.omegaB0_B') ' rad/s']);
                disp('RHtoW      = [');
                disp(this.RHtoW);
                disp(']');
            end
                        
        end
        %% Print values
        function print(this)
            disp(['J: [' num2str(this.J') ']']);
            disp(['omegaB0_W: [' num2str(this.omegaB0_B') ']']);
            disp('RB0toW: [');
            disp(this.RB0toW);
            disp(']');
        end
        %% Check if equal
        function isEqual = equals(this, rb, tol)
            inertiaRatiosEqual = norm(this.J-rb.J) < tol;
            angularVelocityEqual = max(abs(this.omegaB0_B-rb.omegaB0_B)) < tol;
            initialOrientationEqual = max(max(abs(this.RB0toW-rb.RB0toW))) < tol;
            isEqual = inertiaRatiosEqual & angularVelocityEqual & initialOrientationEqual;
        end
        %% Plot the two intersecting ellipsoids and the polhode
        function [xE,yE,zE,xH,yH,zH] = plotPolhode(this,hAlpha,eAlpha,suppressPlot)
            
            % Get energy, angular momentum and inertia ratios
            Ek = this.Ek; h = this.h; J1 = this.J(1); J2 = this.J(2); J3 = this.J(3);
            
            % Get energy ellipsoid E in the body frame
            aE = sqrt(2*Ek/J1); bE = sqrt(2*Ek/J2); cE = sqrt(2*Ek/J3);
            [xE,yE,zE] = ellipsoid(0,0,0,aE,bE,cE,100);
            
            % Get angular momentum ellipsoid H in the body frame
            aH = h/J1; bH = h/J2; cH = h/J3;
            [xH,yH,zH] = ellipsoid(0,0,0,aH,bH,cH,100);
            
            % Get polhode for full period
            omegaB_B = this.predictOmega(linspace(0,4*this.T,1000));
            
            % Uncomment for misaligning
%             RBtoG = Exp([pi/4 0 pi/4]);
%             omegaB_B = RBtoG*omegaB_B;
            
            % Set the default transparencies of the ellipsoids 
            if nargin < 3
                hAlpha = 1; eAlpha = 1;
            end
            
            % Plot
            if nargin < 4 || ~suppressPlot
%                 surf(xH,yH,zH,'FaceAlpha',hAlpha,'EdgeAlpha',0.1,'FaceColor',[0.9 0.9 0.9]);
                hold on; axis equal; view(135,30); grid on;
%                 surf(xE,yE,zE,'FaceAlpha',eAlpha,'EdgeAlpha',0.1,'FaceColor',[0 0.8 1]); 
                plot3(omegaB_B(1,:),omegaB_B(2,:),omegaB_B(3,:),'-k','LineWidth',1); hold on;
%                 plot3(omegaB_B(1,1),omegaB_B(2,1),omegaB_B(3,1),'*g','MarkerSize',15,'LineWidth',4);
                plot3([0 max([aE aH])+0.25],[0 0],[0 0],'-r','LineWidth',2);   % x-axis
                plot3([0 0],[0 max([bE bH])+0.25],[0 0],'-g','LineWidth',2);   % y-axis
                plot3([0 0],[0 0],[0 max([cE cH])+0.25],'-b','LineWidth',2);   % z-axis
                xlabel('\omega_1'); ylabel('\omega_2'); zlabel('\omega_3');                 
                % Uncomment for projections
                set(gca,'XTickLabel','');
                set(gca,'YTickLabel','');
                set(gca,'ZTickLabel','');
                xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([-1.5 1.5]);
                plot3(omegaB_B(1,:),omegaB_B(2,:),-1.5*ones(1,length(omegaB_B(3,:))),'-m','LineWidth',1); 
                plot3(omegaB_B(1,:),-1.5*ones(1,length(omegaB_B(2,:))),omegaB_B(3,:),'-m','LineWidth',1); 
                plot3(-1.5*ones(1,length(omegaB_B(1,:))),omegaB_B(2,:),omegaB_B(3,:),'-m','LineWidth',1);
%                 xB_G = (max([aE aH])+0.25)*RBtoG(:,1); yB_G = (max([bE bH])+0.25)*RBtoG(:,2); 
%                 zB_G = (max([cE cH])+0.25)*RBtoG(:,3); 
%                 plot3([0 xB_G(1)],[0 xB_G(2)],[0 xB_G(3)],'r','LineWidth',2);   % x-axis
%                 plot3([0 yB_G(1)],[0 yB_G(2)],[0 yB_G(3)],'g','LineWidth',2);   % y-axis
%                 plot3([0 zB_G(1)],[0 zB_G(2)],[0 zB_G(3)],'b','LineWidth',2);   % z-axis
            end
            
        end
        %% Plot the ellipsoid, invariable plane, polhode, and herpelhode from time t=0 to time t
        % publish the plot to figure(figNum)
        function plotEllipsoid(this, t, figNum)
            
            disp('in plotEllipsoid 1');
            % Setup the time interval to plot over
            dt = 0.1;
            tP = 0:dt:t;
            nT = length(tP);
            
            % Get the angular velocities and orientations
            omegaB_B = this.predictOmega(tP);   omegaB_W = zeros(3,nT);
            RBtoW = this.predictOrientation(tP, omegaB_B);            
            disp('in plotEllipsoid 2');
            for i = 1:nT
                omegaB_W(:,i) = RBtoW(:,:,i) * omegaB_B(:,i);
            end
                        
            disp('in plotEllipsoid 3');
            
            % Get the normal of the invariable plane, which is coincident with the angular momentum
            h_W = this.getAngularMomentum();
            nIP_W = h_W/norm(h_W);
            
            % Get the distance and vector of the of the body omega axes origin from the invariable plane
            omegaParallel = 2*this.Ek/norm(h_W);
            O = omegaParallel*nIP_W;
            
            % Get the semi-major axes of the ellipsoid, and a generic ellipsoid mesh
            J = this.inertiaRatios.matrix(); J1 = J(1,1); J2 = J(2,2); J3 = J(3,3);
            aE = sqrt(2*this.Ek/J1); bE = sqrt(2*this.Ek/J2); cE = sqrt(2*this.Ek/J3);
            nE = 100;
            [xEg,yEg,zEg] = ellipsoid(0,0,0,aE,bE,cE,nE);
            
            
            % Rotate the ellipsoid to its current position and translate its center to O
            xE = zeros(nE,nE); yE = zeros(nE,nE); zE = zeros(nE,nE);
            for i = 1:nE+1
                for j = 1:nE+1
                    omegaRot = RBtoW(:,:,end) * [xEg(i,j) yEg(i,j) zEg(i,j)]';
                    xE(i,j) = omegaRot(1); 
                    yE(i,j) = omegaRot(2);
                    zE(i,j) = omegaRot(3);
                end
            end
            
            % Rotate and translate the polhode into its current position
            omegaB_Bc = RBtoW(:,:,end) * omegaB_B;
            
            % Create 4 points for plotting the invariable plane. Base x-axis on the intial contact 
            % point of omega and the invariable plane, P0.
            P0 = omegaB_W(:,1);
            xIP = (O-P0)/norm(O-P0); yIP = cross(-O,xIP)/norm(cross(-O,xIP)); zIP = -O/norm(-O);
            pt1 = O + (cE+0.1)*(xIP+yIP);   pt2 = O + (cE+0.1)*(-xIP+yIP); 
            pt3 = O + (cE+0.1)*(-xIP-yIP);  pt4 = O + (cE+0.1)*(xIP-yIP);
            
            % Get a scaled version of the body axes at the current timestep
            axis1 = aE*RBtoW(:,1,end); axis2 = bE*RBtoW(:,2,end); axis3 = cE*RBtoW(:,3,end);
            
            % So far, everything has been calculated in coordinates with O at the origin. For
            % visual purposes, it is desired that the invariable plane be the x-y plane and the
            % coordinates of O be [0 0 omegaParallel]. The block of code below performs the required
            % rotations and translations to put everything into the invariable plane coordinates,
            % defined by xIP, yIP, zIP.
            RWtoIP = [xIP yIP zIP];
            O = [0 0 omegaParallel]';
            for i = 1:nE+1
                for j = 1:nE+1
                    omegaRot = RWtoIP(:,:,end) * [xE(i,j) yE(i,j) zE(i,j)]';
                    xE(i,j) = omegaRot(1) + O(1); 
                    yE(i,j) = omegaRot(2) + O(2);
                    zE(i,j) = omegaRot(3) + O(3);
                end
            end
            omegaB_Bc = RWtoIP * omegaB_Bc + repmat(O,1,nT);
            pt1 = RWtoIP*pt1 + O; pt2 = RWtoIP*pt2 + O; pt3 = RWtoIP*pt3 + O; pt4 = RWtoIP*pt4 + O;
            axis1 = RWtoIP*axis1; axis2 = RWtoIP*axis2; axis3 = RWtoIP*axis3;
            omegaB_W = RWtoIP * omegaB_W + repmat(O,1,nT);
            
            % Plot the invariable plane
            figure(figNum); clf;
            fill3([pt1(1) pt2(1) pt3(1) pt4(1)],[pt1(2) pt2(2) pt3(2) pt4(2)], ...;
                [pt1(3) pt2(3) pt3(3) pt4(3)],'c','FaceAlpha',0.15,'FaceColor',[0.8 0.8 0.8], ...
                'EdgeAlpha',0.05);
            hold on;
            
            % Plot the herpelhode and omegaParallel
            plot3(omegaB_W(1,:),omegaB_W(2,:),omegaB_W(3,:),'k');
            plot3(O(1),O(2),O(3),'or');
            plot3([0 O(1)],[0 O(2)],[0 O(3)],'--k');
            
            % Plot the ellipsoid and the body axes
            surf(xE,yE,zE,'FaceAlpha',0.15,'EdgeAlpha',0.05,'FaceColor',[0 0.8 1]);
            plot3([O(1) O(1)+axis1(1)],[O(2) O(2)+axis1(2)],[O(3) O(3)+axis1(3)],'r','LineWidth',1.5);
            plot3([O(1) O(1)+axis2(1)],[O(2) O(2)+axis2(2)],[O(3) O(3)+axis2(3)],'g','LineWidth',1.5);
            plot3([O(1) O(1)+axis3(1)],[O(2) O(2)+axis3(2)],[O(3) O(3)+axis3(3)],'b','LineWidth',1.5);
            
            % Plot the polhode attached to the body ellipsoid
            plot3(omegaB_Bc(1,:),omegaB_Bc(2,:),omegaB_Bc(3,:),'color',[0.9 0.6 0],'LineWidth',1.5);
            
            grid on; axis equal; title('Rolling of body-fixed energy ellipsoid on the invariable plane');
            
            
        end
    end
        %% Static methods
    methods(Static)
        %% Plot example polhodes for TA, AS1, and AS3 in LE, ME, and HE
        function plotExamplePolhodes()
                        
            % Set the inertia ratios
            Jta = [1.239 1.1905]';      % approx SPHERES fit
            Jas1 = [1.8534 1]';         % from CL_MAJ test
            Jas3 = [1.6839 1.6839]';    % from CL_MIN test
                        
            % Example initial angular velocity and resulting properties for taLE
            omegaB0_BtaLE = [0.939392242898362   0   0.500486277097766]';
            
            % Create the rigid body rotation objects, basing them off the taLE case
            RB0toW = eye(3);
            taLE  = RigidBodyRotation(Jta,RB0toW,omegaB0_BtaLE,'omega0',0);
            taME  = RigidBodyRotation(Jta,RB0toW,[taLE.h^2/(2*Jta(2)) taLE.h 0]','Ekht0',0);
            taHE  = RigidBodyRotation(Jta,RB0toW,[taLE.h^2/(2*Jta(2))+0.005 taLE.h 0]','Ekht0',0);
            as1LE = RigidBodyRotation(Jas1,RB0toW,[taLE.Ek taLE.h 0]','Ekht0',0);
            as1ME = RigidBodyRotation(Jas1,RB0toW,[taLE.h^2/(2*Jas1(2)) taLE.h 0]','Ekht0',0);
            as3ME = RigidBodyRotation(Jas3,RB0toW,[taLE.h^2/(2*Jas3(2)) taLE.h 0]','Ekht0',0);
            as3HE = RigidBodyRotation(Jas3,RB0toW,[taLE.h^2/(2*Jas3(2))+0.15 taLE.h 0]','Ekht0',0);
            
            % Plot in a 3x3 grid with rows TA, AS1, AS3 and columns LE, ME, HE. Use transparency
            % of one ellipsoid when necessary.
            theta = linspace(0,2*pi); x = cos(theta); y = sin(theta); % for AS circles
            subplot = @(m,n,p) subtightplot (m, n, p, [0.075 0.01], [0.05 0.05], [0.001 0.001]);      
            subplot(3,3,1); taLE.plotPolhode();
            subplot(3,3,2); taME.plotPolhode();
            subplot(3,3,3); taHE.plotPolhode();
            subplot(3,3,4); as1LE.plotPolhode();
            subplot(3,3,5); as1ME.plotPolhode();
              plot3(zeros(size(x)),x*as1ME.h/as1ME.J(2),y*as1ME.h/as1ME.J(3), ...
                  'LineWidth',2,'Color',[0.9 0.9 0.9]);
            subplot(3,3,8); as3ME.plotPolhode(); 
              plot3(x*as3ME.h/as3ME.J(2),y*as3ME.h/as3ME.J(2),zeros(size(x)), ...
                  'LineWidth',2,'Color',[0 0.8 1]);
            subplot(3,3,9); as3HE.plotPolhode();
            
        end
    end
    
end

