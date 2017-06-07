classdef InertiaConstraints < handle
    %INERTIACONSTRAINTS Contains six conic constraints on inertia ratios from fitting to omegaB_B.
    %   Conics fit to the rigid body angular velocity in planes [1 2], [2 3], and [1 3] impose
    %   constraints on possible values of J = [J1 J2]'. This class creates parametric versions of
    %   these conic constraints for optimization purposes.
    
    properties
        omegaB_B;               % the rigid body angular velocity data in the body frame
        energyState;            % string, 'LE' when h^2 > 2*Ek*J2, 'HE' when h^2 < 2*Ek*J2
        inertiaSymmetry;        % either tri-axial 'TA', or axis-symmetric 'AS1' or 'AS3'
        conicFits;              % a 3x1 vector of conic fits in planes [1 2], [2 3], and [1 3]
        omegaConstraints;       % a 3x1 vector of conic equality constraints from fit ang vels
        energyConstraint;       % a conic inequality constraint from the high / low energy level
        linConstraintMat = [-1 1; 0 -1; 1 -1]; % matrix A in the inequality constraint A*J <= b
        linConstraintVec = [0 -1 1]';          % vector b in the inequality constraint A*J <= b
    end
    
    methods
        %% Construct from body angular velocities and an indication of 
        function this = InertiaConstraints(omegaB_B,energyState,inertiaSymmetry)
            
            % Store angular velocities and energy state, and get the number of samples
            this.omegaB_B = omegaB_B;
            this.energyState = energyState;      
            this.inertiaSymmetry = inertiaSymmetry;
            nT = size(omegaB_B,2);
            J3 = 1;
            
            % Get the appropriate conic fits in each plane ([1 2], [2 3], and [1 3]) and 
            % normalize the K's. See logbook #3 pp 85 and theory_inertiaRatiosConicConstraints
            [Kfit1,~,~,~] = Conic.fitConicCanonical([omegaB_B(1,:); omegaB_B(2,:)]); % fit ell
            [Kfit2,~,~,~] = Conic.fitConicCanonical([omegaB_B(2,:); omegaB_B(3,:)]); % fit ell
            [~,Kfit3,~,~] = Conic.fitConicCanonical([omegaB_B(1,:); omegaB_B(3,:)]); % fit hyp
            Kfit1 = Kfit1/Kfit1(end); Kfit2 = Kfit2/Kfit2(end); Kfit3 = Kfit3/Kfit3(end);
            this.conicFits = [Conic(Kfit1); Conic(Kfit2); Conic(Kfit3)];
            
            if strcmp(this.inertiaSymmetry,'TA')                % Tri-Axial: conic constraints
                
                omegaBarSq = 1/nT*sum(omegaB_B.^2,2);
                
                Kalpha1 = [ Kfit1(1)*omegaBarSq(1) + 1;         % A (J1^2)
                            0;                                  % B (J1*J2)
                            Kfit1(1)*omegaBarSq(2);             % C (J2^2)
                           -Kfit1(1)*J3*omegaBarSq(1) - J3;     % D (J1)
                           -Kfit1(1)*J3*omegaBarSq(2);          % E (J2)
                            0 ];                                % F ( )
                Kbeta1 = [  Kfit1(3)*omegaBarSq(1);             % A (J1^2)
                            0;                                  % B (J1*J2)
                            Kfit1(3)*omegaBarSq(2) + 1          % C (J2^2)
                           -Kfit1(3)*J3*omegaBarSq(1);          % D (J1)
                           -Kfit1(3)*J3*omegaBarSq(2) - J3;     % E (J2)
                            0 ];                                % F ( )
                Kalpha2 = [ 0;                                  % A (J1^2)
                           -Kfit2(1)*omegaBarSq(2) - 1;         % B (J1*J2)
                            Kfit2(1)*omegaBarSq(2) + 1;         % C (J2^2)
                           -Kfit2(1)*J3*omegaBarSq(3);          % D (J1)
                            0;                                  % E (J2)
                            Kfit2(1)*J3^2*omegaBarSq(3) ];      % F ( )
                Kbeta2 = [  0;                                  % A (J1^2)
                           -Kfit2(3)*omegaBarSq(2);             % B (J1*J2)
                            Kfit2(3)*omegaBarSq(2);             % C (J2^2)
                           -Kfit2(3)*J3*omegaBarSq(3) - J3;     % D (J1)
                            0;                                  % E (J2)
                            Kfit2(3)*J3^2*omegaBarSq(3) + J3^2];% F ( )
                Kalpha3 = [ Kfit3(1)*omegaBarSq(1) + 1;         % A (J1^2)
                           -Kfit3(1)*omegaBarSq(1) - 1;         % B (J1*J2)
                            0;                                  % C (J2^2)
                            0;                                  % D (J1)
                           -Kfit3(1)*J3*omegaBarSq(3);          % E (J2)
                            Kfit3(1)*J3^2*omegaBarSq(3)];       % F ( )
                Kbeta3 = [  Kfit3(3)*omegaBarSq(1);             % A (J1^2)
                           -Kfit3(3)*omegaBarSq(1);             % B (J1*J2)
                            0;                                  % C (J2^2)
                            0;                                  % D (J1)
                           -Kfit3(3)*J3*omegaBarSq(3) - J3;     % E (J2)
                            Kfit3(3)*J3^2*omegaBarSq(3) + J3^2];% F ( )
                this.omegaConstraints = [ Conic(Kalpha1) ; %Conic(Kbeta1); ...
                                          Conic(Kalpha2) ; %Conic(Kbeta2); ...
                                          Conic(Kalpha3) ]; %Conic(Kbeta3)  ];
                            
            end
            
            % Find the conic which divides the low energy cases from the high energy cases
            % See logbook #3 pp 79
            Ae = -sum(omegaB_B(1,:).^2);    Be = sum(omegaB_B(1,:).^2);
            Ce = 0;                         De = 0;
            Ee = J3*sum(omegaB_B(3,:).^2);  Fe = -J3^2*sum(omegaB_B(3,:).^2);
            Ke = [Ae Be Ce De Ee Fe]'/Fe;
            this.energyConstraint = Conic(Ke);
            
            % TODO: determine intersections of equality and inequality constraints
            % and corresponding parametric limits
            
        end
        %% Get the inertia ratios from their parameterized form, from the equality constraint 
        %  num (1-3); if p is a 1xN matrix, the returned J will be a 2xN matrix of inertia ratios
        function J = ratios(this,p,constraintNum)
            
            if strcmp(this.inertiaSymmetry,'TA')        % Tri-axial: conic omega constraints
                J = this.omegaConstraints(constraintNum).conic(p);
            elseif strcmp(this.inertiaSymmetry,'AS1')   % Axis-symmetric 1: linear omega constraints
                J = [p; ones(1,length(p))];             % J1 = p, J2 = 1
            elseif strcmp(this.inertiaSymmetry,'AS3')   % Axis-symmetric 3: linear omega constraints
                J = [p; p];                             % J1 = J2 = p
            end
        end
        %% Find the inertia ratios on the equality constraints num (1-3) that are closest to J
        function [pOpt,Jopt] = findClosestPoint(this,J,constraintNum)
            
            % Find the value closest to the equality constraint
            if strcmp(this.inertiaSymmetry,'TA')        % Tri-axial: conic omega constraints
                [pOpt,Jopt] = this.omegaConstraints(constraintNum).findClosestPoint(J);
            elseif strcmp(this.inertiaSymmetry,'AS1')   % Axis-symmetric 1: linear omega constraints
                pOpt = J(1);
                Jopt = [pOpt 1]';
            elseif strcmp(this.inertiaSymmetry,'AS3')   % Axis-symmetric 3: linear omega constraints
                pOpt = (J(1)+J(2))/2;                   % see logbook #3 pp 121
                Jopt = [pOpt pOpt]';
            end
            
        end
        %% Plot all of the constraints on inertia ratios
        function plot(this)
            
            % Create a range of J1 to plot over
            J1 = 0:0.01:4; N = length(J1);
            
            % Create line styles for linear constraints
            lnStyles = {'-k','--k','-.k','.k'};
            
            % Create shaded area for inequality constraints
%             fill([1 2 4 4],[1 1 3 4],[0.9 0.9 0.9],'EdgeAlpha',0); hold on;
            
            % Plot the linear constraints
            for i = 1:size(this.linConstraintMat,1)
                J2 = (repmat(this.linConstraintVec(i),1,N) - this.linConstraintMat(i,1)*J1) ...
                        / this.linConstraintMat(i,2);
                plot(J1,J2,lnStyles{i});
            end
            
            % Plot the energy constraint
            this.energyConstraint.plot('actual','--k');
            
            % Plot the omega equality constraints
            if strcmp(this.inertiaSymmetry,'TA')
                this.omegaConstraints(1).plot('actual','-r');
                this.omegaConstraints(2).plot('actual','--g');
                this.omegaConstraints(3).plot('actual','-.b');
%                 this.omegaConstraints(4).plot('actual','--m');
%                 this.omegaConstraints(5).plot('actual','-b');
%                 this.omegaConstraints(6).plot('actual','--c');            
            elseif strcmp(this.inertiaSymmetry,'AS1')
                plot([1 2],[1 1],'r','LineWidth',2); hold on;
            elseif strcmp(this.inertiaSymmetry,'AS3')
                plot([1 15],[1 15],'r','LineWidth',2); hold on;
            end
            
            
        end
    end
    
end

