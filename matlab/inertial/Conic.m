classdef Conic < handle
    %CONIC Represents a 2D conic.
    %   The conic can be represented parametrically in its canonical frame C (i.e. centered at 
    %   the origin and aligned with the x-y axes), or parametrically in its actual frame A (i.e.
    %   centered and aligned as dictated by K).
    
    properties
        K;                  % vector [A B C D E F]' in Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
                            % in actual frame
        Kc;                 % vector [A B C D E F]' in Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
                            % in canonical frame
        PAtoC;              % the 2D pose change [RCtoA tAtoC_A; 0 0 1] from the canonical 
                            % (parametric) frame C to the actual frame C
        aSemi; bSemi;       % semi-major axes values of the conic
        ecc;                % eccentricity fo the conic
        type;               % conic type, 'ell', 'hyp', or 'par'
        dir;                % conic alignment 'xy' for ell, 'xy' or 'yx' for hyp or par
    end
    
    methods
        %% Construct from K = [A B C D E F]' in Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
        function this = Conic(K)
            this.K = K;
            [this.Kc,this.PAtoC,this.aSemi,this.bSemi,this.ecc,this.type,this.dir] = this.findCanonical();
        end
        %% Find the transformation [RCtoA tAtoC_A; 0 0 1] from the canonical frame C to the actual
        %  frame A, the semi-major axes, and the conic type and direction
        function [Kc,PAtoC,aSemi,bSemi,ecc,type,dir] = findCanonical(this) 
                        
            % Get the coefficients in the actual frame A
            [aA,bA,cA,dA,eA,fA] = this.getCoeffs();
            
            % Find the center of the conic in actual coordinates (see logbook #3 pp 88)
            Amat = [aA bA/2; bA/2 cA]; bvec = [dA eA]';
            tAtoC_A = (-2*Amat)\bvec;
            
            % Find the rotation matrix between the canonical and actual frames
            theta = 1/2*atan(bA/(aA-cA));
            RCtoA = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            
            % Form transformation PAtoC
            PAtoC = [RCtoA tAtoC_A; 0 0 1];
            
            % Get coefficients for the de-translated and de-rotated conic
            aC = aA*cos(theta)^2 + bA*cos(theta)*sin(theta) + cA*sin(theta)^2;
            bC = 0;
            cC = aA*sin(theta)^2 - bA*sin(theta)*cos(theta) + cA*cos(theta)^2;
            dC = 0; 
            eC = 0;
            fC = fA - (aA*tAtoC_A(1)^2 + bA*tAtoC_A(1)*tAtoC_A(2) + cA*tAtoC_A(2)^2);
            Kc = [aC bC cC dC eC fC]';
            if sign(fC) ~= 0
                Kc = -Kc/fC;
                aC = Kc(1); bC = Kc(2); cC = Kc(3); dC = Kc(4); eC = Kc(5); fC = Kc(6);
            end
                        
            % Get the semi-major axes, conic type, and direction
            aSemi = sqrt(abs( -fC/aC )); bSemi = sqrt(abs( -fC/cC ));
            if aC*cC > 0
                type = 'ell';
                dir = 'xy';
                ecc = sqrt(1 - min([-fC/aC -fC/cC])/max([-fC/aC -fC/cC]));
            elseif aC*cC < 0
                type = 'hyp';
                ecc = sqrt(1 + abs(min([-fC/aC -fC/cC]))/abs(max([-fC/aC -fC/cC])));
                if aC > 0 && cC < 0
                    dir = 'xy';
                elseif aC < 0 && cC > 0
                    dir = 'yx';
                end
            elseif aC*cC == 0
                type = 'par';       % TODO: more classification of parabola
            end
                
        end
        %% Evaluate the canonical (parametric) conic at a 1xN vector of parameters, returning a 
        %  2xN matrix of (x,y) points; twoSided is an optional boolean indicating whether to 
        %  provide both sides of the hyperbola.
        function pts = conicCanonical(this,p,twoSided)
            if strcmp(this.type,'ell')
                pts(1,:) = this.aSemi*cos(p);           % x
                pts(2,:) = this.bSemi*sin(p);           % y
            elseif strcmp(this.type,'hyp')
                if strcmp(this.dir,'xy')
                    pts = [ this.aSemi*cosh(p) ;        % x
                            this.bSemi*sinh(p) ];       % y
                    if nargin > 2 && twoSided
                        pts = [ pts(1,:) -this.aSemi*cosh(p) ;
                                pts(2,:)  this.bSemi*sinh(p) ];
                    end
                elseif strcmp(this.dir,'yx')
                    pts = [ this.aSemi*sinh(p) ;        % x
                            this.bSemi*cosh(p) ];       % y
                    if nargin > 2 && twoSided
                        pts = [ pts(1,:)  this.aSemi*sinh(p) ;
                                pts(2,:) -this.bSemi*cosh(p) ];
                    end
                end
            elseif strcmp(this.type,'par')
                % TODO: suppor parabolas
            end
        end
        %% Evaluate the actual parametric conic at a 1xN vecotor of paramters, returning a 2xN
        %  matrix of (x,y) points; twoSided is an optional boolean indicating whether to 
        %  provide both sides of the hyperbola.
        function pts = conic(this,p,twoSided)
            
            % Default twoSided boolean
            if nargin < 3
                twoSided = 0;
            end
            
            % Get points in canonical frame
            pts_C = this.conicCanonical(p,twoSided);
            
            % Get rotation and translation
            RCtoA = this.PAtoC(1:2,1:2);    tAtoC_A = this.PAtoC(1:2,3);
            pts =  RCtoA*pts_C + repmat(tAtoC_A,1,size(pts_C,2));
            
        end
        %% Find the closest point on the conic to the given (x,y) point in actual frame A
        function [pOpt,ptOpt] = findClosestPoint(this,pt)
            
            % Initialize
            p0 = 0;
            
            % Perform unconstrained minimization to find the optimal parameter
            opts = optimoptions(@fminunc,'Algorithm','quasi-newton','TolFun',1e-10,'display','off');
            pOpt = fminunc(@(p) norm(this.conic(p) - pt)^2,p0,opts);
            
            % Find the optimal point
            ptOpt = this.conic(pOpt);
            
        end
        %% Plot the conic, with either frame = 'actual' or 'canonical', with the specified 
        %   linespec (e.g. '-k'). Suppress print when that boolean is true. Return two sided
        %   hyperbola when that boolean is true.
        function pts = plot(this,frame,lineSpec,supressPrint,twoSided)
            
            % Default for twoSided boolean
            if nargin < 5
                twoSided = 0;
            end
            
            % Get a range of parameters to plot over
            if strcmp(this.type,'ell')
                p = 0:0.01:2*pi;           
            elseif strcmp(this.type,'hyp')
                p = -pi:0.01:pi;    
            elseif strcmp(this.type,'par')
                % TODO: handle parabola case
            end
            nP = length(p);
            
            % Get and plot points
            if strcmp(frame,'actual')
                pts = this.conic(p,twoSided);
            elseif strcmp(frame,'canonical')
                pts = this.conicCanonical(p,twoSided);
            end
            
            if nargin < 4 || ~supressPrint
                plot(pts(1,1:nP),pts(2,1:nP),lineSpec,'LineWidth',2); hold on;
                if nargin > 4 && twoSided
                    plot(pts(1,nP+1:end),pts(2,nP+1:end),lineSpec,'LineWidth',2);
                end
            end
            
        end
        %% Get the coefficients of the conic
        function [A,B,C,D,E,F] = getCoeffs(this)
            A = this.K(1); B = this.K(2); C = this.K(3); 
            D = this.K(4); E = this.K(5); F = this.K(6);
        end
        
    end
    
    %% Static methods
    methods(Static)
        %% Fit conic sections (an ellipse and a hyperbola) given data, a 2xN matrix of (x,y) points 
        %  in the canonical frame (i.e. un-rotated and un-translated)
        function [Kell,Khyp,rSqEll,rSqHyp] = fitConicCanonical(data)
            
            %   Solve the problem argmin_z z'*D'*D*z, subject to z'*C*z = k
            %   where: z = [a b f]' (coefficient vector), D = [x1^2 y1^2 1; ... ; xN^2 yN^2 1]
            %   (data matrix), C = [0 0.5 0; 0.5 0 0; 0 0 0] (constraint matrix), k = 1 for ellipse,
            %   and k = -1 for hyperbola (makes ab = 1 for ellipse, ab = -1 for hyperbola).
            %   See logbook #3 pp 67-68

            % Constraint matrix
            C = [0 0.5 0; 0.5 0 0; 0 0 0];

            % Get the data and data size
            x = data(1,:)'; y = data(2,:)'; N = length(x);

            % Data matrix
            D = [x.^2 y.^2 ones(N,1)];

            % Optimize by setting the derivative of the Lagrangian to zero
            % gradient(L)_z = D'*D*z - lambda*C*z = [0 0 0]'
            % Solve by solving the the general eigenvalue problem
            [eigVec,eigVal] = eig(D'*D,C);

            % Get the eigenvalues eigenvectors of the ellipse and hyperbola
            nZer = find(diag(eigVal) >= -1e-6 & diag(eigVal) <= 1e-6);
            nPos = find(diag(eigVal) >   1e-6 & ~isinf(diag(eigVal)));
            nNeg = find(diag(eigVal) <  -1e-6 & ~isinf(diag(eigVal)));
            if size(nZer,1) == 0                % noisy data, no perfect fit
                eigVecEll = eigVec(:,nPos);
                eigVecHyp = eigVec(:,nNeg);
                lambdaEll = eigVal(nPos,nPos);
                lambdaHyp = eigVal(nNeg,nNeg);
            elseif size(nZer,1) ~= 0            % noisless data, contains a perfect fit
                if size(nPos,1) == 0            % no positive noninfinite eigenvalues
                    eigVecEll = eigVec(:,nZer);
                    eigVecHyp = eigVec(:,nNeg);
                    lambdaEll = eigVal(nZer,nZer);
                    lambdaHyp = eigVal(nNeg,nNeg);
                elseif size(nNeg,1) == 0        % no negative noninfinite eigenvalues
                    eigVecEll = eigVec(:,nPos);
                    eigVecHyp = eigVec(:,nZer);
                    lambdaEll = eigVal(nPos,nPos);
                    lambdaHyp = eigVal(nZer,nZer);
                end
            end

            % Scale the resultant eigenvectors [a b f]' so that f = -1
            ell = -eigVecEll/eigVecEll(3);
            hyp = -eigVecHyp/eigVecHyp(3);

            % The residual squared r'*r is given by the absolute value of the eigenvalues
            rSqEll = abs(lambdaEll);
            rSqHyp = abs(lambdaHyp);
            
            % Rearrange the fits ell and hyp, which will be in the form [A 0 C 0 0 F]'
            Kell = [ell(1) 0 ell(2) 0 0 ell(3)]';
            Khyp = [hyp(1) 0 hyp(2) 0 0 hyp(3)]';

        end
        %% Fit the conic section with the lowest residual squared to a 2xN matrix of (x,y) points 
        %  in the canonical frame (i.e. un-rotated and un-translated)
        function [Kopt,rSqOpt,type] = fitBestConicCanonical(data)
            
            % Get both elliptic and hyperbolic solutions
            [Kell,Khyp,rSqEll,rSqHyp] = Conic.fitConicCanonical(data);
            
            % Select the best one
            if rSqEll < rSqHyp
                Kopt = Kell;
                rSqOpt = rSqEll;
                type = 'ell';
            else
                Kopt = Khyp;
                rSqOpt = rSqHyp;
                type = 'hyp';
            end
            
        end
    end
    
end

