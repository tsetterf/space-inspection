classdef PrincipalAxesOpt < handle
    %PRINCIPALAXESOPT Finds optimal principal axes alignment from a polhode in an arbitrary frame.
    %   Given an angular velocity polhode in an arbitrary geometric frame G, this class seeks the
    %   optimal rotation matrix from the body frame B to the geometric frame G, RBtoG. It finds
    %   by finding the rotation which creates the best-fitting conic polhode projections: a
    %   hyperbola on the cn-dn plane and ellipses on the sn-cn and sn-dn planes.
    
    properties
        omegaB_G;               % body angular velocities in geometric frame G   
        omegaB_B;               % body angular velocities in body frame B after optimization
        conicFits1;             % a 3x1 vector of init conic fits in planes [1 2], [2 3], and [1 3]
        conicFits2;             % a 3x1 vector of final conic fits in planes [1 2], [2 3], and [1 3]
        REtoG;                  % the optimal rotation from elliptical frame E to geometric frame G
        RBtoG;                  % the optimal rotation from body frame B to geometric frame G
        costOpt;                % cost of the optimization
        energyState;            % string, 'LE' when h^2 > 2*Ek*J2, 'HE' when h^2 < 2*Ek*J2
        planeAxes = [ [1; 2] [2; 3] [1; 3] ];   % axis combinations forming planes
        inertiaSymmetry;        % either tri-axial 'TA', or axis-symmetric 'AS1' or 'AS3'
    end
    
    methods
        %% Construct from body angular velocities in geometric frame G
        function this = PrincipalAxesOpt(omegaB_G)
            this.omegaB_G = omegaB_G;
        end
        %% Get the cost of the conic fits given an so(3) alignment vector theta
        function Cost = cost(this,theta)
            
            % Get the proposed rotation matrix using the SO(3) exponential map
            REtoGp = Exp(theta);
            
            % Find the best conic fits in the E planes [1 2], [2 3], and [1 3]
            [~,rSqOpt,~] = this.getBestConicFits(REtoGp);
            
            % Cost is now the sum of the deviations from the conic equation a*x^2 + b*y^2 + f = 0,
            % where each plane has the best conic fit possible
            Cost = sum(rSqOpt);
            
        end
        %% Constrain the so(3) vector theta to the R^3 ball <= pi 
        function [c,ceq] = nonlinearConstraints(this,theta)
            
            % Inequality constraints of the form c <= 0 (i.e. c = ||theta||^2 - pi^2 <= 0)
            c = norm(theta)^2 - pi^2;
            
            % No equality constraints
            ceq = [];
            
        end
        %% Perform the optimization to find the elliptic frame E, and then process the results 
        %  to further find RGtoB and the energy state.
        function [RBtoG,costOpt] = optimize(this)
            
            % Initial guess for so(3) vector theta
            theta0 = [0 0 0]';
            
            % Optimize the rotation into the ellipse frame REtoG
            opts = optimset('TolFun',1e-13,'MaxFunEvals',5e3,'display','off', ...
                'Algorithm','active-set');  % interior-point sometimes performs poorly
            [thetaOpt,costOpt,~,~,~,~,~] = fmincon(@(theta) this.cost(theta), ...
                theta0,[],[],[],[],[],[],@(theta) this.nonlinearConstraints(theta),opts);

            disp('Principal axes optimization complete');
            disp(['Cost   : ' num2str(costOpt)]);
            disp(['theta* : ' num2str(thetaOpt')]);
            
            % Get the optimal rotation
            this.REtoG = Exp(thetaOpt);
            
            % Get the best initial conic fits where they can be any combination of ell and hyp
            [Kopt1,~,~] = this.getBestConicFits(this.REtoG);
            this.conicFits1{1} = Conic(Kopt1(:,1)); this.conicFits1{2} = Conic(Kopt1(:,2));
            this.conicFits1{3} = Conic(Kopt1(:,3));
            
            % Get the final conic fits where there must be two ell and one hyp
            [Kopt2,~,~] = this.getBestConicFitsConstrained(this.REtoG);
            this.conicFits2{1} = Conic(Kopt2(:,1)); this.conicFits2{2} = Conic(Kopt2(:,2));
            this.conicFits2{3} = Conic(Kopt2(:,3));
                        
            % Get the body frame and energy state
            [RBtoG,energyState] = this.getBodyFrame(this.REtoG);
            this.RBtoG = RBtoG;
            this.costOpt = costOpt;
            this.energyState = energyState;
            
            % Get body aligned angular velocities
            this.omegaB_B = RBtoG' * this.omegaB_G;
            
        end
        %% Given an optimal rotation REtoG = [xE_G yE_G zE_G], find the cn and dn-axes in the 
        %  elliptic frame E.
        function [cnHat_E,dnHat_E] = getHypAxes(this,REtoG)
            
            % Rotate all data from geometric frame to the proposed elliptical frame E
            omegaB_E = REtoG' * this.omegaB_G;
                     
            % Find the plane with the best hyperbolic fit
            types = {this.conicFits2{1}.type this.conicFits2{2}.type this.conicFits2{3}.type};
            indHyp = find(strcmp(types,'hyp'));
            Khyp = this.conicFits2{indHyp}.Kc;
            
            % Get all non-hyperbolic planes with low eccentricity fits
            inds = (1:3)'; 
            ecc = [this.conicFits2{1}.ecc this.conicFits2{2}.ecc this.conicFits2{3}.ecc]';
            circPlanes = inds(inds~=indHyp & ecc <= 0.35)';
            
            % Check if hyperbola opens normal to any of the circular planes
            isAxisSymmetric = 0;
            planeAxesStack = [this.planeAxes(1,:); zeros(1,3); this.planeAxes(2,:)];
            for j = 1:length(circPlanes)
                
                % Get the axis normal to the circular plane
                circNormAx = inds(inds~=this.planeAxes(1,circPlanes(j)) ...
                                    & inds~=this.planeAxes(2,circPlanes(j)));
                % Get the index in the conic fit Kc corresponding to the hyperbola's coefficient
                indKc = find(planeAxesStack(:,indHyp)==circNormAx);
                % Check if hyperbola opens normal to circular plane
                if -this.conicFits2{indHyp}.Kc(indKc)/this.conicFits2{indHyp}.Kc(6) > 0
                    isAxisSymmetric = 1;
                    eccCirc = ecc(circPlanes(j));
                    break;
                end
                
            end
             
            % Behave differently if polhode is axis-symmetric
            if isAxisSymmetric
                                
                % Axis-symmetric; whether it is 'AS1' or 'AS3' is determined with the energy
                this.inertiaSymmetry = 'AS'; 
                disp(['Treating as axis-symmetric, eccentricity = ' num2str(eccCirc)]);
                
                % The equivalent to dnHat_E is the circle normal (hyperbolic axis) in the direction
                % that makes omegaDn always positive
                dnHat_E = sign(mean(omegaB_E(circNormAx,:)))* (inds==circNormAx);
                
                % Choose cnHat_E so that omegaB0_E is maximum in the cn(cos) direction
                omegaB0_E = omegaB_E(:,1);  omegaB0_E(circNormAx) = -Inf;
                [~,indCn] = max(omegaB0_E);
                cnHat_E = inds==indCn;
                
            else
                            
                % Tri-axial
                this.inertiaSymmetry = 'TA'; 
                disp(['Treating as tri-axial, min eccentricity = ' num2str(min(ecc))]);
                
                % Get the normal vectors of asymptotes to the hyperbolic fit (logbook #3 pp 76)
                n1 = [-sqrt(-Khyp(1)/Khyp(3)) 1]';
                n2 = [ sqrt(-Khyp(1)/Khyp(3)) 1]';

                % Assign each angular velocity vector in frame E to a quadrant (Q1, Q2, Q3, or Q4), and
                % get the quadrant in which the dn-axis is located
                n1Dot = ( n1' * omegaB_E(this.planeAxes(:,indHyp),:) )';     % n1^T * x
                n2Dot = ( n2' * omegaB_E(this.planeAxes(:,indHyp),:) )';     % n2^T * x
                numInQuad = sum( [ (n1Dot < 0 & n2Dot > 0) (n1Dot > 0 & n2Dot > 0) ...
                                   (n1Dot > 0 & n2Dot < 0) (n1Dot < 0 & n2Dot < 0) ]);
                [~,dnSec] = max(numInQuad);

                % Get the axis index of the dn-axis, and then determine the index of the sn and cn-axes
                % Make the dn-axis sign so as to make all omegas +ve along the dn-axis
                if dnSec == 1
                    dnHat_E = inds==this.planeAxes(1,indHyp);  
                    cnHat_E = inds==this.planeAxes(2,indHyp);
                elseif dnSec == 2
                    dnHat_E =   inds==this.planeAxes(2,indHyp);  
                    cnHat_E = inds==this.planeAxes(1,indHyp);
                elseif dnSec == 3
                    dnHat_E = -(inds==this.planeAxes(1,indHyp)); 
                    cnHat_E = inds==this.planeAxes(2,indHyp);
                elseif dnSec == 4
                    dnHat_E = -(inds==this.planeAxes(2,indHyp)); 
                    cnHat_E = inds==this.planeAxes(1,indHyp);
                end
                % Make the cn-axis sign so as to make the intial omega along the cn-axis +ve
                cnHat_E = sign(omegaB_E(:,1)'*cnHat_E)*cnHat_E;
            
            end
            
        end
        %% Get the best conic fits in each plane given an optimal rotation REtoG
        function [Kopt,rSqOpt,types] = getBestConicFits(this,REtoG)
            
            % Rotate all data from geometric frame to the proposed elliptical frame E
            omegaB_E = REtoG' * this.omegaB_G;
            
            % Get the best fits in each plane
            Kopt = zeros(6,3); rSqOpt = zeros(1,3); types = [];
            [Kopt(:,1),rSqOpt(1),types{1}] = ...
                Conic.fitBestConicCanonical(omegaB_E(this.planeAxes(:,1),:));
            [Kopt(:,2),rSqOpt(2),types{2}] = ...
                Conic.fitBestConicCanonical(omegaB_E(this.planeAxes(:,2),:));
            [Kopt(:,3),rSqOpt(3),types{3}] = ...
                Conic.fitBestConicCanonical(omegaB_E(this.planeAxes(:,3),:));
            
        end
        %% Get the best conic fits in each plane given an optimal rotation REtoG and constraining 
        %  the set of fits to have two ellipsoids and one hyperbola
        function [Kopt,rSqOpt,types] = getBestConicFitsConstrained(this,REtoG)
            
            % Rotate all data from geometric frame to the proposed elliptical frame E
            omegaB_E = REtoG' * this.omegaB_G;
            
            % Get all best fits
            Kell = zeros(6,3); rSqEll = zeros(1,3); Khyp = zeros(6,3); rSqHyp = zeros(1,3);
            [Kell(:,1),Khyp(:,1),rSqEll(1),rSqHyp(1)] = ...
                    Conic.fitConicCanonical(omegaB_E(this.planeAxes(:,1),:));
            [Kell(:,2),Khyp(:,2),rSqEll(2),rSqHyp(2)] = ...
                    Conic.fitConicCanonical(omegaB_E(this.planeAxes(:,2),:));
            [Kell(:,3),Khyp(:,3),rSqEll(3),rSqHyp(3)] = ...
                    Conic.fitConicCanonical(omegaB_E(this.planeAxes(:,3),:));
            
            % Keep track of number of elliptical fits and number of hyperbolic fits
            numEll = 0; numHyp = 0;
            
            % Concatenate options and sort to find the lowest values error values
            Ks = [Kell Khyp]; rSq = [rSqEll rSqHyp]; 
            planes = [1 2 3 1 2 3];
            typesAll = {'ell','ell','ell','hyp','hyp','hyp'};
            [~,rSqRank] = sort(rSq); 
            typesAll = typesAll(rSqRank); planes = planes(rSqRank); Ks = Ks(:,rSqRank);
            planesChosen = zeros(1,3); Kopt = zeros(6,3); rSqOpt = zeros(1,3); 
                        
            % Loop through the different options from best to worst, finding the best combination 
            % of chosen fits
            for i = 1:6
                
                if (strcmp(typesAll{i},'ell') && numEll < 2 && sum(planesChosen==planes(i))==0) ...
                || (strcmp(typesAll{i},'hyp') && numHyp < 1 && sum(planesChosen==planes(i))==0)                    
                    planesChosen(i)   = planes(i);
                    types{planes(i)}  = typesAll{i};
                    Kopt(:,planes(i)) = Ks(:,i);
                    rSqOpt(i)         = rSq(i); 
                    
                    if strcmp(typesAll{i},'ell')
                        numEll = numEll + 1;
                    else
                        numHyp = numHyp + 1;
                    end                    
                end
                
            end
            
        end
        %% Obtain the rotation matrix RGtoB and energy state given an optimal rotation REtoG
        function [RBtoG,energyState] = getBodyFrame(this,REtoG)
           
            % Assign the axes from REtoG = [xE_G yE_G zE_G] to RGtoB = [xB_G yB_G zB_G]
            [cnHat_E,dnHat_E] = this.getHypAxes(REtoG);
            
            % Get the unit vectors for the cn, and dn-axes in the G frame
            cnHat_G = REtoG * cnHat_E;
            dnHat_G = REtoG * dnHat_E;
                        
            % Find the energy state (LE) or (HE) based on the direction of evolution of the polhode 
            % around the dn-axis. See logbook #3 pp. 10 & 19. Use the (LE) / (HE) determination to 
            % form a proper RGtoB
            
            % Get the average projection of omega(i) x omega(i+1) on the dn-axis
            nT = size(this.omegaB_G,2); muProjDn = 0;
            for i = 1:nT-1
                muProjDn = muProjDn + ...
                        1/(nT-1)*dnHat_G' * (skew(this.omegaB_G(:,i)) * this.omegaB_G(:,i+1));
            end

            % If the projection follows the right hand rule on average, it is the (LE) case; 
            % otherwise it is the (HE) case. Specify RGtoBE using the appropriate (LE)/(HE) 
            % elliptic functions relationships. The sn-axis = cn_E x dn_E (LE) or dn_E x cn_E (HE)
            % If axis-symmetric, the energy state will determine whether it is the AS1 or AS3 case
            if muProjDn >= 0
                energyState = 'LE';   
                snHat_E = skew(cnHat_E)*dnHat_E;
                snHat_G = REtoG * snHat_E;                
                RBtoG = [dnHat_G snHat_G cnHat_G]; 
                if strcmp(this.inertiaSymmetry,'AS')
                    this.inertiaSymmetry = [this.inertiaSymmetry '1'];    % AS1
                end
            else
                energyState = 'HE';
                snHat_E = skew(dnHat_E)*cnHat_E;
                snHat_G = REtoG * snHat_E;       
                RBtoG = [cnHat_G snHat_G dnHat_G]; 
                if strcmp(this.inertiaSymmetry,'AS')
                    this.inertiaSymmetry = [this.inertiaSymmetry '3'];    % AS3
                end
            end
            
        end
        %% Plot the fit, plotting the polhode in the estimated body frame in 3D, including the
        %  the axes of frame E, the axes of estimated frame B, and the conic fits on each plane.
        %  Call after optimize()
        function plot(this,suppressB,suppressConicFits1)
           
            % Plot polhode
            plot3(this.omegaB_B(1,:),this.omegaB_B(2,:),this.omegaB_B(3,:),'Color',[0.5 0 0.5], ...
                    'LineWidth',1);
            hold on; grid on; axis equal;
                
            % Plot axes of frame E
            REtoB = this.RBtoG' * this.REtoG;
            sE = max(max(this.omegaB_G)) * 0.75;         % axes scale for E
            plot3([0 REtoB(1,1)*sE],[0 REtoB(2,1)*sE],[0 REtoB(3,1)*sE],...
                    '--r','LineWidth',2);
            plot3([0 REtoB(1,2)*sE],[0 REtoB(2,2)*sE],[0 REtoB(3,2)*sE],...
                    '--g','LineWidth',2);
            plot3([0 REtoB(1,3)*sE],[0 REtoB(2,3)*sE],[0 REtoB(3,3)*sE],...
                    '--b','LineWidth',2);
            
            sB = max(max(this.omegaB_G)) * 0.5;         % axes scale for B
            if nargin < 2 || ~suppressB
                % Plot axes of frame B
                plot3([0 sB],[0 0],[0 0],'-.r','LineWidth',4);
                plot3([0 0],[0 sB],[0 0],'-.g','LineWidth',4);
                plot3([0 0],[0 0],[0 sB],'-.b','LineWidth',4);
            end
            
            % Threshold fits for points too far away from from angular velocities
            omegaPlotMax = 1.5 * sqrt( max( this.omegaB_B(1,:).^2 + this.omegaB_B(2,:).^2 + ...
                                                  this.omegaB_B(3,:).^2 ) );
            
            if nargin < 3 || ~suppressConicFits1
                % Get the points for each initial conic fit
                pts1Plane12 = this.conicFits1{1}.plot('actual','-k',1,1);
                pts1Plane23 = this.conicFits1{2}.plot('actual','-k',1,1);
                pts1Plane13 = this.conicFits1{3}.plot('actual','-k',1,1);            

                % Cut off pts too far away from angular velocities                
                pts1Plane12 = pts1Plane12(:,pts1Plane12(1,:).^2+pts1Plane12(2,:).^2<omegaPlotMax^2);
                pts1Plane23 = pts1Plane23(:,pts1Plane23(1,:).^2+pts1Plane23(2,:).^2<omegaPlotMax^2);
                pts1Plane13 = pts1Plane13(:,pts1Plane13(1,:).^2+pts1Plane13(2,:).^2<omegaPlotMax^2);

                % Change into 3D coordinates in the E frame
                pts1Plane12 = [pts1Plane12(1,:); pts1Plane12(2,:); zeros(1,size(pts1Plane12,2))];
                pts1Plane23 = [zeros(1,size(pts1Plane23,2)); pts1Plane23(1,:); pts1Plane23(2,:)];
                pts1Plane13 = [pts1Plane13(1,:); zeros(1,size(pts1Plane13,2)); pts1Plane13(2,:)];

                % Rotate into estimated B frame for plotting
                pts1Plane12 = REtoB*pts1Plane12;
                pts1Plane23 = REtoB*pts1Plane23;
                pts1Plane13 = REtoB*pts1Plane13;
                
                % Force projections to the wall                
                pts1Plane12 = this.projectToWall(pts1Plane12,-omegaPlotMax);
                pts1Plane23 = this.projectToWall(pts1Plane23,-omegaPlotMax);
                pts1Plane13 = this.projectToWall(pts1Plane13,-omegaPlotMax);

                % Plot the conic fits in the correct planes
                plot3(pts1Plane12(1,:),pts1Plane12(2,:),pts1Plane12(3,:),'.m','MarkerSize',3, ...
                    'Color',[0.5 0.5 0.5]);
                plot3(pts1Plane23(1,:),pts1Plane23(2,:),pts1Plane23(3,:),'.m','MarkerSize',3, ...
                    'Color',[0.5 0.5 0.5]);
                plot3(pts1Plane13(1,:),pts1Plane13(2,:),pts1Plane13(3,:),'.m','MarkerSize',3, ...
                    'Color',[0.5 0.5 0.5]);
            end
            
            if nargin < 2 || ~suppressB
                % Get and plot points from final conic fit
                pts2Plane12 = this.conicFits2{1}.plot('actual','-k',1,1);
                pts2Plane23 = this.conicFits2{2}.plot('actual','-k',1,1);
                pts2Plane13 = this.conicFits2{3}.plot('actual','-k',1,1);
                pts2Plane12 = pts2Plane12(:,pts2Plane12(1,:).^2 + pts2Plane12(2,:).^2 < omegaPlotMax^2);
                pts2Plane23 = pts2Plane23(:,pts2Plane23(1,:).^2 + pts2Plane23(2,:).^2 < omegaPlotMax^2);
                pts2Plane13 = pts2Plane13(:,pts2Plane13(1,:).^2 + pts2Plane13(2,:).^2 < omegaPlotMax^2);
                pts2Plane12 = [pts2Plane12(1,:); pts2Plane12(2,:); zeros(1,size(pts2Plane12,2))];
                pts2Plane23 = [zeros(1,size(pts2Plane23,2)); pts2Plane23(1,:); pts2Plane23(2,:)];
                pts2Plane13 = [pts2Plane13(1,:); zeros(1,size(pts2Plane13,2)); pts2Plane13(2,:)];
                pts2Plane12 = REtoB*pts2Plane12; 
                pts2Plane23 = REtoB*pts2Plane23;
                pts2Plane13 = REtoB*pts2Plane13;
                pts2Plane12 = this.projectToWall(pts2Plane12,-omegaPlotMax);
                pts2Plane23 = this.projectToWall(pts2Plane23,-omegaPlotMax);
                pts2Plane13 = this.projectToWall(pts2Plane13,-omegaPlotMax);
                plot3(pts2Plane12(1,:),pts2Plane12(2,:),pts2Plane12(3,:),'.k','MarkerSize',3);
                plot3(pts2Plane23(1,:),pts2Plane23(2,:),pts2Plane23(3,:),'.k','MarkerSize',3);
                plot3(pts2Plane13(1,:),pts2Plane13(2,:),pts2Plane13(3,:),'.k','MarkerSize',3);
                xlim([-omegaPlotMax omegaPlotMax]); ylim([-omegaPlotMax omegaPlotMax]); 
                zlim([-omegaPlotMax omegaPlotMax]);
            end
            view(140,30);
            
        end
        %% Take the dimension with zero data and change it to the value of dist
        function pts = projectToWall(this, pts, dist)
            
            [~,zeroDim] = min(sum(abs(pts),2));
            pts(zeroDim,:) = dist;
            
        end
    end
    
end

