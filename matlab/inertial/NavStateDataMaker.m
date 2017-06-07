classdef NavStateDataMaker < handle
    %NAVSTATEDATAMAKER Takes the final isam results from TP_1510 and creates the best navState_data.txt
    
    properties
        resultsDir;         % results directory for TP_1510
        globMetData;        % global metrology data to align to
        isamData;           % final isam data to be realigned
        navStateData;       % aligned navigation state data for saving as csv
        PWtoBisam;          % original isam poses
        vB_Wisam;           % velocity in isam
        PWtoBalign;         % aligned isam poses
        vB_Walign;          % aligned isam velocities
        PWtoBgm;            % global metrology poses
        t;                  % timesteps
        nT;                 % number of timesteps 
        col = [[0 0   1]' [1 0.6 0]' [0.5 0 0.5]' [0 1 0]' [0 1 1]']; % colors
    end
    
    methods
        %% Constructor taking the TP_1510 results data directory
        function this = NavStateDataMaker(resultsDir)
            
            this.resultsDir     = resultsDir;
            this.globMetData    = csvread([resultsDir '\groundTruth_data.txt']);
            this.isamData       = csvread([resultsDir '\gtsamIncrementalFinal_data.txt']);
            this.navStateData   = zeros(size(this.isamData));
            this.t              = this.isamData(:,1)/1000;
            this.nT             = size(this.isamData,1);
            this.PWtoBisam      = zeros(4,4,this.nT);
            this.vB_Wisam       = zeros(3,this.nT);
            this.PWtoBgm        = zeros(4,4,this.nT);
            
            
            % Get the poses
            for i = 1:this.nT           
                this.PWtoBisam(:,:,i)   = reshape(this.isamData(i,2:17), 4, 4)';
                this.vB_Wisam(:,i)      = this.isamData(i,18:20)';
                this.PWtoBgm(:,:,i)     = reshape(this.globMetData(i,2:17), 4, 4)';                
                
            end
            
            % Initialize
            this.PWtoBalign = this.PWtoBisam;
            this.vB_Walign = this.vB_Wisam;
            
            % Align the data
            this.align();
            
            % Plot
            figure(1); clf;
            this.plot();
            
            % Format the nav state data
            for i = 1:this.nT
                this.navStateData(i,:) = [this.isamData(i,1) ...
                        reshape(this.PWtoBalign(:,:,i)',1,16) this.vB_Walign(:,i)' ...
                                    this.isamData(i,21:end)];
            end
            dlmwrite('navState_data.txt',this.navStateData,'precision',6);
            
        end
        %% Align the data by finding the pose change that minimizes the rotation error
        function align(this)
            
            % Minimize rotation error by finding RVtoW = Exp(theta) to best minimize rot error
            theta0 = zeros(3,1);
            [thetaOpt,costThetaOpt,~,~,~,~,~] = fmincon(@(theta) this.rotCost(theta), theta0);
            RVtoW = Exp(thetaOpt);
            disp(['Optimal rotation cost: ' num2str(costThetaOpt)]);
            disp(['Optimal rotation (' num2str(rad2deg(norm(thetaOpt))) ' deg):']);
            disp(RVtoW);
            
            % Align data from using the found rotation PWtoB = [RVtoW*RBtoV RVtoW*tVtoB_V; 0 0 0 1]
            for i = 1:this.nT
                this.PWtoBalign(:,:,i) = ...
                    [RVtoW*this.PWtoBisam(1:3,1:3,i) RVtoW*this.PWtoBisam(1:3,4,i); 0 0 0 1];
                this.vB_Walign(:,i) = RVtoW * this.vB_Wisam(:,i);
            end
            
            % Translate to the same starting translation as global metrology
            tOpt = this.PWtoBgm(1:3,4,1) - this.PWtoBalign(1:3,4,1);
            for i = 1:this.nT
                this.PWtoBalign(1:3,4,i) = this.PWtoBalign(1:3,4,i) + tOpt;
            end
            
            %{
            % Minimize translation error of rotated data
            t0 = zeros(3,1);            
            [tOpt,costTOpt,~,~,~,~,~] = fmincon(@(t) this.transCost(t), t0);
            disp(['Optimal translation cost: ' num2str(costTOpt)]);
            disp('Optimal translation:');
            disp(tOpt);           
            
            % Align data translationally using the found translation
            for i = 1:this.nT
                this.PWtoBalign(1:3,4,i) = this.PWtoBalign(1:3,4,i) + tOpt;
            end
            %}
            
        end
        %% Cost of rotation given angle-axis vector theta
        function Cost = rotCost(this, theta)
           
            Cost = 0;
            for i = 1:this.nT
                RBtoWisam = this.PWtoBisam(1:3,1:3,i);
                RBtoWpred = Exp(theta) * RBtoWisam;
                RBtoWgm   = this.PWtoBgm(1:3,1:3,i);
                Cost = Cost + norm(Log(RBtoWpred'*RBtoWgm))^2;
            end
            
        end
        %% Cost of translation given transltion t
        function Cost = transCost(this, t)
            
            Cost = 0;
            for i = 1:this.nT
                tWtoB_Walign    = this.PWtoBalign(1:3,4,i);
                tWtoB_Wpred     = tWtoB_Walign + t;
                tWtoB_Wgm       = this.PWtoBgm(1:3,4,i);
                Cost = Cost + norm(tWtoB_Wgm - tWtoB_Wpred)^2;
            end
            
        end
        %% Plot the global metrology data, the aligned data, and  the original data
        function plot(this)
            
            % Setup
            subplot(3,2,[1 3 5]);
            Ps(:,:,:,1) = this.PWtoBgm;
            Ps(:,:,:,2) = this.PWtoBalign;
            Ps(:,:,:,3) = this.PWtoBisam;
            items = {'Global Metrology','Nav State','TP 1510 iSAM'};
            
            N  = size(Ps,3);    % maximum number of poses
            F  = size(Ps,4);    % num of frames
            lt = 0.025;         % length of the triad axes
            plines = [];        % holds plot handles for connecting lines
            PItoD = [-1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 1]; % pose change from default to ISS display frame
            axesLimits = [-1.5 1.5 -1.5 1.5 -1.5 1.5]; % [xmin xmax ymin ymax zmin zmax]

            for i = 1:F                     % loop through frames

                numNonZeroPoses = 0;        % keeps track of number of nonZero Poses

                for k = 1:N                 % loop through pose chain

                    % Transform into the the ISS frame to see it from the same angle as in the OVHD videos
                    P = PItoD * Ps(:,:,k,i);

                    if max(max(P)) == 0
                        continue;
                    end

                    % Get triad [ x_k y_k z_k ] and translation
                    numNonZeroPoses = numNonZeroPoses + 1;
                    tri = P(1:3,1:3);
                    t   = P(1:3,4);

                    % Plot marker
                    plot3(t(1),t(2),t(3),'.','MarkerEdgeColor',this.col(:,i),'MarkerFaceColor',this.col(:,i));
                    hold on;

                    % Plot triad
                    plot3([t(1) t(1)+tri(1,1)*lt],[t(2) t(2)+tri(2,1)*lt],[t(3) t(3)+tri(3,1)*lt],'-r');
                    plot3([t(1) t(1)+tri(1,2)*lt],[t(2) t(2)+tri(2,2)*lt],[t(3) t(3)+tri(3,2)*lt],'-g');
                    plot3([t(1) t(1)+tri(1,3)*lt],[t(2) t(2)+tri(2,3)*lt],[t(3) t(3)+tri(3,3)*lt],'-b');

                    if numNonZeroPoses > 1
                        % Plot connecting line
                        tLast = PItoD(1:3,1:3) * Ps(1:3,4,k-1,i);
                        plines(i) = plot3([tLast(1) t(1)],[tLast(2) t(2)],[tLast(3) t(3)],'-', ...
                            'DisplayName',items{i},'Color',this.col(:,i));
                        % Plot projection onto x-y plane
                        plot3([tLast(1) t(1)],[tLast(2) t(2)],[axesLimits(5) axesLimits(5)],'-', ...
                            'Color',[min(this.col(1,i)+0.2,1) min(this.col(2,i)+0.2,1) min(this.col(3,i)+0.2,1)]);
                        % Plot projection onto x-z plane
                        plot3([tLast(1) t(1)],[axesLimits(4) axesLimits(4)],[tLast(3) t(3)],'-', ...
                            'Color',[min(this.col(1,i)+0.2,1) min(this.col(2,i)+0.2,1) min(this.col(3,i)+0.2,1)]);
                        % Plot projection onto y-z plane
                        plot3([axesLimits(1) axesLimits(1)],[tLast(2) t(2)],[tLast(3) t(3)],'-', ...
                            'Color',[min(this.col(1,i)+0.2,1) min(this.col(2,i)+0.2,1) min(this.col(3,i)+0.2,1)]);            
                    else
                        % Plot start position marker
                        plot3(t(1),t(2),t(3),'o','MarkerEdgeColor',[0 1 0], ...
                            'MarkerFaceColor',[0 1 0],'MarkerSize',10);
                    end   
                end

                % Plot current position marker
                plot3(t(1),t(2),t(3),'o','MarkerEdgeColor',this.col(:,i), ...
                            'MarkerFaceColor',this.col(:,i),'MarkerSize',10);
            end


            % Setup axes and grid
            axis equal;
            axis(axesLimits);
            grid on;
            set(gca,'xtick',axesLimits(1):0.25:axesLimits(2));
            set(gca,'ytick',axesLimits(3):0.25:axesLimits(4));
            set(gca,'ztick',axesLimits(5):0.25:axesLimits(6));
            xlabel('-{}^Wx [m]'); ylabel('{}^Wy [m]'); zlabel('-{}^Wz [m]');
            view(30,30);
            leg = legend(plines);
            set(leg,'Location','NorthWest');
            title('Poses');
            
            errPos = zeros(size(Ps,3),2); errRot = zeros(size(Ps,3),2);
            for j = 1:this.nT
                errPos(j,1) = norm( Ps(1:3,4,j,2) - Ps(1:3,4,j,1) );
                errPos(j,2) = norm( Ps(1:3,4,j,3) - Ps(1:3,4,j,1) );
                errRot(j,1) = rad2deg( norm( Log( Ps(1:3,1:3,j,1)' * Ps(1:3,1:3,j,2) ) ) );   
                errRot(j,2) = rad2deg( norm( Log( Ps(1:3,1:3,j,1)' * Ps(1:3,1:3,j,3) ) ) );              
            end
                      
            disp(['Average remaining rotational error: ' num2str(mean(errRot(:,1))) 'deg']);
            
            subplot(3,2,2);
            plot(this.t,errPos(:,1),'Color',this.col(:,2)); title('Position Error'); hold on;
            plot(this.t,errPos(:,2),'Color',this.col(:,3));
            ylabel('Error [m]'); grid on;
            subplot(3,2,4);
            plot(this.t,errRot(:,1),'Color',this.col(:,2)); title('Orientation Error');  hold on;
            plot(this.t,errRot(:,2),'Color',this.col(:,3)); 
            ylabel('Error [deg]'); grid on; xlabel('Time [s]');
            
            
        end
        
        
        
    end
    
end

