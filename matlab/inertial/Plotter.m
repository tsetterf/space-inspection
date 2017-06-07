classdef Plotter < handle
    %PLOTTER Plots the results from TP_1511_DynamicPoseSLAM (uses vertigo time)
    
    properties
        % Results
        r;
        % Colors        
        col = [[0 0   1]' [1 0.6 0]' [0.5 0 0.5]' [0 1 0]' [0 1 1]'];
    end
    
    methods
        %% Constructor using a Results object
        function this = Plotter(results)
            
            % Store the results
            this.r = results;
            
        end
        
        %% Plot pose chains
        function plotPoseChains( this, items, stateNum, axesLimits )
        %  Plots pose chains as if viewed from the ISS OVHD-PORT-AFT corner
        %  inputs:  Ps  =   [4x4xNxF] double, ...) where N is the length of the longest pose chain
        %                   and F is the number of frames. For unequal length pose chains, leave 
        %                   leading or trailing matrices of zeros(4).
        %           items = choose from some or all of {'PWtoBgm','PWtoBtgm','PWtoGgm','PGtoBgm',
        %                   'PWtoBbsam','PWtoGbsam','PGtoBbsam','PWtoBisam','PWtoGisam','PGtoBisam'}
        %           axesLimits = [xmin xmax ymin ymax zmin zmax] for plotting
        
        % Get index of state number
        ind = find(this.r.stateNums==stateNum);
        
        % Setup
        Ps = [];
        for i = 1:length(items)
            if strcmp(items{i},'PWtoBisam') || strcmp(items{i},'PWtoGisam')
                eval(['Ps(:,:,:,' num2str(i) ') = this.r.' items{i} '(:,:,:,ind);']); 
            else               
                eval(['Ps(:,:,:,' num2str(i) ') = this.r.' items{i} ';']);
            end
        end
                
        N  = length(this.r.stateNums(this.r.stateNums<=stateNum));   % maximum number of poses
        F  = size(Ps,4);    % num of frames
        lt = 0.025;         % length of the triad axes
        plines = [];        % holds plot handles for connecting lines
        PItoD = [-1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 1]; % pose change from default to ISS display frame
        if nargin < 4
            axesLimits = [-1 1 -1 1 -1 1]; % [xmin xmax ymin ymax zmin zmax]
        end
                
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
                    
            % Plot target centroid position
            plot3(this.r.PWtoBtgm(1,4,N),this.r.PWtoBtgm(2,4,N),this.r.PWtoBtgm(3,4,N),'o', ...
                'MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',10);
            
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
        set(leg,'Location','North');
%         title('Poses');
%         drawnow;

        end
        
        %% Plot errors given the pose item symbol of the truth (e.g. 'PWtoBvio') and the one to
        %  compare (e.g. 'PWtoBisam'), or pass the actual Ptruth and Ptest data with times
        function plotPosError( this, itemTruth, itemTest, stateNum, times, lineSpec )
            
            % Get index of state number
            ind = find(this.r.stateNums==stateNum);
            
            if ischar(itemTruth) && ischar(itemTest)
                eval(['Ptruth = this.r.' itemTruth ';']);
                if eval(['size(this.r.' itemTest ',4) == 1'])
                    eval(['Ptest  = this.r.' itemTest  '(:,:,:);']);
                else
                    eval(['Ptest  = this.r.' itemTest  '(:,:,:,ind);']);
                end
            else
                Ptruth = itemTruth;
                Ptest  = itemTest;
            end
            
            if nargin < 5
                times  = this.r.tInspState(this.r.stateNums <= stateNum);% - this.r.tInspState(1);
            else
                times = times;% - times(1);
            end
            
            % Get the comparisons between the two
            errPos = zeros(length(times),1);
            for j = 1:length(times)
                errPos(j) = norm( Ptest(1:3,4,j) - Ptruth(1:3,4,j) );
            end
            plot(times,errPos,lineSpec); title('Position Error');
            xlim([times(1) times(end)]);
            ylabel('Error [m]'); grid on; hold on;
            
            disp(['Avg pos err: ' num2str(mean(errPos)) ' m']);
            
        end
        
        %% Plot errors given the pose item symbol of the truth (e.g. 'PWtoBvio') and the one to
        %  compare (e.g. 'PWtoBisam') or pass the actual Ptruth and Ptest data with times
        function plotRotError( this, itemTruth, itemTest, stateNum, times, lineSpec )
                        
            % Get index of state number
            ind = find(this.r.stateNums==stateNum);
            
            if ischar(itemTruth) && ischar(itemTest)
                eval(['Ptruth = this.r.' itemTruth ';']);
                if eval(['size(this.r.' itemTest ',4) == 1'])
                    eval(['Ptest  = this.r.' itemTest  '(:,:,:);']);
                else
                    eval(['Ptest  = this.r.' itemTest  '(:,:,:,ind);']);
                end
            else
                Ptruth = itemTruth;
                Ptest  = itemTest;             
            end
                        
            if nargin < 5
                times  = this.r.tInspState(this.r.stateNums <= stateNum); %- this.r.tInspState(1);
            else
                times = times; % - times(1);
            end
            
            % Get the comparisons between the two
            errRot = zeros(length(times),1);
            for j = 1:length(times)
                errRot(j) = rad2deg(norm( Log(Ptruth(1:3,1:3,j)' * Ptest(1:3,1:3,j)) ));
            end
            plot(times,errRot,lineSpec); %title('Orientation Error');
            xlim([times(1) times(end)]);
            ylabel('Error [deg]'); grid on; hold on;
                        
            disp(['Avg rot err: ' num2str(mean(errRot)) ' deg']);
        end
        
        %% Plot errors given the item symbol of the truth (e.g. 'vio') and the one to
        %  compare (e.g. 'isam') or pass the actual vTruth and vTest data with times
        function plotVelError( this, itemTruth, itemTest, stateNum, times, lineSpec )
            
            % Get index of state number
            ind = find(this.r.stateNums==stateNum);
            
            if ischar(itemTruth) && ischar(itemTest)
                eval(['vTruth = this.r.vB_W' itemTruth ';']);
                if eval(['size(this.r.vB_W' itemTest ',3) == 1'])
                    eval(['vTest  = this.r.vB_W' itemTest  ';']);
                else
                    eval(['vTest  = this.r.vB_W' itemTest  '(:,:,ind);']);
                end
            else
                vTruth = itemTruth;
                vTest  = itemTest;
            end
            
            if nargin < 5
                times  = this.r.tInspState(this.r.stateNums <= stateNum);% - this.r.tInspState(1);
            else
                times = times;% - times(1);
            end
            
            % Get the comparisons between the two
            errVel = zeros(length(times),1);
            for j = 1:length(times)
                errVel(j) = norm( vTest(:,j) - vTruth(:,j) );
            end
            plot(times,errVel,lineSpec); title('Velocity Error');
            xlim([times(1) times(end)]); 
            ylabel('Error [m/s]'); grid on; hold on;
            disp(['Avg vel err: ' num2str(mean(errVel)) ' m/s' ]);
        end
        %% Plot errors given the item symbol of the truth (e.g. 'vio') given omegaTruth and 
        %  omegaTest data with times
        function plotAngVelError( this, itemTruth, itemTest, stateNum, times, lineSpec )
            
            omegaTruth = itemTruth;
            omegaTest  = itemTest;
            
            % Get the comparisons between the two
            errOmega = zeros(length(times),1);
            for j = 1:length(times)
                errOmega(j) = rad2deg( norm( omegaTest(:,j) - omegaTruth(:,j) ) );
            end
            plot(times,errOmega,lineSpec); %title('Angular Velocity Error');
            xlim([times(1) times(end)]); ylim([0 max(errOmega)*1.25]);
            ylabel('Error [deg/s]'); grid on; hold on;
            disp(['Avg ang vel err: ' num2str(mean(errOmega)) ' deg/s' ]);
        end
        %% Get the state number associated with an image number
        function stateNum = getStateNum( this, imageNum )
            stateImageNums = [this.r.stateNums this.r.imageNums];
            [~,iNearest] = min(abs(stateImageNums(:,2)-imageNum));
            stateNum = this.r.stateNums(iNearest);
        end
        
        %% Plot image given a state number
        function plotStateImage( this, stateNum )            
            this.plotImage(this.r.imageNums(stateNum+1));            
        end
             
        %% Plot image given an image number
        function plotImage( this, imageNum )
            
            motion   = im2double(imread([this.r.resultsDir '\MotionImage' ...
                                            num2str(imageNum) '.bmp']));                     
            imshow(motion); axis equal;
            
        end
        
        %% Plot map of the SURF features
        function plotMap( this, stateNum, RGitoW, tGtoBt_G, ptScale )
                        
            
            if nargin < 5
                ptScale = 1;
            end
                        
            % Get all of the data up to the stateNum
            featMap         = this.r.featureMap(this.r.featureMap(:,1) <= stateNum,2:4);
            featMapAll      = this.r.featureMapAll(this.r.featureMapAll(:,1) <= stateNum,2:4);
            firstInliers    = this.r.firstInliers(:,2:4);
            firstInliers    = repmat( - this.r.tGtoBt_Gisam(:,stateNum+1) ...
                                + this.r.PGtoBisam(1:3,4,1,stateNum+1), 1, size(firstInliers,1)) ...
                                + this.r.PGtoBisam(1:3,1:3,1,stateNum+1) * firstInliers';
            firstInliers    = firstInliers';
            
            featMapCols     = linspace(0,1,size(featMap,1));
            
            % Rotate into world frame
            featMap      = ( RGitoW * featMap')';
            featMapAll   = ( RGitoW * featMapAll')';
            firstInliers = ( RGitoW * firstInliers')';
            tGtoBt_W     = RGitoW * tGtoBt_G;
                        
            % Get axes limits
            if this.r.isSimulated
                plotLims = 0.2;
            else
                plotLims = 0.15;
            end
            
            % Get geometric axis triad directions
            xGi_W = RGitoW(:,1) * 0.05; 
            yGi_W = RGitoW(:,2) * 0.05; 
            zGi_W = RGitoW(:,3) * 0.05;
            
            if ptScale > 0.2
                showPointCloud(featMapAll,[0.5 0.5 0.5],'MarkerSize',0.5*ptScale); hold on;
            end
            showPointCloud(featMap,featMapCols,'MarkerSize',30*ptScale); hold on;           
            showPointCloud(firstInliers,[0.5 0 0.5],'MarkerSize',max(10,20*ptScale));
            plot3([-tGtoBt_W(1) 0],[-tGtoBt_W(2) 0],[-tGtoBt_W(3) 0],'-k','LineWidth',max(1,2*ptScale));
            plot3([-tGtoBt_W(1) -tGtoBt_W(1)+xGi_W(1)],[-tGtoBt_W(2) -tGtoBt_W(2)+xGi_W(2)], ...
                [-tGtoBt_W(3) -tGtoBt_W(3)+xGi_W(3)],'-.r','LineWidth',max(1,2*ptScale));
            plot3([-tGtoBt_W(1) -tGtoBt_W(1)+yGi_W(1)],[-tGtoBt_W(2) -tGtoBt_W(2)+yGi_W(2)], ...
                [-tGtoBt_W(3) -tGtoBt_W(3)+yGi_W(3)],'-.g','LineWidth',max(1,2*ptScale));
            plot3([-tGtoBt_W(1) -tGtoBt_W(1)+zGi_W(1)],[-tGtoBt_W(2) -tGtoBt_W(2)+zGi_W(2)], ...
                [-tGtoBt_W(3) -tGtoBt_W(3)+zGi_W(3)],'-.b','LineWidth',max(1,2*ptScale));
            scatter3(0,0,0,'b');
            axis equal; 
%             axis([xmin xmax ymin ymax zmin zmax]);
            axis([-plotLims plotLims -plotLims plotLims -plotLims plotLims]);
%             title('Map');
            
            xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
            
        end
        
        %% Plot map of the SURF features given the current isam estimate
        function plotMapCurrent( this, stateNum, bodyFrameView, ptScale )
            
            if nargin < 3
                bodyFrameView = 0;
            end
            
            if nargin < 4                
                ptScale = 1;
            end
            
            % Get index of this state number
            ind     = find(this.r.stateNums == stateNum);
            
            % Get current estimates of required rotations and translations
            RB0toG0     = this.r.PGtoBisam(1:3,1:3,1,ind);
            RGitoW      = this.r.PWtoGisam(1:3,1:3,ind,ind);
            RBitoW      = this.r.PWtoBisam(1:3,1:3,ind,ind);
            tGtoB_G     = this.r.PGtoBisam(1:3,4,1,ind);
            tGtoBt_G    = this.r.tGtoBt_Gisam(:,ind);
            tGtoBt      = RGitoW * tGtoBt_G;
                        
            % Get the position of the first inliers in the world frame
            firstInliers = this.r.firstInliers(:,2:4);          % B0 frame
            firstInliers = repmat( -tGtoBt_G + tGtoB_G, 1, size(firstInliers,1)) ...
                                + RB0toG0 * firstInliers';      % Gc frame
            firstInliers = (RGitoW * firstInliers)';            % W frame
            
            % Loop through data and put it into the world frame according to current estimates
            featMap = []; featMapAll = []; featMapCols = [];
            for i = 1:ind
                
                tBtoP_B     = this.r.stereoFrames{i};
                tBtoP_Ball  = this.r.stereoFramesAll{i};
                Nin         = size(tBtoP_B,2);
                Nall        = size(tBtoP_Ball,2);
                RBitoGi     = this.r.PGtoBisam(1:3,1:3,i,ind);
                tGitoBi_Gi  = this.r.PGtoBisam(1:3,4,i,ind);
                
                tBttoP_Gc   = repmat(-tGtoBt_G + tGitoBi_Gi,1,Nin) + RBitoGi * tBtoP_B;
                tBttoP_GcAll= repmat(-tGtoBt_G + tGitoBi_Gi,1,Nall) + RBitoGi * tBtoP_Ball;
                
                featMap     = [featMap    ; (RGitoW * tBttoP_Gc)'   ];
                featMapAll  = [featMapAll ; (RGitoW * tBttoP_GcAll)'];
                featMapCols = [featMapCols repmat(i,1,Nin)];
                
            end
            
            % Get current inliers for larger display
            RBitoGi     = this.r.PGtoBisam(1:3,1:3,ind,ind);
            tGitoBi_Gi  = this.r.PGtoBisam(1:3,4,ind,ind); 
            tBtoP_B     = this.r.stereoFrames{ind};         
            Nin = size(tBtoP_B,2);
            tBttoP_Gc   = repmat(-tGtoBt_G + tGitoBi_Gi,1,Nin) + RBitoGi * tBtoP_B;
            currentInliers =  (RGitoW * tBttoP_Gc)';
            
            
            % Get geometric axis triad directions
            xGi = RGitoW(:,1) * 0.05;  yGi = RGitoW(:,2) * 0.05;  zGi = RGitoW(:,3) * 0.05;
            
            % Rotate to body frame view if desired
            if bodyFrameView
                firstInliers = (RBitoW' * firstInliers')';
                featMap      = (RBitoW' * featMap')';
                featMapAll   = (RBitoW' * featMapAll')';
                currentInliers =  (RBitoW' * currentInliers')';               
                tGtoBt       = RBitoW' * tGtoBt;
                xGi          = RBitoW' * xGi;
                yGi          = RBitoW' * yGi;
                zGi          = RBitoW' * zGi;
            end
                        
            % Get axes limits
            if this.r.isSimulated
                plotLims = 0.2;
            else
                plotLims = 0.15;
            end
            
            
            if ptScale > 0.2
                showPointCloud(featMapAll,[0.5 0.5 0.5],'MarkerSize',0.5*ptScale); hold on;
            end
            showPointCloud(featMap,featMapCols,'MarkerSize',30*ptScale); hold on;           
            showPointCloud(firstInliers,[0.5 0 0.5],'MarkerSize',max(10,20*ptScale));
            showPointCloud(currentInliers,[0 1 0],'MarkerSize',max(10,15*ptScale));
            plot3([-tGtoBt(1) 0],[-tGtoBt(2) 0],[-tGtoBt(3) 0],'-k','LineWidth',max(1,2*ptScale));
            plot3([-tGtoBt(1) -tGtoBt(1)+xGi(1)],[-tGtoBt(2) -tGtoBt(2)+xGi(2)], ...
                [-tGtoBt(3) -tGtoBt(3)+xGi(3)],'-.r','LineWidth',max(1,2*ptScale));
            plot3([-tGtoBt(1) -tGtoBt(1)+yGi(1)],[-tGtoBt(2) -tGtoBt(2)+yGi(2)], ...
                [-tGtoBt(3) -tGtoBt(3)+yGi(3)],'-.g','LineWidth',max(1,2*ptScale));
            plot3([-tGtoBt(1) -tGtoBt(1)+zGi(1)],[-tGtoBt(2) -tGtoBt(2)+zGi(2)], ...
                [-tGtoBt(3) -tGtoBt(3)+zGi(3)],'-.b','LineWidth',max(1,2*ptScale));
            scatter3(0,0,0,'b');
            axis equal; 
%             axis([xmin xmax ymin ymax zmin zmax]);
            axis([-plotLims plotLims -plotLims plotLims -plotLims plotLims]);
%             title('Map');
            
            xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
            
        end
        
        %% Plot all projections of a map
        function plotMapOrtho( this, stateNum )
            
%             subplot = @(m,n,p) subtightplot (m, n, p, [0.075 0.01], [0.05 0.01], [0.01 0.01]);
            subInds = {3,6,9,[1 2 4 5 7 8]};
            subViews = {[90 0],[0,0],[0,90],[-60,30]};
            ptScales = [0.1 0.1 0.1 1];
            for i = 1:4
                sp = subplot(3,3,subInds{i});
                sB = 0.2;
                if this.r.isSimulated
                    RBtoBhat = this.r.RBttoGisam(:,:,stateNum+1)'*this.r.RBttoGtrue;
                    tBhattoB_Bhat = -this.r.RBttoGisam(:,:,stateNum+1)'*this.r.tGtoBt_Gisam(:,stateNum+1) ...
                        + this.r.RBttoGisam(:,:,stateNum+1)'*this.r.tGtoBt_Gtrue;
%                     tBhattoB_Bhat = -this.r.tGtoBt_Gisam(:,stateNum+1) ...
%                         + this.r.tGtoBt_Gtrue;
                    tB = tBhattoB_Bhat;
                    xB_Bhat = RBtoBhat(:,1); yB_Bhat = RBtoBhat(:,2); zB_Bhat = RBtoBhat(:,3);
%                     plot3([0 sB*xB_Bhat(1)],[0 sB*xB_Bhat(2)],[0 sB*xB_Bhat(3)],'-r','LineWidth',max(1,4*ptScales(i))); hold on; grid on;
%                     plot3([0 sB*yB_Bhat(1)],[0 sB*yB_Bhat(2)],[0 sB*yB_Bhat(3)],'-g','LineWidth',max(1,4*ptScales(i)));
%                     plot3([0 sB*zB_Bhat(1)],[0 sB*zB_Bhat(2)],[0 sB*zB_Bhat(3)],'-b','LineWidth',max(1,4*ptScales(i)));
                    plot3([tB(1) tB(1)+sB*xB_Bhat(1)],[tB(2) tB(2)+sB*xB_Bhat(2)], ...
                        [tB(3) tB(3)+sB*xB_Bhat(3)],'-r','LineWidth',max(1,4*ptScales(i))); hold on; grid on;
                    plot3([tB(1) tB(1)+sB*yB_Bhat(1)],[tB(2) tB(2)+sB*yB_Bhat(2)], ...
                        [tB(3) tB(3)+sB*yB_Bhat(3)],'-g','LineWidth',max(1,4*ptScales(i)));
                    plot3([tB(1) tB(1)+sB*zB_Bhat(1)],[tB(2) tB(2)+sB*zB_Bhat(2)], ...
                        [tB(3) tB(3)+sB*zB_Bhat(3)],'-b','LineWidth',max(1,4*ptScales(i)));
                end
                plot3([0 sB],[0 0],[0 0],'--r','LineWidth',max(1,4*ptScales(i))); hold on; grid on;
                plot3([0 0],[0 sB],[0 0],'--g','LineWidth',max(1,4*ptScales(i)));
                plot3([0 0],[0 0],[0 sB],'--b','LineWidth',max(1,4*ptScales(i)));
                this.plotMap(stateNum, this.r.RBttoGisam(:,:,stateNum+1)', ...
                                            this.r.tGtoBt_Gisam(:,stateNum+1),ptScales(i));
                v = subViews{i};
                view(v(1),v(2));
                
                % Increase height and width
                if i < 4
                    p = get(sp,'position');
                    p(3:4) = p(3:4)*1.5;    
                    set(sp, 'position', p);
                end
                
            end
            if this.r.isSimulated
                legend('True x_B','True y_B','True x_B','Estimated x_B','Estimated y_B', ...
                    'Estimated z_B','All Features','Inliers','First Frame Inliers', ...
                    'Estimated t_{GtoB}','Estimated x_G', 'Estimated y_G','Estimated z_G', ...
                    'Location','East');
            else
                legend('Estimated x_B','Estimated y_B', ...
                    'Estimated z_B','All Features','Inliers','First Frame Inliers', ...
                    'Estimated t_{GtoB}','Estimated x_G',...
                    'Estimated y_G','Estimated z_G','Location','East');
            end
            
        end
        
        %% Plot the angular velocities of the body in the geometric frame in 3D
        function plotPolhode( this, type, stateNum )
            
            eval(['PGtoB = this.r.PGtoB' type ';']);

            % Get rotations
            RB0toG0 = PGtoB(1:3,1:3,1,stateNum+1);
            RWtoB0 = this.r.PWtoBgm(1:3,1:3,1)';
            RBt0toW = this.r.PWtoBtgm(1:3,1:3,1);
            
            % Get angular velocities
            [~,ind]       = min(abs(this.r.omegaStateNums-stateNum));
            tOmega        = this.r.tOmegaGisam(:,ind)';
            omegaG_Gisam  = this.r.omegaBt_Gisam(:,:,ind);
            validOmegaInds  = ~isinf(tOmega) & [1 tOmega(2:end) ~= 0] ...
                                & sum(omegaG_Gisam.^2,1) ~= 0;
            omegaG_Gisam  = omegaG_Gisam(:,validOmegaInds);
            omegaG_GgmAll = RB0toG0 * RWtoB0 * RBt0toW * this.r.omegaBt_Btgm; 
            omegaG_Ggm    = omegaG_GgmAll(:,this.r.omegaStateNums<=stateNum); 
            
            % Get axis limits
            xmin = min([omegaG_GgmAll(1,:) this.r.omegaBt_Gisam(1,:) 0]);
            xmax = max([omegaG_GgmAll(1,:) this.r.omegaBt_Gisam(1,:)]);
            ymin = min([omegaG_GgmAll(2,:) this.r.omegaBt_Gisam(2,:) 0]);
            ymax = max([omegaG_GgmAll(2,:) this.r.omegaBt_Gisam(2,:)]);
            zmin = min([omegaG_GgmAll(3,:) this.r.omegaBt_Gisam(3,:) 0]);
            zmax = max([omegaG_GgmAll(3,:) this.r.omegaBt_Gisam(3,:)]);
                       
            plot3(omegaG_Ggm(1,:),omegaG_Ggm(2,:),omegaG_Ggm(3,:),'Color',this.col(:,1)); hold on;
            plot3(omegaG_Gisam(1,:),omegaG_Gisam(2,:),omegaG_Gisam(3,:),'Color',this.col(:,2));
            plot3(omegaG_Ggm(1,1),omegaG_Ggm(2,1),omegaG_Ggm(3,1),'.g','MarkerSize',20);
            plot3(omegaG_Ggm(1,end),omegaG_Ggm(2,end),omegaG_Ggm(3,end),'.','MarkerFaceColor', ...
                this.col(:,2),'MarkerEdgeColor',this.col(:,1),'MarkerSize',20); 
            plot3(omegaG_Gisam(1,1),omegaG_Gisam(2,1),omegaG_Gisam(3,1),'.g','MarkerSize',20);
            plot3(omegaG_Gisam(1,end),omegaG_Gisam(2,end),omegaG_Gisam(3,end),'.','MarkerFaceColor', ...
                this.col(:,2),'MarkerEdgeColor',this.col(:,2),'MarkerSize',20);
            plot3(0,0,0,'.'); 
            axis equal; grid on; axis([xmin xmax ymin ymax zmin zmax]); 
            title('Polhode'); legend('Gyroscopes',type,'Location','NorthWest');
            xlabel('{}^G\omega_1 [rad/s]'); ylabel('{}^G\omega_2 [rad/s]'); zlabel('{}^G\omega_3 [rad/s]'); 
            view(60,30);
            
        end
        
        %% Plot the body polhode together with current inertial properties estimates
        function plotPolhodeCurrent( this, stateNum )
                        
            % Get the inertial properties optimization runs, and remove the really high cost ones
            ipoRuns         = this.r.ipo(~cellfun('isempty',this.r.ipo));
            i = 1;
            while i <= length(ipoRuns)
                if ipoRuns{i}.inertiaRatiosCost >= 1000
                    ipoRuns(i) = [];
                end
                i = i + 1;
            end
            
            % Get the latest ipoRun
            if stateNum < ipoRuns{1}.stateNum
                ipoRun = [];
                RBttoG = eye(3);
                Jopt   = [1 1]';
            else
                ipoRun      = ipoRuns{end};
                for i = 1:length(ipoRuns)
                    if ipoRuns{i}.stateNum >= stateNum
                        ipoRun  = ipoRuns{max(1,i-1)};
                        if ipoRun.isSingleAxis
                            ipoRun = [];
                            RBttoG = eye(3);
                            Jopt   = [1 1]';
                        else
                            RBttoG  = ipoRun.RBttoG;
                            Jopt = ipoRun.inertiaRatiosOpt.getInertiaRatios();
                        end
                        break;
                    end
                end
            end
            
            
            % Get truth up until now, applying in-plane correction for simulated AS1 case
            if this.r.isSimulated
                thetaCorr = Log(this.r.RBttoGtrue'*RBttoG);
                omegaB_BtrueCorr = Exp([-thetaCorr(1) 0 0]) * this.r.omegaBt_Bttrue;
            else
                omegaB_BtrueCorr = this.r.omegaBt_Bttrue;
            end
            omegaB_BtrueCorr = omegaB_BtrueCorr(:,this.r.omegaStateNums <= stateNum);
            plot3(omegaB_BtrueCorr(1,:),omegaB_BtrueCorr(2,:),omegaB_BtrueCorr(3,:),'-c', ...
                'LineWidth',2); hold on; grid on; axis equal;
            
            % Plot the aligned measurements results, with a dot at the end indicating the 
            % current location
            stateDiff = stateNum - this.r.omegaStateNums;
            stateDiff = stateDiff(stateDiff >= 0);
            [~,ind] = min(stateDiff);
            omegaBt_Btmeas = RBttoG' * this.r.omegaBt_Gisam(:, ...
                this.r.omegaStateNums <= stateNum, ind);
            omegaBt_Btmeas(:,all(~any(omegaBt_Btmeas),1)) = []; % remove zero entries
            plot3(omegaBt_Btmeas(1,:),omegaBt_Btmeas(2,:),omegaBt_Btmeas(3,:),'Color',[0.5 0 0.5], ...
                'LineWidth',2);
            plot3(omegaBt_Btmeas(1,1),omegaBt_Btmeas(2,1),omegaBt_Btmeas(3,1),'d','MarkerFaceColor', ...
                'g','MarkerEdgeColor',[0.5 0 0.5],'MarkerSize',10);
            plot3(omegaBt_Btmeas(1,end),omegaBt_Btmeas(2,end),omegaBt_Btmeas(3,end),'o','MarkerFaceColor', ...
                [0.5 0 0.5],'MarkerEdgeColor','r','MarkerSize',10);
            
            % Plot the principal axes optimization, if available
            if ~isempty(ipoRun) && ~ipoRun.isSingleAxis
                ipoRun.principalAxesOpt.plot();                
            end
            
            plot3(omegaB_BtrueCorr(1,1),omegaB_BtrueCorr(2,1),omegaB_BtrueCorr(3,1),'d', ...
                'MarkerFaceColor','g','MarkerEdgeColor','c','MarkerSize',10);
            plot3(omegaB_BtrueCorr(1,end),omegaB_BtrueCorr(2,end),omegaB_BtrueCorr(3,end),'oc', ...
                'MarkerFaceColor','c','MarkerEdgeColor','r','MarkerSize',10);
            
            % Plot the inertia ratios, true and estimated
            lims = get(gca,'zlim');
            if this.r.isSimulated   % avoid dwarfing angular velocities
                JtrueP = this.r.Jtrue*0.3;
                JoptP  = Jopt*0.3;
            else
                JtrueP = this.r.Jtrue;
                JoptP  = Jopt;
            end
            plot3([0 JtrueP(1)],[0 0],[-lims(2) -lims(2)],'-r');
            plot3([0 0],[0 JtrueP(2)],[-lims(2) -lims(2)],'-g');
            plot3(JtrueP(1),0,-lims(2),'or','Markersize',10);
            plot3(0,JtrueP(2),-lims(2),'og','Markersize',10);
            plot3([0 JoptP(1)],[0 0],[-lims(2) -lims(2)],'--r','LineWidth',2);
            plot3([0 0],[0 JoptP(2)],[-lims(2) -lims(2)],'--g','LineWidth',2);
            plot3(JoptP(1),0,-lims(2),'.r','Markersize',30);
            plot3(0,JoptP(2),-lims(2),'.g','Markersize',30);
            if ~isempty(ipoRun) && ~ipoRun.isSingleAxis
                xlim([lims(1) max([lims(2) JtrueP(1)+0.05 JoptP(1)+0.05])]);
                ylim([lims(1) max([lims(2) JtrueP(2)+0.05 JoptP(2)+0.05])]);
            end
            
            legend('True','Measured','Location','NorthWest');
            view(140,30);
            
        end
        
        %% Plot the motion image, the map, and the polhode
        function plotQuad(this, imageNum, suppressImage)
            
            % Get the state number
            stateNum = this.getStateNum(imageNum);
                        
            sp1 = subplot(2,2,[1 2]); 
            set(sp1,'Units','Pixels'); 
            set(sp1,'Position',[0,481,1280,480]);
            if nargin < 3 || ~suppressImage
                this.plotImage(imageNum);
            end
            sp2 = subplot(2,2,3); 
            set(sp2,'Units','Pixels'); 
            set(sp2,'Position',[0,60,640,380]);
            this.plotMap(stateNum, this.r.PWtoGisam(1:3,1:3,stateNum+1,stateNum+1), ...
                                        this.r.tGtoBt_Gisam(:,stateNum+1));
            sp3 = subplot(2,2,4);
            set(sp3,'Units','Pixels'); 
            set(sp3,'Position',[640,60,640,380]);
            this.plotPolhode('isam', stateNum);        
            set(gcf,'units','pixels','Position',[50 50 1280 960]);
            drawnow();
            
        end
        
        %% Plot the motion image, the map, and the inertial property estimation summary
        function plotQuadCurrent(this, imageNum, suppressImage)
            
            % Get the state number
            stateNum = this.getStateNum(imageNum);
                        
            sp1 = subplot(2,2,[1 2]); 
            set(sp1,'Units','Pixels'); 
            set(sp1,'Position',[0,481,1280,480]);
            if nargin < 3 || ~suppressImage
                this.plotImage(imageNum);
            end
            
            sp2 = subplot(2,2,3); 
            set(sp2,'Units','Pixels'); 
            set(sp2,'Position',[0,50,640,380]);
            this.plotMapCurrent(stateNum,1,1);
            view(-70,20);
            title('Map');
            
            sp3 = subplot(2,2,4);
            set(sp3,'Units','Pixels'); 
            set(sp3,'Position',[640,50,640,380]);
            this.plotPolhodeCurrent(stateNum);        
            title('Polhode'); xlabel('\omega_1 [rad/s]'); ylabel('\omega_2 [rad/s]');
            zlabel('\omega_3 [rad/s]');
            set(gcf,'units','pixels','Position',[50 50 1280 960]);
            drawnow();
            
        end
        
        %% Plot the errors in estimation
        function plotErrors(this, imageNum)
            
            addpath('D:\PhD\Thesis\Matlab\breakxaxis');
            
            % Get the state number
            stateNum = this.getStateNum(imageNum);
            
            % Plot error in the inspector position
            subplot(2,2,1);
            if this.r.isSimulated
                this.plotPosError('PWtoBtrue','PWtoBisam',stateNum,this.r.tInspState,'-b');
            else
                this.plotPosError('PWtoBgm','PWtoBisam',stateNum,this.r.tInspState,'-b');
            end
            title('Inspector Position Error');
            
            % Plot error in the inspector orientation
            subplot(2,2,3);
            disp('Rot RBtoW');
            this.plotRotError('PWtoBtrue','PWtoBisam',stateNum,this.r.tInspState,'-b');
            title('Inspector Orientation Error');
            xlabel('Time [s]');
            
            % Plot error in inspector velocity
%             subplot(3,2,5);
%             this.plotVelError('true','isam',stateNum,this.r.tInspState,'-b');
%             title('Inspector Velocity Error');
%             xlabel('Time [s]');
            
            % Angle axis correction in ambiguous plane (simulated), or for flipped axes (TA)
            if strcmp(this.r.ipo{stateNum+1}.principalAxesOpt.inertiaSymmetry,'AS1')
                thetaCorr = Log(this.r.RBttoGtrue'*this.r.ipo{stateNum+1}.RBttoG);
                thetaCorr = [thetaCorr(1) 0 0]';
            else
                thetaCorr = [0 0 0]';
            end   
            disp(['A correction of thetaCorr = ' num2str(thetaCorr') ' will be applied']);
            
            % Plot the error in inertia ratios, principal axis alignment, and geometric center
            % location error
%             subplot(3,2,2);
            Jopt = this.r.ipo{stateNum+1}.inertiaRatiosOpt.getInertiaRatios();
%             if this.r.isSimulated
                errRBttoG = rad2deg(abs(Log(this.r.RBttoGtrue'*this.r.RBttoGisam(:,:,stateNum+1))));
                errtGtoBt = abs(this.r.tGtoBt_Gisam(:,stateNum+1) - this.r.tGtoBt_Gtrue);
% %                 [ax,~,~] = plotyy([1 2],abs(Jopt(1:2) - this.r.Jtrue(1:2)),[3 4 5],errRBttoG,'stem'); grid on;
%                 set(ax(1),'XLim',[0 6],'XTick',1:5,'XTickLabel',{'J_1 Error','J_2 Error', ...
%                         '\delta\theta_1','\delta\theta_2','\delta\theta_3'}); 
%                 set(ax(2),'XLim',[0 6],'Xtick',3:5);
%                 ylabel(ax(1),'Error [J - J_{true}]'); ylabel(ax(2),'Error [deg]');
%                 title('Target Inertia Ratios and Principal Axis Alignment Error');
%             else
%                 errRBttoG = rad2deg(abs(Log(this.r.RBttoGtrue'*(this.r.RBttoGisam(:,:,stateNum+1)*Exp(-thetaCorr)))));
%                 [ax,~,~] = plotyy([1 2],abs(Jopt(1:2) - this.r.Jtrue(1:2)),[3 4 5],errRBttoG,'stem'); grid on;
%                 set(ax(1),'XLim',[0 6],'XTick',1:5,'XTickLabel',{'J_1 Error','J_2 Error', ...
%                         '\delta\theta_1','\delta\theta_2','\delta\theta_3'}); 
%                 set(ax(2),'XLim',[0 6],'Xtick',3:5);
%                 ylabel(ax(1),'Error [J - J_{true}]'); ylabel(ax(2),'Error [deg]');
%                 title('Target Inertia Ratio and Principal Axis Alignment Error');
                
%                 stem(abs(Jopt(1:2) - this.r.Jtrue(1:2))); grid on;
%                 set(gca,'XLim',[0 3],'XTick',1:2,'XTickLabel',{'J_1 Error','J_2 Error'}); 
%                 ylabel('Error [J - J_{true}]');
%                 title('Target Inertia Ratios Error');                
%             end
            
            % Create inertial properties error table
%             if this.r.isSimulated
                disp('Inertial properties error table =====================================');
                disp('{\bf Variable} & {\bf True value} & {\bf Estimated value} & {\bf Error magnitude} \\');
                disp('\hline');
                disp(['    Center of mass ${}^{\underline{G}}\mathbf{t}_{\underline{B}\text{to}\underline{G}}$ (x dir) & ' ...
                    num2str(this.r.tGtoBt_Gtrue(1)) '~m & ' num2str(this.r.tGtoBt_Gisam(1,stateNum+1)) ...
                    '~m & ' num2str(errtGtoBt(1)) '~m \\']);
                disp(['    Center of mass ${}^{\underline{G}}\mathbf{t}_{\underline{B}\text{to}\underline{G}}$ (y dir) & ' ...
                    num2str(this.r.tGtoBt_Gtrue(2)) '~m & ' num2str(this.r.tGtoBt_Gisam(2,stateNum+1)) ...
                    '~m & ' num2str(errtGtoBt(2)) '~m \\']);
                disp(['    Center of mass ${}^{\underline{G}}\mathbf{t}_{\underline{B}\text{to}\underline{G}}$ (z dir) & ' ...
                    num2str(this.r.tGtoBt_Gtrue(3)) '~m & ' num2str(this.r.tGtoBt_Gisam(3,stateNum+1)) ...
                    '~m & ' num2str(errtGtoBt(3)) '~m \\']);
                disp(['    Center of mass ${}^{\underline{G}}\mathbf{t}_{\underline{B}\text{to}\underline{G}}$ (total) & ' ...
                    ' &  & ' num2str(norm(errtGtoBt)) ' \\']);
                disp('\hline');
                disp(['    Principal axes $\delta\theta_1$ &  &  & ' num2str(errRBttoG(1)) '$^\circ$ \\']);
                disp(['    Principal axes $\delta\theta_2$ &  &  & ' num2str(errRBttoG(2)) '$^\circ$ \\']);
                disp(['    Principal axes $\delta\theta_3$ &  &  & ' num2str(errRBttoG(3)) '$^\circ$ \\']);
                disp('\hline');
%                 disp(['    Principal axes $\lVert \boldsymbol{\delta\theta} \rVert$ &  &  & ' num2str(norm(errRBttoG)) ' \\']);
                disp(['    Inertia ratio $J_1$ & ' num2str(this.r.Jtrue(1)) ' & ' ...
                    num2str(Jopt(1)) ' & ' num2str(abs(Jopt(1)-this.r.Jtrue(1))) ' \\']);
                disp(['    Inertia ratio $J_2$ & ' num2str(this.r.Jtrue(2)) ' & ' ...
                    num2str(Jopt(2)) ' & ' num2str(abs(Jopt(2)-this.r.Jtrue(2))) ' \\']);
                disp('\hline');
%             end
            
       
            % Extended time interval and rigid body rotation for prediction after end of data
            nExt = 40;
            tExt = this.r.tInspState(end):mean(gradient(this.r.tInspState)):...
                        nExt*this.r.tInspState(end);
            rbrTrue = RigidBodyRotation(this.r.Jtrue(1:2),this.r.PWtoBttrue(1:3,1:3,1),...
                                            this.r.omegaBt_Bttrue(:,1),'omega0',0);
            
            % Plot the orientation error of the target
            subplot(2,2,2);
            PWtoBtalign = zeros(4,4,this.r.nStates);     % opt aligned meas
            PWtoBtmodel = zeros(4,4,this.r.nStates);     % opt motion model
            RBttoWmodel = this.r.ipo{stateNum+1}.rigidBodyRotationOpt.predictOrientation( ...
                                                                        this.r.tInspState);
            PBtoBtisam  = zeros(4,4,this.r.nStates);    % isam insp body to tgt body
            PBtoBtgm    = zeros(4,4,this.r.nStates);    % global met insp body to tgt body
            dThetaBtoBt = zeros(3,this.r.nStates);

            for i = 1:this.r.nStates
                PWtoBtalign(:,:,i) = this.r.PWtoGisam(:,:,i,stateNum+1) * ...           
                            [this.r.RBttoGisam(:,:,stateNum+1) zeros(3,1); 0 0 0 1];                    
                PWtoBtmodel(:,:,i) = [RBttoWmodel(:,:,i) zeros(3,1); 0 0 0 1];

                % Correct for alignment in ambiguous plane (AS1)
                PWtoBtalign(1:3,1:3,i) = PWtoBtalign(1:3,1:3,i) * Exp(-thetaCorr);
                PWtoBtmodel(1:3,1:3,i) = PWtoBtmodel(1:3,1:3,i) * Exp(-thetaCorr);
                
                % Get insp body to tgt body pose
%                 PBtoBtisam(:,:,i)   = [ this.r.PWtoBisam(1:3,1:3,i,stateNum+1)' * ...
%                     PWtoBtalign(1:3,1:3,i) zeros(3,1); zeros(1,3) 1];
                PBtoBtisam(:,:,i)   = [ this.r.PGtoBisam(1:3,1:3,i,stateNum+1)' * ...
                    this.r.RBttoGisam(:,:,stateNum+1) zeros(3,1); zeros(1,3) 1];
                PBtoBtgm(:,:,i)     = [this.r.PWtoBgm(1:3,1:3,i)' * this.r.PWtoBttrue(1:3,1:3,i) ...
                    zeros(3,1); zeros(1,3) 1];
                dThetaBtoBt(:,i)    = Log( PBtoBtgm(1:3,1:3,i)' * PBtoBtisam(1:3,1:3,i) );

            end            
            disp('Rot RBttoWalign');
            this.plotRotError(this.r.PWtoBttrue,PWtoBtalign,stateNum,this.r.tInspState,'-b'); 
            hold on;
%             this.plotRotError(PBtoBtgm,PBtoBtisam,stateNum,this.r.tInspState,'-k');
%             plot(this.r.tInspState,abs(rad2deg(dThetaBtoBt(1,:))),'r'); hold on;
%             plot(this.r.tInspState,abs(rad2deg(dThetaBtoBt(2,:))),'g');
%             plot(this.r.tInspState,abs(rad2deg(dThetaBtoBt(3,:))),'b');
            title('Target Orientation Error');

            % Model further in time
            if this.r.isSimulated
            	disp('Rot RBttoWmodel');
                this.plotRotError(this.r.PWtoBttrue,PWtoBtmodel,stateNum,this.r.tInspState,'--m');
                legend('Aligned Estimate (In-Plane Corrected)','Modeled (In-Plane Corrected)');
                RBttoWmodelExt = this.r.ipo{stateNum+1}.rigidBodyRotationOpt.predictOrientation(tExt);
                RBttoWtrueExt  = rbrTrue.predictOrientation(tExt);
                PWtoBtmodelExt = zeros(4,4,length(tExt));
                PWtoBttrueExt = zeros(4,4,length(tExt));
                for i = 1:length(tExt)
                    PWtoBtmodelExt(:,:,i) = [RBttoWmodelExt(:,:,i) zeros(3,1); 0 0 0 1];
                    PWtoBttrueExt(:,:,i)  = [RBttoWtrueExt(:,:,i)  zeros(3,1); 0 0 0 1];

                    % Correct for alignment in ambiguous plane (AS1), or flipped axes (TA)  
                    PWtoBtmodelExt(1:3,1:3,i) =  PWtoBtmodelExt(1:3,1:3,i) * Exp(-thetaCorr); 
                end
                
            	disp('Rot RBttoWmodelExt');
                this.plotRotError(PWtoBttrueExt,PWtoBtmodelExt,stateNum,tExt,'--m');
                xlim([0 tExt(end)]);   
                
                breakxaxis([2*this.r.ipo{stateNum+1}.tOmega(end) ...
                    tExt(end)-this.r.ipo{stateNum+1}.tOmega(end)]);  
            end

            % Plot the angular velocity error of the target
            subplot(2,2,4);     
            omegaBt_Btalign = this.r.ipo{stateNum+1}.omegaB_B; 
            omegaBt_Btmodel = this.r.ipo{stateNum+1}.rigidBodyRotationOpt.predictOmega( ...
                                this.r.ipo{stateNum+1}.tOmega);

            % Correct for alignment issues from unobservability or flipped axes
            omegaBt_Btalign =  Exp(thetaCorr) * omegaBt_Btalign;
            omegaBt_Btmodel =  Exp(thetaCorr) * omegaBt_Btmodel; 
            
            disp('Angvel omegaB align');
            this.plotAngVelError(this.r.ipo{stateNum+1}.omegaB_Btrue,omegaBt_Btalign, ...
                                stateNum,this.r.ipo{stateNum+1}.tOmega,'-b');
            hold on;
            xlabel('Time [s]'); title('Target Angular Velocity Error'); 
            
            disp('Angvel omegaB model');
            this.plotAngVelError(this.r.ipo{stateNum+1}.omegaB_Btrue,omegaBt_Btmodel,...
                                stateNum,this.r.ipo{stateNum+1}.tOmega,'--m');
            errAngVelMod = abs(rad2deg(omegaBt_Btmodel - this.r.ipo{stateNum+1}.omegaB_Btrue));
%             plot(this.r.ipo{stateNum+1}.tOmega,errAngVelMod(1,:),'--r');
%             plot(this.r.ipo{stateNum+1}.tOmega,errAngVelMod(2,:),'--g');
%             plot(this.r.ipo{stateNum+1}.tOmega,errAngVelMod(3,:),'--b');
                                
            % Model further in time
            if this.r.isSimulated
                legend('Aligned Estimate (In-Plane Corrected)','Modeled (In-Plane Corrected)');
                omegaBt_BtmodelExt = this.r.ipo{stateNum+1}.rigidBodyRotationOpt.predictOmega(tExt);
                omegaBt_BttrueExt  = rbrTrue.predictOmega(tExt);

                % Correct for alignment issues from unobservability or flipped axes
                omegaBt_BtmodelExt =  Exp(thetaCorr) * omegaBt_BtmodelExt;              
            
                disp('Angvel omegaB modelExt');
                this.plotAngVelError(omegaBt_BttrueExt,omegaBt_BtmodelExt,stateNum,tExt,'--m');
                xlim([0 tExt(end)]);
                
                breakxaxis([2*this.r.ipo{stateNum+1}.tOmega(end) ...
                tExt(end)-this.r.ipo{stateNum+1}.tOmega(end)]);
            else                
                legend('Aligned Estimate','Modeled');
            end

            qBttoWalign = zeros(4,this.r.nStates);
            qBttoWmodel = zeros(4,this.r.nStates);
            qBttoWtrue = zeros(4,this.r.nStates);
            qBttoBisam = zeros(4,this.r.nStates);
            qBttoBgm   = zeros(4,this.r.nStates);
            for i = 1:this.r.nStates
                qBttoWalign(:,i) = rot2quat(PWtoBtalign(1:3,1:3,i));
                qBttoWmodel(:,i) = rot2quat(PWtoBtmodel(1:3,1:3,i));
                qBttoWtrue(:,i) = rot2quat(this.r.PWtoBttrue(1:3,1:3,i));
                qBttoBisam(:,i) = rot2quat(PBtoBtisam(1:3,1:3,i));
                qBttoBgm(:,i) = rot2quat(PBtoBtgm(1:3,1:3,i));
            end
%             qBttoWalign = fix_quats(qBttoWalign);
%             qBttoWmodel = fix_quats(qBttoWmodel);
%             qBttoWtrue = fix_quats(qBttoWtrue);
            figure;
            subplot(4,1,1);
            plot(this.r.tInspState,qBttoWtrue(1,:),'-k'); hold on; grid on;
            plot(this.r.tInspState,qBttoWalign(1,:),'-b');
%             plot(this.r.tInspState,qBttoWmodel(1,:),'-m');
            plot(this.r.tInspState(2:end),qBttoBisam(1,2:end),'-g');
            plot(this.r.tInspState(1:end-1),qBttoBgm(1,1:end-1),'-k');
            subplot(4,1,2);
            plot(this.r.tInspState,qBttoWtrue(2,:),'-k'); hold on; grid on;
            plot(this.r.tInspState,qBttoWalign(2,:),'-b');
%             plot(this.r.tInspState,qBttoWmodel(2,:),'-m');
            plot(this.r.tInspState(2:end),qBttoBisam(2,2:end),'-g');
            plot(this.r.tInspState(1:end-1),qBttoBgm(2,1:end-1),'-k');
            subplot(4,1,3);
            plot(this.r.tInspState,qBttoWtrue(3,:),'-k'); hold on; grid on;
            plot(this.r.tInspState,qBttoWalign(3,:),'-b');
%             plot(this.r.tInspState,qBttoWmodel(3,:),'-m');
            plot(this.r.tInspState(2:end),qBttoBisam(3,2:end),'-g');
            plot(this.r.tInspState(1:end-1),qBttoBgm(3,1:end-1),'-k');
            subplot(4,1,4);
            plot(this.r.tInspState,qBttoWtrue(4,:),'-k'); hold on; grid on;
            plot(this.r.tInspState,qBttoWalign(4,:),'-b');
%             plot(this.r.tInspState,qBttoWmodel(4,:),'-m');
            plot(this.r.tInspState(2:end),qBttoBisam(4,2:end),'-g');
            plot(this.r.tInspState(1:end-1),qBttoBgm(4,1:end-1),'-k');
            
%             figure;
%             subplot(3,1,1);
%             plot(tExt,omegaBt_BttrueExt(1,:)); hold on; grid on;
%             plot(tExt,omegaBt_BtmodelExt(1,:),'-r');
%             subplot(3,1,2);
%             plot(tExt,omegaBt_BttrueExt(2,:)); hold on; grid on;
%             plot(tExt,omegaBt_BtmodelExt(2,:),'-r');
%             subplot(3,1,3);
%             plot(tExt,omegaBt_BttrueExt(3,:)); hold on; grid on;
%             plot(tExt,omegaBt_BtmodelExt(3,:),'-r');
            
%             linkaxes([ax1 ax2],'x');
                      
            
        end
        
        %% Plot a map that includes all the info, plotImages is true by default
        function plotGlory( this, stateNum, plotImages )
                
            if nargin == 2
                plotImages = 1;
            end
                        
            % Determine subplot postiions
            if plotImages
                subDim = [3 3];
                poseSub = [4 5 7 8];
                mapSub = 6;
                polhodeSub = 9;
            else
                subDim = [2 3];
                poseSub = [1 2 4 5];
                mapSub = 3;
                polhodeSub = 6;
            end
                
            subplot = @(m,n,p) subtightplot (m, n, p, [0.075 0.01], [0.05 0.01], [0.01 0.01]);
%                 figure(1); clf;                
            if plotImages
                set(gcf,'units','pixels','Position',[50 50 1024 920]);
                subplot(subDim(1),subDim(2),[1 2 3]);
                this.plotImage(stateNum);
            else
                set(gcf,'units','pixels','Position',[50 50 1400 920]);
            end
            subplot(subDim(1),subDim(2),poseSub);
            this.plotPoseChains({'PWtoBgm','PWtoBisam','PWtoGisam'}, stateNum, ...
                                    [-0.75 0.75 -0.75 0.75 -0.75 0.75]);
            subplot(subDim(1),subDim(2),mapSub);
            this.plotMap(stateNum, this.r.PWtoGisam(1:3,1:3,stateNum+1), this.r.tGtoBt_Gisam);
            subplot(subDim(1),subDim(2),polhodeSub);
            this.plotPolhode('isam', stateNum);
            set(findall(gcf,'type','text'),'FontName','CMU Serif');                
            set(findall(gcf,'type','legend'),'FontName','CMU Serif');
            drawnow;
            
            
        end
        
        %% Plot results summary
        function plotResultsSummary(this, stateNum)
            
            ind = find(this.r.stateNums==stateNum);
        
            % Plot large pose chains
            subplot(5,2,[1 3]);
            if this.r.isSimulated
                axLims = [-0.75 0.25 -0.75 0.75 -0.25 0.75];
                this.plotPoseChains({'PWtoBtrue','PWtoBisam'}, stateNum, axLims);
            else
                axLims = [-0.25 0.75 -0.75 0.75 -1 0.25];
                this.plotPoseChains({'PWtoBgm','PWtoBisam'}, stateNum, axLims);
            end
            title('Inspector Body Frame Pose');
            subplot(6,2,[2 4]);
            this.plotPoseChains({'PWtoGisam'}, stateNum, ...
                [this.r.tWtoBt_Wtrue(1)-0.2 this.r.tWtoBt_Wtrue(1)+0.2 this.r.tWtoBt_Wtrue(2)-0.2 ...
                this.r.tWtoBt_Wtrue(2)+0.2 this.r.tWtoBt_Wtrue(3)-0.2 this.r.tWtoBt_Wtrue(3)+0.2]);
            title('Target Geometric Frame Pose');
            
            % Plot inpsector position
            t = this.r.tInspState(1:ind);
            tWtoB_Wisam = reshape(this.r.PWtoBisam(1:3,4,:,ind),3,size(this.r.PWtoBisam,3));
            tWtoB_Wisam(:,all(~any(tWtoB_Wisam),1)) = [];
            if this.r.isSimulated
                tWtoB_Wtrue = reshape(this.r.PWtoBtrue(1:3,4,1:ind),3,ind);
            else
                tWtoB_Wtrue = reshape(this.r.PWtoBgm(1:3,4,1:ind),3,ind);
            end
            subplot(5,2,[5 6]);
            plot(t,tWtoB_Wtrue(1,:),'-r'); hold on; grid on;
            plot(t,tWtoB_Wtrue(2,:),'-g');
            plot(t,tWtoB_Wtrue(3,:),'-b');
            plot(t,tWtoB_Wisam(1,:),'--r','LineWidth',2);
            plot(t,tWtoB_Wisam(2,:),'--g','LineWidth',2); 
            plot(t,tWtoB_Wisam(3,:),'--b','LineWidth',2);
            if this.r.isSimulated
                legend('True (x dir.)','True (y dir.)','True (z dir.)','Estimated (x dir.)', ...
                    'Estimated (y dir.)','Estimated (z dir.)');
            else
                legend('Global Met. (x dir.)','Global Met. (y dir.)','Global Met. (z dir.)', ...
                    'Estimated (x dir.)', 'Estimated (y dir.)','Estimated (z dir.)');
            end
            ylabel('Position [m]'); title('Inspector Position');
            xlim([0 t(end)]);
            
            % Plot inspector quaternion
            qWtoBisam = zeros(4,ind);
            qWtoBtrue = zeros(4,ind);
            for i = 1:ind
                qWtoBisam(:,i) = rot2quat(this.r.PWtoBisam(1:3,1:3,i,ind)');
                qWtoBtrue(:,i) = rot2quat(this.r.PWtoBtrue(1:3,1:3,i)');
            end
            subplot(5,2,[7 8]);
            plot(t,qWtoBtrue(1,1:ind),'-r'); hold on; grid on;
            plot(t,qWtoBtrue(2,1:ind),'-g');
            plot(t,qWtoBtrue(3,1:ind),'-b');
            plot(t,qWtoBtrue(4,1:ind),'-k');
            plot(t,qWtoBisam(1,1:ind),'--r','LineWidth',2);
            plot(t,qWtoBisam(2,1:ind),'--g','LineWidth',2);
            plot(t,qWtoBisam(3,1:ind),'--b','LineWidth',2);
            plot(t,qWtoBisam(4,1:ind),'--k','LineWidth',2);
            legend('True q_1','True q_2','True q_3','True q_4','Estimated q_1','Estimated q_2', ...
                'Estimated q_3','Estimated q_4');
            ylabel('Quaternion [-]'); title('Inspector Orientation');
            xlim([0 t(end)]);
          
            % Plot target quaternion
            subplot(5,2,[9 10]);
            qWtoBtisam = zeros(4,ind);
            qWtoBttrue = zeros(4,ind);
            for i = 1:ind
                qWtoBttrue(:,i) = rot2quat(this.r.PWtoBttrue(1:3,1:3,i)');
                RBttoWisam = this.r.PWtoGisam(1:3,1:3,i,ind) * this.r.RBttoGisam(:,:,ind);
                qWtoBtisam(:,i) = rot2quat(RBttoWisam');                
            end
            qWtoBtisam = quatfix(qWtoBtisam);
            qWtoBttrue = quatfix(qWtoBttrue);
            % Set to conjugate equivalent negative quaternion if that matches better
            for i = 1:ind
                if norm(qWtoBtisam(:,i)+qWtoBttrue(:,i)) < norm(qWtoBtisam(:,i)-qWtoBttrue(:,i))
                    qWtoBtisam(:,i) = -qWtoBtisam(:,i);
                end
            end
            plot(t,qWtoBttrue(1,1:ind),'-r'); hold on; grid on;
            plot(t,qWtoBttrue(2,1:ind),'-g');
            plot(t,qWtoBttrue(3,1:ind),'-b');
            plot(t,qWtoBttrue(4,1:ind),'-k'); hold on;
            plot(t,qWtoBtisam(1,1:ind),'--r','LineWidth',2);
            plot(t,qWtoBtisam(2,1:ind),'--g','LineWidth',2);
            plot(t,qWtoBtisam(3,1:ind),'--b','LineWidth',2);
            plot(t,qWtoBtisam(4,1:ind),'--k','LineWidth',2);
            ylabel('Quaternion [-]'); title('Target Orientation');
            xlim([0 t(end)]);
            xlabel('Time [s]');
            legend('True q_1','True q_2','True q_3','True q_4','Estimated q_1','Estimated q_2', ...
                'Estimated q_3','Estimated q_4');            
            
        end
        
        %% Make video
        function makeVideo(this)
            
            % Delete old images
            delete('./plotglory/pg*.png');
            
            % Read default style and set to use png format
            sty = hgexport('readstyle','default');
            sty.Format = 'png';
            
            % Make the video
            for i = 2:length(this.imageRange)
                this.plotGlory(this.imageRange(i));
                hgexport(gcf,['./plotglory/pg' num2str(i-1) '.png'],sty);
            end
        end
        
        %% Make quad video
        function makeQuadVideo(this)
            
            % Delete old images
            delete('D:/PhD/Thesis/Matlab/v0.9/plotquad/pq*.png');
            
            % Read default style and set to use png format
            sty = hgexport('readstyle','default');
            sty.Format = 'png';
            
            % Make the video
            for i = this.r.imageNums(2):this.r.imageNums(end)
                disp(['Plotting Quad Chart ' num2str(i) ' of ' num2str(length(this.r.imageNums))]);
                this.plotQuad(i);
                hgexport(gcf,['D:/PhD/Thesis/Matlab/v0.9/plotquad/pq' num2str(i-1) '.png'],sty);
            end
            
        end
        
        %% Make quad video with current isam estimates
        function makeQuadCurrentVideo(this)
            
            % Delete old images
            delete('D:/PhD/Thesis/Matlab/v0.9/plotquad/pq*.png');
            
            % Read default style and set to use png format
            sty = hgexport('readstyle','default');
            sty.Format = 'png';
            sty.Renderer = 'opengl';
            
            % Make the video
            for i = this.r.imageNums(2):this.r.imageNums(end)
                disp(['Plotting Quad Chart ' num2str(i) ' of ' num2str(length(this.r.imageNums))]);
                this.plotQuadCurrent(i);
                hgexport(gcf,['D:/PhD/Thesis/Matlab/v0.9/plotquad/pq' num2str(i-1) '.png'],sty);
            end
            
        end
        
        %% Make a video of the final map, rotating
        function makeMapVideo(this)            
            
            % Delete old images
            delete('D:/PhD/Thesis/Matlab/v0.9/mapvid/map*.png');
            
            % Create plot
            ptScale = 1;
            sB = 0.2;
            if this.r.isSimulated
                RBtoBhat = this.r.RBttoGisam(:,:,end)'*this.r.RBttoGtrue;
                xB_Bhat = RBtoBhat(:,1); yB_Bhat = RBtoBhat(:,2); zB_Bhat = RBtoBhat(:,3);
                plot3([0 sB*xB_Bhat(1)],[0 sB*xB_Bhat(2)],[0 sB*xB_Bhat(3)],'-r','LineWidth',max(1,4*ptScale)); hold on; grid on;
                plot3([0 sB*yB_Bhat(1)],[0 sB*yB_Bhat(2)],[0 sB*yB_Bhat(3)],'-g','LineWidth',max(1,4*ptScale));
                plot3([0 sB*zB_Bhat(1)],[0 sB*zB_Bhat(2)],[0 sB*zB_Bhat(3)],'-b','LineWidth',max(1,4*ptScale));
            end
            plot3([0 sB],[0 0],[0 0],'--r','LineWidth',max(1,4*ptScale)); hold on; grid on;
            plot3([0 0],[0 sB],[0 0],'--g','LineWidth',max(1,4*ptScale));
            plot3([0 0],[0 0],[0 sB],'--b','LineWidth',max(1,4*ptScale));
            this.plotMap(this.r.stateNums(end), this.r.RBttoGisam(:,:,end)', ...
                                        this.r.tGtoBt_Gisam(:,end),ptScale);            
            if this.r.isSimulated
                legend('True x_B','True y_B','True x_B','Estimated x_B','Estimated y_B', ...
                    'Estimated z_B','All Features','Inliers','First Frame Inliers', ...
                    'Estimated t_{GtoB}','Estimated x_G', 'Estimated y_G','Estimated z_G', ...
                    'Location','East');
            else
                legend('Estimated x_B','Estimated y_B', ...
                    'Estimated z_B','All Features','Inliers','First Frame Inliers', ...
                    'Estimated t_{GtoB}','Estimated x_G',...
                    'Estimated y_G','Estimated z_G','Location','East');
            end
            
            % Initial view
            view(-30,30);
            
            % Set position of figure
            set(gcf,'Units','Pixels'); 
            set(gcf,'Position',[50,50,1024,768]);
            
            % Read default style and set to use png format
            sty = hgexport('readstyle','default');
            sty.Format = 'png';
            sty.Renderer = 'opengl';
            
            % Rotate around and export images
            az1 = [-30:2:270];                  
            el1 = repmat(30,1,length(az1));
            el2 = [30:-2:0 zeros(1,25)];
            az2 = repmat(270,1,length(el2));
            az3 = [270:2:360 repmat(360,1,25)];
            el3 = zeros(1,length(az3));
            el4 = [0:2:90 repmat(90,1,25)];
            az4 = repmat(360,1,length(el4));
            az = [az1 az2 az3 az4];
            el = [el1 el2 el3 el4];
            for i = 1:length(az)
                view(az(i),el(i));
                hgexport(gcf,['D:/PhD/Thesis/Matlab/v0.9/mapvid/map' num2str(i) '.png'],sty);
            end
                       
            
        end
        
        %% Plot the execution times
        function plotExecTimes(this)
            
            imageNums = 0:this.r.imageNums(end);
            
            subplot(3,1,1);
            plot(imageNums,this.r.equalizeTime,'-r'); hold on; grid on;
            plot(imageNums,this.r.blobTrackTime,'-g');
            plot(imageNums,this.r.stereoFrameTime,'-c');
            plot(imageNums,this.r.visOdomTime,'-m');
            plot(imageNums,this.r.loopClosureTime,'-b');
            plot(imageNums,this.r.isamTime,'-k');
            legend('Equalize','Range and Bearing Meas','Stereo Frame Acquisition', ...
                'Visual Odometry','Loop Closure Attempts','iSAM2','Location','NorthWest');
            ylabel('C++ Time [s]');
            title('Per-Image Computation Times');
            xlim([0 imageNums(end)]);
            
            % Do a linear regression on iSAM times, without outliers caused by loop closure
            % https://www.mathworks.com/help/matlab/data_analysis/linear-regression.html#f1-15010
            isamTimeImageNums = (this.r.imageNums(1):this.r.imageNums(end))';
            isamTimeImageNums = isamTimeImageNums(this.r.isamTime < 0.01);
            isamTime = this.r.isamTime(this.r.isamTime < 0.01);
            A = [ones(length(isamTimeImageNums),1) isamTimeImageNums];
            b = A\isamTime;    % intercept and slope
            isamTimeLin = b(1) + b(2)*isamTimeImageNums;
            residuals = isamTime - isamTimeLin;
            rsq = 1 - sum(residuals.^2)/( (length(isamTime)-1)*var(isamTime) );
            disp(['Isam time growth is ' num2str(b(2)) ' s / image, with R^2 of ' num2str(rsq)]);
            
            subplot(3,1,2);
            plot(imageNums,this.r.isamTime,'-k'); grid on; hold on;
            plot(isamTimeImageNums,isamTimeLin,'--r');
            ylabel('C++ Time [s]'); legend('iSAM2','iSAM2 Linear Fit','Location','NorthWest');
            xlim([0 imageNums(end)]); ylim([0 0.06]);
            
            subplot(3,1,3);
            ipoRuns = this.r.ipo(~cellfun('isempty',this.r.ipo));
            iroTimes = zeros(1,length(ipoRuns));
            paoTimes = zeros(1,length(ipoRuns));
            ipoImageNums = zeros(1,length(ipoRuns));
            for i = 1:length(ipoRuns)
                paoTimes(i) = ipoRuns{i}.paoExecTime;
                iroTimes(i) = ipoRuns{i}.iroExecTime;
                ipoImageNums(i) = this.r.imageNums(this.r.stateNums==ipoRuns{i}.stateNum);
            end
            plot(ipoImageNums,paoTimes,'-b'); grid on; hold on;
            plot(ipoImageNums,iroTimes,'-g');
            ylabel('MATLAB Time [s]');
            legend('Principal Axes Opt','Inertia Ratios Opt','Location','NorthWest');
            xlim([0 imageNums(end)]);            
            xlabel('Image Number [-]'); 
            
            totTime = sum([this.r.equalizeTime' this.r.blobTrackTime' this.r.stereoFrameTime' ...
                this.r.visOdomTime' this.r.loopClosureTime' this.r.isamTime' ...
                paoTimes(end)+iroTimes(end)]);
                
            disp(['Number of images              & ' num2str(length(imageNums)+1) ' \\']);
            disp('\hline');
            disp(['Equalize                      & ' num2str(sum(this.r.equalizeTime)) ...
                ' s & ' num2str(sum(this.r.equalizeTime)/totTime*100) '\% \\']);
            disp(['Range and bearing measurement & ' num2str(sum(this.r.blobTrackTime)) ...
                ' s & ' num2str(sum(this.r.blobTrackTime)/totTime*100) '\% \\']);
            disp(['Stereo frame acquisition      & ' num2str(sum(this.r.stereoFrameTime)) ...
                ' s & ' num2str(sum(this.r.stereoFrameTime)/totTime*100) '\% \\']);
            disp(['Visual odometry               & ' num2str(sum(this.r.visOdomTime)) ...
                ' s & ' num2str(sum(this.r.visOdomTime)/totTime*100) '\% \\']);
            disp(['Loop closure                  & ' num2str(sum(this.r.loopClosureTime)) ...
                ' s & ' num2str(sum(this.r.loopClosureTime)/totTime*100) '\% \\']);
            disp(['iSAM2                         & ' num2str(sum(this.r.isamTime)) ...
                ' s & ' num2str(sum(this.r.isamTime)/totTime*100) '\% \\']);
            disp(['Final inertial properties     & ' num2str(paoTimes(end)+iroTimes(end)) ...
                ' s & ' num2str((paoTimes(end)+iroTimes(end))/totTime*100) '\% \\']);
            disp('\hline');
            disp(['Total                         & ' num2str(totTime) ' s  & 100\% \\']);
            
            
        end
        
        %% Plot the convergence of the target estimation parameters
        function plotConvergence(this)
            
            % Get the inertial properties estimation runs
            ipoRuns = this.r.ipo(~cellfun('isempty',this.r.ipo));
            ipoImageNums = zeros(1,length(ipoRuns));
            ipoStateNums = zeros(1,length(ipoRuns));
            ipoSingleMulti = zeros(1,length(ipoRuns));  % 0 = single, 2 = multi
            ipoEnergyState = zeros(1,length(ipoRuns));  % 0 = n/a, 1 = low, 2 = high
            ipoInertiaSymm = zeros(1,length(ipoRuns));  % 0 = n/a, 1 = tri-axial, 2 = axis-symm
            for i = 1:length(ipoRuns)                
                ipoImageNums(i) = this.r.imageNums(this.r.stateNums==ipoRuns{i}.stateNum);
                ipoStateNums(i) = ipoRuns{i}.stateNum;
                if ipoRuns{i}.isSingleAxis
                    ipoSingleMulti(i) = 0;
                else
                    ipoSingleMulti(i) = 2;
                    if strcmp(ipoRuns{i}.inertiaRatiosOpt.energyState,'LE')
                        ipoEnergyState(i) = 1;
                    else
                        ipoEnergyState(i) = 2;
                    end
                    if strcmp(ipoRuns{i}.inertiaRatiosOpt.inertiaSymmetry,'TA');
                        ipoInertiaSymm(i) = 1;
                    else
                        ipoInertiaSymm(i) = 2;
                    end
                end
            end
            
            % Plot convergence of the center of mass
            errTransGeom = this.r.tGtoBt_Gisam(:,2:end) ...
                            - repmat(this.r.tGtoBt_Gtrue,1,this.r.nStates-1);
            errTransGeomAbs = abs(errTransGeom);
            totErrTransGeom = sqrt(sum(errTransGeom.^2));
            sigmaTransGeom = zeros(3,this.r.nStates-1);
            for i = 1:this.r.nStates-1
                sigmaTransGeom(:,i) = sqrt(diag(this.r.covtGtoBt_Gisam(:,:,i+1)));
            end
            subplot(3,1,1);
            plot(this.r.imageNums(2:end),errTransGeomAbs(1,:),'-r'); hold on; grid on;
            plot(this.r.imageNums(2:end),errTransGeomAbs(2,:),'-g');
            plot(this.r.imageNums(2:end),errTransGeomAbs(3,:),'-b');
            plot(this.r.imageNums(2:end),totErrTransGeom,'-.k');            
            plot(this.r.imageNums(2:end),sigmaTransGeom(1,:),'--r');
            plot(this.r.imageNums(2:end),sigmaTransGeom(2,:),'--g');
            plot(this.r.imageNums(2:end),sigmaTransGeom(3,:),'--b');
            legend('{}^Gt_{GtoB} (x dir)','{}^Gt_{GtoB} (y dir)', ...
                '{}^Gt_{GtoB} (z dir)','Total Error','1-\sigma Bounds (x dir)', ...
                '1-\sigma Bounds (y dir)','1-\sigma Bounds (z dir)');
            title('Convergence of Center of Mass'); ylabel('Error [m]');
            ylim([-0.01 0.1]); xlim([0 this.r.imageNums(end)]);
            
            % Plot convergence of principal axes
            errPrincAx = zeros(3,length(ipoRuns));
            for i = 1:length(ipoRuns)
                errPrincAx(:,i) = rad2deg( abs( Log(this.r.RBttoGtrue'*ipoRuns{i}.RBttoG) ) );
            end
            totErrPrincAx = sqrt(sum(errPrincAx.^2));
            subplot(3,1,2);
            plot(ipoImageNums,errPrincAx(1,:),'-r'); hold on; grid on;
            plot(ipoImageNums,errPrincAx(2,:),'-g');
            plot(ipoImageNums,errPrincAx(3,:),'-b');
            plot(ipoImageNums,totErrPrincAx,'-.k');
            legend('|\delta\theta_1|','|\delta\theta_2|', ...
                '|\delta\theta_3|','||\delta\theta||');
            title('Convergence of Principal Axes');
            ylim([-1 max(totErrPrincAx)]); ylabel('Error [deg]'); 
            xlim([0 this.r.imageNums(end)]);
            
            % Plot convergence of inertia ratios; if not simulated, only plot this
            errJ = zeros(2,length(ipoRuns));
            for i = 1:length(ipoRuns)
                if ~ipoRuns{i}.isSingleAxis
                    errJ(:,i) = abs( ipoRuns{i}.inertiaRatiosOpt.getInertiaRatios() ...
                                        - this.r.Jtrue(1:2) );
                else
                    errJ(:,i) = [Inf Inf]';
                end
            end
            if ~this.r.isSimulated
                gcf; clf;
                subplot(2,1,1);
            else
                subplot(3,1,3);
            end
            plot(ipoImageNums,errJ(1,:),'-r'); hold on; grid on;
            plot(ipoImageNums,errJ(2,:),'-g');
            legend('J_1','J_2'); title('Convergence of the Inertia Ratios');
%             ylim([-0.01 0.05]); ylabel('Error [-]'); 
            ylim([-0.01 0.2]); ylabel('Error [-]'); 
            xlim([0 this.r.imageNums(end)]);
            
            % Plot whether the rotation was determined to be single or multi-axis, low or high
            % energy, and tri-axial or axis-symmetric
            if ~this.r.isSimulated
                subplot(2,1,2);
            else
                figure();
            end
            pCol = flipud([ zeros(1,length(ipoSingleMulti)+1); [ipoSingleMulti 0] ; ...
                [ipoEnergyState 0]; [ipoInertiaSymm 0] ]);
            pcolor([ipoImageNums ipoImageNums(end)+1],1:4,pCol);
            colormap(gray(3));
            colorbar('southoutside','Ticks',[0 1 2],'TickLabels', ...
                {'Single-Axis, N/A','LE, TA','Multi-axis, HE, AS'});
            xlabel('Image Number [-]'); title('Rotation Classification');
            
        end
        
        
    end
    
end






















