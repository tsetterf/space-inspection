% Plot the C++ results and perform InertialProperties optimization on real or simulated data

addpath('../');

% Select real, simulated visually, or simulated analytically already processed with GTSAM
% dataType = 'real';
dataType = 'simVis';

% Select test case
% testCase = 'TA';
testCase = 'AS1';

% Get data depending on case
disp('Processing results ...');
baseDir = 'D:\PhD\vmshare\data\Results\';
if strcmp(dataType,'real') && strcmp(testCase,'TA')     % real data, tri-axial
%     results = Results('D:\PhD\vmshare\data\TS53T7R3', ...
%         'D:\PhD\vmshare\data\Results\TS53T7R3\TP_1511',2,10,253800,409000);
%         'D:\PhD\vmshare\data\TempResults',2,10,253800,409000);
%         'D:\PhD\vmshare\data\Results\TS53T7R3\TP_1511',2,10,253800,409000);
%         'D:\PhD\vmshare\data\Results\TS53T7R3\TP_1511',2,10,253800,357390);
%     results = Results('D:\PhD\vmshare\data\TS53T7R3', ...
%         'D:\PhD\vmshare\data\TempResults',2,10,253800,357390);
    plotter = Plotter(results);
elseif strcmp(dataType,'real') && strcmp(testCase,'AS1')
%     pmp = PlotMakerPro('D:\PhD\vmshare\data\TS74T3R1', ...
%         'D:\PhD\vmshare\data\TempResults',imgRange,'1511');
    pmp = PlotMakerPro('D:\PhD\vmshare\data\TS53T7R3', ...
        'D:\PhD\vmshare\data\TempResults',imgRange,'1511');
elseif strcmp(dataType,'simVis') && strcmp(testCase,'AS1')
%     disp(['Running simulated visual case AS1 with true J = [' num2str(0.0215/0.0116) ',1,1]']);            
%     results = Results('D:\PhD\vmshare\data\TS00AS1', ...
%         'D:\PhD\vmshare\data\Results\TS00AS1\TP_1511',5,norm([4 1 1]),329840,362400);
%         'D:\PhD\vmshare\data\TempResults',5,4,329840,392050);
    plotter = Plotter(results);
end

% error('I am done!');

%{

% Perform inertial properties optimization at every Nth state
ipoStates = 5:5:results.stateNums(end);
ipo = {};
for i = 1:length(ipoStates)
    disp(['Performing ipo ' num2str(i) ' of ' num2str(length(ipoStates))]);
    ipo{i} = InertialPropertiesOpt(results,ipoStates(i));
    ipo{i}.optimize();
    results.addInertialPropertiesOpt(ipo{i},'isam',ipoStates(i));
end

% Create InertialPropertiesOpt object and optimize at the last timestep
disp('Performing inertial properties optimization...');
inertialPropertiesOpt = InertialPropertiesOpt(results,results.stateNums(end));
tic;
inertialPropertiesOpt.optimize();
ipoTime = toc;
disp(['Inertial properties optimization time: ' num2str(ipoTime) ' s']);
results.addInertialPropertiesOpt(inertialPropertiesOpt,'isam',results.stateNums(end));

%}
% %}

% figure(10); clf;
% plotter.plotConvergence();

% figure(5); clf;
% plotter.plotExecTimes();

disp('Plotting...');

% Plot the results
figure(2); clf;
inertialPropertiesOpt.plot();
if results.isSimulated
    sB = max(max(max(results.omegaBt_Gisam)));      % scale for B axes
    RBttoBthat = results.RBttoGisam(:,:,results.stateNums(end)+1)'*results.RBttoGtrue;
    xBt_Bthat = RBttoBthat(:,1); yBt_Bthat = RBttoBthat(:,2); zBt_Bthat = RBttoBthat(:,3); 
    plot3([0 xBt_Bthat(1)],[0 xBt_Bthat(2)],[0 xBt_Bthat(3)],'-r');    
    plot3([0 yBt_Bthat(1)],[0 yBt_Bthat(2)],[0 yBt_Bthat(3)],'-g'); 
    plot3([0 zBt_Bthat(1)],[0 zBt_Bthat(2)],[0 zBt_Bthat(3)],'-b');
    legend('Measured \omega_B','Estimated x_E','Estimated y_E', ...
    'Estimated z_E','Estimated x_B','Estimated y_B', ...
    'Estimated z_B','Final Conic Fit 1','Final Conic Fit 2','Final Conic Fit 3', ...
    'Optimally Fit \omega_B','True \omega_B','True x_B','True y_B', ...
    'True Bz_B','Location','eastoutside');
else
    legend('Measured \omega_B','Estimated x_E','Estimated y_E', ...
    'Estimated z_E','Estimated x_B','Estimated y_B', ...
    'Estimated z_B','Final Conic Fit 1','Final Conic Fit 2','Final Conic Fit 3', ...
    'Estimated \omega_B','True \omega_B','Location','eastoutside');
end
% error('volunstop');

% Plot C++ results
figure(3); clf;
% plotter.plotGlory(results.stateNums(end-1),0);
% drawnow;
% plotter.plotQuad(results.imageNums(end),1); %results.imageNums(end),1);
% drawnow;

figure(4); clf;
plotter.plotErrors(results.imageNums(end));

figure(6); clf;
plotter.plotMapOrtho(results.stateNums(end));

figure(7); clf;
plotter.plotResultsSummary(results.stateNums(end));

