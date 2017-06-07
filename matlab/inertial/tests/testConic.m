% Test the functionality of the Conic class

addpath('../');

% Example conic (centered off-origin, as given by nonzero D, E)
A = -1; B = 12; C = 2; D = -10; E = -10; F = -1;
K = [A B C D E F]';

% Create conic
conic = Conic(K);

% Fit the canonical conic
[Kc,PAtoC,aSemi,bSemi,type,dir] = conic.findCanonical();
RCtoA = PAtoC(1:2,1:2); tAtoC_A = PAtoC(1:2,3);

% Get axes of the actual, translated and canonical frames
xAx_A = [2 0]'; yAx_A = [0 2]';
xAx_B = xAx_A + tAtoC_A;
yAx_B = yAx_A + tAtoC_A;
xAx_C = RCtoA*xAx_A + tAtoC_A;
yAx_C = RCtoA*yAx_A + tAtoC_A;

% Plot the conic the actual and canonical frames
figure(1); clf;
subplot(1,2,1);
pts_A = conic.plot('actual','-b'); hold on; grid on; axis equal;
plot([0 xAx_A(1)],[0 xAx_A(2)],'-k','LineWidth',2);
plot([tAtoC_A(1) xAx_B(1)],[tAtoC_A(2) xAx_B(2)],'-.k','LineWidth',1);
plot([tAtoC_A(1) xAx_C(1)],[tAtoC_A(2) xAx_C(2)],'--k','LineWidth',2);
plot([0 tAtoC_A(1)],[0 tAtoC_A(2)],'-.r','LineWidth',2);
plot([0 yAx_A(1)],[0 yAx_A(2)],'-k','LineWidth',2);
plot([tAtoC_A(1) yAx_B(1)],[tAtoC_A(2) yAx_B(2)],'-.k','LineWidth',1);
plot([tAtoC_A(1) yAx_C(1)],[tAtoC_A(2) yAx_C(2)],'--k','LineWidth',2);
xlabel('x'); ylabel('y');
title('Conic in Actual Frame A');
legend('Conic in Frame A','Actual Frame A','Translated Frame B','Canonical Frame C','Translation x_0');
subplot(1,2,2);
pts_C = conic.plot('canonical','b'); hold on; grid on; axis equal;
plot([0 xAx_A(1)],[0 xAx_A(2)],'--k','LineWidth',2);
plot([0 yAx_A(1)],[0 yAx_A(2)],'--k','LineWidth',2);
legend('Conic in Frame C','Canonical Frame C');
title('Conic in Canonical Frame C');
xlabel('x'); ylabel('y');

% Find the closest point to a point
figure(2); clf;
pt = [-1 3.5]';
[pOpt,ptOpt] = conic.findClosestPoint(pt);
conic.plot('actual','-b'); hold on; grid on; axis equal;
plot(pt(1),pt(2),'oc');
plot(ptOpt(1),ptOpt(2),'ob');
plot([pt(1) ptOpt(1)],[pt(2) ptOpt(2)],'-r');
title('Finding the Closest Point on a Conic');

% Fit a conic to the canonical form and check to see if it lines up with that specified earlier
% Add a bit of noise to the data (otherwise eigenvalues are zero)
[Kell,Khyp,rSqEll,rSqHyp] = conic.fitConicCanonical(pts_C + normrnd(0,0.01,size(pts_C)));
disp('Comparison of fit canonical hyperbola Khyp (with noise)')
disp('to computed one Kc [Kc Khyp]');
[-Kc/Kc(end) Khyp]
