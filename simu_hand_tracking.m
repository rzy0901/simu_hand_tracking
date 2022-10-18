% Assume zero amplitude loss and phase changed caused by hardwares.
% Assume hand moving with a constant velocity.
% Hand tracking is implemented via intersections of two ellipses, assuming the
% initial position of hand known.
clc; clear; close all;
fc = 250;                       
fs = 1e3;
t = 0:1/fs:0.75;
c = 340;                        % voice velocity (m);
L1 = 2; L2 = 1;
Mic1.x = 0; Mic1.y = L1;        % Microphone location (m);
Mic2.x = 0; Mic2.y = -L2;
Hand.x = 2; Hand.y = 1;         % Hand initial location(m);
Hand.v_x = 1; Hand.v_y = 1;     % Hand velocity (m/s);
Speaker.x = 0; Speaker.y = 1;   % Speaker location (m);   
d1 = sqrt((Hand.x+Hand.v_x*t-Speaker.x).^2+(Hand.y+Hand.v_y*t-Speaker.y).^2) + ...
sqrt((Hand.x+Hand.v_x*t-Mic1.x).^2+(Hand.y+Hand.v_y*t-Mic1.y).^2);
d2 = sqrt((Hand.x+Hand.v_x*t-Speaker.x).^2+(Hand.y+Hand.v_y*t-Speaker.y).^2) + ...
sqrt((Hand.x+Hand.v_x*t-Mic2.x).^2+(Hand.y+Hand.v_y*t-Mic2.y).^2);
% Transmit signal.
% x = cos(2*pi*fc*t);
% Received signal;
y1 = 2*cos(2*pi*fc*t-2*pi*fc*d1/c);
y2 = 2*cos(2*pi*fc*t-2*pi*fc*d2/c);
% Apply lowpass filter to get I/Q.
I1 = lowpass(y1.*cos(2*pi*fc*t),1.5*fc,fs); 
I2 = lowpass(y2.*cos(2*pi*fc*t),1.5*fc,fs); 
Q1 = lowpass(-y1.*sin(2*pi*fc*t),1.5*fc,fs);
Q2 = lowpass(-y2.*sin(2*pi*fc*t),1.5*fc,fs);
% Combine I/Q.
B1 = complex(I1,Q1);
B2 = complex(I2,Q2);
phi1 = angle(B1);
phi2 = angle(B2);
% Use phase change between adjacent time slots to calculate distance change
d1_change = c*(phi1 - [phi1(1) phi1(1:end-1)])/(-2*pi*fc);
d2_change = c*(phi2 - [phi2(1) phi2(1:end-1)])/(-2*pi*fc);
d1_change = medfilt1(d1_change,15); % Remove abnormal peaks.
d2_change = medfilt1(d2_change,15);
temp1=d1 - [d1(1) d1(1:end-1)];
temp2=d2 - [d2(1) d2(1:end-1)];
figure; 
sgtitle(sprintf('Hand moves from (2,1) with speed (%d,%d) m/s',Hand.v_x,Hand.v_y));
subplot(321),plot(I1,Q1);title('Mic1');xlabel('I signal');ylabel('Qsignal'); axis equal;
subplot(322),plot(I2,Q2);title('Mic2');xlabel('I signal');ylabel('Qsignal'); axis equal;
subplot(323),plot(t,phi1); title('Phase change for Mic1'); xlabel('t(s)'); ylabel('Phase (rad)');
subplot(324),plot(t,phi2); title('Phase change for Mic2'); xlabel('t(s)'); ylabel('Phase (rad)');
subplot(325),plot(t,d1_change); title('Distance change per time slot for d1'); xlabel('t(s)'); ylabel('Distance change (m)');
subplot(326),plot(t,d2_change); title('Distance change per time slot for d2'); xlabel('t(s)'); ylabel('Distance change (m)');
figure;
subplot(221),stem(t,temp1); title('Ground Truth d1 change per time slot'); xlabel('t(s)'); ylabel('Distance change (m)');
subplot(222),stem(t,temp2); title('Ground Truth d2 change per time slot'); xlabel('t(s)'); ylabel('Distance change (m)');
subplot(223),stem(t,d1_change); title('Estimated d1 change per time slot'); xlabel('t(s)'); ylabel('Distance change (m)');
subplot(224),stem(t,d2_change); title('Estimated d2 change per time slot'); xlabel('t(s)'); ylabel('Distance change (m)');
%% trajectory estimation (post processing)
% Assume initial position of hand known
D1 = d1(1); D2 = d2(2);
x_true = Hand.x + Hand.v_x*t;
y_true = Hand.y + Hand.v_y*t;
x_estimate = [Hand.x];
y_estimtae = [Hand.y];
for ii = 2:1:length(t)
    D1 = D1 + d1_change(ii);
    D2 = D2 + d2_change(ii);
    % two ellipses' intersection
    syms x y;
    eqns = [sqrt((x-Mic1.x)^2+(y-Mic1.y)^2)+sqrt((x-Speaker.x)^2+(y-Speaker.y)^2)-D1==0, sqrt((x-Mic2.x)^2+(y-Mic2.y)^2)+sqrt((x-Speaker.x)^2+(y-Speaker.y)^2)-D2==0];
    S = solve(eqns,[x y]);
    tempx = double(S.x);
    tempy = double(S.y);
    x_estimate = [x_estimate tempx(tempx>0)]; % Select the solution in the right of y-axis.
    y_estimtae = [y_estimtae tempy(tempx>0)];
    fprintf('%d/%d\n',ii,length(t));
end
figure;scatter(x_estimate,y_estimtae);
hold on;
scatter(x_true,y_true);ylim([1,2]);
legend('estimation','ground truth','Location','northeast');