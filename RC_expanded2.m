
%% Waveform classification
% clear;clc;

rise = csvread('rise.csv');
decline = csvread('decline.csv');

%Fit of the rise edge of memrsitor


[xData, yData] = prepareCurveData( [], rise );

ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Normalize = 'on';
opts.Robust = 'LAR';
opts.StartPoint = [4.52311098647382e-06 0.063538780663872 -4.65651007651062e-08 -1.67406215463767];


% Fit model to data.
[rise_fit, gof] = fit( xData, yData, ft, opts );

% Fit of the decline edge of memrsitor
[xData, yData] = prepareCurveData( [], decline );

ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Normalize = 'on';
opts.Robust = 'LAR';
opts.StartPoint = [5.38119546318644e-08 -1.69411328847156 3.61577644858889e-06 -0.0775012994203096];

% Fit model to data.
[decline_fit, gof] = fit( xData, yData, ft, opts );

x1=-1:0.001:150;
y1=rise_fit(x1); %fitted rise curve
x2=-100:0.001:50;
y2=decline_fit(x2); %fitted decline curve

[M,I1] = min(abs(y1 - 1*10^-5));
[M,I2] = min(abs(y1 - 3.15E-06));
rise = y1(I2:I1);

[M,I1] = min(abs(y2 - 1*10^-5));
[M,I2] = min(abs(y2 - 3.15E-06));
decline = y2(I1:end);
% clear
% rise = csvread('rise.csv');
% decline = csvread('decline.csv');
% 
% %Fit of the rise edge of memrsitor
% [xData, yData] = prepareCurveData( [], rise );
% ft = fittype( 'poly7' );
% opts = fitoptions( 'Method', 'LinearLeastSquares' );
% opts.Robust = 'LAR';
% 
% % Fit model to data.
% [rise_fit, gof] = fit( xData, yData, ft, opts );
% 
% [xData, yData] = prepareCurveData( [], decline );
% 
% % Fit of the decline edge of memrsitor
% ft = fittype( 'poly8' );
% opts = fitoptions( 'Method', 'LinearLeastSquares' );
% opts.Robust = 'LAR';

% Fit model to data.
% [decline_fit, gof] = fit( xData, yData, ft, opts );
% 
% x1=1:0.001:30;
% rise=rise_fit(x1); %fitted rise curve
% x2=1:0.001:32;
% decline=decline_fit(x2); %fitted decline curve

% ----------------------DM_RC Parameters----------------------

N = 5;
ML = 4;
% ----------------------DATASET----------------------
sample = 11;
step = 2000;
Data = zeros(1, 2*step);
p1 = 3*sin(pi*2*(0:sample)/(sample+1));
p2(1:(sample+1)/2) = 3;
p2((sample+1)/2+1:sample+1) = -3;
for i = 1:2*step/sample
    q = unidrnd(2);
    if q == 1
        Data((sample+1)*(i-1)+1:(sample+1)*i) = p1;
        Label((sample+1)*(i-1)+1:(sample+1)*i) = 0;
    else
        Data((sample+1)*(i-1)+1:(sample+1)*i) = p2;
        Label((sample+1)*(i-1)+1:(sample+1)*i) = 1;
    end
end

Data = Data + 6;
 max_len = length(dec2bin(max(Data)));
% ----------------------TRAIN----------------------
% initialize input stream
Input = Data(1:step);

% generate target
Target = Label(1:step);

input = [];
% mask process
Mask = 2*unidrnd(2, N, ML)-3;
% Mask = unidrnd(3, N, ML)/3;
for j = 1:N
    for i = 1:step
        input(j, (i-1)*ML+1:ML*i) = Input(i)*Mask(j, :);
    end
end
max_len = length(dec2bin(max(max(input))));

% memristor output
memout = [];


for i = 1:length(input(:,1))
    G = 3.15E-06;
    % count = 0;
    for j = 1:length(input(1,:))
    % if(count == 8)
    %     count = 0;
    %     G = 3.15E-06;
    % end
    
       
    pulse = dec2bin(input(i,j));
    if(length(pulse) < max_len)
        for nn = 1:max_len-length(pulse)
            G = UpdateMemristor(G,0,rise,decline);
        end
    end
    
    for k = 1: length(pulse)
        if bin2dec(pulse(k)) == 1
 
                G = UpdateMemristor(G,1,rise,decline);

           
        else
             G = UpdateMemristor(G,0,rise,decline);
             
        end
    end
    memout(i,j) = G*10^5;
    % count = count + 1;
    end
end

% states collection
states = [];
for i = 1:step
    a = memout(:, ML*(i-1)+1:ML*i);
    states(:, i) = a(:);
end
X = [ones(1,step); states];

% linear regression
Wout = Target*pinv(X);

% ----------------------TEST----------------------
% initialize input stream
Input = Data(step+1:4000);

% generate target
Target = Label(step+1:4000);

% mask process
input = [];
for j = 1:N
    for i = 1:step
        input(j, (i-1)*ML+1:ML*i) = Input(i)*Mask(j, :);
    end
end

max_len = length(dec2bin(max(max(input))));
% memristor output
memout = [];


for i = 1:length(input(:,1))
    G = 3.15E-06;
    % count = 0;
    for j = 1:length(input(1,:))
    % if(count == 8)
    %     count = 0;
    %     G = 3.15E-06;
    % end
    
       
    pulse = dec2bin(input(i,j));
    if(length(pulse) < max_len)
        for nn = 1:max_len-length(pulse)
            G = UpdateMemristor(G,0,rise,decline);
        end
    end
    for k = 1: length(pulse)
        if bin2dec(pulse(k)) == 1
 
                G = UpdateMemristor(G,1,rise,decline);

           
        else
             G = UpdateMemristor(G,0,rise,decline);
             
        end
    end
    memout(i,j) = G*10^5;
    % count = count + 1;
    end
end

% states collection
for i = 1:step
    a = memout(:, ML*(i-1)+1:ML*i);
    states(:,i) = a(:);
end
X = [ones(1,step);states];

% system output
Out = Wout*X;
NRMSE = sqrt(mean((Out(10:end)-Target(10:end)).^2)./var(Target(10:end)));
% sprintf('%s',['NRMSE:',num2str(NRMSE)])

% ----------------------PLOT----------------------
figure;
subplot(2, 1, 1);
plot(Input, 'linewidth', 1);
hold on;
plot(Input, '.');
axis([30, 230, 2, 10])
ylabel('Input')
set(gca,'FontName', 'Arial', 'FontSize', 11);
subplot(2, 1, 2);
plot(Target, 'linewidth', 2);
hold on;
plot(Out, 'linewidth',1);
axis([30, 230, -0.2, 1.2])
str1 = 'Target';
str2 = 'Output';
lg = legend(str1, str2);
set(lg, 'Orientation', 'horizon');
ylabel('Prediction')
xlabel('Index')

set(gca,'FontName', 'Arial', 'FontSize', 11);
% set(gcf, 'unit', 'normalized', 'position', [0.2, 0.2, 0.6, 0.35]);