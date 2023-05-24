clear
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

[M,I1] = min(abs(y1 - 2*10^-5));
[M,I2] = min(abs(y1 - 3.15E-06));
rise = y1(I2:I1);

[M,I1] = min(abs(y2 - 2*10^-5));
[M,I2] = min(abs(y2 - 3.15E-06));
decline = y2(I1:end);

%dataset

[data,t]=MackeyGlass(10000,17);
%data = load('mackey_glass.txt');
%data = data - min(data);
data = data(1000:end);
data = data*100;
 data = data - 20;
data = round(data);
%data = dec2bin(data);
train_data = data(1001:4000);
train_output = data(1002:4001);

test_data = data(4002:7001);
test_output = data(4003:7002);
% step = 1000;
% data = HenonMap(2*step+1)';
% 
% data = data*50;
% data = data+70; 
% 
% train_data = data(1:step);
% train_output = data(2:step+1);
% test_data = data(step+1:2*step);
% test_output = data(step+2:2*step+1);


N = 25; % number of Memrestor
ML = 15;

% mask process
 Mask = 2*unidrnd(2, N, ML)-3;
 % Mask = unidrnd(3, N, ML)/3;

% Mask = load('mask.mat');
% Mask = Mask.Mask;


input = [];
for i = 1:N
    for j = 1:length(train_data)
        input(i, (j-1)*ML+1:ML*j) = train_data(j)*Mask(i, :);
    end
end

Vmax = 100;
Vmin = 20;
UL = max(max(input));
DL = min(min(input));
input = (input-DL)/(UL-DL)*(Vmax - Vmin)+Vmin;
max_len = length(dec2bin(max(max(input))));
%---------------------train-------------------

W_res = zeros(length(input(:,1)),length(input(1,:)));


for i = 1:length(input(:,1))
    G = 3.15E-06;
    % count = 0;
    for j = 1:length(input(1,:))
    % if(count == 16)
    %     count = 0;
    %     G = 3.15E-06;
    % end


    pulse = dec2bin(round(input(i,j)));

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
    W_res(i,j) = G*10^5;
    % count = count + 1;
    end
end

% W_res = input;
states = [];
for i = 1:length(train_data)
    a = W_res(:, ML*(i-1)+1:ML*i);
    states(:, i) = a(:);
end

X = [ones(1,length(train_data)); states];

% W_out = regress(train_output',X);
%W_out = inv(W_res * W_res' + eye(2500)) * (W_res * bin2dec(train_output)'); 
W_out = train_output'*pinv(X);


% mask process

input = [];

for i = 1:N
    for j = 1:length(train_data)
        input(i, (j-1)*ML+1:ML*j) = test_data(j)*Mask(i, :);
    end
end

UL = max(max(input));
DL = min(min(input));
input = (input-DL)/(UL-DL)*(Vmax - Vmin)+Vmin;
max_len = length(dec2bin(max(max(input))));
%---------------------test-------------------

W = zeros(length(input(:,1)),length(input(1,:)));


for i = 1:length(input(:,1))
    G = 3.15E-06;
    % count = 0;
    for j = 1:length(input(1,:))

   % if(count == 16 )
   %      count = 0;
   %      G = 3.15E-06;
   % end
   pulse = dec2bin(round(input(i,j)));

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
    W(i,j) = G*10^5;
    % count = count + 1;
    end
end

% W = input;
states = [];
for i = 1:length(train_data)
    a = W(:, ML*(i-1)+1:ML*i);
    states(:, i) = a(:);
end

X = [ones(1,length(train_data)); states];

output = W_out * X;
output = output';

% pred = output(1801:end);
% true = test_output(2000-200+1:end);

plot(output(1:200))
hold on
plot(test_output(1:200))
legend('output','test_output')



% plot(pred)
% hold on
% plot(true)
% legend('predictions','true values')
% 
NRMSE = sqrt(mean((output - test_output).^2)./var(test_output))            
        
% for i =1:3:length(output)-3
% output2(i)=(output(i)+output(i+1)+output(i+2))/3;
% output2(i+1)=(output(i)+output(i+1)+output(i+2))/3;
% output2(i+2)=(output(i)+output(i+1)+output(i+2))/3;
% end
% 
% output2 = [];
% for i =1:4:length(output)-4
% output2(i)=(output(i)+output(i+1)+output(i+2)+output(i+4))/4;
% output2(i+1)=(output(i)+output(i+1)+output(i+2)+output(i+4))/4;
% output2(i+2)=(output(i)+output(i+1)+output(i+2)+output(i+4))/4;
% output2(i+3)=(output(i)+output(i+1)+output(i+2)+output(i+4))/4;
% end
