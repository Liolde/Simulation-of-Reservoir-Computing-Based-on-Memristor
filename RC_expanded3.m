
%% Waveform classification
clear;clc;

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
x2=-100:0.001:100;
y2=decline_fit(x2); %fitted decline curve

[M,I1] = min(abs(y1 - 10^-5));
[M,I2] = min(abs(y1 - 3.15E-06));
rise = y1(I2:I1);

[M,I1] = min(abs(y2 - 10^-5));
[M,I2] = min(abs(y2 - 3.15E-06));
decline = y2(I1:end);

%%------------------------spoken digit recognization--------------------
name1 = '-george-';
name2 = '-jackson-';
name3 = '-lucas-';
name4 = '-nicolas-';
name5 = '-theo-';
name6 = '-yweweler-';


N = 30; % number of memrestor
ML = 20;

% bound = 2;

% mask
Mask = 2*randi([0,1],13,ML,N)-1;
% Mask = randi([1,3],13,ML,N)/3;
% train
W_out = [];
X = [];
Target = [];
L= [];

for index = randperm(45)
 for number = randperm(10)
    for nm = randperm(6)
                    name = [];
                    switch nm
                        case 1
                            name = name1;
                        case 2 
                            name = name2;
                        case 3 
                            name = name3;
                        case 4 
                            name = name4;
                        case 5
                            name = name5;
                        case 6
                            name = name6;
                    end
                 
             str = 'recordings\' + "" + num2str(number-1) + "" + name + "" + num2str(index-1) + "" + '.wav';
             [x,fs] = audioread(str);
             c = mfcc2(x,fs)*10;
             % 
             
             [x,fs] = audioread(str);

             % f = LyonPassiveEar(x,8000,250);
            
             l = length(c(1,:));

             t = zeros(10,l);
             t(nm,1:end) = 1;
             Target = [Target, t];

             % f = f*10^7;

             %Mask process
             input = [];
                 for i=1:N
                     for j=1:l
                         input(i, ML*(j-1)+1:ML*j) = abs(c(:,j))'*Mask(:,:,i);            
                     end
                 end
              max_len = length(dec2bin(max(max(input))));
             
             W_res = zeros(length(input(:,1)),length(input(1,:)));

            for i = 1:length(input(:,1))
                G = 3.15E-06;
                % count = 0;
                for j = 1:length(input(1,:))
                % if(count == bound)
                %     count = 0;
                %     G = 3.15E-06;
                % end
                pulse = dec2bin(round(input(i,j)));
                    if(length(pulse) < max_len)
                    for nn = 1:max_len-length(pulse)
                        G = UpdateMemristor(G,0,rise,decline);
                    end
                    end
               % P = 0; % Number of the positive pulse

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


            states = [];
            for i = 1:l
                a = W_res(:, ML*(i-1)+1:ML*i);
                states(:, i) = a(:);
            end
            X = [X states];
sprintf('%s',['number:',num2str(number-1),',name:',name,',index:',num2str(index)])
        end
    end
end

            
W_out = Target*X'*pinv(X*X');         
             
%---------------------test------------------             


Target = [];
X = [];
Label = [];
L = [];

           for index = 45+randperm(5)
                  for number = randperm(10)
                for nm = randperm(6)
                    name = [];
                    switch nm
                        case 1
                            name = name1;
                        case 2 
                            name = name2;
                        case 3 
                            name = name3;
                        case 4 
                            name = name4;
                        case 5
                            name = name5;
                        case 6
                            name = name6; 
                    end


             str = 'recordings\' + "" + num2str(number-1) + "" + name + "" + num2str(index-1) + "" + '.wav';
             % data = audioread(str);


             [x,fs] = audioread(str);
             c = mfcc2(x,fs)*10;
             % f = LyonPassiveEar(x,8000,250);
             
             l = length(c(1,:));
             L = [L, l];

             t = zeros(10,l);
             t2 = zeros(10,1);
             t(nm,1:end) = 1;
             t2(nm) = 1;
             Target = [Target, t];
             Label = [Label, t2];

             % f = f*10^7;
             %Mask process

           
             input = [];
             for i=1:N
                 for j=1:l
                     input(i, ML*(j-1)+1:ML*j) = abs(c(:,j))'*Mask(:,:,i);            
                 end
             end

             max_len = length(dec2bin(max(max(input))));
             W_res = zeros(length(input(:,1)),length(input(1,:)));

            for i = 1:length(input(:,1))
                 G = 3.15E-06;
                 % count = 0;
                for j = 1:length(input(1,:))
                % if(count == bound)
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
            %             if P ~= 0
            %                 G = UpdateMemristor(G,0,rise,decline);
            %                 P = P - 1;
            %             end

                           G = UpdateMemristor(G,1,rise,decline);

                   
                    else  
                         G = UpdateMemristor(G,0,rise,decline);

                    end
                end
                W_res(i,j) = G*10^5;
                % count = count + 1;
                end
            end


            states = [];
            for i = 1:l
                a = W_res(:, ML*(i-1)+1:ML*i);
                states(:, i) = a(:);
            end
            X = [X states];
sprintf('%s',['number:',num2str(number-1),',name:',name,',index:',num2str(index)])
        end
    end
end

Y = W_out*X;

M = [];
Out = zeros(10,length(Label(1,:)));
total = length(Label);

for i=1:total
    if i == 1
        st = 1;
    else
        st = (sum(L(1:i-1)))+1;
    end  
    ed = sum(L(1:i));
    M(:,i) = mean(Y(:,st:ed),2);
    [~,id] = max(M(:,i));
    Out(id,i) = 1;
end

count = 0;
for i=1:length(Out)
    [~,id1] = max(Label(:,i));
    [~,id2] = max(Out(:,i));
    if id1 == id2
        count = count +1;
    end
end