function [fitresult, gof] = createFit1(decline)
%CREATEFIT1(DECLINE)
%  创建一个拟合。
%
%  要进行 '无标题拟合 1' 拟合的数据:
%      Y 输出: decline
%  输出:
%      fitresult: 表示拟合的拟合对象。
%      gof: 带有拟合优度信息的结构体。
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 05-May-2023 23:00:39 自动生成


%% 拟合: '无标题拟合 1'。
[xData, yData] = prepareCurveData( [], decline );

% 设置 fittype 和选项。
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Normalize = 'on';
opts.Robust = 'LAR';
opts.StartPoint = [5.38119546318644e-08 -1.69411328847156 3.61577644858889e-06 -0.0775012994203096];

% 对数据进行模型拟合。
[fitresult, gof] = fit( xData, yData, ft, opts );

% 绘制数据拟合图。
figure( 'Name', '无标题拟合 1' );
h = plot( fitresult, xData, yData );
legend( h, 'decline', '无标题拟合 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% 为坐标区加标签
ylabel( 'decline', 'Interpreter', 'none' );
grid on


