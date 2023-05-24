function [fitresult, gof] = createFit1(rise)
%CREATEFIT1(RISE)
%  创建一个拟合。
%
%  要进行 '无标题拟合 1' 拟合的数据:
%      Y 输出: rise
%  输出:
%      fitresult: 表示拟合的load拟合对象。
%      gof: 带有拟合优度信息的结构体。
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 05-May-2023 22:45:01 自动生成


%% 拟合: '无标题拟合 1'。
[xData, yData] = prepareCurveData( [], rise );

% 设置 fittype 和选项。
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Normalize = 'on';
opts.Robust = 'LAR';
opts.StartPoint = [4.52311098647382e-06 0.063538780663872 -4.65651007651062e-08 -1.67406215463767];

% 对数据进行模型拟合。
[fitresult, gof] = fit( xData, yData, ft, opts );

% 绘制数据拟合图。
figure( 'Name', '无标题拟合 1' );
h = plot( fitresult, xData, yData );
legend( h, 'rise', '无标题拟合 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% 为坐标区加标签
ylabel( 'rise', 'Interpreter', 'none' );
grid on


