function [value, isterminal, direction] = myEventFunction(t, y)
    %弹目距离小于1——10m中断
    global rr
    value = ((y(14)-y(7))^2+(y(15)-y(8))^2+(y(16)-y(9))^2)^0.5-1; % 触发事件的条件，例如 y(1) - 1
    rr=[rr;value+1];
    isterminal = 1; % 1: 中断求解，0: 不中断
    direction = -1; % 0: 不限制方向, 1: 正向, -1: 反向
end

