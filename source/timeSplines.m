function [f] = timeSplines(refinement,index,derivative)
T = linspace(0,1,2^refinement+1);
T = [-T(2), T, T(end) + T(2)];

if derivative == 0
    f = @(x) -(x - T(index)).^3 .* (x - T(index+2)).^3 .* (x >= T(index)) .* (x <= T(index+2));
    if index == 1
        f = @(x)  1e3 * - (x - T(index+2)).^3 .* (x >= T(index)) .* (x <= T(index+2));
    end
elseif derivative == 1
    f = @(x)  3 * (x - T(index)).^2 .* (T(index+2)-x).^2 .* ...
        (T(index) + T(index+2) - 2 * x) .* (x >= T(index)) .* (x <= T(index+2));
    if index == 1
        f = @(x) 1e3 *  - 3*(x - T(index+2)).^2 .* (x >= T(index)) .* (x <= T(index+2));
    end
elseif derivative == 2
    f = @(x)  1e6 * 6 * (T(index) -x) .* (x - T(index+2)) .* ...
        (T(index)^2 + 3 * T(index) * T(index+2) - 5 * T(index) * x + T(index+2)^2 - 5 * T(index+2) * x ...
        + 5 * x.^2).* (x >= T(index)) .* (x <= T(index+2));
    if index == 1
        f = @(x) 1e3 * - 6*(x - T(index+2)) .* (x >= T(index)) .* (x <= T(index+2));
    end
else
    error('Not implemented')
end


end

