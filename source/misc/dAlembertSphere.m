% Returns the solution of the wave equation alongside a cut through a
% sphere
% Input:
% r   - Distance to the center
% t   - Time
% c   - Wave speed
% u0  - Initial dilation
% Output: 
% sol - Solution
% Note: Produces NaN for r = 0
function sol = dAlembertSphere(r,t,c, u0)
if r == 0
    sol = NaN;
else
    sol = ( (r+c*t) .* u0(r + c*t) + (r - c*t).*u0(r-c*t))./(2.*r) ;
end
end
