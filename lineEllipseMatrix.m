function [ x1, x2, y1, y2 ] = lineEllipseMatrix(a, b, h, k, alpha, p, q )
%  LineEllipse  Intersections of generalized ellipse and lines in Cartesian plane
%
%  [x1, x2, y1, y2] = LineEllipse(a, b, h, k, alpha, p, q ), where
%     ellipse is defined by 
%        (x*cos(alpha) + y * sin(alpha))^2 / a^2 + (x*sin(alpha) + y *
%        cos(alpha))^2 / b^2 = 1, with alpha (in radian) as the clockwise 
%        axis of rotation for the ellipse, and
%     line is defined by 
%         x/p + y/q = 1, with p and q x- and y- intercepts, respectively
%  The two points of intersection are (x1,x2, y1, y2) 
%  POST CONDITION:
%  If (line does not intersect) NaN is returned
%  If a single point of intersection, (x1, y1) = (x2, y2)
%  If two points of intersetion, two distinct (x1,y1) and (x2, y2) returned
        
% Copyright: Paurakh Rajbhandary, Stanford University
%%
% Parameter of line
m = -q./p;
c = q;

% Parameter for ellipse
A = ( cos(alpha)^2 / a^2  + sin(alpha)^2 / b^2 );
B = - 2 * cos (alpha) * sin(alpha) * (1/a^2 - 1/b^2);
C = ( sin(alpha)^2 / a^2  + cos(alpha)^2 / b^2 );

% Handle case with infinite slope
if (q == inf)
    disp('vertical line');
    x1 = p;
    x2 = p;
    M = C;
    N = B*p-(2*C*k+B*h);
    O = A*p^2 - p*(2*A*h + k*B) + (A*h^2+B*h*k+C*k^2-1);
    determinant = (N^2 - 4* M * O);
    y1  = (-N + sqrt(determinant))/ (2*M);
    y2  = (-N - sqrt(determinant))/ (2*M);
else    
    M = A + B*m + C* m.^2;
    N = B*c + 2*C*m.*c - 2*A*h - k*B - m* (2*C*k+ B* h);
    O = C*c.^2 - c* (2*C*k + B*h) + A*h^2 + B*h*k + C*k^2 -1 ;
    
    determinant = (N .^ 2 - 4 * M .* O);
    
    x1  = (- N + sqrt(determinant)) ./ (2*M);
    x2  = (- N - sqrt(determinant)) ./ (2*M);
    
    y1 = m.*x1 + c;
    y2 = m.*x2 + c;
end
% cond1 = ( abs(  ((x1-h)*cos(alpha) + (y1-k)*sin(alpha))^2 /a^2 + ((x1-h)*sin(alpha) - (y1-k)*cos(alpha))^2/b^2    - 1) <= (1e-5));
% 
% cond2 =  ( abs(((x2-h)*cos(alpha) + (y2-k)*sin(alpha))^2/a^2 + ((x2-h)*sin(alpha) - (y2-k)*cos(alpha))^2/b^2      - 1) <= (1e-5)) ;
% if ((imag(x1) ~=0) || (imag(x2) ~=0) || (imag(y1) ~=0) || (imag(y2) ~=0))
%     x1 = NaN;
%     x2 = NaN;
%     y1 = NaN;
%     y2 = NaN;
% elseif (cond1 == 1 && cond2 == 0)
%     x2 = x1;
%     y2 = y1;
% elseif (cond1 == 0 && cond2 == 1)
%     x1 = x2;
%     y1 = y2;
% elseif (cond1 == 0 && cond2 == 0)
%     x1 = NaN;
%     x2 = NaN;
%     y1 = NaN;
%     y2 = NaN;
% end

end


