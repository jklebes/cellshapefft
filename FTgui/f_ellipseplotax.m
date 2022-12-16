function [hEllipse,hAxes] = f_ellipseplotax(a,b,x0,y0,phi,col,ldw,axh)
%  Isabelle Bonnet
%  Begin: 2009/10/12 
%  Draw an ellipse with given parameters.
%  
%   Input parameters:
%       a           Value of the HALF-major axis
%       b           Value of the HALF-minor axis
%       x0          Abscissa of the center of the ellipse
%       y0          Ordinate of the center of the ellipse
%       phi         Angle between X-axis and the MAJOR axis - units: RADIANS
%       lineStyle   Definition of the plotted line style
%
%   Output:
%       h    Handle of the ellipse
%
%  Usage:
%       f_ellipseplot(maxi,mini,a,b,phi);
%       or  
%       hEllipse = f_ellipseplot(5,3,1,-2,pi/4);
%       set(h,'LineWidth',2);
%   A slight change added by Melina Durande regarding the handling of the
%   axes and the gestion of the colors and line width
%


if (nargin < 2)|(nargin > 8),
   error('Please see help for INPUT DATA.'); 
return;
end

angle = [-0.003:0.01:2*pi];

% Parametric equation of the ellipse
%----------------------------------------
 x = a*cos(angle);
 y = b*sin(angle);


% Coordinate transform 
%----------------------------------------
 X = cos(phi)*x - sin(phi)*y;
 Y = sin(phi)*x + cos(phi)*y;
 X = X + x0;
 Y = Y + y0;

       
% Plot the ellipse
%----------------------------------------
 hEllipse = plot(X,Y,'Color',col,'parent',axh);
 set(hEllipse,'LineWidth',ldw);
 axis equal;
 
%%% Add lines representing the major and minor axes of an
% ellipse on the current plot.
 
 % INPUT
%   axes      [major minor] axes of the ellipse

% Major Axis coordinates
x_maj = a*cos(phi);
y_maj = a*sin(phi);

% Draw the line
Major = line(x_maj*[-1 1]+x0,y_maj*[-1 1]+ y0,'Color',col,'parent',axh);
set(Major,'LineWidth',ldw);

% Minor Axis coordinates
x_min = b*cos(phi+pi/2);
y_min = b*sin(phi+pi/2);

% Draw the line
Minor = line(x_min*[-1 1]+ x0,y_min*[-1 1]+ y0,'Color',col,'parent',axh);
set(Minor,'LineWidth',ldw);
% Return the line element handles so it can be customized.
hAxes = [Major Minor];




 
 
 
 
 
 
 
 
 
 
 
 
 