function visualization_ellipse(param, abphi, Posi, regl, out_folder, c)
% affich_size(obj.param,results,col)
% A function that plots maps of cell size, registers the figures
% as .png and .fig
%
% + col the colour for the plot of the cell size
% for each averaged time
display(['Representation of size map, time: ',num2str(c)]); % for the user to keep track
I = imadjust(param.reader.readSpecificImage(c));      
line_array = zeros(regl, 4);
for win =1:regl % for each subimage
    xo =  Posi(win,1)+floor(param.tile_size/2); % find positions of the center of the subimages
    yo =  Posi(win,2)+floor(param.tile_size/2);

    a=abphi(1,win);
    b=abphi(2,win);
    phi=abphi(3,win);
    %this seems to take too long, I will just draw a line instead
    %f_ellipseplotax(a,b,xo,yo,phi,param.col,0.5,Figure.CurrentAxes);
    
    amp = .5*log( a/b); %equivalent 

    x1 = xo+cos(phi+pi/2)*amp*param.scale; 
    y1 = yo+sin(phi+pi/2)*amp*param.scale;
    x2 = xo-cos(phi+pi/2)*amp*param.scale;
    y2 = yo-sin(phi+pi/2)*amp*param.scale;

    line_array(win,:) = [x1 y1 x2 y2];
end
RGB = insertShape(I, 'line', line_array, 'LineWidth', 5, 'Color', param.col);
imwrite(RGB, [out_folder filesep 'img_' num2str(c, '%04d') '.png']);

% set(gca,'XTickLabel','','YTickLabel','');
% X3 = [param.siz(1)-250 param.siz(1)-100];
% Y3 = [param.siz(2)-250 param.siz(2)-250];
% plot(X3,Y3,'.b-','LineWidth',1);                                   % scaling bar
% text(param.siz(1)-250,param.siz(2)-300,'0.5','Color','black','FontSize',12);
% img = getframe(gcf);
% imwrite(img.cdata, [obj.param.pathout  ['Size_map' num2str(param.timestep) '_'] num2str(c) '.png']); % save image png
% %savefig(Figure,[obj.param.pathout   ['Size_map' num2str(param.timestep) '_'] num2str(c) '.fig'])     % save image fig
% close(Figure);


end

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

% parametric equation of the ellipse
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
end