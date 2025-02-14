function visualization_strain(param, SangS, quality, tileCoords, nTiles, out_folder, c)
% affich_result(obj.param,results,col)
% A function that plots maps of cell deformations, registers the figures
% as .png and .fig
%
%--------------------------------------------------------------------------

display(['Representation of deformation map, time: ',num2str(c)]); % for the user to keep track
%                 I = imadjust(imcomplement(obj.param.reader.readSpecificImage(obj.param.time_points(c))));
I = imadjust(param.reader.readSpecificImage(c));
line_array = zeros(nTiles, 4);
colors = zeros(nTiles, 3);
c_map= interp1(0:.5:1,[1 0 0; 1 1 0; 0 .9 .2],linspace(0,1,256))*255; %bright red to yellow to green
for win = 1:nTiles % for each subimage
    xo =  tileCoords(win,1)+floor(param.tileSize/2); % find positions of the center of the subimages
    yo =  tileCoords(win,2)+floor(param.tileSize/2);

    amp = SangS(1,win);    % extract amplitude
    ang = SangS(2,win); % extract angle of deformation

    x1 = xo-cos(ang+pi/2)*amp*param.scale; % define the extremities for the plot of deformation
    y1 = yo-sin(ang+pi/2)*amp*param.scale;
    x2 = xo+cos(ang+pi/2)*amp*param.scale;
    y2 = yo+sin(ang+pi/2)*amp*param.scale;
if isfinite([x1 y1 x2 y2])
    line_array(win,:) = [x1 y1 x2 y2];
    % TODO could do outside the loop
    q = quality(win);
    colors(win,:) = c_map(cast(q*255, 'uint8')+1,:);
end
end
RGB = insertShape(I, 'line', line_array, 'LineWidth', 15, 'Color', colors);
imwrite(RGB, [out_folder filesep 'img_' num2str(c, '%04d') '.png']);
end
