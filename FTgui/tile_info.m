function [Posi,regl,numX,numY] = tile_info(Param)  
%returns 
%Posi a list of locations of subimages, 
%reg1 number of possible submiages
%numX, numY number of subimages in x and y directions

%renamed
%rec -> overlap
%pas2 -> subwindow_size
spacing = (1-Param.overlap)*Param.subwindow_size;  % number of shared pixels between
                                                   % subimages
%already been done, why not take from Param?
%Anamesi = cleanfolder(Param.name{1}); %and not used subsequently
%Param.siz = size_image(Param.name{1},Anamesi(1).name);

% The positions needed
% possible starting positions of subimages

AllX = 1:spacing:(Param.siz(2)- Param.subwindow_size);
AllY = 1:spacing:(Param.siz(1)- Param.subwindow_size);
% AllX = 1:rec:(Param.siz(2)- Param.pas2+1);
% AllY = 1:rec:(Param.siz(1)- Param.pas2+1);
%AllX = 1:rec:(Param.siz(2)- rec+1);
%AllY = 1:rec:(Param.siz(1)-rec+1);


% The number of positions in X and Y directions

numX = numel(AllX);
numY = numel(AllY);

for ii = 1:numX
   Posi((ii-1)*numY+1:ii*numY,:) = [floor(AllX(ii))*ones(numY,1),floor(AllY')]; 
end

regl = size(Posi,1);               % register number of subimages at each time

end