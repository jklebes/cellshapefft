function tile_slopes=basis_from_heightmap(folder, voxel, n_tiles, tilecoords, tiledims, n_timepoints)
% read heightmap files, for each tile at each timepoint
% determine a plane best approximating
numdims = size(tiledims);
if numdims(2)>1
    tiledim_x=tiledims(1);
    tiledim_y=tiledims(2);
else
    tiledim_x = tiledims;
    tiledim_y = tiledims;
end
% the tile heightmap, return its defining vectors
file_names = dir(folder); %files in dir, hopeuflly all heightmap files
file_names = file_names(3:end); %remove files '.', '..'
%empty array
tile_slopes=zeros(n_timepoints, n_tiles,2);
%for each timepoint
for c =1:n_timepoints
    %read heightmap to array
    heightmap_file = [folder, file_names(c).name];
    heightmap = read(Tiff(heightmap_file,"r"));  
    %maybe smooth
    heightmap = imgaussfilt(double(heightmap),3);
    %for each tile
    for tile =1:n_tiles
        tilecoord = tilecoords(tile,:);
        tile_heightmap = heightmap(tilecoord(2):tilecoord(2)+tiledim_x-1, tilecoord(1):tilecoord(1)+tiledim_y-1);
        %x,y,z coords of points in heightmap as 1D lists
        xs = reshape(repmat(tilecoord(2):tilecoord(2)+tiledim_x-1, tiledim_y, 1),[],1);
        ys =  reshape(repmat((tilecoord(1):tilecoord(1)+tiledim_y-1)', 1, tiledim_x),[],1);
        zs = reshape(tile_heightmap, [], 1);
        %fit a plane to heightmap
        planefit = fit([xs,ys,],zs, 'poly11'); %returns an sfit fit object
        %holding coeffs in f=p00+p10*x+p01*y, save second two coeffs
        tile_slopes(c,tile,:)=[planefit.p10, planefit.p01];
    end
end