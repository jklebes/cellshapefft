function  affich_size(Param,Results,col)


    for c = 1:Results.ci                            % for each averaged time
    display(['Representation of size map, time: ',num2str(c)]); % for the user to keep track
    Figure = figure(20);                            % Open figure
    imshow(ones(Param.siz(1),Param.siz(2)));        % Initialize figure with zeros, same size as the original
    hold on                                         % wait
         for win =1:Results.regl % for each subimage
            xo =  Results.Posi(win,1)+Param.pas2/2; % find positions of the center of the subimages
            yo =  Results.Posi(win,2)+Param.pas2/2;
            
            a = Results.im_regav.a(c,win);    % extract half major axis
            b = Results.im_regav.b(c,win); % extract half minor axis
            phi = Results.im_regav.phi(c,win)+pi/2; % extract angle

           
            f_ellipseplotax(a,b,xo,yo,phi,col,0.5,Figure.CurrentAxes);   
        end
    set(gca,'XTickLabel','','YTickLabel','');
    img = getframe(gcf);
    imwrite(img.cdata, [Param.pathout  ['Size_map_averagetime' num2str(Param.timestep) '_averagespace_' num2str(Param.pas2) '.tif']],'compression','none','WriteMode','append'); % save image png
    savefig(Figure,[Param.pathout   ['Size_map_averagetime' num2str(Param.timestep) '_averagespace_' num2str(Param.pas2) '_'] num2str(c) '.fig'])     % save image fig
    close(Figure);

    
    
            
    Figure = figure('units','normalized');      
    name = cleanfolder([Param.name{1}]);
    a = choose_image(Param.name{1},name,c);
    a = im2double(a);
    imshow(a,'InitialMagnification',400);        % Initialize figure with zeros, same size as the original

    hold on                                         % wait
     for win =1:Results.regl % for each subimage
        xo =  Results.Posi(win,1)+Param.pas2/2; % find positions of the center of the subimages
        yo =  Results.Posi(win,2)+Param.pas2/2;

        a = Results.im_regav.a(c,win);    % extract half major axis
        b = Results.im_regav.b(c,win); % extract half minor axis
        phi = Results.im_regav.phi(c,win)+pi/2; % extract angle


        f_ellipseplotax(a,b,xo,yo,phi,col,0.5,Figure.CurrentAxes);   
    end
    set(gca,'XTickLabel','','YTickLabel','');
    img = getframe(gcf);
    imwrite(img.cdata, [Param.pathout  ['Size_map_onim_averagetime' num2str(Param.timestep) '_averagespace_' num2str(Param.pas2) '.tif']],'compression','none','WriteMode','append'); % save image png
    close(Figure);


    
    
    
    end
end
