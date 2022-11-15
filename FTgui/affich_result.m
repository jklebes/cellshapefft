function affich_result(Param,Results,col)

  
  [Udef,Vdef,Xdef,Ydef] = format_matrix_plot(Param,Results,0,0);
  Anamesi = cleanfolder(Param.name{1});

    for c = 1:Results.ci                            % for each averaged time
  
        display(['Representation of deformation map, time: ',num2str(c)]); % for the user to keep track

        Figure = figure('units','normalized');                            % Open figure
        F = imshow(ones(Param.siz(1),Param.siz(2)));        % Initialize figure with zeros, same size as the original
        hold on                                         % wait
        quiver(Xdef(:,:,c),Ydef(:,:,c),(Udef(:,:,c))/2,(Vdef(:,:,c))/2,'color','r','Autoscale','off','MaxHeadSize',0,'linewidth',1);
        quiver(Xdef(:,:,c),Ydef(:,:,c),-(Udef(:,:,c))/2,-(Vdef(:,:,c))/2,'color','r','Autoscale','off','MaxHeadSize',0,'linewidth',1);

        set(gca,'XTickLabel','','YTickLabel','');
        plot((floor(Param.siz(2)*0.9)-Param.scale):(floor(Param.siz(2)*0.9)-Param.scale+Param.scale-1),repmat(0+Param.scale,Param.scale,1),'linewidth',1,'color','black')
        text((floor(Param.siz(2)*0.9)-Param.scale),Param.scale-30,'100%','color','black')
        set(gcf,'color','white');
        img = getframe(gcf);
        imwrite(img.cdata, [Param.pathout  ['Deformation_map_averagedtime' num2str(Param.timestep) '_averaged_space_' num2str(Param.pas2) '.tif' ]],'compression','none','WriteMode','append'); % save image png
        close(Figure);

        
        Figure = figure('units','normalized');      
        name = cleanfolder([Param.name{1}]);
        a = choose_image(Param.name{1},name,c);
        a = im2double(a);
        imshow(a,'InitialMagnification',400);        % Initialize figure with zeros, same size as the original
 
        hold on                                         % wait
        quiver(Xdef(:,:,c),Ydef(:,:,c),(Udef(:,:,c))/2,(Vdef(:,:,c))/2,'color','y','Autoscale','off','MaxHeadSize',0,'linewidth',1);
        quiver(Xdef(:,:,c),Ydef(:,:,c),-(Udef(:,:,c))/2,-(Vdef(:,:,c))/2,'color','y','Autoscale','off','MaxHeadSize',0,'linewidth',1);

        set(gca,'XTickLabel','','YTickLabel','');
        plot((floor(Param.siz(2)*0.9)-Param.scale):(floor(Param.siz(2)*0.9)-Param.scale+Param.scale-1),repmat(0+Param.scale,Param.scale,1),'linewidth',1,'color','black')
        text((floor(Param.siz(2)*0.9)-Param.scale),Param.scale-30,'100%','color','black')
        set(gcf,'color','white');
        img = getframe(gcf);
        imwrite(img.cdata, [Param.pathout  ['Deformation_map_onim_averagedtime' num2str(Param.timestep) '_averaged_space_' num2str(Param.pas2) '.tif' ]],'compression','none','WriteMode','append'); % save image png
        close(Figure);


    end
end








 