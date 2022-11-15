function [Udef,Vdef,Xdef,Ydef] = format_matrix_plot(Param,Results,Row,Col)
    
    Ys = unique((Results.Posi(:,2,1)));
    Reg_xy = zeros(numel(Ys),numel(unique((Results.Posi(:,1,1)))));

    % Le tableau Reg_xy contient les indices de chaque régions
    % correspondantes aux positions x,y,x en colonne
    
    for ii = 1:numel(Ys)
        ind = find(Results.Posi(:,2,1) == Ys(ii));
        Reg_xy(ii,:) = ind;
    end
    
%%
    Udef = zeros(size(Reg_xy,1),size(Reg_xy,1),Results.ci);
    Vdef = zeros(size(Reg_xy,1),size(Reg_xy,1),Results.ci);
    Xdef = zeros(size(Reg_xy,1),size(Reg_xy,1),Results.ci);
    Ydef = zeros(size(Reg_xy,1),size(Reg_xy,1),Results.ci);
    %%
    for ml = 1:Results.ci
        for win =1:Results.regl % for each subimage
            [alk,blk] = find(Reg_xy == win);
          xo =  Results.Posi(win,1)+Param.pas2/2-1; % find positions of the center of the subimages
          yo =  Results.Posi(win,2)+Param.pas2/2-1;
          [I,J] = ind2sub([Results.numY,Results.numX],win);

            if sum(yo == Row)>0 && sum(Col(yo==Row)==xo)>0
                amp= NaN;
                ang=NaN;
            else
                amp = Results.im_regav.S(ml,I,J);    % extract amplitude 
                ang = Results.im_regav.angS(ml,I,J); % extract angle of deformation

                x1 = cos(ang+pi/2)*amp*Param.scale; % define the extremities for the plot of deformation
                y1 = sin(ang+pi/2)*amp*Param.scale;
                
                 if round(x1) <0
                     Vdef(alk,blk,ml) = -(round(y1)); % define the extremities for the plot of deformation
                     Udef(alk,blk,ml) = -(round(x1));
                 else
                    Vdef(alk,blk,ml) = (round(y1)); % define the extremities for the plot of deformation
                    Udef(alk,blk,ml) = (round(x1));
                 end
                Xdef(alk,blk,ml) = xo;
                Ydef(alk,blk,ml) = yo;
                
            end

        end
    end

end