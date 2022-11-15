function Msft = inertia_matp_sigma(Param,SFT)
%% Definition de M 

    Msft = zeros(size(SFT,1),2,2);
    
%% Filter definition
    filtre = gaussianFilter(ceil(2*Param.sigma),Param.sigma);
%% Find center of the image
    ic = Param.pas2*1/2+1;
    jc = Param.pas2*1/2+1;              

%%  Define pixels to remove each time  
    [I,J] = ndgrid(1:Param.pas2,1:Param.pas2);
    M = double((I-ic).^2+(J-jc).^2<=(Param.cut)^2); 

%% Boolean matrix used after and strel
     se=strel('disk',Param.strel);                % smooth with a strel


%% Loop on all regions

    for im = 1:size(SFT,1)
        ff = squeeze(SFT(im,:,:));
        th_spec = filter2(filtre,ff);      % th_spec gaussian filter of the spectrum
        th_spec(M>0) = 0;                  % Remove pixels in the center
        reshaped = reshape(th_spec,[Param.pas2*Param.pas2 ,1]);  % reshape the spectrum to order
        [~,index]= sort(reshaped,'descend');                     % sort in intensity and keep indexes
        nbel = Param.propor * Param.pas2*Param.pas2;             % nbel number of points determined by Param.propor
   
        if nbel > Param.pas2*Param.pas2            % if nbel is to large, take the whole picture
            nbel = Param.pas2*Param.pas2;
        else
        end
        
        th_spec(index(1:floor(nbel))) = 1;           % put to one the points kept by proportion
        th_spec(index(floor(nbel)+1:end)) = 0;       % put to zero the others
        th_spec(M>0) = 1;                            % get back the pixels in the center
        
        B = false(Param.pas2);     

        B(th_spec>0)=1;                              % put rights where th_spec = 1
        th_spec=imclose(B,se);                       % close image
        
        %% Inertia matrix computation
        
        norm = sum(sum((th_spec).^2));               % computes norm 
        [Yg,Xg]=meshgrid(1:(Param.pas2),1:(Param.pas2)); % create meshgrid of correct size
        xc = Param.pas2/2;                               % position of the center 
        yc = Param.pas2/2;
        Xo = Xg-xc;                           
        Yo = Yg-yc;
        Y = (Yg-yc).^2; % Lignes
        X = (Xg-xc).^2; % Colonnes
        YY = sum(sum((th_spec).^2.*Y))/norm;         % compute the 3 different values in the inertia matrix
        XX = sum(sum((th_spec).^2.*X))/norm;
        XY = sum(sum((th_spec).^2.*Xo.*Yo))/norm;
        
        %% Registering of inertia matrix
        
        Msft(im,:,:) = [XX,XY;...
            XY,YY];                           % keep matrix of inertia



    end


     
 
end


 