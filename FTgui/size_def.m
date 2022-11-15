function [A,B,Phi,Yu,Xu]= size_def(Param,SFT) 
%% Definition de M 

    A = zeros(size(SFT,1),1);
    B = zeros(size(SFT,1),1);
    Yu = zeros(size(SFT,1),Param.nbpoints);
    Xu = zeros(size(SFT,1),Param.nbpoints);
    Phi = zeros(size(SFT,1),1);
    
%% Filter definition
    filtre = gaussianFilter(ceil(2*Param.sigma),Param.sigma);
%% Find center of the image
    ic = Param.pas2*1/2+1;
    jc = Param.pas2*1/2+1;              

%%  Define pixels to remove each time  
    [I,J] = ndgrid(1:Param.pas2,1:Param.pas2);
    M = double((I-ic).^2+(J-jc).^2<=(Param.cut)^2); 


%% Loop on all regions

    for im = 1:size(SFT,1)
        ff = squeeze(SFT(im,:,:));
        th_spec = filter2(filtre,ff);      % th_spec gaussian filter of the spectrum
        th_spec(M>0) = 0;   
        
        % Search local maximas
        [yu,xu] = localMaximum_h(th_spec,1,0,Param.nbpoints);   
        if numel(unique(yu)) < 3 || numel(unique(xu)) < 3
            ai = 0;
            bi  = 0;
            phi = 0;
        else
        [~,~,ai,bi,phi,~]=ellipsefit(xu,yu);  
        end
        
        A(im) = Param.pas2/2*1/bi;  % On divise par deux car on obtient le diam?tre de la 
    % cellule car distance entre les cot?s qui ressort en Fourier
        B(im) = Param.pas2/2*1/ai;
        Xu(im,:) =xu';
        Yu(im,:) =yu';
        Phi(im) = phi;
    end

  
    
end