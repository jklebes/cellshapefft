function Results = def_analysis(Param,Results)

     for c = 1:Results.ci                                   % for each averaged time
        display(['Computation of cell deformation, time: ',num2str(c)]); % for the user to keep track
       
        %% We will group calculus
        SFT = Results.im_regav.spect(c,:,:,:);
        SFT = squeeze(SFT);
        ST = sum(sum(SFT,2),3);
        %indexx = find();
        %% Inertia matrix 
        
        Msft = inertia_matp_sigma(Param,SFT((ST~=0 | isnan(ST)~=0),:,:));
        
        %% Computation of cell deformation 
        
        [Si,angSi] = S_angS(Msft);
        AlS = zeros(Results.regl,1);
        AlS((ST~=0)) = Si;
        Results.im_regav.S(c,:,:) = reshape(AlS,1,Results.numY,Results.numX);
        AlangS = Results.im_regav.angS(c,:);
        AlangS(1,(ST~=0)) = angSi;
        Results.im_regav.angS(c,:,:) = reshape(AlangS,1,Results.numY,Results.numX);
     end

    if Param.register == 0 && Param.regsize == 0
           Results.im_regav.spect = [];       % delete the spectrum to keep space
    else
    end

end