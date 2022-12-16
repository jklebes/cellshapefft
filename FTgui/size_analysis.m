function Results = size_analysis(Param,Results)
     for c = 1:Results.ci                                   % for each averaged time
        display(['Computation of cell size, time: ',num2str(c)]); % for the user to keep track

        SFT = Results.im_regav.spect(c,:,:,:);
        SFT = squeeze(SFT);
        ST = sum(sum(SFT,2),3);

        %%
        
        [A,B,Phi,~,~]= size_def(Param,SFT((ST~=0 | isnan(ST)~=0),:,:));
        
        %% Registering in correct format
        
        interA = zeros(Results.regl,1);
        interB = zeros(Results.regl,1);
        interphi = zeros(Results.regl,1);
        interA((ST~=0 | isnan(ST)~=0)) = A;
        interB((ST~=0 | isnan(ST)~=0)) = B;
        interphi((ST~=0 | isnan(ST)~=0)) = Phi;
        
        
        Results.im_regav.a(c,:,:) = reshape(interA,1,Results.numY,Results.numX);
        Results.im_regav.b(c,:,:) = reshape(interB,1,Results.numY,Results.numX);
        Results.im_regav.phi(c,:,:) = reshape(interphi,1,Results.numY,Results.numX);
        
     end
end
