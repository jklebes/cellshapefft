function Results = def_analysis_ellips(Results)

     for c = 1:Results.ci                                   % for each averaged time
        display(['Computation of cell deformation, time: ',num2str(c)]); % for the user to keep track
       Results.im_regav.S(c,:,:) = log(Results.im_regav.a(c,:,:)./Results.im_regav.b(c,:,:));
       Results.im_regav.angS(c,:,:) =  Results.im_regav.phi(c,:,:);
      
     end

end