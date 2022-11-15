function [Results] = spectrum_analysis(Param,Results)
    

    h =  waitbar(0,'Please wait, computation of FT');
    
%     % Number of slice we will have
%     
%     Len = numel(Param.tsart:Param.tleng-Param.timestep+1);
%     
%     % Initialization of Results
%     
%     Results.im_regav = struct(repmat(struct('c',struct(repmat(struct('spect',zeros(Param.pas2),'S',[]...
%         ,'angS',[],'a',[],'b',[],'phi',[]),max(Results.regl),1))),(Param.tleng-Param.tsart+1)-Param.timestep+1,1));
%     
%     Results.im_regav = struct(repmat(struct('c',struct(repmat(struct('spect',zeros(Param.pas2),'S',[]...
%     ,'angS',[],'a',[],'b',[],'phi',[]),max(Results.regl),1))),Len,1));
%     struct('spect',zeros(Len,Param.pas2),'S',zeros(Len,Param.pas2),'angS',[],'a',[],'b',[],'phi',[])
    
 %   Results.ci =numel(Param.tsart:Param.tleng-Param.timestep+1);                             % ci will be the final number of slices
    Results.im_regav.spect = avfft_forall(path,wins,Param,Positions);
    close(h);
end