function FTs = avfft(path,wins,Param,Positions)
% Initialization
   nbel = size(path,1);% Nb of different experiments

   S = size_image(Param.name{1},cleanfolder(path{1}));
   as = zeros(nbel,S(1),S(2));
   asubs = zeros(size(wins,2),nbel,Param.pas2,Param.pas2);
   FTs = zeros(size(wins,2),Param.pas2,Param.pas2);
   
    x = Positions(wins,1);
    y = Positions(wins,2);
% Fill in small subimages + fft computation
for t = Param.tsart:Param.tsart+Param.timestep-1
    for ii = 1: nbel    
        
        % Clean folders to extract names to treat 
        name = cleanfolder([Param.name{ii}]);
        as(ii,:,:) = choose_image(path{ii},name,t);
        for kk = 1:size(wins,2)                        %  Nb of subimages
            asubs(kk,ii,:,:) = as(ii,y(kk):y(kk)+Param.pas2-1,x(kk):x(kk)+Param.pas2-1);
            FTs(kk,:,:) = FTs(kk,:,:)+reshape(fft_adir(reshape(asubs(kk,ii,:,:),Param.pas2,Param.pas2)),1,Param.pas2,Param.pas2);
        end
    end
end

FTs = FTs/(nbel+Param.timestep);
end