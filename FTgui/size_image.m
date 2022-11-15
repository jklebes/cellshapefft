function siz = size_image(path,name)   
   infos = imfinfo([path,  name(1).name]); 
   infos = infos(1);

   siz = [infos.Height,infos.Width];
end