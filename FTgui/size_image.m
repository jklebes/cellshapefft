function siz = size_image(path,name)   
   infos = imfinfo([path,  name]); 
   siz = [infos.Height,infos.Width];
end