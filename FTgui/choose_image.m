function a = choose_image(path,name,t)

    if numel(extractfield(name,'name'))==1
        ibl = extractfield(name,'name');

        if numel(strfind([path,  ibl{1}],'.png'))>0
            a = imread([path,  ibl{1}]);
            a = im2double(a);
             if size(a,3)>1
                a= rgb2gray(a);
            else
            end

        else
            a = imread([path,  ibl{1}],t);
            a = im2double(a);
            if size(a,3)>1
                a= rgb2gray(a);
            else
            end

        end

    else
         a = imread([path,  name(t).name]);
         a = im2double(a);         
            if size(a,3)>1
                a= rgb2gray(a);
            else
            end
    end
end