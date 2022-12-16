function cleanednames = cleanfolder(folder)
    int1 = dir(folder);
    int1([int1.isdir]) = [];
    A = strfind({int1.name},'.DS_Store');
    A2 = strfind({int1.name},'.db');
    cleanednames = int1(cellfun('isempty',A2));