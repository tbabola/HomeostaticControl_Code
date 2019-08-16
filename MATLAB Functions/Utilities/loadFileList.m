function files = loadFileList(dname)
    list = dir(dname);
    for i = 1:size(list,1)
        files{i,1} = [list(i).folder '\' list(i).name];
    end 

end