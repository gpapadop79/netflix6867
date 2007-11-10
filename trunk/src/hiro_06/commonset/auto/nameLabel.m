fid=fopen('movie_titles.txt');
names={};

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end    
    [id,year,name]=strread(tline,'%d%d%s',3,'delimiter',',');
    names{id}=name;
end

fclose(fid);