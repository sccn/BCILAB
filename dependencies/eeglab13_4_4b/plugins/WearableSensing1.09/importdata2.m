function S = importdata2(file_str)
%IMPORTDATA2   Read a CSV file and extract its information.
%   S = importdata2(file_str) reads the csv file with the path file_str and
%   returns parameter headers text as a cell column vector S.textdata, column headers
%   as a cell row vector S.colheaders and data matrix as S.data.
%
%Author: Umut Orhan
%Copyright: QUASAR 2011
%Modified by Steven Pillen
%Copyright: Wearable Sensing 2015

fid=fopen(file_str,'r');
line_cnt=1;
istextline=1;
comment = 0;
S.data = [];
while(istextline)
    tline = fgets(fid);
    istextline=any(isletter(tline));
    if any(strfind(tline,',Comment'))
        comment = 1;
        %tline = strsplit(tline, ',');
        %S.data = [S.data; {tline}];
        if comment == 1
            
            if strcmp(tline(1:4), 'Time') % Time is the indicator of the montage
                %S.colheaders = strsplit(tline, ',');
                tline = strrep(tline, ',', ' ');
                tline = textscan(tline, '%s');
                S.colheaders = tline{1}';
                S.colheaders = S.colheaders(1:end-1);
                %S.colheaders = S.colheaders';
                %headers{line_cnt}=tline;
                line_cnt=line_cnt+1;
                break
            end
        end
    end
    
    
    headers{line_cnt}=tline;
    line_cnt=line_cnt+1;
    
end

if comment == 1
    
    %buffer = tline;
    %     buffer = fgets(fid);
    %     while length(buffer) > 1
    %         tline = buffer;
    %         if class(tline) == 'char'
    %             tline = strsplit(tline, ',');
    %             tline = cellfun(@str2double, tline(1:end-1));
    %         end
    %         S.data = [S.data, tline];
    %         buffer = fgets(fid);
    %     end
    
    Filedata = fread(fid); %loads the binary data into the workspace
    Filedata = char(Filedata)'; %translates the binary into chars
    %Filedata = strsplit(Filedata, char(10)); %splits the chars into rows
   
    Filedata = textscan(Filedata, '%s', 'Delimiter',char(10));
    Filedata = Filedata{1};
    S.data = zeros(length(S.colheaders),length(Filedata));
    for line = 1:(length(Filedata)-1);
        %neoline =  strsplit(Filedata{line}, ',')';
        
        currentline = textscan(Filedata{line}, '%s', 'Delimiter', ',');
        S.data(:,line) = str2double(currentline{1}(1:length(S.colheaders)));
   %     S.data = [S.data, neoline];
    end
    S.data = S.data';
    S.textdata=headers';
    fclose(fid);
    %S.data = S.data'
else
    fclose(fid);
    S.data=csvread(file_str,line_cnt-2);
    
    
    
    
    
    
    colheader_locs=findstr(headers{line_cnt-2},',');
    
    if(length(colheader_locs)+1==size(S.data,2))
        S.colheaders=cell(1,length(colheader_locs)+1);
        extended_colheader_locs=[0 colheader_locs length(headers{line_cnt-2})-1];
        for(ii=1:length(colheader_locs)+1)
            S.colheaders{ii}=headers{line_cnt-2}(extended_colheader_locs(ii)+1:extended_colheader_locs(ii+1)-1);
        end
        header_end=line_cnt-3;
    else
        header_end=line_cnt-2;
    end
    
    S.textdata=headers(1:header_end)';
end
