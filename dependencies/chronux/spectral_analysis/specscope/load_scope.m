function load_scope(files,refresh_rate,sample_frequency)
close all
shot_size=round(sample_frequency/refresh_rate);
ring=zeros(shot_size*20,1)';

fig = figure('Position',[500,500,800,600],...
        'NumberTitle','off',...
        'Name','Scope',...
        'doublebuffer','on',...
        'HandleVisibility','on',...
	'Renderer', 'openGL');

plot_ref=plot(zeros(10,1));


for i=1:length(files)
    %files{i};
    fid=fopen(files{i},'r');
    fseek(fid,0,1);
    e_o_f=ftell(fid);
    fseek(fid,0,-1);
    
    
    while (ftell(fid)<e_o_f)
        tic
        [data,c]=fread(fid,shot_size,'short');
        data=(data/10000)';
        %data=rand(1,length(data));
        ring=[ring data];
        ring(1:length(data))=[];
        start_display(ring,plot_ref);
        stop_time=toc;
        while stop_time<(1/refresh_rate)
            stop_time=toc;
            
        end
    end
    fclose(fid)
    
end
    
    
