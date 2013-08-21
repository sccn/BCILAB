function rtf(plot_frq,flag_save)


close all
evalin('base','stop=0;');



%=========SET THE BASIC FIGURE=================
fig = figure('Position',[500,500,800,600],...
        'NumberTitle','off',...
        'Name','Scope',...
        'doublebuffer','on',...
        'HandleVisibility','on',...
        'KeyPressFcn', @keypress, ...
        'Renderer', 'openGL');
%=============================================


%=============OPEN THE DEVICE FOR RECORD======
sample_frequency = 44100;
samples_per_frame = 1024;
%plot_frq=10;
record_time=600;
samples_to_acquire = record_time * sample_frequency;


%PREPARE    THE    DEVICE
ai = analoginput('winsound');   
chan = addchannel( ai, 1 );
set( ai, 'SampleRate', sample_frequency )
set( ai, 'SamplesPerTrigger', samples_to_acquire )
set(ai, 'StopFcn', @stop_dev)

sample_frequency = get( ai, 'SampleRate' );


%SETTING CALL BACK FUNCTIONS:
%The first for capture   the 
%second     for      display
set(ai, 'SamplesAcquiredFcnCount',samples_per_frame);
set(ai, 'SamplesAcquiredFcn',@flag);
set(ai, 'TimerPeriod',(1/plot_frq));
set(ai, 'TimerFcn',@disply);


%=============SAVE THE CONFIGURATION======
plot_ref=plot(zeros(10,1));
fid=-1;


%SAVE THE CURRENT PARAMETERS:
name_of_file=sprintf('%s-%d','real-anal',(round(sample_frequency/samples_per_frame)));
remark={1,...
        zeros(samples_per_frame*20,1)',...
        0,...
        plot_ref,...
        plot_frq,...
        cputime,...
        flag_save,...
        -1,...
        name_of_file
        };

set(ai, 'UserData',remark)




%=============START TO RECORD================
fprintf ('To stop the program set <stop=1> or press q in the figure window\n');
start (ai)





%=============================================
%=========THE MAIN PROGRAM====================
%=============================================
%                   *
%                  ***
%                 *****
%                  ***
%                  ***
%                  ***
%                  ***
%                  ***
%                 *****
%                  ***
%                   *
%=============================================
%==========CALLBACK FUNCTIONS=================
%=============================================



%=========Keypress callback===========
function keypress(src, e)
  keypressed=get(gcf,'CurrentCharacter');

  % ignore raw control, shift, alt keys
  if keypressed
    % Quit
    if strcmp( keypressed, 'q')
        evalin('base','stop=1;');
    end
  end
return


%============FLAG FUNCTION===================
%This function activated when we capture
%certain    amount  of          samples
function flag(obj,event)

  % CHECK FOR STOP SIGNAL
  if  evalin('base','stop')
      stop(obj)
  end
 
% GET THE OLD DATA
 remark=get(obj,'UserData');
 flag_write=remark{1};  %Do I have to 
 buffer=remark{3};      %What is the old picture
 flag_save=remark{7};   %Are we in saving mode?
 fid=remark{8};         %What file descriptor to save
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %IN CASE - DELETE/SAVE THE OLD DATA
 
 if flag_write>20
     
     
     % IN CASE WE HAVE TO SAVE - CLOSE THE OLD FILE AND MAKE A NEW
     if flag_save>0
        fclose(fid);
        name_of_data=sprintf('%s-%d.dat','dat',(round(cputime*1000)));
        fid=fopen(name_of_data,'w');
    end
    
      %DELETE OLD DATA
     flag_write=1;
     buffer=[];
     remark{1}=flag_write; % SET THE POSITION OF THE READING SHIFT
 end
  
   
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % TAKE THE NEW DATA
 
 samples_per_frame=get(obj,'SamplesAcquiredFcnCount');
 data=(getdata(obj,samples_per_frame))';    
 
 % IN CASE - WRITE THE DATA
 if flag_save>0
     if fid==-1
         name_of_data=sprintf('%s-%d.dat','dat',(round(cputime*1000)));
         fid=fopen(name_of_data,'w');
     end
     fwrite(fid,(data*10000),'short');
     remark{8}=fid;
 end
 
 
 % Add to buffer
 buffer=[buffer data];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 remark{3}=buffer;
 set(obj,'UserData',remark); 
 
  
return


function stop_dev(obj,event)

     remark=get(obj,'UserData');
     if (remark{8}>0)                 %FID>0 == There is open file
         fclose (remark{8});
     end
     close all
     
     fprintf('\n\nThanks for using Erlich Real-Time scope\n');
     save (remark{9},'remark');
     
     delete(obj)
     clear obj
return











function disply(obj,event)


 

  sample_frequency=get(obj,'SampleRate');
  remark=get(obj,'UserData');
  refresh_frq=remark{5};
  read_shift=remark{1};
  
  ring=remark{2};
  buffer=remark{3};
  
 
  end_shift=min((read_shift+round(sample_frequency/refresh_frq)),length(buffer));
  new_data=buffer(read_shift:end_shift);
  ring=[ring new_data];
  ring(1:length(new_data))=[];
  remark{1}=end_shift;
  remark{2}=ring;
  
  
  start_display(ring,remark{4})
  
  set(obj,'UserData',remark);
return
 
 
