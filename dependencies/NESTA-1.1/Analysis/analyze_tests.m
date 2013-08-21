% run this script once for analysis, and once for synthesis,
% changing the value of the "ANALYSIS" variable appropriately:
ANALYSIS = true;

% Stephen Becker, spring '09 srbecker@caltech.edu

if ANALYSIS, ID = 3; else ID = 2; end

date = datestr(now,29);
% date = '2009-03-25';
% date = '2009-05-08';    % enter this manually if necessary
FILE = sprintf('radar_id%02d_%s',ID,date);

load(FILE)
N = param.N;
%% Plots, showing exact and recovered on the same plot
% (didn't put these in the paper)
figure(1); clf;
% plot the exact answer
pgram = @(x) periodogram(x,[],N,param.FS);
[Pxx,FreqList] = pgram(x);
Pxx = db(Pxx,'power');
plot( FreqList, Pxx,'linewidth',5 );
ylabel('dB');xlabel('Frequency (Hz)');
% plot the recovered answer
[Pxx,FreqList] = pgram(xk);
Pxx = db(Pxx,'power');
hold on
plot( FreqList, Pxx, '-r','linewidth',2 );

legend('Exact','Recovered');

% add some lines
%     y = [-160,-120];
y = [ min( Pxx ), max( Pxx ) ];
f1 = .8345*1e9;
line( f1*[1,1], y , 'color','k');
line( param.f2*[1,1], y , 'color','g');
for ff = param.f3
    line( ff*[1,1], y , 'color','g');
end

ylim( [-95,20]);

%% Plots for the paper
pgram = @(x) periodogram(x,[],N,param.FS);
% -- the exact signal
[Pxx,FreqList] = pgram(x);Pxx = db(Pxx,'power');
tit = 'Spectrum estimate of the Nyquist-sampled signal';
% -- recovered
[Pxx,FreqList] = pgram(xk);Pxx = db(Pxx,'power');
if ANALYSIS
    tit = 'Spectrum estimate of the signal recovered via analysis';
else
    tit = 'Spectrum estimate of the signal recovered via synthesis';
end


fsize = 12;
figure(1); clf;
plot( FreqList, Pxx,'linewidth',2 );
ylabel('dB','fontsize',fsize);xlabel('Frequency (Hz)','fontsize',fsize);
hold on
y = [ min( Pxx ), max( Pxx ) ];
f1 = .8345*1e9;
line( f1*[1,1], y , 'color','k','linestyle','-.');
line( param.f2*[1,1], y , 'color','k','linestyle','--');
for ff = param.f3
    line( ff*[1,1], y , 'color','k','linestyle',':');
end
title(tit,...
    'fontsize',fsize);
ylim( [-105,20]);
legend('Spectrum','Exact frequency of signal 1',...
    'Exact frequency of Doppler Radar',...
    'Exact frequencies of frequency hopping pulse' );
