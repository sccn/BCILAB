function vis_benchmark(timing)
% display results of benchmarking (timing)

fprintf('\n\n------------------------------\n');
fnames = fieldnames(timing);
ttot = 0; % total sum of timing
for k=1:length(fnames);
    tsec = timing.(fnames{k});
    if isnan(tsec), continue; end;
    fprintf('%s\t: %0.5g\n',fnames{k},tsec);
    ttot = ttot+tsec;
end
fprintf('Total\t: %0.5g\n',ttot);
fprintf('------------------------------\n\n');