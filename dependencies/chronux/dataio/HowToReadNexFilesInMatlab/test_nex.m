[nvar, names, types] = nex_info('test.nex');

[n, ts] = nex_ts('test.nex', 'neuron1');

[adfreq, n, ts, fn, d] = nex_cont('test.nex', 'continuous1');

[adfreq, n, ts, nf, w] = nex_wf('test.nex', 'wave1');

[n, nm, nl, ts, names, m] = nex_marker('test.nex', 'marker1');