function [parameters,states] = bci_Construct

quicklog('bcilab-log.txt','passing on to bci_ConstructReal');

[parameters,states] = bci_ConstructReal;

quicklog('bcilab-log.txt','finished successfully');
