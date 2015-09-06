function disp(A)

%@SSMODEL/DISP Display of state space model.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if isempty(A.name), name = '(unnamed)';
else name = A.name; end
fprintf(1, '\t%s\n\t%s\n\n', name, repmat('=', 1, length(name)));
fprintf(1, '\tp = %d\n', A.p);
fprintf(1, '\tm = %d\n', size(A.Z, 2));
fprintf(1, '\tr = %d\n', size(A.R, 2));
fprintf(1, '\tn = %d\n\n\n', A.n);

if ~isempty(A.Hinfo), info = [{A.Hinfo} A.info];
else info = A.info; end
if ~isempty(info)
    fprintf(1, '\tComponents (%d):\n\t==================\n\n', length(info));
    for i = 1 : length(info)
        fprintf(1, '\t[%s]\n', info{i}.type);
        for j = fieldnames(info{i})'
            field   = j{1};
            if ~strcmp(field, 'type')
                %%%%%%%% TODO: Add printing support for cell arrays
                fprintf(1, '\t\t%-16s', field);
                if isnumeric(info{i}.(field)), fprintf(1, '= %s', mat2str(info{i}.(field)));
                elseif ischar(info{i}.(field)), fprintf(1, '= %s', info{i}.(field));
                end
                fprintf(1, '\n');
            end
        end
    end
    fprintf(1, '\n\n');
end

printmatname('H', size(A.H, 1), size(A.H, 2), A.H.n); disp(A.H, true);
printmatname('Z', size(A.Z, 1), size(A.Z, 2), A.Z.n); disp(A.Z, true);
printmatname('T', size(A.T, 1), size(A.T, 2), A.T.n); disp(A.T, true);
printmatname('R', size(A.R, 1), size(A.R, 2), A.R.n); disp(A.R, true);
printmatname('Q', size(A.Q, 1), size(A.Q, 2), A.Q.n); disp(A.Q, true);
if ~issta(A.c) || ~isconst(A.c) || any(A.c.mat ~= 0), printmatname('c', size(A.c, 1), size(A.c, 2), A.c.n); disp(A.c, true); end
if ~isconst(A.a1) || any(A.a1.mat ~= 0), printmatname('a1', size(A.a1, 1), size(A.a1, 2), A.a1.n); disp(A.a1, true); end
printmatname('P1', size(A.P1, 1), size(A.P1, 2), A.P1.n); disp(A.P1, true);

fprintf(1, '\tAdjacency matrix (%d functions):\n', length(A.func));
fprintf(1, '\t-----------------------------------------------------------------------------\n');
fprintf(1, '\tFunction                H Z T R Q c a1 P1 Hd Zd Td Rd Qd cd Hng Qng Znl Tnl\n');
for i = 1 : length(A.func)
    fline   = '\t%-22s  %1d %1d %1d %1d %1d %1d %1d  %1d  %1d  %1d  %1d  %1d  %1d  %1d  %1d   %1d   %1d   %1d\n';
    fprintf(1, fline, func2str(A.func{i}), A.A(i, :));
end
fprintf(1, '\n\n');

fprintf(1, '\tParameters (%d)\n', A.psi.w);
fprintf(1, '\t------------------------------------\n');
fprintf(1, '\tName                Value\n');
psi = get(A.psi);
for i = 1 : A.psi.w, fprintf(1, '\t%-20s%-16.4g\n', A.psi.name{i}, psi(i)); end
fprintf(1, '\n');

function printmatname(name, p, m, n)
if n == 1, fprintf(1, '\t%s matrix (%d, %d)\n\t-------------------\n\n', name, p, m);
else fprintf(1, '\t%s matrix (%d, %d, %d)\n\t----------------------\n\n', name, p, m, n);
end

