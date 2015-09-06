clear('pi');
y           = load('data/M3_21.dat')';
airline     = ssm_airline;
airline     = estimate(y, airline, 0.1);
theta       = -airline.param(1);
Theta12     = (-airline.param(2))^(1/12);
zetavar     = airline.param(3);

param0  = [theta+Theta12 -theta*Theta12 Theta12 Theta12 zetavar];
result4 = zeros(31, 6);
n       = 1;
% 4-5-1
for i = 1 : 6
    genair          = ssm_genair(4, 5, i);
    [genair logL]   = estimate(y, genair, param0, [], 'fmin', 'simplex', 'disp', 'off');
    result4(n, :)   = [roots([1 -genair.param(1:2)])' genair.param(3:4) reallog(genair.param(5))/2 -2*logL-13*reallog(2*pi)+2*4];
    n               = n + 1;
end
% 4-4-2
for i = 1 : 5
    for j = i+1 : 6
        genair          = ssm_genair(4, 4, [i j]);
        [genair logL]   = estimate(y, genair, param0, [], 'fmin', 'simplex', 'disp', 'off');
        result4(n, :)   = [roots([1 -genair.param(1:2)])' genair.param(3:4) reallog(genair.param(5))/2 -2*logL-13*reallog(2*pi)+2*4];
        n               = n + 1;
    end
end
% 4-3-3
for i = 2 : 5
    for j = i+1 : 6
        genair          = ssm_genair(4, 3, [1 i j]);
        [genair logL]   = estimate(y, genair, param0, [], 'fmin', 'simplex', 'disp', 'off');
        result4(n, :)   = [roots([1 -genair.param(1:2)])' genair.param(3:4) reallog(genair.param(5))/2 -2*logL-13*reallog(2*pi)+2*4];
        n               = n + 1;
    end
end
fline   = '%.5f\t\t%.5f\t\t%.5f\t\t%.5f\t\t%.5f\t\t%.5f\n';
fprintf(1, fline, result4');

param0  = [theta Theta12 Theta12 zetavar];
result3 = zeros(41, 5);
n       = 1;
% 3-5-1
for i = 1 : 6
    genair          = ssm_genair(3, 5, i);
    [genair logL]   = estimate(y, genair, param0, [], 'fmin', 'simplex', 'disp', 'off');
    result3(n, :)   = [genair.param(1:3) reallog(genair.param(4))/2 -2*logL-13*reallog(2*pi)+2*3];
    n               = n + 1;
end
% 3-4-2
for i = 1 : 5
    for j = i+1 : 6
        genair          = ssm_genair(3, 4, [i j]);
        [genair logL]   = estimate(y, genair, param0, [], 'fmin', 'simplex', 'disp', 'off');
        result3(n, :)   = [genair.param(1:3) reallog(genair.param(4))/2 -2*logL-13*reallog(2*pi)+2*3];
        n               = n + 1;
    end
end
% 3-3-3
for i = 1 : 4
    for j = i+1 : 5
        for k = j+1 : 6
            genair          = ssm_genair(3, 3, [i j k]);
            [genair logL]   = estimate(y, genair, param0, [], 'fmin', 'simplex', 'disp', 'off');
            result3(n, :)   = [genair.param(1:3) reallog(genair.param(4))/2 -2*logL-13*reallog(2*pi)+2*3];
            n               = n + 1;
        end
    end
end
fline   = '%.5f\t\t%.5f\t\t%.5f\t\t%.5f\t\t%.5f\n';
fprintf(1, fline, result3');
