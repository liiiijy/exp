%
% A file that makes a surface plot of many included problems.
%

% This file is part of the 'Expint'-package, 
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.5 $ $Date: 2005/10/22 02:50:14 $


scheme = 'lawson4';

pr{1} = nls;
[t, U, UF] = expglm(pr{1}, [0 10], 0.01, scheme, [0:0.1:10]);
surfplot(UF, t, pr{1});


pr{2} = kursiv;
[t, U, UF] = expglm(pr{2}, [0 100], 0.1, scheme, [0:0.2:100]);
surfplot(UF, t, pr{2});

pr{3} = burgers;
[t, U, UF] = expglm(pr{3}, [0 20], 0.01, scheme, [0:0.1:20]);
surfplot(UF, t, pr{3});


pr{4} = kdv;
[t, U, UF] = expglm(pr{4}, [0 2*pi/625], 2*pi/625/1000, ...
    scheme, [0:2*pi/625/100:2*pi/625]);
surfplot(UF, t, pr{4});

pr{6} = hochostint;
[t, U, UF] = expglm(pr{6}, [0 1], 0.01, scheme, [0:0.01:1]);
surfplot(UF, t, pr{6});

pr{7} = allencahn;
[t, U, UF] = expglm(pr{7}, [0 10], 0.01, scheme, [0:0.01:10]);
surfplot(UF, t, pr{7});

pr{8} = sinegordon;
[t, U, UF] = expglm(pr{8}, [0 40], 0.01, scheme, [0:0.1:10]);
surfplot(UF, t, pr{8});

