%% Constant values common to much numeric code.

%% A very large number.
-define(M_HUGE, 1.7976931348623157e308).

%% A very small number.
-define(M_TINY, 2.2250738585072014e-308).

%% The largest 'Int' /x/ such that 2**(/x/-1) is approximately representable
%% as a 'Double'.
-define(M_MAX_EXP, 1024). 

%% sqrt(2)
-define(M_SQRT_2, 1.4142135623730950488016887242096980785696718753769480731766).

%% sqrt(2 * pi)
-define(M_SQRT_2_PI, 2.5066282746310005024157652848110452530069867406099383166299).

%% 2 / sqrt(pi)
-define(M_2_SQRT_PI, 1.1283791670955125738961589031215451716881012586579977136881).

%% 1 / sqrt(2)
-define(M_1_SQRT_2, 0.7071067811865475244008443621048490392848359376884740365883).

%% The smallest 'float' such that 1 + X = X.
-define(M_EPSILON, 1.1102230246251566637e-16).

%% log( sqrt(2 * pi) )
-define(M_LN_SQRT_2_PI, 0.9189385332046727417803297364056176398613974736377834128171).

%% The biggest 'small integer' (60 bits).
%% http://www.erlang.org/doc/efficiency_guide/advanced.html
-define(MAX64, 576460752303423488).
