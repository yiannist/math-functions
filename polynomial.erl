-module(polynomial).

-export([
          chebyshev/2
        , chebyshev_broucke/2  
        ]).

%% Chebyshev polynomials.
%%
%% A Chebyshev polynomial of the first kind is defined by the
%% following recurrence:
%%
%% > t 0 _ = 1
%% > t 1 x = x
%% > t n x = 2 * x * t (n-1) x - t (n-2) x

%% Evaluate a Chebyshev polynomial of the first kind. Uses
%% Clenshaw's algorithm.
-spec chebyshev(float(), [float()]) -> float().
chebyshev(X, A) ->
    X2 = 2 * X,
    Step = fun (K, {B0, B1}) ->
                   {K + X2 * B0 - B1, B0}
           end,
    Fini = fun ({B0, B1}) ->
                   hd(A) + X * B0 - B1
           end,
    Fini(lists:foldr(Step, {0, 0}, tl(A))).

%% Evaluate a Chebyshev polynomial of the first kind. Uses Broucke's
%% ECHEB algorithm, and his convention for coefficient handling, and so
%% gives different results than 'chebyshev' for the same inputs.
-spec chebyshev_broucke(float(), [float()]) -> float().
chebyshev_broucke(X, A) ->
    X2 = 2 * X,
    Step = fun (K, {B0, B1, _}) ->
                   {K + X2 * B0 - B1, B0, B1}
           end,
    Fini = fun ({B0, _, B2}) ->
                   (B0 - B2) * 0.5
           end,
    Fini(lists:foldr(Step, {0, 0, 0}, A)).


%% References:
%%
%% * Broucke, R. (1973) Algorithm 446: Ten subroutines for the
%%   manipulation of Chebyshev series. /Communications of the ACM/
%%   16(4):254&#8211;256.  <http://doi.acm.org/10.1145/362003.362037>
%%
%% * Clenshaw, C.W. (1962) Chebyshev series for mathematical
%%   functions. /National Physical Laboratory Mathematical Tables 5/,
%%   Her Majesty's Stationery Office, London.
