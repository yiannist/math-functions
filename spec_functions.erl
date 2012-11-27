-module(spec_functions).

-export([
        %% * Gamma function
          log_gamma/1
        , log_gamma_l/1
        , incomplete_gamma/2
        , inv_incomplete_gamma/2
        %% * Beta function
        , log_beta/2
        , incomplete_beta/3
        , incomplete_beta_/4
        , inv_incomplete_beta/3
        %% * Logarithm
        , log1p/1
        , log2/1
        %% * Factorial
        , factorial/1
        , log_factorial/1
        , stirling_error/1
        %% * Combinatorics
        , choose/2
        %% * Extras
        , bd0/2 
]).

-export_type([
               float64/0
             ]).

-include("constants.hrl").

%% Special functions and factorials.

-type float64() :: float()
                 | 'pos_inf'
                 | 'neg_inf'.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gamma function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Adapted from http://people.sc.fsu.edu/~burkardt/f_src/asa245/asa245.html
%%
%% Compute the logarithm of the gamma function Gamma(x). Uses
%% Algorithm AS 245 by Macleod.
%%
%% Gives an accuracy of 10-12 significant decimal digits, except
%% for small regions around x = 1 and x = 2, where the function
%% goes to zero.  For greater accuracy, use 'logGammaL'.
%%
%% Returns 'pos_inf' if the input is outside of the range
%% (0 < x <= 1e305).
-spec log_gamma(float()) -> float64().
log_gamma(X) ->
    %% A bunch of constants:
    Y = math:log(X),
    Alr2pi = 0.918938533204673,
    K = X * (Y - 1) - 0.5 * Y + Alr2pi,
    X1 = 1 / X,
    X2 = X1 * X1,
    R1_0 =  -2.66685511495,   R1_1 =  -24.4387534237,    R1_2 = -21.9698958928,
    R1_3 =  11.1667541262,    R1_4 =    3.13060547623,   R1_5 =   0.607771387771,
    R1_6 =  11.9400905721,    R1_7 =   31.4690115749,    R1_8 =  15.2346874070,

    R2_0 = -78.3359299449,    R2_1 = -142.046296688,     R2_2 = 137.519416416,
    R2_3 =  78.6994924154,    R2_4 =    4.16438922228,   R2_5 =  47.0668766060,
    R2_6 = 313.399215894,     R2_7 =  263.505074721,     R2_8 =  43.3400022514,

    R3_0 =  -2.12159572323e5, R3_1 =    2.30661510616e5, R3_2 =   2.74647644705e4,
    R3_3 =  -4.02621119975e4, R3_4 =   -2.29660729780e3, R3_5 =  -1.16328495004e5,
    R3_6 =  -1.46025937511e5, R3_7 =   -2.42357409629e4, R3_8 =  -5.70691009324e2,

    R4_0 = 0.279195317918525,  R4_1 = 0.4917317610505968,
    R4_2 = 0.0692910599291889, R4_3 = 3.350343815022304,
    R4_4 = 6.012459259764103,

    %% Use a helper instead of guards! ;-)
    LogGamma_wrap =
        fun (X00) ->
                {A, B, C} =
                    case X00 < 0.5 of
                        true  -> {-Y, X00 + 1, X00};
                        false -> {0,  X00,     X00 - 1}
                    end,
                LogGamma = fun (X0) when X0 =< 0    -> 'pos_inf';
                               (X0) when X0 < 1.5   ->
                                   A + C * ((((R1_4 * B + R1_3) * B + R1_2) * B
                                     + R1_1) * B + R1_0) / ((((B + R1_8) * B
                                     + R1_7) * B + R1_6) * B + R1_5);
                               (X0) when X0 < 4     ->
                                   (X - 2) * ((((R2_4 * X + R2_3) * X + R2_2) * X
                                     + R2_1) * X + R2_0) / ((((X + R2_8) * X
                                     + R2_7) * X + R2_6) * X + R2_5);
                               (X0) when X0 < 12    ->
                                   ((((R3_4 * X + R3_3) * X + R3_2) * X + R3_1) * X
                                     + R3_0) / ((((X + R3_8) * X + R3_7) * X
                                     + R3_6) * X + R3_5);
                               (X0) when X0 > 5.1e5 -> K;
                               (_)                ->
                                   K + X1 * ((R4_2 * X2 + R4_1) * X2 + R4_0) /
                                     ((X2 + R4_4) * X2 + R4_3)
                           end,
                LogGamma(X00)
        end,
    LogGamma_wrap(X).

%% Compute the logarithm of the gamma function, G(x). Uses a
%% Lanczos approximation.
%%
%% This function is slower than 'log_gamma', but gives 14 or more
%% significant decimal digits of accuracy, except around x = 1 and
%% x = 2, where the function goes to zero.
%%
%% Returns M_HUGE; if the input is outside of the range
%% (0 < x <= 1e305).
-spec log_gamma_l(float()) -> float64().
log_gamma_l(X) when X =< 0 -> 'pos_inf';
log_gamma_l(X) ->
    A0  = 0.9999999999995183,
    X65 = X + 6.5,
    Go  = fun (K, {L, T}) -> {L + K / T, T - 1} end,
    A   =  [ 0.1659470187408462e-06, 0.9934937113930748e-05
           , -0.1385710331296526,    12.50734324009056
           , -176.6150291498386,     771.3234287757674
           , -1259.139216722289,     676.5203681218835
           ],
    Fini = fun ({L, _}) ->
                   math:log(L + A0) + math:log(?M_SQRT_2_PI) - X65 + (X - 0.5) *
                       math:log(X65)
           end,
    Fini(lists:foldl(Go, {0, X + 7}, A)).

%% Compute the log gamma correction factor for X > 10.  This
%% correction factor is suitable for an alternate (but less
%% numerically accurate) definition of 'log_gamma':
%%
%% > lgg x = 0.5 * log(2*pi) + (x-0.5) * log x - x + logGammaCorrection x
-spec log_gamma_correction(float()) -> float().
log_gamma_correction(X) when X < 10 ->
    exit({?MODULE, log_gamma_correction, "NaN"});
log_gamma_correction(X) ->
    Big = 94906265.62425156,
    T = 10 / X,
    Coeffs = [  0.1666389480451863247205729650822e+0
             , -0.1384948176067563840732986059135e-4
             ,  0.9810825646924729426157171547487e-8
             , -0.1809129475572494194263306266719e-10
             ,  0.6221098041892605227126015543416e-13
             , -0.3399615005417721944303330599666e-15
             ,  0.2683181998482698748957538846666e-17
             ],
    case X < Big of
        true  -> polynomial:chebyshev_broucke(T * T * 2 - 1, Coeffs) / X;
        false -> 1 / (X * 12)
    end.

%% Compute the normalized lower incomplete gamma function
%% g(s, x). Normalization means that g(s, inf)=1.
%% Uses Algorithm AS 239 by Shea.
-spec incomplete_gamma(float(), float()) -> float64().
incomplete_gamma(P0, X0) ->
    Norm = fun (A) -> 0.5 * math:erfc( -A / ?M_SQRT_2 ) end,
    Limit = -88,
    Tolerance = 1.0e-14,
    Overflow = 1.0e37,
    CF = fun () ->
                 A = 1 - P0,
                 B = A + X0 + 1,
                 P3 = X0 + 1,
                 P4 = X0 * B,
                 cont_frac(Tolerance, Overflow, A, B, 0, 1.0, X0, P3, P4, P3 / P4)
         end,
    IncompleteGamma = fun (P, X) when ((X < 0) or (P =< 0)) -> 'pos_inf';
                          (_, 0) -> 0;
                          (P, X) when P >= 1000 ->
                              Norm(3 * math:sqrt(P) *
                                   (math:pow(X / P, 1 / 3) + 1 / (9 * P) - 1));
                          (_, X) when X > 1.0e8 -> 1;
                          (P, X) when ((X =< 1) or (X < P)) ->
                              A = P * math:log(X) - X - log_gamma(P + 1),
                              G = A + math:log(pearson(Tolerance, X, P, 1.0, 1.0)),
                              case G > Limit of
                                  true  -> math:exp(G);
                                  false -> 0
                              end;
                          (P, X) ->
                              G = P*math:log(X) - X - log_gamma(P) + math:log(CF()),
                              case G > Limit of
                                  true  -> 1 - math:exp(G);
                                  false -> 1
                              end
                      end,
    IncompleteGamma(P0, X0).

-spec pearson(float(), float(), float(), float(), float()) -> float().
pearson(Tolerance, X, A, C, G) ->
    A2 = A + 1,
    C2 = C * X / A2,
    G2 = G + C2,
    case C2 =< Tolerance of
        true  -> G2;
        false -> pearson(Tolerance, X, A2, C2, G2)
    end.

-spec cont_frac(float(), float(), float(), float(), integer(), float()
                , float(), float(), float(), float()) -> float().
cont_frac(Tolerance, Overflow, A, B, C, P1, P2, P3, P4, G) ->
    A2 = A + 1,
    B2 = B + 2,
    C2 = C + 1,
    AN = A2 * C2,
    P5 = B2 * P3 - AN * P1,
    P6 = B2 * P4 - AN * P2,
    RN = P5 / P6,
    F  = fun (N) ->  case abs(P5) > Overflow of
                         true  -> N / Overflow;
                         false -> N
                     end
         end,
    case abs(G - RN) =< min(Tolerance, Tolerance * RN) of
        true  -> G;
        false -> cont_frac(Tolerance, Overflow, A2, B2, C2,
                           F(P3), F(P4), F(P5), F(P6), RN)
    end.


%% Adapted from Numerical Recipes Â§6.2.1

%%  Inverse incomplete gamma function. It's approximately inverse of
%%  'incomplete_gamma' for the same s. So following equality
%%  approximately holds:
%%
%%  > invIncompleteGamma s . incompleteGamma s = id
%%
%%  For inv_incomplete_gamma s p, s must be positive and p must be
%%  in [0,1] range.
%%STUB: Might notice some deviation from Haskell results (~1e-11).
-spec inv_incomplete_gamma(float(), float()) -> float64().
inv_incomplete_gamma(A, P) ->
    %% Calculate initial guess for root
    Guess =
        fun () when A > 1 ->
                Lg = case P < 0.5 of
                         true  -> P;
                         false -> 1 - P
                     end,
                T = math:sqrt(-2 * math:log(Lg)),
                X1 = (2.30753 + T * 0.27061) / (1 + T * (0.99229
                      + T * 0.04481)) - T,
                X2 = case P < 0.5 of
                         true  -> -X1;
                         false -> X1
                     end,
                max(1.0e-3, (A * math:pow((1 - 1/(9 * A) - X2 /
                      (3 * math:sqrt(A))), 3)));
            %% For a <= 1 use following approximations:
            %%   γ(a,1) ≈ 0.253a + 0.12a²
            %%
            %%   γ(a,x) ≈ γ(a,1)·x^a                               x <  1
            %%   γ(a,x) ≈ γ(a,1) + (1 - γ(a,1))(1 - exp(1 - x))    x >= 1
            () ->
                T = 1 - A * (0.253 + A * 0.12),
                case P < T of
                    true  -> math:pow(P / T, 1 / A);
                    false -> 1 - math:log( 1 - (P - T) / (1 - T) )
                end
        end,
    Fun = fun () when A =< 0 ->
                  exit({?MODULE, inv_incomplete_gamma, "'a' must be positive!"});
              () when ((P < 0) or (P > 1)) ->
                  exit({?MODULE, inv_incomplete_gamma, "'p' must be in [0,1]"});
              () when P == 0 -> 0;
              () when P == 1 -> 'pos_inf';
              () -> inv_incomplete_gamma__go({A, P}, 0, Guess())
          end,
    Fun().

%% Solve equation γ(a,x) = p using Halley method
-spec inv_incomplete_gamma__go({float(), float()}, integer(), float())
                              -> float().
inv_incomplete_gamma__go(_, I, X) when I >= 12 -> X;
inv_incomplete_gamma__go({A, P}, I, X) ->
    %% Value of γ(a,x) - p
    F  = incomplete_gamma(A, X) - P,
    %% Constants
    EPS  = 1.0e-8,
    A1   = A - 1,
    GLN  = log_gamma(A),
    %% dγ(a,x)/dx
    F2 = case A > 1 of
             true  -> LnA1 = math:log(A1),
                      AFAC = math:exp(A1 * (LnA1 - 1) - GLN),
                      AFAC * math:exp( -(X - A1) + A1 * (math:log(X) - LnA1));
             false -> math:exp( -X + A1 * math:log(X) - GLN)
         end,
    U  = F / F2,
    %% Halley correction to Newton-Rapson step
    Corr = U * (A1 / X - 1),
    DX = U / (1 - 0.5 * min(1.0, Corr)),
    %% New approximation to x
    X2 = case X < DX of
             true  -> 0.5 * X; %% Do not go below 0
             false -> X - DX
         end,
    case abs(DX) < EPS * X2 of
        true  -> X2;
        false -> inv_incomplete_gamma__go({A, P}, I + 1, X2)
    end.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beta function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute the natural logarithm of the beta function.
-spec log_beta(float(), float()) -> float64().
log_beta(A, B) ->
    P = min(A, B),
    Q = max(A, B),
    F = fun () when P < 0   -> exit({?MODULE, log_beta, "NaN"});
            () when P == 0  -> 'pos_inf';
            () when P >= 10 -> PQ  = P + Q,
                               PPQ = P / PQ,
                               C   = log_gamma_correction(Q) - log_gamma_correction(PQ),
                               math:log(Q) * (-0.5) + ?M_LN_SQRT_2_PI
                                   + log_gamma_correction(P) + C
                                   + (P-0.5) * math:log(PPQ) + Q * log1p(-PPQ);
            () when Q >= 10 -> PQ  = P + Q,
                               PPQ = P / PQ,
                               C = log_gamma_correction(Q) - log_gamma_correction(PQ),
                               log_gamma(P) + C + P - P * math:log(PQ)
                                   + (Q - 0.5) * log1p(-PPQ);
            ()              -> log_gamma(P) + log_gamma(Q) - log_gamma(P + Q)
        end,
    F().

%% Regularized incomplete beta function. Uses algorithm AS63 by
%% Majumder and Bhattachrjee.
-spec incomplete_beta(float(), float(), float()) -> float().
incomplete_beta(P, Q, X) ->
    incomplete_beta_(log_beta(P, Q), P, Q, X).

%% Regularized incomplete beta function. Same as 'incomplete_beta'
%% but also takes logarithm of beta function as parameter.
-spec incomplete_beta_(float(), float(), float(), float()) -> float().
incomplete_beta_(_, P, Q, _) when ((P =< 0) or (Q =< 0)) ->
    exit({?MODULE, incomplete_beta_, "p <= 0 || q <= 0"});
incomplete_beta_(_, _, _, X) when ((X < 0) or (X > 1)) ->
    exit({?MODULE, incomplete_beta_, "x out of [0,1] range"});
incomplete_beta_(_, _, _, X) when ((X == 0) or (X == 1)) ->
    X;
incomplete_beta_(Beta, P, Q, X) when P >= (P + Q) * X ->
    incomplete_beta_worker(Beta, P, Q, X);
incomplete_beta_(Beta, P, Q, X) ->
    1 - incomplete_beta_worker(Beta, Q, P, 1 - X).

%% Worker for incomplete beta function. It is separate function to
%% avoid confusion with parameter during parameter swapping
-spec incomplete_beta_worker(float(), float(), float(), float()) -> float().
incomplete_beta_worker(Beta, P, Q, X) ->
    incomplete_beta_worker__go({Beta, P, Q, X}, P + Q,
                               trunc(Q + (1 - X) * (P + Q)), 1, 1.0, 1.0).

%% Loop.
-spec incomplete_beta_worker__go({float(), float(), float(), float()},
                                 float(), integer(), integer(), float(), float())
                                -> float().
incomplete_beta_worker__go({Beta, P, Q, X}, PSQ, NS, AI, Term, Betain) ->
    %% New values
    EPS = 1.0e-15,
    CX = 1 - X,
    Fact = case NS > 0 of
               true  -> (Q - AI) * X / CX;
               false -> case NS == 0 of
                            true  -> (Q - AI) * X;
                            false -> PSQ * X
                        end
           end,
    Term2 = Term * Fact / (P + AI),
    Betain2 = Betain + Term2,
    %% Iterations are complete
    DB = abs(Term2),
    Done = (DB =< EPS) and (DB =< EPS * Betain2),
    PSQ2 = case NS < 0 of
               true  -> PSQ + 1;
               false -> PSQ
           end,
    case Done of
        true  ->
            Betain2 * math:exp(P*math:log(X) + (Q-1) * math:log(CX) - Beta) / P;
        false ->
            incomplete_beta_worker__go({Beta, P, Q, X}, PSQ2, NS-1
                                       , AI+1, Term2, Betain2)
    end.

%% Compute inverse of regularized incomplete beta function. Uses
%% initial approximation from AS109 and Halley method to solve equation.
-spec inv_incomplete_beta(float(), float(), float()) -> float().
inv_incomplete_beta(P, Q, _) when ((P =< 0) or (Q =< 0)) ->
    exit({?MODULE, inv_incomplete_beta, "p <=0 || q <=0"});
inv_incomplete_beta(_, _, A) when ((A < 0) or (A > 1)) ->
    exit({?MODULE, inv_incomplete_beta, "bad a"});
inv_incomplete_beta(_, _, A) when ((A == 0) or (A == 1)) -> A;
inv_incomplete_beta(P, Q, A) when A > 0.5 ->
    1 - inv_incomplete_beta_worker(log_beta(P, Q), Q, P, 1 - A);
inv_incomplete_beta(P, Q, A) ->
    inv_incomplete_beta_worker(log_beta(P, Q), P, Q, A).

-spec inv_incomplete_beta_worker(float(), float(), float(), float()) -> float().
inv_incomplete_beta_worker(Beta, P, Q, A) ->
    %% Calculate initial guess
    R = math:sqrt( -math:log(A * A) ),
    Y = R - ( 2.30753 + 0.27061 * R ) / ( 1.0 + ( 0.99229 + 0.04481 * R ) * R ),
    T = 1 / (9 * Q),
    T2 = 2 * Q * math:pow(1 - T + Y * math:sqrt(T), 3),
    T3 = (4*P + 2*Q - 2) / T2,
    Guess = fun () when ((P > 1) and (Q > 1)) ->
                    RR = (Y*Y - 3) / 6,
                    SS = 1 / (2*P - 1),
                    TT = 1 / (2*Q - 1),
                    HH = 2 / (SS + TT),
                    WW = Y * math:sqrt(HH + RR) / HH - (HH - SS)
                        * (RR + 5/6 - 2 / (3 * HH)),
                    P / (P + Q * math:exp(2 * WW));
                () when T2 =< 0 ->
                    1 - math:exp( (math:log((1 - A) * Q) + Beta) / Q );
                () when T3 =< 1 ->
                    math:exp( (math:log(A * P) + Beta) / P );
                () -> 1 - 2 / (T3 + 1)
            end,
    inv_incomplete_beta_worker__go({Beta, P, Q, A}, 0, Guess()).

-spec inv_incomplete_beta_worker__go({float(), float(), float(), float()},
                                     integer(), float()) -> float().
inv_incomplete_beta_worker__go({Beta, P, Q, A}, I, X) ->
    Fun = fun () when ((X == 0) or (X == 1) or (I >= 10)) -> X;
              () -> P1 = P - 1,
                    Q1 = Q - 1,
                    F  = incomplete_beta_(Beta, P, Q, X) - A,
                    F2  = math:exp(P1 * math:log(X) + Q1 * math:log(1-X) - Beta),
                    U  = F / F2,
                    DX = U / (1 - 0.5 * min(1, (U * (P1 / X - Q1 / (1 - X))))),
                    Z  = X - DX,
                    X2 = fun () when Z < 0 -> X / 2;
                             () when Z > 1 -> (X + 1) / 2;
                             () -> Z
                         end,
                    case abs(DX) =< 16 * ?M_EPSILON * X of
                        true  -> X;
                        false ->
                            inv_incomplete_beta_worker__go({Beta, P, Q, A}
                                                           , I + 1, X2())
                    end
          end,
    Fun().

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Logarithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute the natural logarithm of 1 + x.  This is accurate even
%% for values of x near zero, where use of log(1+x) would lose
%% precision.
-spec log1p(float()) -> float64().
log1p(0)             -> 0;
log1p(-1)            -> 'neg_inf';
log1p(X) when X < -1 -> exit({?MODULE, log1p, "NaN"});
log1p(X)             ->
    X2 = abs(X),
    Coeffs = [ 0.10378693562743769800686267719098e+1,
              -0.13364301504908918098766041553133e+0,
               0.19408249135520563357926199374750e-1,
              -0.30107551127535777690376537776592e-2,
               0.48694614797154850090456366509137e-3,
              -0.81054881893175356066809943008622e-4,
               0.13778847799559524782938251496059e-4,
              -0.23802210894358970251369992914935e-5,
               0.41640416213865183476391859901989e-6,
              -0.73595828378075994984266837031998e-7,
               0.13117611876241674949152294345011e-7,
              -0.23546709317742425136696092330175e-8,
               0.42522773276034997775638052962567e-9,
              -0.77190894134840796826108107493300e-10,
               0.14075746481359069909215356472191e-10,
              -0.25769072058024680627537078627584e-11,
               0.47342406666294421849154395005938e-12,
              -0.87249012674742641745301263292675e-13,
               0.16124614902740551465739833119115e-13,
              -0.29875652015665773006710792416815e-14,
               0.55480701209082887983041321697279e-15,
              -0.10324619158271569595141333961932e-15
             ],
    Fun = fun () when X2 < ?M_EPSILON * 0.5       -> X;
              () when ((X >= 0) and (X < 1.0e-8)) ->
                  X * (1 - X * 0.5);
              () when X2 < 0.375                  ->
                  X * (1 - X * polynomial:chebyshev_broucke(X / 0.375, Coeffs));
              () -> math:log(1 + X)
          end,
    Fun().

%% Compute the logarithm in base 2 of the given value.
-spec log2(integer()) -> integer().
log2(V0) when V0 =< 0 -> exit({?MODULE, log2, "invalid input"});
log2(V0) ->
    log2__go(5, 0, V0).

-spec log2__go(integer(), integer(), integer()) -> integer().
log2__go(-1, R, _) -> R;
log2__go( I, R, V) ->
    BV = [16#2, 16#c, 16#f0, 16#ff00, 16#ffff0000, 16#ffffffff00000000],
    SV = [1,2,4,8,16,32],
    case (V band lists:nth(I+1, BV)) =/= 0 of
        true  -> SI = lists:nth(I+1, SV),
                 log2__go(I-1, R bor SI, V bsr SI);
        false -> log2__go(I-1, R, V)
    end.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Factorial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute the factorial function n!.  Returns ?M_HUGE if the
%% input is above 170 (above which the result cannot be represented by
%% a 64-bit 'float').
-spec factorial(integer()) -> float64().
factorial(N) when N  < 0   -> exit({?MODULE, factorial, "negative input"});
factorial(N) when N =< 1   -> 1.0;
factorial(N) when N =< 170 -> Prod = fun (X, Acc) -> X * Acc end,
                              lists:foldl(Prod, 1.0, lists:seq(2, N));
factorial(_)               -> 'pos_inf'.

%% Compute the natural logarithm of the factorial function. Gives
%% 16 decimal digits of precision.
-spec log_factorial(integer()) -> float().
log_factorial(N) when N =< 14 -> math:log(factorial(N));
log_factorial(N) ->
    X = N + 1,
    Y = 1 / (X * X),
    Z = ((-(5.95238095238e-4 * Y) + 7.936500793651e-4) * Y -
             2.7777777777778e-3) * Y + 8.3333333333333e-2,
    (X - 0.5) * math:log(X) - X + 9.1893853320467e-1 + Z / X.

%% Calculate the error term of the Stirling approximation. This is
%% only defined for non-negative values.
%%
%% > stirlingError n = log(n!) - log(sqrt(2*pi*n)*(n/e)^n)
-spec stirling_error(float()) -> float().
stirling_error(N) ->
    NN = N * N,
    S0 = 0.083333333333333333333,        %% 1/12
    S1 = 0.00277777777777777777778,      %% 1/360
    S2 = 0.00079365079365079365079365,   %% 1/1260
    S3 = 0.000595238095238095238095238,  %% 1/1680
    S4 = 0.0008417508417508417508417508, %% 1/1188
    Sfe = [ 0.0,
            0.1534264097200273452913848,   0.0810614667953272582196702,
            0.0548141210519176538961390,   0.0413406959554092940938221,
            0.03316287351993628748511048,  0.02767792568499833914878929,
            0.02374616365629749597132920,  0.02079067210376509311152277,
            0.01848845053267318523077934,  0.01664469118982119216319487,
            0.01513497322191737887351255,  0.01387612882307074799874573,
            0.01281046524292022692424986,  0.01189670994589177009505572,
            0.01110455975820691732662991,  0.010411265261972096497478567,
            0.009799416126158803298389475, 0.009255462182712732917728637,
            0.008768700134139385462952823, 0.008330563433362871256469318,
            0.007934114564314020547248100, 0.007573675487951840794972024,
            0.007244554301320383179543912, 0.006942840107209529865664152,
            0.006665247032707682442354394, 0.006408994188004207068439631,
            0.006171712263039457647532867, 0.005951370112758847735624416,
            0.005746216513010115682023589, 0.005554733551962801371038690 ],
    Fun = fun () when N =< 15  ->
                  ProperFraction = fun (N0) -> (2 * N0) == (2 * trunc(N0)) end,
                  case ProperFraction(N) of
                      true  -> lists:nth(2 * trunc(N), Sfe);
                      false -> log_gamma(N+1.0) - (N+0.5) * math:log(N)
                                   + N - ?M_LN_SQRT_2_PI
                  end;
              () when N  > 500 -> (S0-S1/NN)/N;
              () when N  > 80  -> (S0-(S1-S2/NN)/NN)/N;
              () when N  > 35  -> (S0-(S1-(S2-S3/NN)/NN)/NN)/N;
              ()               -> (S0-(S1-(S2-(S3-S4/NN)/NN)/NN)/NN)/N
          end,
    Fun().

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combinatorics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Quickly compute the natural logarithm of n `choose` k, with
%% no checking.
-spec log_choose_fast(float(), number()) -> float().
log_choose_fast(N, K) ->
    -math:log(N + 1) - log_beta(N - K + 1, K + 1).

%% Compute the binomial coefficient n `choose` k. For
%% values of k > 30, this uses an approximation for performance
%% reasons. The approximation is accurate to 12 decimal places in the
%% worst case.
%%
%% Example:
%%
%% 1> choose(7, 3).
%% 35
-spec choose(integer(), integer()) -> float().
choose(N, K) when K > N -> 0;
choose(N, K) ->
    K2 = min(K, N-K),
    NK = N - K2,
    Approx = math:exp(log_choose_fast(N, K2)),
    Fun = fun () when K2 < 50         -> Go = fun (I, A) -> A * (NK+I) / I end,
                                         lists:foldl(Go, 1, lists:seq(1, K2));
              () when Approx < ?MAX64 -> round(Approx); %STUB: Is ?MAX64 correct?
              ()                      -> Approx
          end,
    Fun().


%% From Numeric.SpecFunctions.Extra:

%% Evaluate the deviance term: x log(x/np) + np - x.
%% STUB: stripped isInfinite guards!
-spec bd0(float(), float()) -> float().
bd0(_, 0)  -> exit({?MODULE, bd0, "NaN"});
bd0(X, NP) ->
    X_NP = X - NP,
    case abs(X_NP) >= 0.1 * (X + NP) of
        true  -> X * math:log(X/NP) - X_NP;
        false -> V = X_NP / (X + NP),
                 S0 = X_NP * V,
                 EJ0 = 2 * X * V,
                 VV = V * V,
                 bd0__go(VV, 1, EJ0 * VV, S0)
    end.

-spec bd0__go(float(), float(), float(), float()) -> float().                     
bd0__go(VV, J, EJ, S) ->
    S2 = S + EJ / (2 * J + 1), 
    case S2 == S of
        true  -> S2; %% FIXME: Comparing floats for equality!
        false -> bd0__go(VV, J + 1, (EJ * VV), S2)
    end.                  


%% References:
%%
%% * Lanczos, C. (1964) A precision approximation of the gamma
%%   function.  /SIAM Journal on Numerical Analysis B/
%%   1:86&#8211;96. <http://www.jstor.org/stable/2949767>
%%
%% * Loader, C. (2000) Fast and Accurate Computation of Binomial
%%   Probabilities. <http://projects.scipy.org/scipy/raw-attachment/ticket/620/loader2000Fast.pdf>
%%
%% * Macleod, A.J. (1989) Algorithm AS 245: A robust and reliable
%%   algorithm for the logarithm of the gamma function.
%%   /Journal of the Royal Statistical Society, Series C (Applied Statistics)/
%%   38(2):397&#8211;402. <http://www.jstor.org/stable/2348078>
%%
%% * Shea, B. (1988) Algorithm AS 239: Chi-squared and incomplete
%%   gamma integral. /Applied Statistics/
%%   37(3):466&#8211;473. <http://www.jstor.org/stable/2347328>
%%
%% * Majumder, K.L., Bhattacharjee, G.P. (1973) Algorithm AS 63: The
%%   Incomplete Beta Integral. /Journal of the Royal Statistical
%%   Society. Series C (Applied Statistics)/ Vol. 22, No. 3 (1973),
%%   pp. 409-411. <http://www.jstor.org/pss/2346797>
%%
%% * Majumder, K.L., Bhattacharjee, G.P. (1973) Algorithm AS 64:
%%   Inverse of the Incomplete Beta Function Ratio. /Journal of the
%%   Royal Statistical Society. Series C (Applied Statistics)/
%%   Vol. 22, No. 3 (1973), pp. 411-414
%%   <http://www.jstor.org/pss/2346798>
%%
%% * Cran, G.W., Martin, K.J., Thomas, G.E. (1977) Remark AS R19
%%   and Algorithm AS 109: A Remark on Algorithms: AS 63: The
%%   Incomplete Beta Integral AS 64: Inverse of the Incomplete Beta
%%   Function Ratio. /Journal of the Royal Statistical Society. Series
%%   C (Applied Statistics)/ Vol. 26, No. 1 (1977), pp. 111-114
%%   <http://www.jstor.org/pss/2346887>
