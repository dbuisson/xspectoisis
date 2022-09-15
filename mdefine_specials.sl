
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define the special functions allowed by the XSPEC mdefine syntax.
%
% Many of these are already ISIS built-ins so are not defined here:
%
% exp
% sin
% cos
% tan
% sinh
% cosh
% tanh
% sqrt
% abs
% asin
% acos
% atan
% asinh
% acosh
% atanh
% sign
% mean
% atan2
%
% There are some convention differences between XSPEC mdefine and ISIS builtin functions:
% Logarithm bases are corrected by the parser.
% (Vector) smin/smax and (binary) min/max are converted by the parser.
%
% Note that erf/erfc require the statistics ("stats") module to be available.
%
% Not yet implemented:
%
% dim (vector expression dimension) might be hard. Hopefully rarely needed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Trig functions in degrees

define sind(x){return sin(x/57.29577951308232);};
define cosd(x){return cos(x/57.29577951308232);};
define tand(x){return tan(x/57.29577951308232);};

define sinhd(x){return sinh(x/57.29577951308232);};
define coshd(x){return cosh(x/57.29577951308232);};
define tanhd(x){return tanh(x/57.29577951308232);};

%% Gamma function
% Lanczos approximation
define gamma(x)
{
    variable p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    variable g = 7;
    
    variable a, t, lowfac, where_small;
    
    %% Use reflection for x<0.5:
    where_small = x < 0.5;
    lowfac = ( PI / sin(PI * x) )^where_small;
    
    x     += (1-2*x)*where_small;
    
    x -= 1;
    a = p[0];
    t = x + g + 0.5;
    
    variable i;
    for (i = 1; i<=8; i++)
    {
        a += p[i] / (x + i);
    };
    
    return lowfac * (2.5066282746310002 * (t ^ (x + 0.5)) * exp(-t) * a)^(1-2*where_small);
};

%% Step functions
define heaviside(x){return sign(x)/2.+0.5;};
define boxcar(x){return (sign(x)-sign(x-1))/2.;};

%% Error function
% Using non-standard ERF/ERFC, which seem to be what XSPEC uses.
require("stats");
define erf (z){ return   2*normal_cdf(1.4142135623730951*z); };
define erfc(z){ return 2-2*normal_cdf(1.4142135623730951*z); };
%define erf (z){ return   2*normal_cdf(1.4142135623730951*z)-1  ; };
%define erfc(z){ return 1-2*normal_cdf(1.4142135623730951*z)    ; };
