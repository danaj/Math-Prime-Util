#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util
   qw/moebius mertens euler_phi jordan_totient divisor_sum exp_mangoldt
      chebyshev_theta chebyshev_psi carmichael_lambda znorder liouville
      znprimroot znlog kronecker legendre_phi gcd lcm is_power valuation
      invmod vecsum binomial
     /;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my @moeb_vals = (qw/ 1 -1 -1 0 -1 1 -1 0 0 1 -1 0 -1 1 1 0 -1 0 -1 0 /);
my %mertens = (
        1 =>    1,
        2 =>    0,
        3 =>   -1,
        4 =>   -1,
        5 =>   -2,
       10 =>   -1,
      100 =>    1,
     1000 =>    2,
    10000 =>  -23,
        8 =>   -2,
       16 =>   -1,
       32 =>   -4,
       64 =>   -1,
      128 =>   -2,
      256 =>   -1,
      512 =>   -4,
     1024 =>   -4,
     2048 =>    7,
     4096 =>  -19,
     8192 =>   22,
);
my %big_mertens = (
   100000 =>  -48,
  1000000 =>  212,
 10000000 => 1037,
);
if (!$extra && !Math::Prime::Util::prime_get_config->{'xs'}) {
  delete $big_mertens{10000000};
}
if ($extra && $use64) {
  %big_mertens = ( %big_mertens,
          2 =>  0,      # A087987, mertens at primorials
          6 => -1,
         30 => -3,
        210 => -1,
       2310 => -1,
      30030 => 16,
     510510 => -25,
    9699690 => 278,
  223092870 => 3516,

    6433477 => 900,     # 30^2
  109851909 => -4096,   # A084235, 2^12

      2**14 =>  -32,    # A084236
      2**15 =>   26,
      2**16 =>   14,
      2**17 =>  -20,
      2**18 =>   24,
      2**19 => -125,
      2**20 =>  257,
      2**21 => -362,
      2**22 =>  228,
      2**23 =>  -10,

     10**8  => 1928,
#     10**9  => -222,
#  1*10**10 => -33722,  # From Deleglise and Rivat
#  2*10**10 => -48723,  # Too slow with current method
  );
}

my %totients = (
     123456 => 41088,
     123457 => 123456,
  123456789 => 82260072,
);
my @A000010 = (0,1,1,2,2,4,2,6,4,6,4,10,4,12,6,8,8,16,6,18,8,12,10,22,8,20,12,18,12,28,8,30,16,20,16,24,12,36,18,24,16,40,12,42,20,24,22,46,16,42,20,32,24,52,18,40,24,36,28,58,16,60,30,36,32,48,20,66,32,44);
#@totients{0..$#A000010} = @A000010;

my @A001615 = (1,3,4,6,6,12,8,12,12,18,12,24,14,24,24,24,18,36,20,36,32,36,24,48,30,42,36,48,30,72,32,48,48,54,48,72,38,60,56,72,42,96,44,72,72,72,48,96,56,90,72,84,54,108,72,96,80,90,60,144,62,96,96,96,84,144,68,108,96);

my %jordan_totients = (
  # A000010
  1 => [1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4, 12, 6, 8, 8, 16, 6, 18, 8, 12, 10, 22, 8, 20, 12, 18, 12, 28, 8, 30, 16, 20, 16, 24, 12, 36, 18, 24, 16, 40, 12, 42, 20, 24, 22, 46, 16, 42, 20, 32, 24, 52, 18, 40, 24, 36, 28, 58, 16, 60, 30, 36, 32, 48, 20, 66, 32, 44],
  # A007434
  2 => [1, 3, 8, 12, 24, 24, 48, 48, 72, 72, 120, 96, 168, 144, 192, 192, 288, 216, 360, 288, 384, 360, 528, 384, 600, 504, 648, 576, 840, 576, 960, 768, 960, 864, 1152, 864, 1368, 1080, 1344, 1152, 1680, 1152, 1848, 1440, 1728, 1584, 2208, 1536],
  # A059376
  3 => [1, 7, 26, 56, 124, 182, 342, 448, 702, 868, 1330, 1456, 2196, 2394, 3224, 3584, 4912, 4914, 6858, 6944, 8892, 9310, 12166, 11648, 15500, 15372, 18954, 19152, 24388, 22568, 29790, 28672, 34580, 34384, 42408, 39312, 50652, 48006, 57096],
  # A059377
  4 => [1, 15, 80, 240, 624, 1200, 2400, 3840, 6480, 9360, 14640, 19200, 28560, 36000, 49920, 61440, 83520, 97200, 130320, 149760, 192000, 219600, 279840, 307200, 390000, 428400, 524880, 576000, 707280, 748800, 923520, 983040, 1171200],
  # A059378
  5 => [1, 31, 242, 992, 3124, 7502, 16806, 31744, 58806, 96844, 161050, 240064, 371292, 520986, 756008, 1015808, 1419856, 1822986, 2476098, 3099008, 4067052, 4992550, 6436342, 7682048, 9762500, 11510052, 14289858, 16671552, 20511148, 23436248, 28629150, 32505856, 38974100, 44015536, 52501944, 58335552, 69343956, 76759038, 89852664, 99168256, 115856200, 126078612, 147008442, 159761600, 183709944, 199526602, 229345006, 245825536, 282458442, 302637500, 343605152, 368321664],
  # A069091
  6 => [1, 63, 728, 4032, 15624, 45864, 117648, 258048, 530712, 984312, 1771560, 2935296, 4826808, 7411824, 11374272, 16515072, 24137568, 33434856, 47045880, 62995968, 85647744, 111608280, 148035888, 187858944, 244125000, 304088904, 386889048],
  # A069092
  7 => [1, 127, 2186, 16256, 78124, 277622, 823542, 2080768, 4780782, 9921748, 19487170, 35535616, 62748516, 104589834, 170779064, 266338304, 410338672, 607159314, 893871738, 1269983744, 1800262812, 2474870590, 3404825446],
);

my %sigmak = (
  # A0000005
  0 => [1,2,2,3,2,4,2,4,3,4,2,6,2,4,4,5,2,6,2,6,4,4,2,8,3,4,4,6,2,8,2,6,4,4,4,9,2,4,4,8,2,8,2,6,6,4,2,10,3,6,4,6,2,8,4,8,4,4,2,12,2,4,6,7,4,8,2,6,4,8,2,12,2,4,6,6,4,8,2,10,5,4,2,12,4,4,4,8,2,12,4,6,4,4,4,12,2,6,6,9,2,8,2,8],
  # A000203
  1 => [1, 3, 4, 7, 6, 12, 8, 15, 13, 18, 12, 28, 14, 24, 24, 31, 18, 39, 20, 42, 32, 36, 24, 60, 31, 42, 40, 56, 30, 72, 32, 63, 48, 54, 48, 91, 38, 60, 56, 90, 42, 96, 44, 84, 78, 72, 48, 124, 57, 93, 72, 98, 54, 120, 72, 120, 80, 90, 60, 168, 62, 96, 104, 127, 84, 144, 68, 126, 96, 144],
  # A001157
  2 => [1, 5, 10, 21, 26, 50, 50, 85, 91, 130, 122, 210, 170, 250, 260, 341, 290, 455, 362, 546, 500, 610, 530, 850, 651, 850, 820, 1050, 842, 1300, 962, 1365, 1220, 1450, 1300, 1911, 1370, 1810, 1700, 2210, 1682, 2500, 1850, 2562, 2366, 2650, 2210, 3410, 2451, 3255],
  # A001158
  3 => [1, 9, 28, 73, 126, 252, 344, 585, 757, 1134, 1332, 2044, 2198, 3096, 3528, 4681, 4914, 6813, 6860, 9198, 9632, 11988, 12168, 16380, 15751, 19782, 20440, 25112, 24390, 31752, 29792, 37449, 37296, 44226, 43344, 55261, 50654, 61740, 61544],
);

my @tau4 = (1,4,4,10,4,16,4,20,10,16,4,40,4,16,16,35,4,40,4,40,16,16,4,80,10,16,20,40,4,64,4,56,16,16,16,100,4,16,16,80,4,64,4,40,40,16,4,140,10,40,16,40,4,80,16,80,16,16,4,160,4,16,40,84,16,64,4,40,16,64,4,200,4,16,40,40,16);

my %mangoldt = (
-13 => 1,
  0 => 1,
  1 => 1,
  2 => 2,
  3 => 3,
  4 => 2,
  5 => 5,
  6 => 1,
  7 => 7,
  8 => 2,
  9 => 3,
 10 => 1,
 11 => 11,
 25 => 5,
 27 => 3,
 399981 => 1,
 399982 => 1,
 399983 => 399983,
 823543 => 7,
 83521 => 17,
 130321 => 19,
);

my %chebyshev1 = (
       0 =>       0,
       1 =>       0,
       2 =>       0.693147180559945,
       3 =>       1.79175946922805,
       4 =>       1.79175946922805,
       5 =>       3.40119738166216,
     243 =>     226.593507136467,
  123456 =>  123034.091739914,
);
my %chebyshev2 = (
       0 =>       0,
       1 =>       0,
       2 =>       0.693147180559945,
       3 =>       1.79175946922805,
       4 =>       2.484906649788,
       5 =>       4.0943445622221,
     243 =>     245.274469978683,
  123456 =>  123435.148054491
);
if ($extra) {
  $chebyshev1{1234567} = 1233272.80087825;
  $chebyshev2{1234567} = 1234515.17962833;
}
if (!$usexs && !$extra) {
  delete $chebyshev1{$_} for grep { $_ > 50000 } keys %chebyshev1;
  delete $chebyshev2{$_} for grep { $_ > 50000 } keys %chebyshev2;
}

my @A002322 = (0,1,1,2,2,4,2,6,2,6,4,10,2,12,6,4,4,16,6,18,4,6,10,22,2,20,12,18,6,28,4,30,8,10,16,12,6,36,18,12,4,40,6,42,10,12,22,46,4,42,20,16,12,52,18,20,6,18,28,58,4,60,30,6,16,12,10,66,16,22,12,70,6,72,36,20,18,30,12,78,4,54,40,82,6,16,42,28,10,88,12,12,22,30,46,36,8,96,42,30,20,100,16,102,12,12,52,106,18,108,20,36,12,112,18,44,28,12,58,48,4,110,60,40,30,100,6,126,32,42,12,130,10,18,66,36,16,136,22,138,12,46,70,60,12,28,72,42,36,148,20,150,18,48,30,60,12,156,78,52,8,66,54,162,40,20,82,166,6,156,16,18,42,172,28,60,20,58,88,178,12,180,12,60,22,36,30,80,46,18,36,190,16,192,96,12,42,196,30,198,20);

my @mult_orders = (
  [1, 35, 1],
  [2, 35, 12],
  [4, 35, 6],
  [6, 35, 2],
  [7, 35],
  #[2,1000000000000031,81788975100],
  [1, 1, 1],
  [0, 0],
  [1, 0],
  [25, 0],
  [1, 1, 1],
  [19, 1, 1],
  [1, 19, 1],
  [2, 19, 18],
  [3, 20, 4],
  [57,1000000003,189618],
  [67,999999749,30612237],
  [22,999991815,69844],
  [10,2147475467,31448382],
  [141,2147475467,1655178],
  [7410,2147475467,39409],
  [31407,2147475467,266],
);

my %primroots = (
   -11 => 2,
     0 => undef,
     1 => 0,
     2 => 1,
     3 => 2,
     4 => 3,
     5 => 2,
     6 => 5,
     7 => 3,
     8 => undef,
     9 => 2,
    10 => 3,          # 3 is the smallest root.  Pari gives the other root 7.
      1729 => undef,  # Pari goes into an infinite loop.
   5109721 =>  94,
  17551561 =>  97,
  90441961 => 113,
1407827621 =>   2,
1520874431 =>  17,
1685283601 => 164,
 100000001 => undef,  # Without an early exit, this will essentially hang.
);
if ($use64) {
  $primroots{2232881419280027} = 6;         # factor divide goes to FP
  $primroots{14123555781055773271} = 6;     # bmodpow hits RT 71548
  $primroots{89637484042681} = 335;         # smallest root is large
}

my @kroneckers = (
  [ 109981, 737777,  1],
  [ 737779, 121080, -1],
  [-737779, 121080,  1],
  [ 737779,-121080, -1],
  [-737779,-121080, -1],
  [12345,331,-1],
  [1001,9907,-1],
  [19,45,1],
  [8,21,-1],
  [5,21,1],
  [5,1237,-1],
  [10, 49, 1],
  [123,4567,-1],
  [3,18,0], [3,-18,0],
  [-2, 0, 0],  [-1, 0, 1],  [ 0, 0, 0],  [ 1, 0, 1],  [ 2, 0, 0],
  [-2, 1, 1],  [-1, 1, 1],  [ 0, 1, 1],  [ 1, 1, 1],  [ 2, 1, 1],
  [-2,-1,-1],  [-1,-1,-1],  [ 0,-1, 1],  [ 1,-1, 1],  [ 2,-1, 1],
  # Some cases trying to make sure we're not turning UVs into IVs
  [ 3686556869,  428192857,  1],
  [-1453096827,  364435739, -1],
  [ 3527710253, -306243569, 1],
  [-1843526669, -332265377, 1],
  [  321781679, 4095783323, -1],
  [  454249403,  -79475159, -1],
);
if ($use64) {
  push @kroneckers, [17483840153492293897, 455592493, 1];
  push @kroneckers, [-1402663995299718225, 391125073, 1];
  push @kroneckers, [16715440823750591903, -534621209, -1];
  push @kroneckers, [13106964391619451641,16744199040925208803, 1];
  push @kroneckers, [11172354269896048081,10442187294190042188,-1];
  push @kroneckers, [-5694706465843977004,9365273357682496999,-1];
}

my @valuations = (
  [-4,2, 2],
  [0,0, 0],
  [1,0, 0],
  [96552,6, 3],
  [1879048192,2, 28],
);

my @legendre_sums = (
  [ 0,  92372, 0],
  [ 5,  15, 1],
  [ 89, 4, 21 ],
  [ 46, 4, 11 ],
  [ 47, 4, 12 ],
  [ 48, 4, 12 ],
  [ 52, 4, 12 ],
  [ 53, 4, 13 ],
  [10000, 5, 2077],
  [526, 7, 95],
  [588, 6, 111],
  [100000, 5, 20779],
  [5882, 6, 1128],
  [100000, 7, 18053],
  [10000, 8, 1711],
  [1000000, 168, 78331],
  [800000, 213, 63739],
);

my @gcds = (
  [ [], 0],
  [ [8], 8],
  [ [9,9], 9],
  [ [0,0], 0],
  [ [1, 0, 0], 1],
  [ [0, 0, 1], 1],
  [ [17,19], 1 ],
  [ [54,24], 6 ],
  [ [42,56], 14],
  [ [ 9,28], 1 ],
  [ [48,180], 12],
  [ [2705353758,2540073744,3512215098,2214052398], 18],
  [ [2301535282,3609610580,3261189640], 106],
  [ [694966514,510402262,195075284,609944479], 181],
  [ [294950648,651855678,263274296,493043500,581345426], 58 ],
  [ [-30,-90,90], 30],
  [ [-3,-9,-18], 3],
);
my @lcms = (
  [ [], 0],
  [ [8], 8],
  [ [9,9], 9],
  [ [0,0], 0],
  [ [1, 0, 0], 0],
  [ [0, 0, 1], 0],
  [ [17,19], 323 ],
  [ [54,24], 216 ],
  [ [42,56], 168],
  [ [ 9,28], 252 ],
  [ [48,180], 720],
  [ [36,45], 180],
  [ [-36,45], 180],
  [ [-36,-45], 180],
  [ [30,15,5], 30],
  [ [2,3,4,5], 60],
  [ [30245, 114552], 3464625240],
  [ [11926,78001,2211], 2790719778],
  [ [1426,26195,3289,8346], 4254749070],
);
if ($use64) {
  push @gcds, [ [12848174105599691600,15386870946739346600,11876770906605497900], 700];
  push @gcds, [ [9785375481451202685,17905669244643674637,11069209430356622337], 117];
  push @lcms, [ [26505798,9658520,967043,18285904], 15399063829732542960];
  push @lcms, [ [267220708,143775143,261076], 15015659316963449908];
}

my @znlogs = (
 [ [5,2,1019], 10],
 [ [2,4,17], undef],
 [ [7,3,8], undef],
 [ [3,3,8], 1],
 [ [10,2,101], 25],
 [ [2,55,101], 73],         # 2 = 55^73 mod 101
 [ [228,2,383], 110],
 [ [3061666278, 499998, 3332205179], 22],
);
if ($usexs) {
  push @znlogs, [ [5678,5,10007], 8620];  # 5678 = 5^8620 mod 10007
  push @znlogs, [ [5675,5,10000019], 2003974];  # 5675 = 5^2003974 mod 10000019
}

my %powers = (
  0 => [-2, -1, 0, 1, 2, 3, 5, 6, 7, 10, 11, 12, 13, 14, 15, 17, 18, 19],
  2 => [4, 9, 25, 36, 49],
  3 => [8, 27, 125, 343, 17576],
  4 => [16, 38416],
  9 => [19683, 1000000000],
);
if ($use64) {
  push @{$powers{0}}, 9908918038843197151;
  push @{$powers{2}}, 18446743927680663841;
  push @{$powers{3}}, 2250923753991375;
  push @{$powers{4}}, 1150530828529256001;
  push @{$powers{9}}, 118587876497;
}

my @invmods = (
 [ 0, 0, undef],
 [ 1, 0, undef],
 [ 0, 1, undef],
 [ 1, 1, 0],
 [ 45, 59, 21],
 [  42,  2017, 1969],
 [  42, -2017, 1969],
 [ -42,  2017, 48],
 [ -42, -2017, 48],
 [ 14, 28474, undef],
);
if ($use64) {
 push @invmods, [ 13, 9223372036854775808, 5675921253449092805 ];
 push @invmods, [ 14, 18446744073709551615, 17129119497016012214 ];
}

my @vecsums = (
  [ 0 ],
  [ -1, -1 ],
  [ 0, 1,-1 ],
  [ 0, -1,1 ],
  [ 0, -1,1 ],
  [ 0, -2147483648,2147483648 ],
  [ 0, "-4294967296","4294967296" ],
  [ 0, "-9223372036854775808","9223372036854775808" ],
  [ "18446744073709551615", "18446744073709551615","-18446744073709551615","18446744073709551615" ],
  [ "55340232221128654848", "18446744073709551616","18446744073709551616","18446744073709551616" ],
);

my @binomials = (
 [ 0,0, 1 ],
 [ 0,1, 0 ],
 [ 1,0, 1 ],
 [ 1,1, 1 ],
 [ 1,2, 0 ],
 [ 13,13, 1 ],
 [ 13,14, 0 ],
 [ 35,16, 4059928950 ],             # We can do this natively even in 32-bit
 [ 40,19, "131282408400" ],         # We can do this in 64-bit
 [ 67,31, "11923179284862717872" ], # ...and this
 [ 228,12, "30689926618143230620" ],# But the result of this is too big.
 [ 177,78, "3314450882216440395106465322941753788648564665022000" ],
 [ -10,5, -2002 ],
 [ -11,22, 64512240 ],
 [ -12,23, -286097760 ],
 [ -12,-23, 0 ],
 [  12,-23, 0 ],
 [  12,-12, 0 ],
 [ -12,0, 1 ],
 [  0,-1, 0 ],
);

# These are slow with XS, and *really* slow with PP.
if (!$usexs) {
  %big_mertens = map { $_ => $big_mertens{$_} }
                 grep { $_ < 100000000 }
                 keys %big_mertens;
}

my @liouville_pos = (qw/24 51 94 183 294 629 1488 3684 8006 8510 32539 57240
   103138 238565 444456 820134 1185666 3960407 4429677 13719505 29191963
   57736144 134185856 262306569 324235872 563441153 1686170713 2489885844/);
my @liouville_neg = (qw/23 47 113 163 378 942 1669 2808 8029 9819 23863 39712
   87352 210421 363671 562894 1839723 3504755 7456642 14807115 22469612
   49080461 132842464 146060791 279256445 802149183 1243577750 3639860654/);
if ($use64) {
  push @liouville_pos, (qw/1260238066729040 10095256575169232896/);
  push @liouville_neg, (qw/1807253903626380 12063177829788352512/);
}


plan tests => 0 + 1
                + 1 # Small Moebius
                + 3*scalar(keys %mertens)
                + 1*scalar(keys %big_mertens)
                + 2 # Small Phi
                + 8 + scalar(keys %totients)
                + 1 # Small Carmichael Lambda
                + scalar(@kroneckers)
                + scalar(@gcds)
                + scalar(@lcms)
                + scalar(@mult_orders)
                + scalar(@znlogs)
                + scalar(@legendre_sums)
                + scalar(@valuations)
                + 3 + scalar(@invmods)
                + scalar(@vecsums)
                + scalar(@binomials)
                + scalar(keys %powers)
                + scalar(keys %primroots) + 2
                + scalar(keys %jordan_totients)
                + 2  # Dedekind psi calculated two ways
                + 2  # Calculate J5 two different ways
                + 2 * $use64 # Jordan totient example
                + 1 + 2*scalar(keys %sigmak) + 3
                + scalar(keys %mangoldt)
                + scalar(keys %chebyshev1)
                + scalar(keys %chebyshev2)
                + scalar(@liouville_pos) + scalar(@liouville_neg);

ok(!eval { moebius(0); }, "moebius(0)");

{
  my @moebius = map { moebius($_) } (1 .. scalar @moeb_vals);
  is_deeply( \@moebius, \@moeb_vals, "moebius 1 .. " . scalar @moeb_vals );
}

while (my($n, $mertens) = each (%mertens)) {
  my $M = 0;
  $M += moebius($_) for (1 .. $n);
  is( $M, $mertens, "sum(moebius(k) for k=1..$n) == $mertens" );
  # Calculate using ranged moebius
  $M = 0;
  $M += $_ for moebius(1,$n);
  is( $M, $mertens, "sum(moebius(1..$n) == $mertens" );
  # Now with mertens function
  is( mertens($n), $mertens, "mertens($n) == $mertens" );
}
while (my($n, $mertens) = each (%big_mertens)) {
  is( mertens($n), $mertens, "mertens($n)" );
}

{
  my @phi = map { euler_phi($_) } (0 .. $#A000010);
  is_deeply( \@phi, \@A000010, "euler_phi 0 .. $#A000010" );
}
{
  my @phi = euler_phi(0, $#A000010);
  is_deeply( \@phi, \@A000010, "euler_phi with range: 0, $#A000010" );
}
while (my($n, $phi) = each (%totients)) {
  is( euler_phi($n), $phi, "euler_phi($n) == $phi" );
}
is_deeply( [euler_phi(0,0)], [0],     "euler_phi(0,0)" );
is_deeply( [euler_phi(1,0)], [],      "euler_phi with end < start" );
is_deeply( [euler_phi(0,1)], [0,1],   "euler_phi 0-1" );
is_deeply( [euler_phi(1,2)], [1,1],   "euler_phi 1-2" );
is_deeply( [euler_phi(1,3)], [1,1,2], "euler_phi 1-3" );
is_deeply( [euler_phi(2,3)], [1,2],   "euler_phi 2-3" );
is_deeply( [euler_phi(10,20)], [4,10,4,12,6,8,8,16,6,18,8], "euler_phi 10-20" );
is_deeply( [euler_phi(1513,1537)],
   [qw/1408 756 800 756 1440 440 1260 576 936 760 1522 504 1200 648
       1016 760 1380 384 1530 764 864 696 1224 512 1456/],
           "euler_phi(1513,1537)" );

###### Jordan Totient
while (my($k, $tref) = each (%jordan_totients)) {
  my @tlist = map { jordan_totient($k, $_) } 1 .. scalar @$tref;
  is_deeply( \@tlist, $tref, "Jordan's Totient J_$k" );
}

{
  my @psi_viaj;
  my @psi_viamobius;
  foreach my $n (1 .. scalar @A001615) {
    push @psi_viaj, int(jordan_totient(2, $n) / jordan_totient(1, $n));
    push @psi_viamobius, int($n * divisor_sum( $n, sub { moebius($_[0])**2 / $_[0] } ) + 0.5);
  }
  is_deeply( \@psi_viaj, \@A001615, "Dedekind psi(n) = J_2(n)/J_1(n)" );
  is_deeply( \@psi_viamobius, \@A001615, "Dedekind psi(n) = divisor_sum(n, moebius(d)^2 / d)" );
}

{
  my $J5 = $jordan_totients{5};
  my @J5_jordan = map { jordan_totient(5, $_) } 1 .. scalar @$J5;
  is_deeply( \@J5_jordan, $J5, "Jordan totient 5, using jordan_totient");
  my @J5_moebius = map { my $n = $_; divisor_sum($n, sub { my $d=shift; $d**5 * moebius($n/$d); }) } 1 .. scalar @$J5;
  is_deeply( \@J5_moebius, $J5, "Jordan totient 5, using divisor sum" );
}

if ($use64) {
  is( jordan_totient(4, 12345), 22902026746060800, "J_4(12345)" );
  # Apostal page 48, 17a.
  is( divisor_sum( 12345, sub { jordan_totient(4,$_[0]) } ),
      # was int(12345 ** 4), but Perl 5.8.2 gets it wrong.
      int(12345*12345*12345*12345),
      "n=12345, k=4  :   n**k = divisor_sum(n, jordan_totient(k, d))" );
}

###### Divisor sum
while (my($k, $sigmaref) = each (%sigmak)) {
  my @slist;
  foreach my $n (1 .. scalar @$sigmaref) {
    push @slist, divisor_sum( $n, sub { int($_[0] ** $k) } );
  }
  is_deeply( \@slist, $sigmaref, "Sum of divisors to the ${k}th power: Sigma_$k" );
  @slist = ();
  foreach my $n (1 .. scalar @$sigmaref) {
    push @slist, divisor_sum( $n, $k );
  }
  is_deeply( \@slist, $sigmaref, "Sigma_$k using integer instead of sub" );
}
# k=1 standard sum -- much faster
{
  my @slist = map { divisor_sum($_) } 1 .. scalar @{$sigmak{1}};
  is_deeply(\@slist, $sigmak{1}, "divisor_sum(n)");
}
# tau two ways
{
  my $len = scalar @{$sigmak{0}};
  my @slist1 = map { divisor_sum($_, sub {1}) } 1 .. $len;
  my @slist2 = map { divisor_sum($_, 0      ) } 1 .. $len;
  is_deeply( \@slist1, $sigmak{0}, "tau as divisor_sum(n, sub {1})" );
  is_deeply( \@slist2, $sigmak{0}, "tau as divisor_sum(n, 0)" );
}

{
  # tau_4 A007426
  my @t;
  foreach my $n (1 .. scalar @tau4) {
    push @t, divisor_sum($n, sub { divisor_sum($_[0],sub { divisor_sum($_[0],0) }) });
  }
  is_deeply( \@t, \@tau4, "Tau4 (A007426), nested divisor sums" );
}

###### Exponential of von Mangoldt
while (my($n, $em) = each (%mangoldt)) {
  is( exp_mangoldt($n), $em, "exp_mangoldt($n) == $em" );
}

###### first Chebyshev function
while (my($n, $c1) = each (%chebyshev1)) {
  cmp_closeto( chebyshev_theta($n), $c1, 1e-9*abs($n), "chebyshev_theta($n)" );
}
###### second Chebyshev function
while (my($n, $c2) = each (%chebyshev2)) {
  cmp_closeto( chebyshev_psi($n), $c2, 1e-9*abs($n), "chebyshev_psi($n)" );
}

###### Carmichael Lambda
{
  my @lambda = map { carmichael_lambda($_) } (0 .. $#A002322);
  is_deeply( \@lambda, \@A002322, "carmichael_lambda with range: 0, $#A000010" );
}
###### kronecker
foreach my $karg (@kroneckers) {
  my($a, $n, $exp) = @$karg;
  my $k = kronecker($a, $n);
  is( $k, $exp, "kronecker($a, $n) = $exp" );
}
###### gcd
foreach my $garg (@gcds) {
  my($aref, $exp) = @$garg;
  my $gcd = gcd(@$aref);
  is( $gcd, $exp, "gcd(".join(",",@$aref).") = $exp" );
}
###### lcm
foreach my $garg (@lcms) {
  my($aref, $exp) = @$garg;
  my $lcm = lcm(@$aref);
  is( $lcm, $exp, "lcm(".join(",",@$aref).") = $exp" );
}
###### znorder
foreach my $moarg (@mult_orders) {
  my ($a, $n, $exp) = @$moarg;
  my $zn = znorder($a, $n);
  is( $zn, $exp, "znorder($a, $n) = " . ((defined $exp) ? $exp : "<undef>") );
}
###### znprimroot
while (my($n, $root) = each (%primroots)) {
  is( znprimroot($n), $root, "znprimroot($n) == " . ((defined $root) ? $root : "<undef>") );
}
is( znprimroot("-100000898"), 31, "znprimroot(\"-100000898\") == 31" );
is( znprimroot("+100000898"), 31, "znprimroot(\"+100000898\") == 31" );
###### znlog
foreach my $arg (@znlogs) {
  my($aref, $exp) = @$arg;
  my ($a, $g, $p) = @$aref;
  my $k = znlog($a,$g,$p);
  is( $k, $exp, "znlog($a,$g,$p) = " . ((defined $exp) ? $exp : "<undef>") );
}
###### liouville
foreach my $i (@liouville_pos) {
  is( liouville($i),  1, "liouville($i) = 1" );
}
foreach my $i (@liouville_neg) {
  is( liouville($i), -1, "liouville($i) = -1" );
}
###### Legendre phi
foreach my $r (@legendre_sums) {
  my($x, $a, $exp) = @$r;
  is( legendre_phi($x, $a), $exp, "legendre_phi($x,$a) = $exp" );
}

###### is_power
while (my($e, $vals) = each (%powers)) {
  my @fail;
  foreach my $val (@$vals) {
    push @fail, $val unless is_power($val) == $e;
  }
  ok( @fail == 0, "is_power returns $e for " . join(",",@fail) );
}

###### valuation
foreach my $r (@valuations) {
  my($n, $k, $exp) = @$r;
  is( valuation($n, $k), $exp, "valuation($n,$k) = $exp" );
}
###### invmod
ok(!eval { invmod(undef,11); }, "invmod(undef,11)");
ok(!eval { invmod(11,undef); }, "invmod(11,undef)");
ok(!eval { invmod('nan',11); }, "invmod('nan',11)");

foreach my $r (@invmods) {
  my($a, $n, $exp) = @$r;
  is( invmod($a,$n), $exp, "invmod($a,$n) = ".((defined $exp)?$exp:"<undef>") );
}
###### vecsum
foreach my $r (@vecsums) {
  my($exp, @vals) = @$r;
  is( vecsum(@vals), $exp, "vecsum(@vals) = $exp" );
}
###### binomial
foreach my $r (@binomials) {
  my($n, $k, $exp) = @$r;
  is( binomial($n,$k), $exp, "binomial($n,$k)) = $exp" );
}

sub cmp_closeto {
  my $got = shift;
  my $expect = shift;
  my $tolerance = shift;
  my $message = shift;
  cmp_ok( abs($got - $expect), '<=', $tolerance, $message );
}
