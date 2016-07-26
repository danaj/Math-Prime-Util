#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util
   qw/moebius mertens euler_phi jordan_totient divisor_sum exp_mangoldt
      chebyshev_theta chebyshev_psi carmichael_lambda znorder liouville
      znprimroot znlog kronecker legendre_phi gcd lcm is_power valuation
      binomial gcdext chinese vecmin vecmax factorial
      hammingweight sqrtint rootint logint is_square_free
      is_carmichael is_quasi_carmichael
      is_primitive_root
      ramanujan_sum hclassno ramanujan_tau
     /;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
my $usegmp= Math::Prime::Util::prime_get_config->{'gmp'};
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
my %isf = (
  0 => 0,
  1 => 1,
  2 => 1,
  3 => 1,
  4 => 0,
  5 => 1,
  6 => 1,
  7 => 1,
  8 => 0,
  9 => 0,
  10 => 1,
  11 => 1,
  12 => 0,
  13 => 1,
  14 => 1,
  15 => 1,
  16 => 0,
  758096738 => 0,
  434420340 => 0,
  870589313 => 0,
  752518565 => 1,
  695486396 => 0,
  723570005 => 1,
  602721315 => 0,
  418431087 => 0,
  506916483 => 1,
  617459403 => 1,
);

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

my @tau4 = (1,4,4,10,4,16,4,20,10,16,4,40,4,16,16,35,4,40,4,40,16,16,4,80,10,16,20,40,4,64,4,56,16,16,16,100);
push @tau4, (4,16,16,80,4,64,4,40,40,16,4,140,10,40,16,40,4,80,16,80,16,16,4,160,4,16,40,84,16,64,4,40,16,64,4,200,4,16,40,40,16) if $extra;

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

my %hclassno = (
      -3 => 0,
       0 => -1,
       1 => 0,
       2 => 0,
       3 => 4,
       4 => 6,
       7 => 12,
       8 => 12,
      11 => 12,
      12 => 16,
      20 => 24,
      23 => 36,
      39 => 48,
      47 => 60,
      71 => 84,
     163 => 12,
     427 => 24,
     907 => 36,
    1555 => 48,
    6307 => 96,
   20563 => 156,
   30067 => 168,
   31243 => 192,
   34483 => 180,
    4031 => 1008,
);

my %rtau = (
       0 => 0,
       1 => 1,
       2 => -24,
       3 => 252,
       4 => -1472,
       5 => 4830,
      53 => -1596055698,
     106 => 38305336752,
     243 => 13400796651732,
   16089 => "12655813883111729342208",
);

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
if ($use64) {
  push @mult_orders, [2, 2405286912458753, 1073741824];  # Pari #1031
}

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
  $primroots{9223372036854775837} = 5;      # Pari #905
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

my @popcounts = (
  [0, 0],
  [1, 1],
  [2, 1],
  [3, 2],
  [452398, 12],
  [-452398, 12],
  [4294967295, 32],
  ["777777777777777714523989234823498234098249108234236", 83],
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

my @gcdexts = (
  [ [0,  0], [0, 0, 0] ],
  [ [0, 28], [0, 1,28] ],
  [ [ 28,0], [ 1,0,28] ],
  [ [0,-28], [0,-1,28] ],
  [ [-28,0], [-1,0,28] ],
  [ [ 3706259912, 1223661804], [ 123862139,-375156991, 4] ],
  [ [ 3706259912,-1223661804], [ 123862139, 375156991, 4] ],
  [ [-3706259912, 1223661804], [-123862139,-375156991, 4] ],
  [ [-3706259912,-1223661804], [-123862139, 375156991, 4] ],
  [ [22,242], [1, 0, 22] ],
  [ [2731583792,3028241442], [-187089956, 168761937, 2] ],
  [ [42272720,12439910], [-21984, 74705, 70] ],
);
if ($use64) {
  push @gcdexts, [ [10139483024654235947,8030280778952246347], [-2715309548282941287,3428502169395958570,1] ];
}

my @crts = (
  [ [], 0 ],
  [ [[4,5]], 4 ],
  [ [[77,11]], 0 ],
  [ [[0,5],[0,6]], 0 ],
  [ [[14,5],[0,6]], 24 ],
  [ [[10,11],[4,22],[9,19]], undef ],
  [ [[77,13],[79,17]], 181 ],
  [ [[2,3],[3,5],[2,7]], 23 ],
  [ [[10,11],[4,12],[12,13]], 1000 ],
  [ [[42,127],[24,128]], 2328 ],             # Some tests from Mod::Int
  [ [[32,126],[23,129]], 410 ],
  [ [[2328,16256],[410,5418]], 28450328 ],
  [ [[1,10],[11,100]], 11 ],
  [ [[11,100],[22,100]], undef ],
  [ [[1753051086,3243410059],[2609156951,2439462460]], "6553408220202087311"],
  [ [ ["6325451203932218304","2750166238021308"],
      ["5611464489438299732","94116455416164094"] ],
    "1433171050835863115088946517796" ],
  [ [ ["1762568892212871168","8554171181844660224"],
      ["2462425671659520000","2016911328009584640"] ],
    "188079320578009823963731127992320" ],
  [ [ ["856686401696104448","11943471150311931904"],
      ["6316031051955372032","13290002569363587072"] ],
    "943247297188055114646647659888640" ],
  [ [[-3105579549,3743000622],[-1097075646,1219365911]], "2754322117681955433"],
  [ [ ["-925543788386357567","243569243147991"],
      ["-1256802905822510829","28763455974459440"] ],
    "837055903505897549759994093811" ],
  [ [ ["-2155972909982577461","8509855219791386062"],
      ["-5396280069505638574","6935743629860450393"] ],
    "12941173114744545542549046204020289525" ],
);

my @znlogs = (
 [ [5,2,1019], 10],
 [ [2,4,17], undef],
 [ [7,3,8], undef],
 [ [7,17,36], undef],       # No solution (Pari #1463)
 [ [1,8,9], 0],
 [ [3,3,8], 1],
 [ [10,2,101], 25],
 [ [2,55,101], 73],         # 2 = 55^73 mod 101
 [ [5,2,401], 48],          # 5 = 2^48 mod 401  (Pari #1285)
 [ [228,2,383], 110],
 [ [3061666278, 499998, 3332205179], 22],
 [ [5678,5,10007], 8620],   # 5678 = 5^8620 mod 10007
 [ [7531,6,8101], 6689],    # 7531 = 6^6689 mod 8101
 # Some odd cases.  Pari pre-2.6 and post 2.6 have issues with them.
 [ [0,30,100], 2],          # 0 = 30^2 mod 100
 [ [1,1,101], 0],           # 1 = 1^0 mod 101
 [ [8,2,102], 3],           # 8 = 2^3 mod 102
 [ [18,18,102], 1],         # 18 = 18^1 mod 102
);
if ($usexs || $extra) {
  push @znlogs, [[5675,5,10000019], 2003974];  # 5675 = 5^2003974 mod 10000019
  push @znlogs, [[18478760,5,314138927], 34034873];
}
if ($usexs && $use64) {
  # Nice case for PH
  push @znlogs, [[32712908945642193,5,71245073933756341], 5945146967010377];
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
my @negpowers = (0,0,0,3,0,5,3,7,0,9,5);

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
 [ -23,-26, -2300 ],     # Kronenburg extension
 [ -12,-23, -705432 ],   # same
 [  12,-23, 0 ],
 [  12,-12, 0 ],
 [ -12,0, 1 ],
 [  0,-1, 0 ],
);

my @roots = (
  [25,  3, 15625],
  [13,  4, 28561],
  [13,  5, 371293],
  [25,  6, 244140625],
  [ 7,  7, 823543],
  [13,  8, 815730721],
  [ 7,  9, 40353607],
  [13, 10, "137858491849"],
  [21, 11, "350277500542221"],
  [25, 12, "59604644775390625"],
  [ 7, 13, "96889010407"],
  [ 7, 14, "678223072849"],
  [13, 16, "665416609183179841"],
  [13, 17, "8650415919381337933"],
  [ 7, 18, "1628413597910449"],
  [ 6, 19, "609359740010496"],
  [ 3, 21, "10460353203"],
  [ 3, 23, "94143178827"],
  [ 3, 25, "847288609443"],
  [ 3, 29, "68630377364883"],
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
                + 1*scalar(keys %isf)
                + 2 # is_carmichael
                + 4 # is_quasi_carmichael
                + 2 # Small Phi
                + 9 + scalar(keys %totients)
                + 1 # Small Carmichael Lambda
                + scalar(@kroneckers)
                + scalar(@gcds)
                + scalar(@lcms)
                + scalar(@gcdexts)
                + scalar(@crts)
                + scalar(@mult_orders)
                + scalar(@znlogs)
                + scalar(@legendre_sums)
                + scalar(@valuations)
                + scalar(@popcounts)
                + 4  # sqrtint
                + 6  # rootint
                + 5  # logint
                + 2 + scalar(@binomials)
                + 6 + scalar(keys %powers) + scalar(@negpowers)
                + scalar(keys %primroots) + 1
                + scalar(keys %primroots) + 2  # is_primitive_root
                + scalar(keys %jordan_totients)
                + 2  # Dedekind psi calculated two ways
                + 2  # Calculate J5 two different ways
                + 2 * $use64 # Jordan totient example
                + 1 + 2*scalar(keys %sigmak) + 3
                + scalar(keys %mangoldt)
                + scalar(keys %chebyshev1)
                + scalar(keys %chebyshev2)
                + 3  # Ramanujan sum
                + scalar(keys %hclassno)
                + scalar(keys %rtau)
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
while (my($n, $isf) = each (%isf)) {
  is( is_square_free($n), $isf, "is_square_free($n)" );
}
{
  is_deeply( [grep { is_carmichael($_) } 1 .. 20000],
             [561,1105,1729,2465,2821,6601,8911,10585,15841],
             "Carmichael numbers to 20000" );
  SKIP: {
    skip "Skipping large Carmichael", 1 unless $usegmp;
    ok( is_carmichael("341627175004511735787409078802107169251"), "Large Carmichael" );
  }
}
{
  is_deeply( [grep { is_quasi_carmichael($_) } 1 .. 400],
             [35,77,143,165,187,209,221,231,247,273,299,323,357,391,399],
             "Quasi-Carmichael numbers to 400" );
  is( scalar(grep { is_quasi_carmichael($_) } 1 .. 5000),
             95,
             "95 Quasi-Carmichael numbers under 5000" );
  is(is_quasi_carmichael(5092583), 1, "5092583 is a Quasi-Carmichael number with 1 base");
  is(is_quasi_carmichael(777923), 7, "777923 is a Quasi-Carmichael number with 7 bases");
}

{
  my @phi = map { euler_phi($_) } (0 .. $#A000010);
  is_deeply( \@phi, \@A000010, "euler_phi 0 .. $#A000010" );
}
{
  my @phi = euler_phi(0, $#A000010);
  is_deeply( \@phi, \@A000010, "euler_phi with range: 0, $#A000010" );
}
{
  my $s = 0;
  $s += $_ for euler_phi(1, 240);
  is($s, 17544, "sum of totients to 240");
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
  my @tlist = map { jordan_totient(0+$k, $_) } 1 .. scalar @$tref;
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
  is( exp_mangoldt(0+$n), $em, "exp_mangoldt($n) == $em" );
}

###### first Chebyshev function
while (my($n, $c1) = each (%chebyshev1)) {
  cmp_closeto( chebyshev_theta(0+$n), $c1, 1e-9*abs($n), "chebyshev_theta($n)" );
}
###### second Chebyshev function
while (my($n, $c2) = each (%chebyshev2)) {
  cmp_closeto( chebyshev_psi(0+$n), $c2, 1e-9*abs($n), "chebyshev_psi($n)" );
}

###### Ramanujan Sum
{
  is( ramanujan_sum(0, 34), 0, "Ramanujan Sum  c_0(34) = 0" );
  is( ramanujan_sum(34, 0), 0, "Ramanujan Sum  c_34(0)" );
  # A 30x30 grid of c_k(n)
  my @expect = (qw/1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 -1 2 -1 -1 2 -1 -1 2 -1 -1 2 -1 -1 2 -1 -1 2 -1 -1 2 -1 -1 2 -1 -1 2 -1 -1 2 0 -2 0 2 0 -2 0 2 0 -2 0 2 0 -2 0 2 0 -2 0 2 0 -2 0 2 0 -2 0 2 0 -2 -1 -1 -1 -1 4 -1 -1 -1 -1 4 -1 -1 -1 -1 4 -1 -1 -1 -1 4 -1 -1 -1 -1 4 -1 -1 -1 -1 4 1 -1 -2 -1 1 2 1 -1 -2 -1 1 2 1 -1 -2 -1 1 2 1 -1 -2 -1 1 2 1 -1 -2 -1 1 2 -1 -1 -1 -1 -1 -1 6 -1 -1 -1 -1 -1 -1 6 -1 -1 -1 -1 -1 -1 6 -1 -1 -1 -1 -1 -1 6 -1 -1 0 0 0 -4 0 0 0 4 0 0 0 -4 0 0 0 4 0 0 0 -4 0 0 0 4 0 0 0 -4 0 0 0 0 -3 0 0 -3 0 0 6 0 0 -3 0 0 -3 0 0 6 0 0 -3 0 0 -3 0 0 6 0 0 -3 1 -1 1 -1 -4 -1 1 -1 1 4 1 -1 1 -1 -4 -1 1 -1 1 4 1 -1 1 -1 -4 -1 1 -1 1 4 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 10 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 10 -1 -1 -1 -1 -1 -1 -1 -1 0 2 0 -2 0 -4 0 -2 0 2 0 4 0 2 0 -2 0 -4 0 -2 0 2 0 4 0 2 0 -2 0 -4 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 12 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 12 -1 -1 -1 -1 1 -1 1 -1 1 -1 -6 -1 1 -1 1 -1 1 6 1 -1 1 -1 1 -1 -6 -1 1 -1 1 -1 1 6 1 -1 1 1 -2 1 -4 -2 1 1 -2 -4 1 -2 1 1 8 1 1 -2 1 -4 -2 1 1 -2 -4 1 -2 1 1 8 0 0 0 0 0 0 0 -8 0 0 0 0 0 0 0 8 0 0 0 0 0 0 0 -8 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 16 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 3 0 0 -3 0 0 -6 0 0 -3 0 0 3 0 0 6 0 0 3 0 0 -3 0 0 -6 0 0 -3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 18 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 2 0 -2 0 2 0 -2 0 -8 0 -2 0 2 0 -2 0 2 0 8 0 2 0 -2 0 2 0 -2 0 -8 1 1 -2 1 1 -2 -6 1 -2 1 1 -2 1 -6 -2 1 1 -2 1 1 12 1 1 -2 1 1 -2 -6 1 -2 1 -1 1 -1 1 -1 1 -1 1 -1 -10 -1 1 -1 1 -1 1 -1 1 -1 1 10 1 -1 1 -1 1 -1 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 22 -1 -1 -1 -1 -1 -1 -1 0 0 0 4 0 0 0 -4 0 0 0 -8 0 0 0 -4 0 0 0 4 0 0 0 8 0 0 0 4 0 0 0 0 0 0 -5 0 0 0 0 -5 0 0 0 0 -5 0 0 0 0 -5 0 0 0 0 20 0 0 0 0 -5 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 -12 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 12 1 -1 1 -1 0 0 0 0 0 0 0 0 -9 0 0 0 0 0 0 0 0 -9 0 0 0 0 0 0 0 0 18 0 0 0 0 2 0 -2 0 2 0 -2 0 2 0 -2 0 -12 0 -2 0 2 0 -2 0 2 0 -2 0 2 0 12 0 2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 28 -1 -1 1 2 1 4 -2 -1 1 2 -4 -1 -2 -1 1 -8 1 -1 -2 -1 -4 2 1 -1 -2 4 1 2 1 -1 8/);
  my @got;
  for my $k (1..30) {
    for my $n (1..30) {
      push @got, ramanujan_sum($k, $n);
    }
  }
  is_deeply( \@got, \@expect, "Ramanujan sum c_{1..30}(1..30)" );
}

###### Hurwitz Class Number
while (my($n, $h) = each (%hclassno)) {
  is( hclassno(0 + $n), $h, "H($n) = $h" );
}

###### Ramanujan Tau
while (my($n, $tau) = each (%rtau)) {
  is( ramanujan_tau(0 + $n), $tau, "Ramanujan Tau($n) = $tau" );
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
###### gcdext
foreach my $garg (@gcdexts) {
  my($aref, $eref) = @$garg;
  my($x,$y) = @$aref;
  is_deeply( [gcdext($x,$y)], $eref, "gcdext($x,$y) = [@$eref]" );
}
###### chinese
foreach my $carg (@crts) {
  my($aref, $exp) = @$carg;
  my $crt = chinese(@$aref);
  is( $crt, $exp, "crt(".join(",",map { "[@$_]" } @$aref).") = " . ((defined $exp) ? $exp : "<undef>") );
}
###### znorder
foreach my $moarg (@mult_orders) {
  my ($a, $n, $exp) = @$moarg;
  my $zn = znorder($a, $n);
  is( $zn, $exp, "znorder($a, $n) = " . ((defined $exp) ? $exp : "<undef>") );
}
###### znprimroot
while (my($n, $root) = each (%primroots)) {
  is( znprimroot(0+$n), $root, "znprimroot($n) == " . ((defined $root) ? $root : "<undef>") );
}
is( znprimroot("-100000898"), 31, "znprimroot(\"-100000898\") == 31" );
# I don't think we should rely on this parsing correctly.
#is( znprimroot("+100000898"), 31, "znprimroot(\"+100000898\") == 31" );

###### is_primitive_root
while (my($n, $root) = each (%primroots)) {
  if (defined $root) {
    is( is_primitive_root(0+$root,0+$n), 1, "$root is a primitive root mod $n" );
  } else {
    is( is_primitive_root(2,0+$n), 0, "2 is not a primitive root mod $n" );
  }
}
is(is_primitive_root(19,191), 1, "19 is a primitive root mod 191");
is(is_primitive_root(13,191), 0, "13 is not a primitive root mod 191");

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
foreach my $e (0 .. $#negpowers) {
  is( is_power(-7 ** $e), $negpowers[$e], "is_power(-7^$e ) = $negpowers[$e]" );
}
{
  my($ispow, $root);
  $ispow = is_power(24, 2, \$root);
  is( $ispow, 0, "24 isn't a perfect square...");
  is( $root, undef, "...and the root wasn't set");
  $ispow = is_power( "1000093002883029791", 3, \$root);
  is( $ispow, 1, "1000031^3 is a perfect cube...");
  is( $root, 1000031, "...and the root was set");
  $ispow = is_power( 36**5 , 0, \$root);
  is( $ispow, 10, "36^5 is a 10th power...");
  is( $root, 6, "...and the root is 6");
}
###### valuation
foreach my $r (@valuations) {
  my($n, $k, $exp) = @$r;
  is( valuation($n, $k), $exp, "valuation($n,$k) = $exp" );
}
###### hammingweight
foreach my $r (@popcounts) {
  my($n, $exp) = @$r;
  is( hammingweight($n), $exp, "hammingweight($n) = $exp" );
}
###### sqrtint
is_deeply( [map { sqrtint($_) } 0..100], [map { int(sqrt($_)) } 0..100], "sqrtint 0 .. 100" );
is( sqrtint(1524155677489), 1234567, "sqrtint(1234567^2) = 1234567" );
is( sqrtint(1524158146623), 1234567, "sqrtint(1234568^2-1) = 1234567" );
is( sqrtint(1524155677488), 1234566, "sqrtint(1234567^2-1) = 1234566" );
###### rootint
# TODO: croak if n < 0 or k < 1
is(rootint(928342398,1), 928342398, "rootint(928342398,1) returns 928342398");
is(rootint(88875,3), 44, "rootint(88875,3) returns 44");
is(rootint("266667176579895999",3), 643659, "integer third root of 266667176579895999 is 643659");
{
  my(@got, @expected);
  for my $arr (@roots) {
    my($b, $k, $n) = @$arr;
    push @expected, [$b,$n];
    my $rk;
    my $r = rootint($n,$k,\$rk);
    push @got, [$r,$rk];
  }
  is_deeply( \@got, \@expected, "rootint on perfect powers where log fails" );
}
is( rootint("43091031920942300256108314560009772304748698124094750326895058640841523270081624169128280918534127523222564290447104831706207227117677890695945149868732770531628297914633063561406978145215542597509491443634033203125",23), 2147483645, "integer 23rd root of a large 23rd power" );
is( rootint("43091031920942300256108314560009772304748698124094750326895058640841523270081624169128280918534127523222564290447104831706207227117677890695945149868732770531628297914633063561406978145215542597509491443634033203124",23), 2147483644, "integer 23rd root of almost a large 23rd power" );
###### logint
is_deeply( [map { logint($_,2) } 1..200], [map { int(log($_)/log(2)+1e-10) } 1..200], "logint base 2: 0 .. 200" );
is_deeply( [map { logint($_,3) } 1..200], [map { int(log($_)/log(3)+1e-10) } 1..200], "logint base 3: 0 .. 200" );
is_deeply( [map { logint($_,5) } 1..200], [map { int(log($_)/log(5)+1e-10) } 1..200], "logint base 3: 0 .. 200" );
{
  my $be;
  is( logint(19284098234,16,\$be), 8, "logint(19284098234,16) = 8" );
  is( $be, 16**8, "power is 16^8" );
}
###### binomial
foreach my $r (@binomials) {
  my($n, $k, $exp) = @$r;
  is( binomial($n,$k), $exp, "binomial($n,$k)) = $exp" );
}
is_deeply( [map { binomial(10, $_) } -15 .. 15],
           [qw/0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 10 45 120 210 252 210 120 45 10 1 0 0 0 0 0/],
           "binomial(10,n) for n in -15 .. 15" );
is_deeply( [map { binomial(-10, $_) } -15 .. 15],
           [qw/-2002 715 -220 55 -10 1 0 0 0 0 0 0 0 0 0 1 -10 55 -220 715 -2002 5005 -11440 24310 -48620 92378 -167960 293930 -497420 817190 -1307504/],
           "binomial(-10,n) for n in -15 .. 15" );

sub cmp_closeto {
  my $got = shift;
  my $expect = shift;
  my $tolerance = shift;
  my $message = shift;
  cmp_ok( abs($got - $expect), '<=', $tolerance, $message );
}
