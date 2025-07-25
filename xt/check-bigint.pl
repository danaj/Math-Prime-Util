#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/:all/;
use Time::HiRes qw(gettimeofday tv_interval);
use feature 'say';
use Math::GMPz;
use Math::GMP;
use Math::BigInt;
$| = 1;  # fast pipes

my $class = 'Math::GMPz';
prime_set_config(bigint=>$class);

my($rnat, $rbig);
my @S;

my $r;

# There is no particular reason for these numbers, they are just convenient
# bigints that are comfortably larger than 2^64 but not too much larger.

my $x = $class->new("1180591620717411303424");              # powint(2,70);
my $y = $class->new("2503155504993241601315571986085849");  # powint(3,70);
my $z = $class->new("1104427674243920646305299201");        # powint(7,32);

print "Operations should produce native ints (UV/IV) if possible.\n";
print "They should produce bigint results of the requested class if large.\n";
print "The intention is to check return types, not the functions themselves.\n";

################################################################################

print "\nbasic  ";

$rnat = sqrtint(77777777);
$rbig = sqrtint($x * $x + 31);
die "sqrtint" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 8819 && $rbig == $x;
print ".";

$rnat = addint(15, 31);
$rbig = addint($x, 31);
die "addint" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 46 && $rbig == $x+31;
print ".";

$rnat = subint(15, 31);
$rbig = subint($x, 31);
die "subint" unless !ref($rnat) && ref($rbig) eq $class && $rnat == -16 && $rbig == $x-31;
print ".";

$rnat = add1int(15);
$rbig = add1int($x);
die "add1int" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 16 && $rbig == $x+1;
print ".";

$rnat = sub1int(15);
$rbig = sub1int($x);
die "sub1int" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 14 && $rbig == $x-1;
print ".";

$rnat = mulint(15, 31);
$rbig = mulint($x, 31);
die "mulint" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 465 && $rbig == $x*31;
print ".";

$rnat = divint(310, 15);
$rbig = divint($x, 31);
die "divint" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 20 && $rbig == $class->new("38083600668303590433");
print ".";

$rnat = cdivint(310, 15);
$rbig = cdivint($x, 31);
die "divint" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 21 && $rbig == $class->new("38083600668303590434");
print ".";

$rnat = modint(310, 15);
die "modint" unless !ref($rnat) && $rnat == 10;
$rnat = modint($x, 31);
die "modint" unless !ref($rnat) && $rnat == 1;
$rbig = modint($y,$x);
die "modint" unless ref($rbig) eq $class && $rbig == $class->new("297440975187654227929");;
print ".";

$rnat = powint(310, 3);
$rbig = powint(317, 31);
die "powint" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 29791000 && $rbig == $class->new("341064979770708619981328007297127930890658843588985670937159050598503087466133");
print ".";

$rnat = absint(-10001);
$rbig = absint(-$x);
die "absint" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 10001 && $rbig == $x;
print ".";

$rnat = negint(10001);
$rbig = negint($x);
die "negint" unless !ref($rnat) && ref($rbig) eq $class && $rnat == -10001 && $rbig == -$x;
print ".";

$rnat = lshiftint(10001);
$rbig = lshiftint($x);
die "lshiftint" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 20002 && $rbig == 2*$x;
print ".";

$rnat = rshiftint($x,66);
$rbig = rshiftint($x,3);
die "rshiftint" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 16 && $rbig == $class->new("147573952589676412928");
print ".";

$rnat = rashiftint(-$y,82);
$rbig = rashiftint(-$y,20);
die "rashiftint" unless !ref($rnat) && ref($rbig) eq $class && $rnat == -517640426 && $rbig == $class->new("-2387195115082971192660877215");
print ".";

$rnat = logint($y, 4);
die "logint" unless !ref($rnat) && $rnat == 55;

$rnat = rootint($y,4);
$rbig = rootint($y*$y*$y-1,3);
die "rootint" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 223677323 && $rbig == $y-1;
print ".";

################################################################################

print "\nmodint ";

$rnat = addmod(15, 31, 17);
$rbig = addmod($x, $y, $z);
die "addmod" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 12 && $rbig == $class->new("867780633942779001401200");
print ".";

$rnat = submod(15, 31, 17);
$rbig = submod($x, $y, $z);
die "submod" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 1 && $rbig == $class->new("1103562254793219302126504849");
print ".";

$rnat = mulmod(15, 31, 17);
$rbig = mulmod($x, $y, $z);
die "mulmod" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 6 && $rbig == $class->new("749232057443429182743012873");
print ".";

$rnat = powmod(15, 31, 17);
$rbig = powmod($x, $y, $z);
die "powmod" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 8 && $rbig == $class->new("573707816919163066477349622");
print ".";

$rnat = divmod(15, 31, 17);
$rbig = divmod($x, $y, $z);
die "divmod" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 12 && $rbig == $class->new("814177443435415951196839160");
print ".";

$rnat = muladdmod(15, 31, 2, 17);
$rbig = muladdmod($x, $y, 2, $z);
die "muladdmod" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 8 && $rbig == $class->new("749232057443429182743012875");
print ".";

$rnat = mulsubmod(15, 31, 2, 17);
$rbig = mulsubmod($x, $y, 2, $z);
die "mulsubmod" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 4 && $rbig == $class->new("749232057443429182743012871");
print ".";

$rnat = factorialmod(15, 17);
$rbig = factorialmod(1000, $z+1);
die "factorialmod" unless !ref($rnat) && $rnat == 1;
die "factorialmod" unless ref($rbig) eq $class && $rbig == $class->new("13829182657193587281465406");
print ".";

$rnat = sqrtmod(310, 17);
$rbig = sqrtmod($y, $z);
die "sqrtmod" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 2 && $rbig == $class->new("50031545098999707");
print ".";

# TODO divrem, tdivrem, fdivrem, cdivrem

$rnat = lucasumod(4,-3,2379,377);
$rbig = lucasumod(4,-3,2379,$x);
die "lucasumod" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 68 && $rbig == $class->new("705649089465838257763");
print ".";

$rnat = lucasvmod(4,-3,2379,377);
$rbig = lucasvmod(4,-3,2379,$x);
die "lucasumod" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 105 && $rbig == $class->new("751898682006794398852");
print ".";

{
  my($r1,$r2);
  ($r1,$r2) = lucasuvmod(31,-17,88891112, 570);
  die "lucasuvmod" unless !ref($r1) && !ref($r2);
  die "lucasuvmod" unless $r1 == 535 && $r2 == 257;
  ($r1,$r2) = lucasuvmod(31,-17,88891112, $x);
  die "lucasuvmod" unless ref($r1) eq $class && ref($r2) eq $class;
  die "lucasuvmod" unless $r1 == $class->new("585481810409681695659") && $r2 == $class->new("1102837250424259606255");
  print ".";
}

################################################################################

################################################################################

# TODO
#  next_prime  prev_prime

# urandomb urandomm random_nbig_prime random_ndigit_prime
# random_strong_prime random_safe_prime random_prime
# random_maurer_prime random_shawe_taylor_prime

# consecutive_integer_lcm
# partitions
# gcd lcm gcdext
# chinese ramanujan_tau
# exp_mangoldt
# jordan_totient
# carmichael_lambda
# binomial
# stirling
# primorial pn_primorial
# permtonum
# subfactorial falling_factorial rising_factorial
# pisano_period powersum fromdigits
# lucsasu lucasv lucasuv

# powerful_count powerfree_count prime_power_count perfect_power_count
# nth_powerfree nth_perfect_power nth_perfect_power_approx
# next_perfect_power prev_perfect_power
# bernfrac harmfrac

# ... more ...
# many _approx, _lower, _upper


################################################################################

print "\nsemipr ";

$r = is_semiprime(165);
die "is_semiprime $r" unless $r == 0 && !ref($r);
print ".";
$r = is_semiprime(166);
die "is_semiprime $r" unless $r == 1 && !ref($r);
print ".";
$r = is_semiprime($class->new("27273137616939507011"));
die "is_semiprime $r" unless $r == 1 && !ref($r);
print ".";

$r = random_semiprime(32);
die "random_semiprime 32" unless is_semiprime($r) && !ref($r);
print ".";
$r = random_semiprime(65);
die "random_semiprime 65" unless is_semiprime($r) && ref($r) eq $class;
print ".";

$r = random_unrestricted_semiprime(32);
die "random_unrestricted_semiprime 32" unless is_semiprime($r) && !ref($r);
print ".";
$r = random_unrestricted_semiprime(65);
die "random_unrestricted_semiprime 65" unless is_semiprime($r) && ref($r) eq $class;
print ".";

$r = semi_primes(1000000000,1000000010);
die "semiprimes small" unless $r->[0] == 1000000006 && !ref($r->[0]);
print ".";

$r = semi_primes($x,$x+10);
die "semiprimes large" unless $r->[0] == $x+9 && ref($r->[0]) eq $class;
print ".";

$r = semiprime_count(1000000000);
die "semiprime_count" unless $r == 160788536 && !ref($r);
print ".";

$r = semiprime_count_approx($x);
die "semiprime_count_approx" unless ref($r) eq $class;
print ".";

$r = nth_semiprime(1000000);
die "nth_semiprime" unless $r == 5109839 && !ref($r);
print ".";

$r = nth_semiprime_approx($x);
die "nth_semiprime_approx" unless ref($r) eq $class;
print ".";

forsemiprimes { push @S, $_; } 3000000000,3000000005;
die "forsemiprimes small" unless vecequal(\@S,[3000000001,3000000002,3000000005]) && vecall { !ref($_) } @S;
print ".";

@S=();
forsemiprimes { push @S, $_; } $x,$x+10;
die "forsemiprimes large" unless vecequal(\@S,["1180591620717411303433"]) && vecall { ref($_) eq $class } @S;
print ".";

################################################################################


print "\nPASS\n";
