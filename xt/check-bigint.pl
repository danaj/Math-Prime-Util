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
die "sqrtint native" unless checkv($rnat,8819);
die "sqrtint bigint" unless checkv($rbig,$x);
print ".";

$rnat = addint(15, 31);
$rbig = addint($x, 31);
die "addint" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 46 && $rbig == $x+31;
die "addint native" unless checkv($rnat,46);
die "addint bigint" unless checkv($rbig,$x+31);
print ".";

$rnat = subint(15, 31);
$rbig = subint($x, 31);
die "subint native" unless checkv($rnat,-16);
die "subint bigint" unless checkv($rbig,$x-31);
print ".";

$rnat = add1int(15);
$rbig = add1int($x);
die "add1int native" unless checkv($rnat,16);
die "add1int bigint" unless checkv($rbig,$x+1);
print ".";

$rnat = sub1int(15);
$rbig = sub1int($x);
die "sub1int" unless !ref($rnat) && ref($rbig) eq $class && $rnat == 14 && $rbig == $x-1;
die "sub1int native" unless checkv($rnat,14);
die "sub1int bigint" unless checkv($rbig,$x-1);
print ".";

$rnat = mulint(15, 31);
$rbig = mulint($x, 31);
die "mulint native" unless checkv($rnat,465);
die "mulint bigint" unless checkv($rbig,$x*31);
print ".";

$rnat = divint(310, 15);
$rbig = divint($x, 31);
die "divint native" unless checkv($rnat,20);
die "divint bigint" unless checkv($rbig,"38083600668303590433");
print ".";

$rnat = cdivint(310, 15);
$rbig = cdivint($x, 31);
die "cdivint native" unless checkv($rnat,21);
die "cdivint bigint" unless checkv($rbig,"38083600668303590434");
print ".";

$rnat = modint(310, 15);
die "modint native" unless checkv($rnat,10);
$rnat = modint($x, 31);
die "modint native" unless checkv($rnat,1);
$rbig = modint($y,$x);
die "modint bigint" unless checkv($rbig,"297440975187654227929");
print ".";

$rnat = powint(310, 3);
$rbig = powint(317, 31);
die "powint native" unless checkv($rnat,29791000);
die "powint bigint" unless checkv($rbig,"341064979770708619981328007297127930890658843588985670937159050598503087466133");
print ".";

$rnat = absint(-10001);
$rbig = absint(-$x);
die "absint native" unless checkv($rnat,10001);
die "absint bigint" unless checkv($rbig,$x);
print ".";

$rnat = negint(10001);
$rbig = negint($x);
die "negint" unless !ref($rnat) && ref($rbig) eq $class && $rnat == -10001 && $rbig == -$x;
die "negint native" unless checkv($rnat,-10001);
die "negint bigint" unless checkv($rbig,-$x);
print ".";

$rnat = lshiftint(10001);
$rbig = lshiftint($x);
die "lshiftint native" unless checkv($rnat,20002);
die "lshiftint bigint" unless checkv($rbig,2*$x);
print ".";

$rnat = rshiftint($x,66);
$rbig = rshiftint($x,3);
die "rshiftint native" unless checkv($rnat,16);
die "rshiftint bigint" unless checkv($rbig,"147573952589676412928");
print ".";

$rnat = rashiftint(-$y,82);
$rbig = rashiftint(-$y,20);
die "rashiftint native" unless checkv($rnat,-517640426);
die "rashiftint bigint" unless checkv($rbig,"-2387195115082971192660877215");
print ".";

$rnat = logint($y, 4);
die "logint" unless checkv($rnat,55);

$rnat = rootint($y,4);
$rbig = rootint($y*$y*$y-1,3);
die "rootint native" unless checkv($rnat,223677323);
die "rootint bigint" unless checkv($rbig,$y-1);
print ".";

################################################################################

print "\nmodint ";

$rnat = addmod(15, 31, 17);
$rbig = addmod($x, $y, $z);
die "addmod native" unless checkv($rnat,12);
die "addmod bigint" unless checkv($rbig,"867780633942779001401200");
print ".";

$rnat = submod(15, 31, 17);
$rbig = submod($x, $y, $z);
die "submod native" unless checkv($rnat,1);
die "submod bigint" unless checkv($rbig,"1103562254793219302126504849");
print ".";

$rnat = mulmod(15, 31, 17);
$rbig = mulmod($x, $y, $z);
die "mulmod native" unless checkv($rnat,6);
die "mulmod bigint" unless checkv($rbig,"749232057443429182743012873");
print ".";

$rnat = powmod(15, 31, 17);
$rbig = powmod($x, $y, $z);
die "powmod native" unless checkv($rnat,8);
die "powmod bigint" unless checkv($rbig,"573707816919163066477349622");
print ".";

$rnat = divmod(15, 31, 17);
$rbig = divmod($x, $y, $z);
die "divmod native" unless checkv($rnat,12);
die "divmod bigint" unless checkv($rbig,"814177443435415951196839160");
print ".";

$rnat = muladdmod(15, 31, 2, 17);
$rbig = muladdmod($x, $y, 2, $z);
die "muladdmod native" unless checkv($rnat,8);
die "muladdmod bigint" unless checkv($rbig,"749232057443429182743012875");
print ".";

$rnat = mulsubmod(15, 31, 2, 17);
$rbig = mulsubmod($x, $y, 2, $z);
die "mulsubmod native" unless checkv($rnat,4);
die "mulsubmod bigint" unless checkv($rbig,"749232057443429182743012871");
print ".";

$rnat = factorialmod(15, 17);
$rbig = factorialmod(1000, $z+1);
die "factorialmod native" unless checkv($rnat,1);
die "factorialmod bigint" unless checkv($rbig,"13829182657193587281465406");
print ".";

$rnat = sqrtmod(310, 17);
$rbig = sqrtmod($y, $z);
die "sqrtmod native" unless checkv($rnat,2);
die "sqrtmod bigint" unless checkv($rbig,"50031545098999707");
print ".";

# divrem (Euclidean): remainder always >= 0
{
  my($q,$r);
  ($q,$r) = divrem(310, 17);
  die "divrem native" unless checkv($q,18) && checkv($r,4);
  ($q,$r) = divrem(-310, 17);
  die "divrem native neg" unless checkv($q,-19) && checkv($r,13);
  ($q,$r) = divrem($y, $x);
  die "divrem bigint" unless checkv($q,2120255184830) && checkv($r,"297440975187654227929");
  ($q,$r) = divrem(-$y, $x);
  die "divrem bigint neg" unless checkv($q,-2120255184831) && checkv($r,"883150645529757075495");
  print ".";
}

# tdivrem (truncated toward zero)
{
  my($q,$r);
  ($q,$r) = tdivrem(310, 17);
  die "tdivrem native" unless checkv($q,18) && checkv($r,4);
  ($q,$r) = tdivrem(-310, 17);
  die "tdivrem native neg" unless checkv($q,-18) && checkv($r,-4);
  ($q,$r) = tdivrem($y, $x);
  die "tdivrem bigint" unless checkv($q,2120255184830) && checkv($r,"297440975187654227929");
  ($q,$r) = tdivrem(-$y, $x);
  die "tdivrem bigint neg" unless checkv($q,-2120255184830) && checkv($r,"-297440975187654227929");
  print ".";
}

# fdivrem (floored)
{
  my($q,$r);
  ($q,$r) = fdivrem(310, 17);
  die "fdivrem native" unless checkv($q,18) && checkv($r,4);
  ($q,$r) = fdivrem(-310, 17);
  die "fdivrem native neg" unless checkv($q,-19) && checkv($r,13);
  ($q,$r) = fdivrem($y, $x);
  die "fdivrem bigint" unless checkv($q,2120255184830) && checkv($r,"297440975187654227929");
  ($q,$r) = fdivrem(-$y, $x);
  die "fdivrem bigint neg" unless checkv($q,-2120255184831) && checkv($r,"883150645529757075495");
  print ".";
}

# cdivrem (ceiling)
{
  my($q,$r);
  ($q,$r) = cdivrem(310, 17);
  die "cdivrem native" unless checkv($q,19) && checkv($r,-13);
  ($q,$r) = cdivrem(-310, 17);
  die "cdivrem native neg" unless checkv($q,-18) && checkv($r,-4);
  ($q,$r) = cdivrem($y, $x);
  die "cdivrem bigint" unless checkv($q,2120255184831) && checkv($r,"-883150645529757075495");
  ($q,$r) = cdivrem(-$y, $x);
  die "cdivrem bigint neg" unless checkv($q,-2120255184830) && checkv($r,"-297440975187654227929");
  print ".";
}

$rnat = lucasumod(4,-3,2379,377);
$rbig = lucasumod(4,-3,2379,$x);
die "lucasumod native" unless checkv($rnat,68);
die "lucasumod bigint" unless checkv($rbig,"705649089465838257763");
print ".";

$rnat = lucasvmod(4,-3,2379,377);
$rbig = lucasvmod(4,-3,2379,$x);
die "lucasvmod native" unless checkv($rnat,105);
die "lucasvmod bigint" unless checkv($rbig,"751898682006794398852");
print ".";

{
  my($r1,$r2);
  ($r1,$r2) = lucasuvmod(31,-17,88891112, 570);
  die "lucasuvmod" unless checkv($r1,535) && checkv($r2,257);
  ($r1,$r2) = lucasuvmod(31,-17,88891112, $x);
  die "lucasuvmod" unless checkv($r1,"585481810409681695659") && checkv($r2,"1102837250424259606255");
  print ".";
}

# binomialmod
# negmod
# invmod
# rootmod
# allsqrtmod
# allrootmod

################################################################################

print "\nnthry  ";

# gcd: small result from small inputs, big result from big inputs
$rnat = gcd(15, 35);
$rbig = gcd(mulint($z,6), mulint($z,10));
die "gcd native" unless checkv($rnat, 5);
die "gcd bigint" unless checkv($rbig, "2208855348487841292610598402");
print ".";

# lcm
$rnat = lcm(15, 35);
$rbig = lcm(mulint($z,6), mulint($z,10));
die "lcm native" unless checkv($rnat, 105);
die "lcm bigint" unless checkv($rbig, "33132830227317619389158976030");
print ".";

# gcdext
{
  my @ge = gcdext(15, 35);
  die "gcdext native" unless checkv($ge[0],-2) && checkv($ge[1],1) && checkv($ge[2],5);
  @ge = gcdext(mulint($z,6), mulint($z,10));
  die "gcdext bigint" unless checkv($ge[0],2) && checkv($ge[1],-1) && checkv($ge[2],"2208855348487841292610598402");
  print ".";
}

# euler_phi:  phi(100) = 40,  phi(7^32) = 7^31 * 6
$rnat = euler_phi(100);
$rbig = euler_phi($z);
die "euler_phi native" unless checkv($rnat, 40);
die "euler_phi bigint" unless checkv($rbig, "946652292209074839690256458");
print ".";

# carmichael_lambda:  lambda(100) = 20,  lambda(7^32) = 7^31 * 6
$rnat = carmichael_lambda(100);
$rbig = carmichael_lambda($z);
die "carmichael_lambda native" unless checkv($rnat, 20);
die "carmichael_lambda bigint" unless checkv($rbig, "946652292209074839690256458");
print ".";

# jordan_totient:  J_2(100) = 7200,  J_2(7^32) = 7^62 * 48
$rnat = jordan_totient(2, 100);
$rbig = jordan_totient(2, $z);
die "jordan_totient native" unless checkv($rnat, 7200);
die "jordan_totient bigint" unless checkv($rbig, "1194867416459594155237786640878013212168766002414274352");
print ".";

# kronecker:  always returns native -1, 0, or 1
$rnat = kronecker(15, 35);
die "kronecker" unless checkv($rnat, 0);
$rnat = kronecker($x, $z);
die "kronecker bigint args" unless checkv($rnat, 1);
print ".";

# znorder:  znorder(2,101) = 100
$rnat = znorder(2, 101);
die "znorder native" unless checkv($rnat, 100);
# znorder with bigint prime modulus: znorder(2, next_prime(2^70))
$rbig = znorder(2, $class->new("1180591620717411303449"));
die "znorder bigint" unless checkv($rbig, "295147905179352825862");
print ".";

# znprimroot:  znprimroot(101) = 2
$rnat = znprimroot(101);
die "znprimroot native" unless checkv($rnat, 2);
$rnat = znprimroot($class->new("1180591620717411303449"));
die "znprimroot bigint" unless checkv($rnat, 3);
print ".";

# exp_mangoldt:  prime power -> prime,  otherwise -> 1
$rnat = exp_mangoldt(128);  # 2^7 -> 2
die "exp_mangoldt native" unless checkv($rnat, 2);
$rnat = exp_mangoldt($z);   # 7^32 -> 7
die "exp_mangoldt bigint" unless checkv($rnat, 7);
$rnat = exp_mangoldt(30);   # not a prime power -> 1
die "exp_mangoldt 30" unless checkv($rnat, 1);
print ".";

# consecutive_integer_lcm:  result can exceed UV
$rnat = consecutive_integer_lcm(20);
$rbig = consecutive_integer_lcm(100);
die "consecutive_integer_lcm native" unless checkv($rnat, 232792560);
die "consecutive_integer_lcm bigint" unless checkv($rbig, "69720375229712477164533808935312303556800");
print ".";

# ramanujan_tau
$rnat = ramanujan_tau(20);
$rbig = ramanujan_tau(1000);
die "ramanujan_tau native" unless checkv($rnat, -7109760);
die "ramanujan_tau bigint" unless checkv($rbig, "-30328412970240000");
print ".";

# partitions
$rnat = partitions(30);
$rbig = partitions(1000);
die "partitions native" unless checkv($rnat, 5604);
die "partitions bigint" unless checkv($rbig, "24061467864032622473692149727991");
print ".";

# chinese (CRT)
$rnat = chinese([3,7],[5,11]);
$rbig = chinese([3,$x],[5,$z]);
die "chinese native" unless checkv($rnat, 38);
die "chinese bigint" unless checkv($rbig, "274780012208942403819037834081377580846524923907");
print ".";

# binomial
$rnat = binomial(20, 10);
$rbig = binomial(100, 50);
die "binomial native" unless checkv($rnat, 184756);
die "binomial bigint" unless checkv($rbig, "100891344545564193334812497256");
$rbig = binomial($x, 3);
die "binomial bigint n" unless checkv($rbig, "274250759553534340358464632138771002190616713273449955064807424");
print ".";

# stirling (type 1 and type 2)
$rnat = stirling(10, 3, 1);
$rbig = stirling(30, 3, 1);
die "stirling s1 native" unless checkv($rnat, -1172700);
die "stirling s1 bigint" unless checkv($rbig, "-62262192842035613491057459200000");
$rnat = stirling(10, 3, 2);
$rbig = stirling(50, 3, 2);
die "stirling s2 native" unless checkv($rnat, 9330);
die "stirling s2 bigint" unless checkv($rbig, "119649664052358811373730");
print ".";

# permtonum
$rnat = permtonum([2,1,0]);
$rbig = permtonum([24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0]);
die "permtonum native" unless checkv($rnat, 5);
die "permtonum bigint" unless checkv($rbig, "15511210043330985983999999");
print ".";

# subfactorial
$rnat = subfactorial(10);
$rbig = subfactorial(50);
die "subfactorial native" unless checkv($rnat, 1334961);
die "subfactorial bigint" unless checkv($rbig, "11188719610782480504630258070757734324011354208865721592720336801");
print ".";

# falling_factorial
$rnat = falling_factorial(10, 3);
$rbig = falling_factorial($x, 3);
die "falling_factorial native" unless checkv($rnat, 720);
die "falling_factorial bigint" unless checkv($rbig, "1645504557321206042150787792832626013143700279640699730388844544");
print ".";

# rising_factorial
$rnat = rising_factorial(10, 3);
$rbig = rising_factorial($x, 3);
die "rising_factorial native" unless checkv($rnat, 1320);
die "rising_factorial bigint" unless checkv($rbig, "1645504557321206042159150572282074996821776173992942865953587200");
print ".";

# pisano_period
$rnat = pisano_period(100);
$rbig = pisano_period($x);
die "pisano_period native" unless checkv($rnat, 300);
die "pisano_period bigint" unless checkv($rbig, "1770887431076116955136");
print ".";

# powersum
$rnat = powersum(10, 2);
$rbig = powersum(100, 10);
die "powersum native" unless checkv($rnat, 385);
die "powersum bigint" unless checkv($rbig, "959924142434241924250");
print ".";

# lucasu:  U(1,-1,k) = Fibonacci(k)
$rnat = lucasu(1, -1, 20);
$rbig = lucasu(1, -1, 100);
die "lucasu native" unless checkv($rnat, 6765);
die "lucasu bigint" unless checkv($rbig, "354224848179261915075");
print ".";

# lucasv:  V(1,-1,k) = Lucas(k)
$rnat = lucasv(1, -1, 20);
$rbig = lucasv(1, -1, 100);
die "lucasv native" unless checkv($rnat, 15127);
die "lucasv bigint" unless checkv($rbig, "792070839848372253127");
print ".";

# lucasuv
{
  my($u,$v);
  ($u,$v) = lucasuv(1, -1, 20);
  die "lucasuv native" unless checkv($u, 6765) && checkv($v, 15127);
  ($u,$v) = lucasuv(1, -1, 100);
  die "lucasuv bigint" unless checkv($u, "354224848179261915075") && checkv($v, "792070839848372253127");
  print ".";
}

# TODO
#  next_prime  prev_prime

# urandomb urandomm random_nbig_prime random_ndigit_prime
# random_strong_prime random_safe_prime random_prime
# random_maurer_prime random_shawe_taylor_prime

# primorial pn_primorial
# fromdigits

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


sub checkv {
  my($n,$v)=@_;

  # This is great except it assumes addint always works correctly.
  # That should be true, but this test is trying to verify it.
  #$v = addint($v,0);  # assuming addint works, v is the right type
  #return 0 unless (!ref($n) && !ref($v)) || (ref($v) && ref($n) eq $class);
  #return $n == $v;

  # This assumes cmpint will work correctly.  It doesn't depend on return types.
  if (cmpint($v,~0) <= 0 && cmpint($v,negint(1+(~0 >> 1))) >= 0) {
    return 0 unless !ref($n) && $n == $v;
  } else {
    $v = $class->new($v) unless ref($v);
    return 0 unless ref($n) eq $class && $n == $v;
  }
  1;
}
