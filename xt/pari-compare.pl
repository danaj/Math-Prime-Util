#!/usr/bin/env perl
use strict;
use warnings;
use Math::PariInit qw( primes=10000000 stack=1e8 );
use Math::Pari qw/pari2iv/;
use Math::Prime::Util qw/:all/;
use Data::Dumper;
$|=1;

BEGIN {
  use Config;
  die "Tests have 64-bit assumptions" if $Config{uvsize} < 8;
  die "Tests need double floats" if $Config{nvsize} < 8;
  no Config;
}

my $small = 80_000;
print "Comparing for small inputs: 0 - $small\n";

foreach my $n (0 .. $small) {
  print '.' unless ($n+1) % int($small/80);
  die "isprime($n)" unless Math::Pari::isprime($n) == !!is_prime($n);
  die "is_prob_prime($n)" unless Math::Pari::isprime($n) == !!is_prob_prime($n);
  die "next_prime($n)" unless Math::Pari::nextprime($n+1) == next_prime($n);
  die "prev_prime($n)" unless Math::Pari::precprime($n-1) == prev_prime($n);
  next if $n == 0;

  my($pn,$pc) = @{Math::Pari::factorint($n)};
  my @f1 = map { [ pari2iv($pn->[$_]), pari2iv($pc->[$_])] } 0 .. $#$pn;
  array_compare( \@f1, [factor_exp($n)], "factor_exp($n)" );
  @f1 = map { ($_->[0]) x $_->[1] } @f1;
  array_compare( \@f1, [factor($n)], "factor($n)" );

  array_compare( [map { pari2iv($_) } @{Math::Pari::divisors($n)}], [divisors($n)], "divisors($n)" );

  die "omega($n)" unless Math::Pari::omega($n) == factor_exp($n);
  die "bigomega($n)" unless Math::Pari::bigomega($n) == factor($n);
  die "numdiv($n)" unless Math::Pari::numdiv($n) == divisors($n);

  for my $k (2,3,9,10) {
    die "valuation($n,$k)" unless Math::Pari::valuation($n,$k) == valuation($n,$k);
  }

  foreach my $k (0..4) {
    die "sigma($n,$k)" unless Math::Pari::sigma($n,$k) == divisor_sum($n,$k);
  }

  die "moebius($n)" unless Math::Pari::moebius($n) == moebius($n);
  die "euler_phi($n)" unless Math::Pari::eulerphi($n) == euler_phi($n);

  my $d = PARI "d";
  die "jordan_totient(2,$n)"
    unless Math::Pari::sumdiv($n,"d","d^2*moebius($n/d)")
        == jordan_totient(2,$n);
  die "jordan_totient(3,$n)"
    unless Math::Pari::sumdiv($n,"d","d^3*moebius($n/d)")
        == jordan_totient(3,$n);

  if ($n > 1) {
    for (1..10) {
      my $k;  do { $k = int(rand(50)) } while !($k % $n);
      die "binomial($n,$k)" unless Math::Pari::binomial($n,$k) == binomial($n,$k);
      my $negn = - ($n >> 1);
      die "binomial($negn,$k)" unless Math::Pari::binomial($negn,$k) == binomial($negn,$k);
    }
  }

  {
    my $d = $n+3;
    my @gmpu  = gcdext($n,$d);
    my $gpari = Math::Pari::bezout($n,$d);
    die "gcdext($n,$d)" unless $gmpu[0] == $gpari->[0] && $gmpu[1] == $gpari->[1] && $gmpu[2] == $gpari->[2];
  }

  die "nth_prime($n)" unless Math::Pari::prime($n) == nth_prime($n);

  # All the pari2iv calls are very time-consuming
  if ($n < 1000) {
    array_compare( [map { pari2iv($_) } @{Math::Pari::primes($n)}], primes(nth_prime($n)), "primes($n)" );
  }

  # Math Pari's forprime is super slow for some reason. Pari/gp isn't this slow.
  if ($n < 1000) {
    my $m = $n+int(rand(10**4));
    PARI "s1=0";
    PARI "forprime(X=$n,$m,s1=s1+X)";
    my $s1 = PARI('s1');

    my $s2 = 0; forprimes { $s2 += $_ } $n,$m;
    die "forprimes($n,$m)  $s1 != $s2" unless $s1 == $s2;
  }
  {
    my $d = PARI "d";
    my @a1; Math::Pari::fordiv($n, $d, sub { push @a1, pari2iv($d)});
    my @a2; fordivisors { push @a2, $_ } $n;
    array_compare( \@a1, \@a2, "fordivisors($n)" );
  }

  { my $m = int(rand($n-1));
    my $invmod = invmod($m, $n);
    if (defined $invmod) {
      die "invmod($m, $n)" unless Math::Pari::lift(PARI "Mod(1/$m,$n)") == $invmod;
    } else {
      eval { PARI "Mod(1/$m,$n)" };
      die "invmod($m, $n) defined in Pari" unless $@ =~ /impossible inverse/
        || ($m == 0 && $@ =~ /division by zero/);
    }
  }

  { my $m = int(rand($n-1));
    my $mn = PARI "Mod($m,$n)";
    my $order = znorder($m, $n);
    if (defined $order) {
      die "znorder($m, $n)" unless Math::Pari::znorder($mn) == $order
    } else {
      eval { Math::Pari::znorder($mn); };
      die "znorder($m, $n) defined in Pari" unless $@ =~ /not an element/;
    }
  }

  # Pari's znprimroot is iffy for non-primes
  if (is_prime($n)) {
    my $g = znprimroot($n);
    die "znprimroot($n)" unless Math::Pari::znprimroot($n) == $g;
    my $a = 1 + int(rand($n-2));
    my $gn = PARI "Mod($g,$n)";
    my $log = znlog($a, $g, $n);
    die "znlog($a, $g, $n) should be defined" unless defined $log;
    die "znlog($a, $g, $n)" unless Math::Pari::znlog($a,$gn) == $log;
  }

  if ($n < 100) {
    foreach my $d (0 .. 9) {
      my $arg = $n + $d/10;
      next if $arg < 0.1;
      my $e1 = -Math::Pari::eint1(-$arg);
      my $e2 = ExponentialIntegral($arg);
      die "ExponentialIntegral($arg)  $e1 != $e2" if abs($e1 - $e2) > $e1*1e-14;
    }
  }
  if ($n > 1) {
    my $arg = $n;
    my $e1 = -Math::Pari::eint1(-log($arg));
    my $e2 = LogarithmicIntegral($arg);
    die "LogarithmicIntegral($arg)  $e1 != $e2" if abs($e1 - $e2) > $e1*1e-14;
  }

  {
    my $s = 50.0/$small;
    if ($s != 1.0) {
      my $zeta1 = Math::Pari::zeta($s) - 1;
      my $zeta2 = RiemannZeta($s);
      die "zeta($s) $zeta1 != $zeta2" if abs($zeta1 - $zeta2) > abs($zeta1) * 1e-14;
    }
  }

  #print "." unless $n % 1250;
}

print "\nkronecker, gcd, and lcm for small values\n";
foreach my $a (-400 .. 400) {
  foreach my $b (-400 .. 400) {
    # Pari 2.1's gcd doesn't work right for 0,-x and -x,0.  Pari 2.2.3 fixed.
    if ($a != 0 && $b != 0) {
      die "gcd($a,$b)" unless Math::Pari::gcd($a,$b) == gcd($a,$b);
    }
    die "kronecker($a,$b)" unless Math::Pari::kronecker($a,$b) == kronecker($a,$b);
    die "lcm($a,$b)" unless Math::Pari::lcm($a,$b) == lcm($a,$b);
  }
  print "." unless (400+$a) % 20;
}

print "\nloop forever with random values\n";

# forcomposites in Pari 2.6, not Math::Pari's 2.1

my $loops = 0;
while (1) {
  my $n;

  {
  do { $n = (int(rand(2**32)) << 32) + int(rand(2**32)) } while $n < $small;
  die "isprime($n)" unless Math::Pari::isprime($n) == !!is_prime($n);
  die "is_prob_prime($n)" unless Math::Pari::isprime($n) == !!is_prob_prime($n);
  die "next_prime($n)" unless Math::Pari::nextprime($n+1) == next_prime($n);
  die "prev_prime($n)" unless Math::Pari::precprime($n-1) == prev_prime($n);

  my($pn,$pc) = @{Math::Pari::factorint($n)};
  my @f1 = map { [ pari2iv($pn->[$_]), pari2iv($pc->[$_])] } 0 .. $#$pn;
  array_compare( \@f1, [factor_exp($n)], "factor_exp($n)" );
  @f1 = map { ($_->[0]) x $_->[1] } @f1;
  array_compare( \@f1, [factor($n)], "factor($n)" );

  array_compare( [map { pari2iv($_) } @{Math::Pari::divisors($n)}], [divisors($n)], "divisors($n)" );

  die "omega($n)" unless Math::Pari::omega($n) == factor_exp($n);
  die "bigomega($n)" unless Math::Pari::bigomega($n) == factor($n);
  die "numdiv($n)" unless Math::Pari::numdiv($n) == divisors($n);

  for my $k (2,3,9,10) {
    die "valuation($n,$k)" unless Math::Pari::valuation($n,$k) == valuation($n,$k);
  }

  foreach my $k (0..4) {
    die "sigma($n,$k)" unless Math::Pari::sigma($n,$k) == divisor_sum($n,$k);
  }

  die "moebius($n)" unless Math::Pari::moebius($n) == moebius($n);
  die "euler_phi($n)" unless Math::Pari::eulerphi($n) == euler_phi($n);

  my $d = PARI "d";
  # TODO: our jordan_totient should auto-bigint
  die "jordan_totient(2,$n)"
    unless Math::Pari::sumdiv($n,"d","d^2*moebius($n/d)")
        == jordan_totient(2,$n);
  die "jordan_totient(3,$n)"
    unless Math::Pari::sumdiv($n,"d","d^3*moebius($n/d)")
        == jordan_totient(3,$n);

  if ($n > 2) {
    for (1..10) {
      my $k;  do { $k = int(rand(10)) } while !($k % $n);
      die "binomial($n,$k)" unless Math::Pari::binomial($n,$k) == binomial($n,$k);
      my $negn = - ($n >> 1);
      die "binomial($negn,$k)" unless Math::Pari::binomial($negn,$k) == binomial($negn,$k);
    }
  }

  # TODO: exp_mangoldt:
  # Lambda(n)={
  #   v=factor(n);
  #   if(matsize(v)[1]!=1,return(0),return(log(v[1,1])));
  # };
  # TODO: chebyshev_theta, chebyshev_psi
  # Chebyshev Psi(x)=sum(n=2,floor(x),Lambda(n));

  # TODO: partitions.  new Pari has this as numbpart.
  #       See OEIS A000041 for some alternate Pari functions

  # TODO: primorial / pn_primorial

  # TODO: carmichael lambda?  Pari doesn't have it.

  { my $m = int(rand($n-1));
    my $invmod = invmod($m, $n);
    if (defined $invmod) {
      die "invmod($m, $n)" unless Math::Pari::lift(PARI "Mod(1/$m,$n)") == $invmod;
    } else {
      eval { PARI "Mod(1/$m,$n)" };
      die "invmod($m, $n) defined in Pari" unless $@ =~ /impossible inverse/
        || ($m == 0 && $@ =~ /division by zero/);
    }
  }

  { my $m = int(rand($n-1));
    my $mn = PARI "Mod($m,$n)";
    my $order = znorder($m, $n);
    if (defined $order) {
      die "znorder($m, $n)" unless Math::Pari::znorder($mn) == $order;
    } else {
      eval { Math::Pari::znorder($mn); };
      die "znorder($m, $n) defined in Pari" unless $@ =~ /not an element/;
    }
  }

  # TODO: znlog with reasonable values

  if ($n > 1) {
    my $arg = $n;
    my $e1 = -Math::Pari::eint1(-log($arg));
    my $e2 = LogarithmicIntegral($arg);
    die "LogarithmicIntegral($arg)  $e1 != $e2" if abs($e1 - $e2) > $e1*1e-12;
  }
  # TODO: RiemannZeta
  }



{ my $a = $small + int(rand(10**6));
  my $b = $a+int(rand(10**4));
  my $x = PARI "x";
  my @a1; Math::Pari::forprime($x,$a,$b,sub { push @a1, pari2iv($x) });
  my @a2; forprimes { push @a2, $_ } $a,$b;
  array_compare( \@a1, \@a2, "forprimes($a,$b)" );
}

# forcomposites in Pari 2.6, not Math::Pari's 2.1

{ my $n = $small + int(rand(10**12));
  my $d = PARI "d";
  my @a1; Math::Pari::fordiv($n, $d, sub { push @a1, pari2iv($d) });
  my @a2; fordivisors { push @a2, $_ } $n;
  array_compare( \@a1, \@a2, "fordivisors($n)" );
}

# Pari's primepi in 2.1-2.5 is strangely lacking

{ my $a = (int(rand(2**32)) << 32) + int(rand(2**32));
  my $b = (int(rand(2**32)) << 32) + int(rand(2**32));
  die "gcd($a,$b)" unless Math::Pari::gcd($a,$b) == gcd($a,$b);

  die "kronecker($a,$b)" unless Math::Pari::kronecker($a,$b) == kronecker($a,$b);
  $a >>= 1 if $a > 2**63;
  die "kronecker(-$a,$b)" unless Math::Pari::kronecker(-$a,$b) == kronecker(-$a,$b);
  $b >>= 1 if $b > 2**63;
  die "kronecker($a,-$b)" unless Math::Pari::kronecker($a,-$b) == kronecker($a,-$b);
  die "kronecker(-$a,-$b)" unless Math::Pari::kronecker(-$a,-$b) == kronecker(-$a,-$b);

  { my @gmpu  = gcdext($a,$b);   my $gpari = Math::Pari::bezout($a,$b);
    die "gcdext($a,$b)" unless $gmpu[0] == $gpari->[0] && $gmpu[1] == $gpari->[1] && $gmpu[2] == $gpari->[2]; }
}
{ my $a = int(rand(2**32));
  my $b = int(rand(2**32));
  die "lcm($a,$b)" unless Math::Pari::lcm($a,$b) == lcm($a,$b);
}

{ my $n = random_prime(10000,~0);
  die "znprimroot($n)" unless Math::Pari::znprimroot($n) == znprimroot($n);
}

$loops++;
print "." unless $loops % 100;
}




use Bytes::Random::Secure qw/random_string_from/;
sub ndigit_rand {
  my($digits, $howmany) = @_;
  die "digits must be > 0" if $digits < 1;
  $howmany = 1 unless defined $howmany;
  my @nums = map { random_string_from("123456789",1) . random_string_from("0123456789",$digits-1) } 1 .. $howmany;
  if (10**$digits > ~0) {  @nums = map { Math::BigInt->new($_) } @nums;  }
  else                  {  @nums = map { int($_) } @nums;                }
  return wantarray ? @nums : $nums[0];
}


sub array_compare {
  my($a1, $a2, $text) = @_;
  #eq_or_diff $a1, $a2, $text;
  die "$text wrong count ",scalar @$a1," ",scalar @$a2 unless @$a1 == @$a2;
  foreach my $i (0 .. $#$a1) {
    if (ref($a1->[$i])) {
      array_compare($a1->[$i],$a2->[$i], "> $text");
    } else {
    #print "a1: ", Dumper($a1), "\na2: ", Dumper($a2), "\n" unless $a1->[$i] == $a2->[$i];
      die "$text entry $i  $a1->[$i] != $a2->[$i]" unless $a1->[$i] == $a2->[$i];
    }
  }
}
