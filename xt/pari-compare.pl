#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util qw/:all/;
use Math::PariInit qw( primes=10000000 stack=1e8 );
use Math::Pari;
use Data::Dumper;
use Test::Differences;
$|=1;

# For these tests we're looking at small inputs.

my $loops = 0;
while (1) {

{ my $n = (int(rand(2**32)) << 32) + int(rand(2**32));
  die "isprime($n)" unless Math::Pari::isprime($n) == !!is_prime($n);
  die "is_prob_prime($n)" unless Math::Pari::isprime($n) == !!is_prob_prime($n);
  die "next_prime($n)" unless Math::Pari::nextprime($n+1) == next_prime($n);
  die "prev_prime($n)" unless Math::Pari::precprime($n-1) == prev_prime($n);

  my($pn,$pc) = @{Math::Pari::factorint($n)};
  my @f1 = map { [ int("$pn->[$_]"), int("$pc->[$_]")] } 0 .. $#$pn;
  array_compare( \@f1, [factor_exp($n)], "factor_exp($n)" );
  @f1 = map { ($_->[0]) x $_->[1] } @f1;
  array_compare( \@f1, [factor($n)], "factor($n)" );

  my @d1 = map { int("$_") } @{Math::Pari::divisors($n)};
  array_compare( \@d1, [divisors($n)], "divisors($n)" );

  die "omega($n)" unless Math::Pari::omega($n) == factor_exp($n);
  die "bigomega($n)" unless Math::Pari::bigomega($n) == factor($n);
  die "numdiv($n)" unless Math::Pari::numdiv($n) == divisors($n);

  die "moebius($n)" unless Math::Pari::moebius($n) == moebius($n);
  die "euler_phi($n)" unless Math::Pari::eulerphi($n) == euler_phi($n);

  my $d = PARI "d";
  # TODO: our jordan_totient should auto-bigint
  die "jordan_totient(2,$n)"
    unless Math::Pari::sumdiv($n,"d","d^2*moebius($n/d)")
        == jordan_totient(2,Math::BigInt->new("$n"));
  die "jordan_totient(3,$n)"
    unless Math::Pari::sumdiv($n,"d","d^3*moebius($n/d)")
        == jordan_totient(3,Math::BigInt->new("$n"));

  # TODO: exp_mangoldt:
  # Lambda(n)={ 
  #   v=factor(n);
  #   if(matsize(v)[1]!=1,return(0),return(log(v[1,1])));
  # };
  # Chebyshev Psi(x)=sum(n=2,floor(x),Lambda(n));

  # TODO: partitions.  new Pari has this as numbpart.
  #       See OEIS A000041 for some alternate Pari functions

# TODO:
#      chebyshev_theta chebyshev_psi
#      divisor_sum
#      carmichael_lambda znorder znprimroot znlog legendre_phi
#      ExponentialIntegral LogarithmicIntegral RiemannZeta RiemannR


}

{ my $n = 1+int(rand(10**5));  # Pari doesn't like 0
  die "nth_prime($n)" unless Math::Pari::prime($n) == nth_prime($n);
}

{ my $n = int(rand(10**4));
  array_compare( [map { int("$_") } @{Math::Pari::primes($n)}], primes(nth_prime($n)), "primes($n)" );
}

{ my $a = int(rand(10**6));
  my $b = $a+int(rand(10**4));
  my $x = PARI "x";
  my @a1; Math::Pari::forprime($x,$a,$b,sub { push @a1, int("$x") });
  my @a2; forprimes { push @a2, $_ } $a,$b;
  array_compare( \@a1, \@a2, "forprimes($a,$b)" );
}

# forcomposites in Pari 2.6, not Math::Pari's 2.1

{ my $n = int(rand(10**12));
  my $d = PARI "d";
  my @a1; Math::Pari::fordiv($n, $d, sub { push @a1, int("$d") });
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
}
{ my $a = int(rand(2**32));
  my $b = int(rand(2**32));
  die "lcm($a,$b)" unless Math::Pari::lcm($a,$b) == lcm($a,$b);
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
