package Math::Prime::Util::PP;
use strict;
use warnings;
use Carp qw/carp croak confess/;

BEGIN {
  $Math::Prime::Util::PP::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::PP::VERSION = '0.08';
}

# The Pure Perl versions of all the Math::Prime::Util routines.
#
# Some of these will be relatively similar in performance, some will be
# horrendously slow in comparison.
#
# TODO: These are currently all terribly simple.

my $_uv_size;
BEGIN {
  use Config;
  $_uv_size =
   (   (defined $Config{'use64bitint'} && $Config{'use64bitint'} eq 'define')
    || (defined $Config{'use64bitall'} && $Config{'use64bitall'} eq 'define')
    || (defined $Config{'longsize'} && $Config{'longsize'} >= 8)
   )
   ? 64
   : 32;
  no Config;
}
sub _maxbits { $_uv_size }

my $_precalc_size = 0;
sub prime_precalc {
  my($n) = @_;
  croak "Input must be a positive integer" unless is_positive_int($n);
  $_precalc_size = $n if $n > $_precalc_size;
}
sub prime_memfree {
  $_precalc_size = 0;
}
sub _get_prime_cache_size { $_precalc_size }
sub _prime_memfreeall { prime_memfree; }


sub is_positive_int {
  ((defined $_[0]) && ($_[0] !~ tr/0123456789//c));
}

my @_primes_small = (
   0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
   193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
   293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,
   409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499);
my @_prime_count_small = (
   0,0,1,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,8,8,8,8,9,9,9,9,9,9,10,10,
   11,11,11,11,11,11,12,12,12,12,13,13,14,14,14,14,15,15,15,15,15,15,
   16,16,16,16,16,16,17,17,18,18,18,18,18,18,19);
my @_prime_next_small = (
   2,2,3,5,5,7,7,11,11,11,11,13,13,17,17,17,17,19,19,23,23,23,23,
   29,29,29,29,29,29,31,31,37,37,37,37,37,37,41,41,41,41,43,43,47,
   47,47,47,53,53,53,53,53,53,59,59,59,59,59,59,61,61,67,67,67,67,67,67,71);

# For wheel-30
my @_prime_indices = (1, 7, 11, 13, 17, 19, 23, 29);
my @_nextwheel30 = (1,7,7,7,7,7,7,11,11,11,11,13,13,17,17,17,17,19,19,23,23,23,23,29,29,29,29,29,29,1);
my @_prevwheel30 = (29,29,1,1,1,1,1,1,7,7,7,7,11,11,13,13,13,13,17,17,19,19,19,19,23,23,23,23,23,23);

sub _is_prime7 {  # n must not be divisible by 2, 3, or 5
  my($x) = @_;
  my $q;
  foreach my $i (7, 11, 13, 17, 19, 23, 29) {
    $q = int($x/$i); return 1 if $q < $i; return 0 if $x == ($q*$i);
  }
  my $i = 31;  # mod-30 loop
  while (1) {
    $q = int($x/$i); return 1 if $q < $i; return 0 if $x == ($q*$i);  $i += 6;
    $q = int($x/$i); return 1 if $q < $i; return 0 if $x == ($q*$i);  $i += 4;
    $q = int($x/$i); return 1 if $q < $i; return 0 if $x == ($q*$i);  $i += 2;
    $q = int($x/$i); return 1 if $q < $i; return 0 if $x == ($q*$i);  $i += 4;
    $q = int($x/$i); return 1 if $q < $i; return 0 if $x == ($q*$i);  $i += 2;
    $q = int($x/$i); return 1 if $q < $i; return 0 if $x == ($q*$i);  $i += 4;
    $q = int($x/$i); return 1 if $q < $i; return 0 if $x == ($q*$i);  $i += 6;
    $q = int($x/$i); return 1 if $q < $i; return 0 if $x == ($q*$i);  $i += 2;
  }
  1;
}

sub is_prime {
  my($n) = @_;
  croak "Input must be an integer" unless defined $_[0];
  return 0 if $n < 2;        # 0 and 1 are composite
  return 2 if ($n == 2) || ($n == 3) || ($n == 5);  # 2, 3, 5 are prime
  # multiples of 2,3,5 are composite
  return 0 if (($n % 2) == 0) || (($n % 3) == 0) || (($n % 5) == 0);

  return is_prob_prime($n) if $n > 10_000_000;

  2*_is_prime7($n);
}

sub _sieve_erat {
  my($end) = @_;
  my $last = int(($end+1)/2);

  my $sieve = '';
  my $n = 3;
  while ( ($n*$n) <= $end ) {
    my $s = $n*$n;
    while ($s <= $end) {
      vec($sieve, $s >> 1, 1) = 1;
      $s += 2*$n;
    }
    do { $n += 2 } while vec($sieve, $n >> 1, 1) != 0;
  }
  # Mark 1 as composite
  vec($sieve, 0, 1) = 1;
  $sieve;
}
# Uses 8x more memory, but about 50% faster
sub _sieve_erat_array {
  my($end) = @_;
  my $last = int(($end+1)/2);

  my @sieve;
  my $n = 3;
  while ( ($n*$n) <= $end ) {
    my $s = $n*$n;
    while ($s <= $end) {
      $sieve[$s>>1] = 1;
      $s += 2*$n;
    }
    do { $n += 2 } while $sieve[$n>>1];
  }
  # Mark 1 as composite
  $sieve[0] = 1;
  \@sieve;
}

sub primes {
  my $optref = (ref $_[0] eq 'HASH')  ?  shift  :  {};
  croak "no parameters to primes" unless scalar @_ > 0;
  croak "too many parameters to primes" unless scalar @_ <= 2;
  my $low = (@_ == 2)  ?  shift  :  2;
  my $high = shift;
  my $sref = [];

  croak "Input must be a positive integer" unless is_positive_int($low)
                                               && is_positive_int($high);

  return $sref if ($low > $high) || ($high < 2);

  # Ignore method options in this code

  push @$sref, 2  if ($low <= 2) && ($high >= 2);
  push @$sref, 3  if ($low <= 3) && ($high >= 3);
  push @$sref, 5  if ($low <= 5) && ($high >= 5);
  $low = 7 if $low < 7;

  if ($low == 7) {
    my $sieve = _sieve_erat_array($high);
    my $n = 7;
    while ($n <= $high) {
      push @$sref, $n  if !$sieve->[$n>>1];
      $n += 2;
    }
  } else {
    my $n = $low;
    my $base = 30 * int($n/30);
    my $in = 0;  $in++ while ($n - $base) > $_prime_indices[$in];
    $n = $base + $_prime_indices[$in];
    while ($n <= $high) {
      push @$sref, $n  if _is_prime7($n);
      if (++$in == 8) {  $base += 30; $in = 0;  }
      $n = $base + $_prime_indices[$in];
    }
  }
  $sref;
}

sub next_prime {
  my($n) = @_;
  croak "Input must be a positive integer" unless is_positive_int($n);
  return 0 if $n >= ((_maxbits == 32) ? 4294967291 : 18446744073709551557);
  return $_prime_next_small[$n] if $n <= $#_prime_next_small;

  my $d = int($n/30);
  my $m = $n - $d*30;
  if ($m == 29) { $d++;  $m = 1;} else { $m = $_nextwheel30[$m]; }
  while (!_is_prime7($d*30+$m)) {
    $m = $_nextwheel30[$m];
    $d++ if $m == 1;
  }
  $d*30 + $m;
}

sub prev_prime {
  my($n) = @_;
  croak "Input must be a positive integer" unless is_positive_int($n);
  if ($n <= 7) {
    return ($n <= 2) ? 0 : ($n <= 3) ? 2 : ($n <= 5) ? 3 : 5;
  }

  $n++ if ($n % 2) == 0;
  do {
    $n -= 2;
  } while ( (($n % 3) == 0) || (!is_prime($n)) );
  $n;

  # This is faster for larger intervals, slower for short ones.
  #my $base = 30 * int($n/30);
  #my $in = 0;  $in++ while ($n - $base) > $_prime_indices[$in];
  #if (--$in < 0) {  $base -= 30; $in = 7;  }
  #$n = $base + $_prime_indices[$in];
  #while (!_is_prime7($n)) {
  #  if (--$in < 0) {  $base -= 30; $in = 7;  }
  #  $n = $base + $_prime_indices[$in];
  #}
  #$n;

  #my $d = int($n/30);
  #my $m = $n - $d*30;
  #do {
  #  $m = $_prevwheel30[$m];
  #  $d-- if $m == 29;
  #} while (!_is_prime7($d*30+$m));
  #$d*30+$m;
}

sub prime_count {
  my($low,$high) = @_;
  if (!defined $high) {
    $high = $low;
    $low = 2;
  }
  croak "Input must be a positive integer" unless is_positive_int($low)
                                               && is_positive_int($high);
  my $count = 0;

  # We should get sieved segments and count them.

  #my $d = $low;   # High can be outside iterator range
  #while ($d <= $high) {
  #  $count++ if is_prime($d);
  #  $d++;
  #}

  $count++ if ($low <= 2) && ($high >= 2);
  $count++ if ($low <= 3) && ($high >= 3);
  $count++ if ($low <= 5) && ($high >= 5);
  $low = 7 if $low < 7;

  if ($low == 7) {
    my $sieve = _sieve_erat_array($high);
    my $n = 7;
    while ($n <= $high) {
      $count++ if !$sieve->[$n>>1];
      $n += 4;
      if ($n <= $high) {
        $count++ if !$sieve->[$n>>1];
        $n += 2;
      }
    }
  } else {
    my $n = $low;
    my $base = 30 * int($n/30);
    my $in = 0;  $in++ while ($n - $base) > $_prime_indices[$in];
    $n = $base + $_prime_indices[$in];
    while ($n <= $high) {
      $count++ if _is_prime7($n);
      if (++$in == 8) {  $base += 30; $in = 0;  }
      $n = $base + $_prime_indices[$in];
    }
  }
  $count;
}

sub prime_count_lower {
  my($x) = @_;
  croak "Input must be a positive integer" unless is_positive_int($x);

  return $_prime_count_small[$x] if $x <= $#_prime_count_small;

  my $flogx = log($x);

  return int( $x / ($flogx - 0.7) ) if $x < 599;

  my $a;
  if    ($x <       2700) { $a = 0.30; }
  elsif ($x <       5500) { $a = 0.90; }
  elsif ($x <      19400) { $a = 1.30; }
  elsif ($x <      32299) { $a = 1.60; }
  elsif ($x <     176000) { $a = 1.80; }
  elsif ($x <     315000) { $a = 2.10; }
  elsif ($x <    1100000) { $a = 2.20; }
  elsif ($x <    4500000) { $a = 2.31; }
  elsif ($x <  233000000) { $a = 2.36; }
  elsif ($x < 5433800000) { $a = 2.32; }
  elsif ($x <60000000000) { $a = 2.15; }
  else                    { $a = 1.80; } # Dusart 1999, page 14

  return int( ($x/$flogx) * (1.0 + 1.0/$flogx + $a/($flogx*$flogx)) );
}

sub prime_count_upper {
  my($x) = @_;
  croak "Input must be a positive integer" unless is_positive_int($x);

  return $_prime_count_small[$x] if $x <= $#_prime_count_small;

  my $flogx = log($x);

  return int( ($x / ($flogx - 1.048)) + 1.0 ) if $x <  1621;
  return int( ($x / ($flogx - 1.071)) + 1.0 ) if $x <  5000;
  return int( ($x / ($flogx - 1.098)) + 1.0 ) if $x < 15900;

  my $a;
  if    ($x <      24000) { $a = 2.30; }
  elsif ($x <      59000) { $a = 2.48; }
  elsif ($x <     350000) { $a = 2.52; }
  elsif ($x <     355991) { $a = 2.54; }
  elsif ($x <     356000) { $a = 2.51; }
  elsif ($x <    3550000) { $a = 2.50; }
  elsif ($x <    3560000) { $a = 2.49; }
  elsif ($x <    5000000) { $a = 2.48; }
  elsif ($x <    8000000) { $a = 2.47; }
  elsif ($x <   13000000) { $a = 2.46; }
  elsif ($x <   18000000) { $a = 2.45; }
  elsif ($x <   31000000) { $a = 2.44; }
  elsif ($x <   41000000) { $a = 2.43; }
  elsif ($x <   48000000) { $a = 2.42; }
  elsif ($x <  119000000) { $a = 2.41; }
  elsif ($x <  182000000) { $a = 2.40; }
  elsif ($x <  192000000) { $a = 2.395; }
  elsif ($x <  213000000) { $a = 2.390; }
  elsif ($x <  271000000) { $a = 2.385; }
  elsif ($x <  322000000) { $a = 2.380; }
  elsif ($x <  400000000) { $a = 2.375; }
  elsif ($x <  510000000) { $a = 2.370; }
  elsif ($x <  682000000) { $a = 2.367; }
  elsif ($x <60000000000) { $a = 2.362; }
  else                    { $a = 2.51; }

  return int( ($x/$flogx) * (1.0 + 1.0/$flogx + $a/($flogx*$flogx)) + 1.0 );
}

sub prime_count_approx {
  my($x) = @_;
  croak "Input must be a positive integer" unless is_positive_int($x);

  return $_prime_count_small[$x] if $x <= $#_prime_count_small;

  #return int( (prime_count_upper($x) + prime_count_lower($x)) / 2);
  return int(RiemannR($x)+0.5);
}

sub nth_prime_lower {
  my($n) = @_;
  croak "Input must be a positive integer" unless is_positive_int($n);

  return $_primes_small[$n] if $n <= $#_primes_small;

  my $flogn  = log($n);
  my $flog2n = log($flogn);

  # Dusart 1999 page 14, for all n >= 2
  my $lower = $n * ($flogn + $flog2n - 1.0 + (($flog2n-2.25)/$flogn));

  if ($lower >= ~0) {
    if (_maxbits == 32) {
      return 4294967291 if $n <= 203280221;
    } else {
      return 18446744073709551557 if $n <= 425656284035217743;
    }
    croak "nth_prime_lower($n) overflow";
  }
  return int($lower);
}

sub nth_prime_upper {
  my($n) = @_;
  croak "Input must be a positive integer" unless is_positive_int($n);

  return $_primes_small[$n] if $n <= $#_primes_small;

  my $flogn  = log($n);
  my $flog2n = log($flogn);

  my $upper;
  if ($n >= 39017) {
    $upper = $n * ( $flogn  +  $flog2n - 0.9484 ); # Dusart 1999 page 14
  } elsif ($n >= 27076) {
    $upper = $n * ( $flogn  +  $flog2n - 1.0 + (($flog2n-1.80)/$flogn) );
  } elsif ($n >= 7022) {
    $upper = $n * ( $flogn  +  0.9385 * $flog2n ); # Robin 1983
  } else {
    $upper = $n * ( $flogn  +  $flog2n );
  }

  if ($upper >= ~0) {
    if (_maxbits == 32) {
      return 4294967291 if $n <= 203280221;
    } else {
      return 18446744073709551557 if $n <= 425656284035217743;
    }
    croak "nth_prime_upper($n) overflow";
  }
  return int($upper + 1.0);
}

sub nth_prime_approx {
  my($n) = @_;
  croak "Input must be a positive integer" unless is_positive_int($n);

  return $_primes_small[$n] if $n <= $#_primes_small;

  my $flogn  = log($n);
  my $flog2n = log($flogn);

  my $approx = $n * ( $flogn + $flog2n - 1
                      + (($flog2n - 2)/$flogn)
                      - ((($flog2n*$flog2n) - 6*$flog2n + 11) / (2*$flogn*$flogn))
                    );

  my $order = $flog2n/$flogn;
  $order = $order*$order*$order * $n;

  if    ($n <        259) { $approx += 10.4 * $order; }
  elsif ($n <        775) { $approx +=  7.52* $order; }
  elsif ($n <       1271) { $approx +=  5.6 * $order; }
  elsif ($n <       2000) { $approx +=  5.2 * $order; }
  elsif ($n <       4000) { $approx +=  4.3 * $order; }
  elsif ($n <      12000) { $approx +=  3.0 * $order; }
  elsif ($n <     150000) { $approx +=  2.1 * $order; }
  elsif ($n <  200000000) { $approx +=  0.0 * $order; }
  else                    { $approx += -0.010 * $order; }

  if ($approx >= ~0) {
    if (_maxbits == 32) {
      return 4294967291 if $n <= 203280221;
    } else {
      return 18446744073709551557 if $n <= 425656284035217743;
    }
    croak "nth_prime_approx($n) overflow";
  }

  return int($approx + 0.5);
}

sub nth_prime {
  my($n) = @_;
  croak "Input must be a positive integer" unless is_positive_int($n);

  return $_primes_small[$n] if $n <= $#_primes_small;

  if (_maxbits == 32) {
    croak "nth_prime($n) overflow" if $n > 203280221;
  } else {
    croak "nth_prime($n) overflow" if $n > 425656284035217743;
  }

  my $prime = 0;
  # This is quite slow, so try to skip some
  if    ($n >= 1000000) { $prime = 15485863; $n -= 1000000; }
  elsif ($n >=  500000) { $prime =  7368787; $n -=  500000; }
  elsif ($n >=  100000) { $prime =  1299709; $n -=  100000; }
  elsif ($n >=   50000) { $prime =   611953; $n -=   50000; }
  elsif ($n >=   10000) { $prime =   104729; $n -=   10000; }
  elsif ($n >=    1000) { $prime =     7919; $n -=    1000; }
  elsif ($n >=     100) { $prime =      541; $n -=     100; }

  while ($n-- > 0) {
    $prime = next_prime($prime);
  }
  $prime;
}

sub _mulmod {
  my($a, $b, $m) = @_;
  my $r = 0;
  while ($b > 0) {
    if ($b & 1) {
      if ($r == 0) {
        $r = $a;
      } else {
        $r = $m - $r;
        $r = ($a >= $r)  ?  $a - $r  :  $m - $r + $a;
      }
    }
    $a = ($a > ($m - $a))  ?  ($a - $m) + $a  :  $a + $a;
    $b >>= 1;
  }
  $r;
}

sub _powmod {
  my($n, $power, $m) = @_;
  my $t = 1;

  while ($power) {
    $t = _mulmod($t, $n, $m) if ($power & 1) != 0;
    $n = _mulmod($n, $n, $m);
    $power >>= 1;
  }
  $t;
}

sub miller_rabin {
  my($n, @bases) = @_;
  croak "Input must be a positive integer" unless is_positive_int($n);
  croak "No bases given to miller_rabin" unless @bases;

  return 0 if ($n == 0) || ($n == 1);
  return 2 if ($n == 2) || ($n == 3);

  # I was using bignum here for a while, but doing "$a ** $d" with a
  # big $d is **ridiculously** slow.  It's thousands of times faster
  # to do it ourselves using _powmod and _mulmod.

  my $s = 0;
  my $d = $n - 1;

  while ( ($d & 1) == 0 ) {
    $s++;
    $d >>= 1;
  }

  foreach my $a (@bases) {
    croak "Base $a is invalid" if $a < 2;
    my $x = _powmod($a, $d, $n);
    next if ($x == 1) || ($x == ($n-1));

    foreach my $r (1 .. $s) {
      $x = _mulmod($x, $x, $n);
      return 0 if $x == 1;
      if ($x == ($n-1)) {
        $a = 0;
        last;
      }
    }
    return 0 if $a != 0;
  }
  1;
}

sub is_prob_prime {
  my($n) = @_;

  return 2 if ($n == 2) || ($n == 3) || ($n == 5) || ($n == 7);
  return 0 if ($n < 2) || (($n % 2) == 0) || (($n % 3) == 0) || (($n % 5) == 0) || (($n % 7) == 0);
  return 2 if $n < 121;

  my @bases;
  if    ($n <          9080191) { @bases = (31, 73); }
  elsif ($n <       4759123141) { @bases = (2, 7, 61); }
  elsif ($n <     105936894253) { @bases = (2, 1005905886, 1340600841); }
  elsif ($n <   31858317218647) { @bases = (2, 642735, 553174392, 3046413974); }
  elsif ($n < 3071837692357849) { @bases = (2, 75088, 642735, 203659041, 3613982119); }
  else                          { @bases = (2, 325, 9375, 28178, 450775, 9780504, 1795265022); }

  # Run our selected set of Miller-Rabin strong probability tests
  my $prob_prime = miller_rabin($n, @bases);
  # We're deterministic for 64-bit numbers
  $prob_prime *= 2 if $n <= ~0; 
  $prob_prime;
}

sub trial_factor {
  my($n) = @_;
  croak "Input must be a positive integer" unless is_positive_int($n);

  return ($n) if $n < 4;

  my @factors;

  while ( ($n & 1) == 0) {
    push @factors, 2;
    $n >>= 1;
  }

  my $limit = int( sqrt($n) + 0.001);
  my $f = 3;
  while ($f <= $limit) {
    if ( ($n % $f) == 0) {
      while ( ($n % $f) == 0) {
        push @factors, $f;
        $n /= $f;
      }
      $limit = int( sqrt($n) + 0.001);
    }
    $f += 2;
  }
  push @factors, $n  if $n > 1;
  @factors;
}

# TODO:
sub factor { trial_factor(@_) }
sub fermat_factor { trial_factor(@_) }
sub holf_factor { trial_factor(@_) }
sub squfof_factor { trial_factor(@_) }
sub pbrent_factor { trial_factor(@_) }
sub prho_factor { trial_factor(@_) }
sub pminus1_factor { trial_factor(@_) }


my $_const_euler = 0.57721566490153286060651209008240243104215933593992;
my $_const_li2 = 1.045163780117492784844588889194613136522615578151;

sub ExponentialIntegral {
  my($x) = @_;
  my $tol = 1e-16;
  my $sum = 0.0;
  my($y, $t);
  my $c = 0.0;

  croak "Invalid input to ExponentialIntegral:  x must be != 0" if $x == 0;

  my $val; # The result from one of the four methods

  if ($x < -1) {
    # Continued fraction
    my $lc = 0;
    my $ld = 1 / (1 - $x);
    $val = $ld * (-exp($x));
    for my $n (1 .. 100000) {
      $lc = 1 / (2*$n + 1 - $x - $n*$n*$lc);
      $ld = 1 / (2*$n + 1 - $x - $n*$n*$ld);
      my $old = $val;
      $val *= $ld/$lc;
      last if abs($val - $old) <= ($tol * abs($val));
    }
  } elsif ($x < 0) {
    # Rational Chebyshev approximation
    my @C6p = ( -148151.02102575750838086,
                 150260.59476436982420737,
                  89904.972007457256553251,
                  15924.175980637303639884,
                   2150.0672908092918123209,
                    116.69552669734461083368,
                      5.0196785185439843791020);
    my @C6q = (  256664.93484897117319268,
                 184340.70063353677359298,
                  52440.529172056355429883,
                   8125.8035174768735759866,
                    750.43163907103936624165,
                     40.205465640027706061433,
                      1.0000000000000000000000);
    my $sumn = $C6p[0]-$x*($C6p[1]-$x*($C6p[2]-$x*($C6p[3]-$x*($C6p[4]-$x*($C6p[5]-$x*$C6p[6])))));
    my $sumd = $C6q[0]-$x*($C6q[1]-$x*($C6q[2]-$x*($C6q[3]-$x*($C6q[4]-$x*($C6q[5]-$x*$C6q[6])))));
    $val = log(-$x) - ($sumn / $sumd);
  } elsif ($x < -log($tol)) {
    # Convergent series
    my $fact_n = 1;
    $y = $_const_euler-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
    $y = log($x)-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
    for my $n (1 .. 200) {
      $fact_n *= $x/$n;
      my $term = $fact_n / $n;
      $y = $term-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
      last if $term < $tol;
    }
    $val = $sum;
  } else {
    # Asymptotic divergent series
    $val = exp($x) / $x;
    $y = 1.0-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
    my $term = 1.0;
    for my $n (1 .. 200) {
      my $last_term = $term;
      $term *= $n/$x;
      last if $term < $tol;
      if ($term < $last_term) {
        $y = $term-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
      } else {
        $y = (-$last_term/3)-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
        last;
      }
    }
    $val *= $sum;
  }
  $val;
}

sub LogarithmicIntegral {
  my($x) = @_;
  return 0 if $x == 0;
  return 0+(-Infinity) if $x == 1;
  return $_const_li2 if $x == 2;
  croak "Invalid input to LogarithmicIntegral:  x must be > 0" if $x <= 0;
  ExponentialIntegral(log($x));
}

# Riemann Zeta function for integers, used for computing Riemann R
my @_Riemann_Zeta_Table = (
  0.6449340668482264364724151666460251892,  # zeta(2)
  0.2020569031595942853997381615114499908,
  0.0823232337111381915160036965411679028,
  0.0369277551433699263313654864570341681,
  0.0173430619844491397145179297909205279,
  0.0083492773819228268397975498497967596,
  0.0040773561979443393786852385086524653,
  0.0020083928260822144178527692324120605,
  0.0009945751278180853371459589003190170,
  0.0004941886041194645587022825264699365,
  0.0002460865533080482986379980477396710,
  0.0001227133475784891467518365263573957,
  0.0000612481350587048292585451051353337,
  0.0000305882363070204935517285106450626,
  0.0000152822594086518717325714876367220,
  0.0000076371976378997622736002935630292,
  0.0000038172932649998398564616446219397,
  0.0000019082127165539389256569577951013,
  0.0000009539620338727961131520386834493,
  0.0000004769329867878064631167196043730,
  0.0000002384505027277329900036481867530,
  0.0000001192199259653110730677887188823,
  0.0000000596081890512594796124402079358,
  0.0000000298035035146522801860637050694,
  0.0000000149015548283650412346585066307,
  0.0000000074507117898354294919810041706,
  0.0000000037253340247884570548192040184,
  0.0000000018626597235130490064039099454,
  0.0000000009313274324196681828717647350,
  0.0000000004656629065033784072989233251,
  0.0000000002328311833676505492001455976,
  0.0000000001164155017270051977592973835,
  0.0000000000582077208790270088924368599,
  0.0000000000291038504449709968692942523,
  0.0000000000145519218910419842359296322,
  0.0000000000072759598350574810145208690,
  0.0000000000036379795473786511902372363,
  0.0000000000018189896503070659475848321,
  0.0000000000009094947840263889282533118,
);

# Compute Riemann Zeta function.  Slow and inaccurate near x = 1, but improves
# very rapidly (x = 5 is quite reasonable).
sub _evaluate_zeta {
  my($x) = @_;
  my $tol = 1e-16;
  my $sum = 0.0;
  my($y, $t);
  my $c = 0.0;

  $y = (1.0/(2.0**$x))-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;

  for my $k (3 .. 100000) {
    my $term = 1.0 / ($k ** $x);
    $y = $term-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
    last if abs($term) < $tol;
  }
  $sum;
}

# Riemann R function
sub RiemannR {
  my($x) = @_;
  my $tol = 1e-16;
  my $sum = 0.0;
  my($y, $t);
  my $c = 0.0;

  croak "Invalid input to ReimannR:  x must be > 0" if $x <= 0;

  $y = 1.0-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
  my $flogx = log($x);
  my $part_term = 1.0;

  for my $k (1 .. 10000) {
    # Small k from table, larger k from function
    my $zeta = ($k <= $#_Riemann_Zeta_Table) ? $_Riemann_Zeta_Table[$k+1-2]
                                             : _evaluate_zeta($k+1);
    $part_term *= $flogx / $k;
    my $term = $part_term / ($k + $k * $zeta);
    $y = $term-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
    last if abs($term) < $tol;
  }
  $sum;
}

1;

__END__


# ABSTRACT: Pure Perl version of Math::Prime::Util

=pod

=head1 NAME

Math::Prime::Util::PP - Pure Perl version of Math::Prime::Util


=head1 VERSION

Version 0.08


=head1 SYNOPSIS

  # TODO

=head1 DESCRIPTION

  # TODO

=head1 LIMITATIONS


=head1 PERFORMANCE


  # TODO

=head1 SEE ALSO

L<Math::Prime::Util>


=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>


=head1 COPYRIGHT

Copyright 2012 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
