package Math::Prime::Util::PP;
use strict;
use warnings;
use Carp qw/carp croak confess/;

BEGIN {
  $Math::Prime::Util::PP::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::PP::VERSION = '0.13';
}

# The Pure Perl versions of all the Math::Prime::Util routines.
#
# Some of these will be relatively similar in performance, some will be
# very slow in comparison.
#
# Most of these are pretty simple.  Also, you really should look at the C
# code for more detailed comments, including references to papers.

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
sub _PP_prime_maxbits { $_uv_size }

# If $n < $_half_word, then $n*$n will be exact.
my $_half_word = (~0 == 18446744073709551615) ? 4294967296 :    # 64-bit
                 (~0 ==           4294967295) ?      65536 :    # 32-bit
                 (~0 ==                   -1) ?   1000**10 :    # bignum
                                                         0 ;    # No idea

my $_precalc_size = 0;
sub prime_precalc {
  my($n) = @_;
  croak "Input must be a positive integer" unless _is_positive_int($n);
  $_precalc_size = $n if $n > $_precalc_size;
}
sub prime_memfree {
  $_precalc_size = 0;
}
sub _get_prime_cache_size { $_precalc_size }
sub _prime_memfreeall { prime_memfree; }


sub _is_positive_int {
  ((defined $_[0]) && ($_[0] !~ tr/0123456789//c));
}

sub _validate_positive_integer {
  my($n, $min, $max) = @_;
  croak "Parameter must be defined" if !defined $n;
  croak "Parameter '$n' must be a positive integer" if $n =~ tr/0123456789//c;
  croak "Parameter '$n' must be >= $min" if defined $min && $n < $min;
  croak "Parameter '$n' must be <= $max" if defined $max && $n > $max;
  if ($n <= ~0) {
    $_[0] = $_[0]->as_number() if ref($_[0]) eq 'Math::BigFloat';
    $_[0] = int($_[0]->bstr) if ref($_[0]) eq 'Math::BigInt';
  } elsif (ref($n) ne 'Math::BigInt') {
    croak "Parameter '$n' outside of integer range" if !defined $bigint::VERSION;
    $_[0] = Math::BigInt->new("$n"); # Make $n a proper bigint object
  }
  # One of these will be true:
  #     1) $n <= max and $n is not a bigint
  #     2) $n  > max and $n is a bigint
  1;
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
  my($n) = @_;

  if ($n < 61*61) {
    foreach my $i (qw/7 11 13 17 19 23 29 31 37 41 43 47 53 59/) {
      return 2 if $i*$i > $n;
      return 0 if !($n % $i);
    }
    return 2;
  }

  return 0 if !($n %  7) || !($n % 11) || !($n % 13) || !($n % 17) ||
              !($n % 19) || !($n % 23) || !($n % 29) || !($n % 31) ||
              !($n % 37) || !($n % 41) || !($n % 43) || !($n % 47) ||
              !($n % 53) || !($n % 59);

  return Math::Prime::Util::is_prob_prime($n) if $n > 10_000_000;

  my $limit = int(sqrt($n));
  my $i = 61;
  while (($i+30) <= $limit) {
    return 0 if !($n % $i);  $i += 6;
    return 0 if !($n % $i);  $i += 4;
    return 0 if !($n % $i);  $i += 2;
    return 0 if !($n % $i);  $i += 4;
    return 0 if !($n % $i);  $i += 2;
    return 0 if !($n % $i);  $i += 4;
    return 0 if !($n % $i);  $i += 6;
    return 0 if !($n % $i);  $i += 2;
  }
  while (1) {
    last if $i > $limit;  return 0 if !($n % $i);  $i += 6;
    last if $i > $limit;  return 0 if !($n % $i);  $i += 4;
    last if $i > $limit;  return 0 if !($n % $i);  $i += 2;
    last if $i > $limit;  return 0 if !($n % $i);  $i += 4;
    last if $i > $limit;  return 0 if !($n % $i);  $i += 2;
    last if $i > $limit;  return 0 if !($n % $i);  $i += 4;
    last if $i > $limit;  return 0 if !($n % $i);  $i += 6;
    last if $i > $limit;  return 0 if !($n % $i);  $i += 2;
  }
  2;
}

sub is_prime {
  my($n) = @_;
  _validate_positive_integer($n);

  return 2 if ($n == 2) || ($n == 3) || ($n == 5);  # 2, 3, 5 are prime
  return 0 if $n < 7;             # everything else below 7 is composite
                                  # multiples of 2,3,5 are composite
  return 0 if !($n % 2) || !($n % 3) || !($n % 5);
  return _is_prime7($n);
}

# Possible sieve storage:
#   1) vec with mod-30 wheel:   8 bits  / 30
#   2) vec with mod-2 wheel :  15 bits  / 30
#   3) str with mod-30 wheel:   8 bytes / 30
#   4) str with mod-2 wheel :  15 bytes / 30
#
# It looks like using vecs is about 2x slower than strs, and the strings also
# let us do some fast operations on the results.  E.g.
#   Count all primes:
#      $count += $$sieveref =~ tr/0//;
#   Loop over primes:
#      foreach my $s (split("0", $$sieveref, -1)) {
#        $n += 2 + 2 * length($s);
#        .. do something with the prime $n
#      }
#
# We're using method 4, though sadly it is memory intensive relative to the
# other methods.  I will point out that it is 30-60x less memory than sieves
# using an array, and the performance of this function is over 10x that
# of naive sieves like found on RosettaCode.

sub _sieve_erat_string {
  my($end) = @_;

  # Prefill with 3 and 5 already marked.
  my $whole = int( ($end>>1) / 15);
  croak "Sieve too large" if $whole > 1_145_324_612;  # ~32 GB string
  my $sieve = "100010010010110" . "011010010010110" x $whole;
  # Make exactly the number of entries requested, never more.
  substr($sieve, ($end>>1)+1) = '';
  my $n = 7;
  while ( ($n*$n) <= $end ) {
    my $s = $n*$n;
    my $filter_s   = $s >> 1;
    my $filter_end = $end >> 1;
    while ($filter_s <= $filter_end) {
      substr($sieve, $filter_s, 1) = "1";
      $filter_s += $n;
    }
    do { $n += 2 } while substr($sieve, $n>>1, 1);
  }
  \$sieve;
}

# TODO: this should be plugged into precalc, memfree, etc. just like the C code
{
  my $primary_size_limit = 15000;
  my $primary_sieve_size = 0;
  my $primary_sieve_ref;
  sub _sieve_erat {
    my($end) = @_;

    return _sieve_erat_string($end) if $end > $primary_size_limit;

    if ($primary_sieve_size == 0) {
      $primary_sieve_size = $primary_size_limit;
      $primary_sieve_ref = _sieve_erat_string($primary_sieve_size);
    }
    my $sieve = substr($$primary_sieve_ref, 0, ($end+1)>>1);
    \$sieve;
  }
}


sub _sieve_segment {
  my($beg,$end) = @_;
  croak "Internal error: segment beg is even" if ($beg % 2) == 0;
  croak "Internal error: segment end is even" if ($end % 2) == 0;
  croak "Internal error: segment end < beg" if $end < $beg;
  croak "Internal error: segment beg should be >= 3" if $beg < 3;
  my $range = int( ($end - $beg) / 2 ) + 1;

  # Prefill with 3 and 5 already marked, and offset to the segment start.
  my $whole = int( ($range+14) / 15);
  my $startp = ($beg % 30) >> 1;
  my $sieve = substr("011010010010110", $startp) . "011010010010110" x $whole;
  # Set 3 and 5 to prime if we're sieving them.
  substr($sieve,0,2) = "00" if $beg == 3;
  substr($sieve,0,1) = "0"  if $beg == 5;
  # Get rid of any extra we added.
  substr($sieve, $range) = '';

  # If the end value is below 7^2, then the pre-sieve is all we needed.
  return \$sieve if $end < 49;

  my $limit = int(sqrt($end)) + 1;
  # For large value of end, it's a huge win to just walk primes.
  my $primesieveref = _sieve_erat($limit);
  my $p = 7-2;
  foreach my $s (split("0", substr($$primesieveref, 3), -1)) {
    $p += 2 + 2 * length($s);
    my $p2 = $p*$p;
    last if $p2 > $end;
    if ($p2 < $beg) {
      $p2 = int($beg / $p) * $p;
      $p2 += $p if $p2 < $beg;
      $p2 += $p if ($p2 % 2) == 0;   # Make sure p2 is odd
    }
    # With large bases and small segments, it's common to find we don't hit
    # the segment at all.  Skip all the setup if we find this now.
    if ($p2 <= $end) {
      # Inner loop marking multiples of p
      # (everything is divided by 2 to keep inner loop simpler)
      my $filter_end = ($end - $beg) >> 1;
      my $filter_p2  = ($p2  - $beg) >> 1;
      while ($filter_p2 <= $filter_end) {
        substr($sieve, $filter_p2, 1) = "1";
        $filter_p2 += $p;
      }
    }
  }
  \$sieve;
}

sub trial_primes {
  my($low,$high) = @_;
  if (!defined $high) {
    $high = $low;
    $low = 2;
  }
  _validate_positive_integer($low);
  _validate_positive_integer($high);

  return if $low > $high;

  my @primes;
  $low-- if $low >= 2;
  my $curprime = next_prime($low);
  while ($curprime <= $high) {
    push @primes, $curprime;
    $curprime = next_prime($curprime);
  }
  return \@primes;
}

sub primes {
  my $optref = (ref $_[0] eq 'HASH')  ?  shift  :  {};
  croak "no parameters to primes" unless scalar @_ > 0;
  croak "too many parameters to primes" unless scalar @_ <= 2;
  my $low = (@_ == 2)  ?  shift  :  2;
  my $high = shift;
  my $sref = [];

  _validate_positive_integer($low);
  _validate_positive_integer($high);

  return $sref if ($low > $high) || ($high < 2);

  # Ignore method options in this code

  # At some point even the pretty-fast pure perl sieve is going to be a
  # dog, and we should move to trials.  This is typical with a small range
  # on a large base.  More thought on the switchover should be done.
  return trial_primes($low, $high) if ref($low)  eq 'Math::BigInt'
                                   || ref($high) eq 'Math::BigInt'
                                   || ($low > 1_000_000_000_000 && ($high-$low) < int($low/1_000_000));

  push @$sref, 2  if ($low <= 2) && ($high >= 2);
  push @$sref, 3  if ($low <= 3) && ($high >= 3);
  push @$sref, 5  if ($low <= 5) && ($high >= 5);
  $low = 7 if $low < 7;
  $low++ if ($low % 2) == 0;
  $high-- if ($high % 2) == 0;
  return $sref if $low > $high;

  if ($low == 7) {
    my $sieveref = _sieve_erat($high);
    my $n = $low - 2;
    foreach my $s (split("0", substr($$sieveref, 3), -1)) {
      $n += 2 + 2 * length($s);
      push @$sref, $n if $n <= $high;
    }
  } else {
    my $sieveref = _sieve_segment($low,$high);
    my $n = $low - 2;
    foreach my $s (split("0", $$sieveref, -1)) {
      $n += 2 + 2 * length($s);
      push @$sref, $n if $n <= $high;
    }
  }
  $sref;
}

sub next_prime {
  my($n) = @_;
  _validate_positive_integer($n);
  if ($n >= ((_PP_prime_maxbits == 32) ? 4294967291 : 18446744073709551557)) {
    return 0 if ref($_[0]) ne 'Math::BigInt';
    $n = $_[0];  # $n is a bigint now
  }
  return $_prime_next_small[$n] if $n <= $#_prime_next_small;

  # Be careful trying to do:
  #     my $d = int($n/30);
  #     my $m = $n - $d*30;
  # See:  int(9999999999999999403 / 30) => 333333333333333312  (off by 1)
  my $m = $n % 30;
  my $d = ($n - $m) / 30;
  if ($m == 29) { $d++;  $m = 1;} else { $m = $_nextwheel30[$m]; }
  while (!_is_prime7($d*30+$m)) {
    $m = $_nextwheel30[$m];
    $d++ if $m == 1;
  }
  $d*30 + $m;
}

sub prev_prime {
  my($n) = @_;
  _validate_positive_integer($n);
  if ($n <= 7) {
    return ($n <= 2) ? 0 : ($n <= 3) ? 2 : ($n <= 5) ? 3 : 5;
  }

  $n++ if ($n % 2) == 0;
  do {
    $n -= 2;
  } while ( (($n % 3) == 0) || (($n % 5) == 0) || (!_is_prime7($n)) );
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

  #my $m = $n % 30;
  #my $d = int( ($n - $m) / 30 );
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
  _validate_positive_integer($low);
  _validate_positive_integer($high);

  my $count = 0;

  $count++ if ($low <= 2) && ($high >= 2);   # Count 2
  $low = 3 if $low < 3;

  $low++ if ($low % 2) == 0;   # Make low go to odd number.
  $high-- if ($high % 2) == 0; # Make high go to odd number.
  return $count if $low > $high;

  if (   ref($low) eq 'Math::BigInt' || ref($high) eq 'Math::BigInt'
      || $high > 16_000_000_000
      || ($high-$low) < int($low/1_000_000) ) {
    # Too big to sieve.
    my $count = 0;
    my $curprime = next_prime($low-1);
    while ($curprime <= $high) {
      $count++;
      $curprime = next_prime($curprime);
    }
    return $count;
  }

  my $sieveref;
  if ($low == 3) {
    $sieveref = _sieve_erat($high);
  } else {
    $sieveref = _sieve_segment($low,$high);
  }

  $count += $$sieveref =~ tr/0//;

  $count;
}


sub nth_prime {
  my($n) = @_;
  _validate_positive_integer($n);

  return $_primes_small[$n] if $n <= $#_primes_small;

  if (!defined $bigint::VERSION) { # This isn't ideal.
    if (_PP_prime_maxbits == 32) {
      croak "nth_prime($n) overflow" if $n > 203280221;
    } else {
      croak "nth_prime($n) overflow" if $n > 425656284035217743;
    }
  }

  my $prime = 0;

  {
    my $count = 1;
    my $start = 3;
    # Make sure incr is an even number.
    my $incr = ($n < 1000) ? 1000 : ($n < 10000) ? 10000 : 100000;
    my $sieveref;
    while (1) {
      $sieveref = _sieve_segment($start, $start+$incr);
      my $segcount = $$sieveref =~ tr/0//;
      last if ($count + $segcount) >= $n;
      $count += $segcount;
      $start += $incr+2;
    }
    # Our count is somewhere in this segment.  Need to look for it.
    $prime = $start - 2;
    while ($count < $n) {
      $prime += 2;
      $count++ if !substr($$sieveref, ($prime-$start)>>1, 1);
    }
  }
  $prime;
}

sub _mulmod {
  my($a, $b, $m) = @_;
  return (($a * $b) % $m) if ($a|$b) < $_half_word;
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

sub _native_powmod {
  my($n, $power, $m) = @_;
  my $t = 1;
  $n = $n % $m;
  while ($power) {
    $t = ($t * $n) % $m if ($power & 1);
    $power >>= 1;
    $n = ($n * $n) % $m if $power;
  }
  $t;
}

sub _powmod {
  my($n, $power, $m) = @_;
  my $t = 1;

  if  ($m < $_half_word) {
    $n %= $m;
    while ($power) {
      $t = ($t * $n) % $m if ($power & 1);
      $power >>= 1;
      $n = ($n * $n) % $m if $power;
    }
  } else {
    while ($power) {
      $t = _mulmod($t, $n, $m) if ($power & 1);
      $power >>= 1;
      $n = _mulmod($n, $n, $m) if $power;
    }
  }
  $t;
}

sub _gcd_ui {
  my($x, $y) = @_;
  if ($y < $x) { ($x, $y) = ($y, $x); }
  while ($y > 0) {
    # y1 <- x0 % y0 ; x1 <- y0
    my $t = $y;
    $y = $x % $y;
    $x = $t;
  }
  $x;
}

sub _is_perfect_power {
  my $n = shift;
  my $log2n = _log2($n);
  $n = Math::BigInt->new("$n") unless ref($n) eq 'Math::BigInt';
  for my $e (@{primes($log2n)}) {
    return 1 if $n->copy()->broot($e)->bpow($e) == $n;
  }
  0;
}

sub _order {
  my($r, $n, $lim) = @_;
  $lim = $r unless defined $lim;

  return 1 if ($n % $r) == 1;
  for (my $j = 2; $j <= $lim; $j++) {
    return $j if _powmod($n, $j, $r) == 1;
  }
  return $lim+1;
}

# same result as:  int($n->blog(2)->floor->bstr)
sub _log2 {
  my $n = shift;
  my $log2n = 0;
  $log2n++ while ($n >>= 1);
  $log2n;
}



sub miller_rabin {
  my($n, @bases) = @_;
  _validate_positive_integer($n);
  croak "No bases given to miller_rabin" unless @bases;

  return 0 if ($n == 0) || ($n == 1);
  return 1 if ($n == 2) || ($n == 3);
  return 0 if !($n % 2);

  if ( ref($n) eq 'Math::BigInt' ) {

    my $s = 0;
    my $nminus1 = $n->copy->bsub(1);
    my $d = $nminus1->copy;
    while ($d->is_even) {
      $s++;
      $d->brsft(1);
    }

    foreach my $a (@bases) {
      croak "Base $a is invalid" if $a < 2;
      my $x = $n->copy->bzero->badd($a)->bmodpow($d,$n);
      next if ($x->is_one) || ($x->bcmp($nminus1) == 0);
      foreach my $r (1 .. $s) {
        $x->bmul($x); $x->bmod($n);
        return 0 if $x->is_one;
        if ($x->bcmp($nminus1) == 0) {
          $a = 0;
          last;
        }
      }
      return 0 if $a != 0;
    }

  } else {

   my $s = 0;
   my $d = $n - 1;
   while ( ($d & 1) == 0 ) {
     $s++;
     $d >>= 1;
   }

   if ($n < $_half_word) {
    foreach my $a (@bases) {
      croak "Base $a is invalid" if $a < 2;
      my $x = _native_powmod($a, $d, $n);
      next if ($x == 1) || ($x == ($n-1));
      foreach my $r (1 .. $s) {
        $x = ($x*$x) % $n;
        return 0 if $x == 1;
        if ($x == ($n-1)) {
          $a = 0;
          last;
        }
      }
      return 0 if $a != 0;
    }
   } else {
    foreach my $a (@bases) {
      croak "Base $a is invalid" if $a < 2;
      my $x = _powmod($a, $d, $n);
      next if ($x == 1) || ($x == ($n-1));

      foreach my $r (1 .. $s) {
        $x = ($x < $_half_word) ? ($x*$x) % $n : _mulmod($x, $x, $n);
        return 0 if $x == 1;
        if ($x == ($n-1)) {
          $a = 0;
          last;
        }
      }
      return 0 if $a != 0;
    }
   }

  }
  1;
}

# Calculate Jacobi symbol (M|N)
sub _jacobi {
  my($n, $m) = @_;
  return 0 if $m <= 0 || ($m % 2) == 0;
  my $j = 1;
  if ($n < 0) {
    $n = -$n;
    $j = -$j if ($m % 4) == 3;
  }
  # Split loop so we can reduce n/m to non-bigints after first iteration.
  if ($n != 0) {
    while (($n % 2) == 0) {
      $n >>= 1;
      $j = -$j if ($m % 8) == 3 || ($m % 8) == 5;
    }
    ($n, $m) = ($m, $n);
    $j = -$j if ($n % 4) == 3 && ($m % 4) == 3;
    $n = $n % $m;
    $n = int($n->bstr) if $n <= ~0 && ref($n) eq 'Math::BigInt';
    $m = int($m->bstr) if $m <= ~0 && ref($m) eq 'Math::BigInt';
  }
  while ($n != 0) {
    while (($n % 2) == 0) {
      $n >>= 1;
      $j = -$j if ($m % 8) == 3 || ($m % 8) == 5;
    }
    ($n, $m) = ($m, $n);
    $j = -$j if ($n % 4) == 3 && ($m % 4) == 3;
    $n = $n % $m;
  }
  return ($m == 1) ? $j : 0;
}

# Find first D in sequence (5,-7,9,-11,13,-15,...) where (D|N) == -1
sub _find_jacobi_d_sequence {
  my($n) = @_;

  # D is typically quite small: 67 max for N < 10^19.  However, it is
  # theoretically possible D could grow unreasonably.  I'm giving up at 4000M.
  my $d = 5;
  my $sign = 1;
  while (1) {
    my $gcd = (ref($n) eq 'Math::BigInt') ? Math::BigInt::bgcd($d, $n)
                                          : _gcd_ui($d, $n);
    return 0 if $gcd > 1 && $gcd != $n;  # Found divisor $d
    my $j = _jacobi($d * $sign, $n);
    last if $j == -1;
    $d += 2;
    croak "Could not find Jacobi sequence for $n" if $d > 4_000_000_000;
    $sign = -$sign;
  }
  return ($sign * $d);
}


sub is_strong_lucas_pseudoprime {
  my($n) = @_;
  _validate_positive_integer($n);

  # We're trying to limit the bignum calculations as much as possible.
  # It's also important to try to match whatever they passed in.  For instance
  # if they use a GMP or Pari object, we must do the same.  Hence instead of:
  #        my $U = Math::BigInt->bone;
  # we do
  #        my $U = $n->copy->bone;
  # so U is the same class as n.  If they passed in a string or a small value,
  # then we just make it up.

  return 1 if $n == 2;
  return 0 if $n < 2 || ($n % 2) == 0;

  # References:
  #     http://www.trnicely.net/misc/bpsw.html
  #     Math::Primality

  # Check for perfect square
  if (ref($n) eq 'Math::BigInt') {
    my $mcheck = int(($n & 127)->bstr);
    if (($mcheck*0x8bc40d7d) & ($mcheck*0xa1e2f5d1) & 0x14020a) {
      # ~82% of non-squares were rejected by the bloom filter
      my $sq = $n->copy->bsqrt->bfloor;
      $sq->bmul($sq);
      return 0 if $sq == $n;
    }
  } else {
    my $mcheck = $n & 127;
    if (($mcheck*0x8bc40d7d) & ($mcheck*0xa1e2f5d1) & 0x14020a) {
      my $sq = int(sqrt($n));
      return 0 if ($sq*$sq) == $n;
    }
  }

  # Determine Selfridge D, P, and Q parameters
  my $D = _find_jacobi_d_sequence($n);
  return 0 if $D == 0;  # We found a divisor in the sequence
  my $P = 1;
  my $Q = int( (1 - $D) / 4 );
  # Verify we've calculated this right
  die "Selfridge error: $D, $P, $Q\n" if ($D != $P*$P - 4*$Q);
  #warn "N: $n  D: $D  P: $P  Q: $Q\n";

  # It's now time to perform the Lucas pseudoprimality test using $D.

  if (ref($n) ne 'Math::BigInt') {
    require Math::BigInt;
    $n = Math::BigInt->new($n);
  }

  my $m = $n->copy->badd(1);
  # Traditional d,s:
  #   my $d=$m->copy; my $s=0; while ($d->is_even) { $s++; $d->brsft(1); }
  #   die "Invalid $m, $d, $s\n" unless $m == $d * 2**$s;
  my $dstr = substr($m->as_bin, 2);
  $dstr =~ s/(0*)$//;
  my $s = length($1);

  my $ZERO = $n->copy->bzero;
  my $U = $ZERO + 1;
  my $V = $ZERO + $P;
  my $U_2m = $U->copy;
  my $V_2m = $V->copy;
  my $Q_m = $ZERO + $Q;
  my $Q_m2 = $Q_m->copy->bmul(2);
  my $Qkd = $Q_m->copy;
  substr($dstr,-1) = '';   #$d->brsft(1);
  #my $i = 1;
  while ($dstr ne '') {    #while (!$d->is_zero) {
    #warn "U=$U  V=$V  Qm=$Q_m  Qm2=$Q_m2\n";
    $U_2m->bmul($V_2m);             $U_2m->bmod($n);
    $V_2m->bmuladd($V_2m, -$Q_m2);  $V_2m->bmod($n);
    #warn "  $i  U2m=$U_2m  V2m=$V_2m\n";  $i++;
    $Q_m->bmul($Q_m);               $Q_m->bmod($n);
    $Q_m2 = $Q_m->copy->bmul(2);    # no mod
    if (substr($dstr,-1)) {   #if ($d->is_odd) {
      my $T1 = $U_2m->copy->bmul($V);
      my $T2 = $U_2m->copy->bmul($U)->bmul($D);
      $U->bmuladd($V_2m, $T1);         # U = U*V_2m + V*U_2m
      $U->badd($n) if $U->is_odd;      # U += n if U & 1
      $U->brsft(1,2);                  # U = floor(U / 2)
      $U->bmod($n);                    # U = U % n

      $V->bmuladd($V_2m, $T2);
      $V->badd($n) if $V->is_odd;
      $V->brsft(1,2);
      $V->bmod($n);

      $Qkd->bmul($Q_m);
      $Qkd->bmod($n);
    }
    substr($dstr,-1) = '';   #$d->brsft(1);
  }
  #warn "l0 U=$U  V=$V\n";
  return 1 if $U->is_zero || $V->is_zero;

  # Compute powers of V
  my $Qkd2 = $Qkd->copy->bmul(2);
  foreach my $r (1 .. $s-1) {
    #warn "l$r U=$U  V=$V  Qkd2=$Qkd2\n";
    $V->bmuladd($V, -$Qkd2);  $V->bmod($n);
    return 1 if $V->is_zero;
    if ($r < ($s-1)) {
      $Qkd->bmul($Qkd);  $Qkd->bmod($n);
      $Qkd2 = $Qkd->copy->bmul(2);
    }
  }
  return 0;
}


my $_poly_bignum;
sub _poly_new {
  my @poly;
  if ($_poly_bignum) {
    foreach my $c (@_) {
      push @poly, (ref $c eq 'Math::BigInt') ? $c->copy : Math::BigInt->new("$c");
    }
  } else {
    push @poly, $_ for (@_);
    push @poly, 0 unless scalar @poly;
  }
  return \@poly;
}

sub _poly_print {
  my($poly) = @_;
  warn "poly has null top degree" if $#$poly > 0 && !$poly->[-1];
  foreach my $d (reverse 1 .. $#$poly) {
    my $coef = $poly->[$d];
    print "", ($coef != 1) ? $coef : "", ($d > 1) ? "x^$d" : "x", " + "
      if $coef;
  }
  my $p0 = $poly->[0] || 0;
  print "$p0\n";
}

sub _poly_mod_mul {
  my($px, $py, $r, $n) = @_;

  my $px_degree = $#$px;
  my $py_degree = $#$py;
  my @res;

  # convolve(px, py) mod (X^r-1,n)
  my @indices_y = grep { $py->[$_] } (0 .. $py_degree);
  for (my $ix = 0; $ix <= $px_degree; $ix++) {
    my $px_at_ix = $px->[$ix];
    next unless $px_at_ix;
    foreach my $iy (@indices_y) {
      my $py_at_iy = $py->[$iy];
      my $rindex = ($ix + $iy) % $r;  # reduce mod X^r-1
      if (!defined $res[$rindex]) {
        $res[$rindex] = $_poly_bignum ? Math::BigInt->bzero : 0
      }
      $res[$rindex] = ($res[$rindex] + ($py_at_iy * $px_at_ix)) % $n;
    }
  }
  # In case we had upper terms go to zero after modulo, reduce the degree.
  pop @res while !$res[-1];
  return \@res;
}

sub _poly_mod_pow {
  my($pn, $power, $r, $mod) = @_;
  my $res = _poly_new(1);
  my $p = $power;

  while ($p) {
    $res = _poly_mod_mul($res, $pn, $r, $mod) if ($p & 1);
    $p >>= 1;
    $pn  = _poly_mod_mul($pn,  $pn, $r, $mod) if $p;
  }
  return $res;
}

sub _test_anr {
  my($a, $n, $r) = @_;
  my $pp = _poly_mod_pow(_poly_new($a, 1), $n, $r, $n);
  $pp->[$n % $r] = (($pp->[$n % $r] || 0) -  1) % $n;  # subtract X^(n%r)
  $pp->[      0] = (($pp->[      0] || 0) - $a) % $n;  # subtract a
  return 0 if scalar grep { $_ } @$pp;
  1;
}

sub is_aks_prime {
  my $n = shift;
  $n = Math::BigInt->new("$n") unless ref($n) eq 'Math::BigInt';

  return 0 if $n < 2;
  return 0 if _is_perfect_power($n);

  # limit = floor( log2(n) * log2(n) ).  o_r(n) must be larger than this
  my $sqrtn = int(Math::BigFloat->new($n)->bsqrt->bfloor->bstr);
  my $log2n = Math::BigFloat->new($n)->blog(2);
  my $log2_squared_n = $log2n * $log2n;
  my $limit = int($log2_squared_n->bfloor->bstr);

  my $r = next_prime($limit);
  foreach my $f (@{primes(0,$r-1)}) {
    return 1 if $f == $n;
    return 0 if !($n % $f);
  }

  while ($r < $n) {
    return 0 if !($n % $r);
    #return 1 if $r >= $sqrtn;
    last if _order($r, $n, $limit) > $limit;
    $r = next_prime($r);
  }

  return 1 if $r >= $n;

  # Since r is a prime, phi(r) = r-1
  my $rlimit = int( Math::BigFloat->new($r)->bsub(1)
                    ->bsqrt->bmul($log2n)->bfloor->bstr);

  $_poly_bignum = 1;
  if ( $n < ( (~0 == 4294967295) ? 65535 : 4294967295 ) ) {
    $_poly_bignum = 0;
    $n = int($n->bstr) if ref($n) eq 'Math::BigInt';
  }

  for (my $a = 1; $a <= $rlimit; $a++) {
    return 0 unless _test_anr($a, $n, $r);
  }

  return 1;
}


sub _basic_factor {
  # MODIFIES INPUT SCALAR
  return ($_[0]) if $_[0] < 4;

  my @factors;
  while ( !($_[0] % 2) ) { push @factors, 2;  $_[0] = int($_[0] / 2); }
  while ( !($_[0] % 3) ) { push @factors, 3;  $_[0] = int($_[0] / 3); }
  while ( !($_[0] % 5) ) { push @factors, 5;  $_[0] = int($_[0] / 5); }

  # Stop using bignum if we can
  $_[0] = int($_[0]->bstr) if ref($_[0]) eq 'Math::BigInt' && $_[0] <= ~0;

  if ( ($_[0] > 1) && _is_prime7($_[0]) ) {
    push @factors, $_[0];
    $_[0] = 1;
  }
  @factors;
}

sub trial_factor {
  my($n) = @_;
  _validate_positive_integer($n);

  my @factors = _basic_factor($n);
  return @factors if $n < 4;

  my $limit = int( sqrt($n) + 0.001);
  my $f = 3;
  while ($f <= $limit) {
    if ( ($n % $f) == 0) {
      while ( ($n % $f) == 0) {
        push @factors, $f;
        $n = int($n/$f);
      }
      $limit = int( sqrt($n) + 0.001);
    }
    $f += 2;
  }
  push @factors, $n  if $n > 1;
  @factors;
}

sub factor {
  my($n) = @_;
  _validate_positive_integer($n);

  return trial_factor($n) if $n < 100000;

  my @factors = _basic_factor($n);
  return @factors if $n < 4;

  # Use 'n = int($n/7)' instead of 'n/=7' to not "upgrade" n to a Math::BigFloat.
  while (($n %  7) == 0) { push @factors,  7;  $n = int($n /  7); }
  while (($n % 11) == 0) { push @factors, 11;  $n = int($n / 11); }
  while (($n % 13) == 0) { push @factors, 13;  $n = int($n / 13); }
  while (($n % 17) == 0) { push @factors, 17;  $n = int($n / 17); }
  while (($n % 19) == 0) { push @factors, 19;  $n = int($n / 19); }
  while (($n % 23) == 0) { push @factors, 23;  $n = int($n / 23); }
  while (($n % 29) == 0) { push @factors, 29;  $n = int($n / 29); }
  if ($n < (31*31)) {
    push @factors, $n  if $n != 1;
    return @factors;
  }

  my @nstack = ($n);
  while (@nstack) {
    $n = pop @nstack;
    # Don't use bignum on $n if it has gotten small enough.
    $n = int($n->bstr) if ref($n) eq 'Math::BigInt' && $n <= ~0;
    #print "Looking at $n with stack ", join(",",@nstack), "\n";
    while ( ($n >= (31*31)) && !_is_prime7($n) ) {
      my @ftry;
      my $holf_rounds = 0;
      if ($n < $_half_word) {
        $holf_rounds = 64*1024;
        #warn "trying holf 64k on $n\n";
        @ftry = holf_factor($n, $holf_rounds);
      }
      if (scalar @ftry < 2) {
        foreach my $add (3, 5, 7, 11, 13) {
          #warn "trying prho 64k {$add} on $n\n" if scalar @ftry < 2;
          @ftry = prho_factor($n, 64*1024, $add) if scalar @ftry < 2;
        }
      }
      if (scalar @ftry < 2) {
        #warn "trying holf 128k on $n\n";
        @ftry = holf_factor($n, 128*1024, $holf_rounds);
        $holf_rounds += 128*1024;
      }
      if (scalar @ftry < 2) {
        #warn "trying prho 128k {17} on $n\n";
        @ftry = prho_factor($n, 128*1024, 17);
      }
      if (scalar @ftry > 1) {
        #print "  split into ", join(",",@ftry), "\n";
        $n = shift @ftry;
        push @nstack, @ftry;
      } else {
        #warn "trial factor $n\n";
        push @factors, trial_factor($n);
        #print "  trial into ", join(",",@factors), "\n";
        $n = 1;
        last;
      }
    }
    push @factors, $n  if $n != 1;
  }
  sort {$a<=>$b} @factors;
}

# TODO:
sub fermat_factor { trial_factor(@_) }
sub squfof_factor { trial_factor(@_) }

sub prho_factor {
  my($n, $rounds, $a) = @_;
  _validate_positive_integer($n);
  $rounds = 4*1024*1024 unless defined $rounds;
  $a = 3 unless defined $a;

  my @factors = _basic_factor($n);
  return @factors if $n < 4;

  my $inloop = 0;
  my $U = 7;
  my $V = 7;

  if ( ref($n) eq 'Math::BigInt' ) {

    $U = $n->copy->bzero->badd($U);
    $V = $n->copy->bzero->badd($V);
    for my $i (1 .. $rounds) {
      # Would use bmuladd here, but old Math::BigInt's barf with scalar $a.
      #$U->bmuladd($U, $a);  $U->bmod($n);
      #$V->bmuladd($V, $a);  $V->bmod($n);
      #$V->bmuladd($V, $a);  $V->bmod($n);
      $U->bmul($U); $U->badd($a); $U->bmod($n);
      $V->bmul($V); $V->badd($a); $V->bmod($n);
      $V->bmul($V); $V->badd($a); $V->bmod($n);
      my $f = Math::BigInt::bgcd( ($U > $V) ? $U-$V : $V-$U,  $n);
      if ($f == $n) {
        last if $inloop++;  # We've been here before
      } elsif ($f != 1) {
        my $f2 = $n->copy->bdiv($f);
        push @factors, $f;
        push @factors, $f2;
        croak "internal error in prho" unless ($f * $f2) == $n;
        return @factors;
      }
    }

  } elsif ($n < $_half_word) {

    for my $i (1 .. $rounds) {
      $U = ($U * $U + $a) % $n;
      $V = ($V * $V + $a) % $n;
      $V = ($V * $V + $a) % $n;
      my $f = _gcd_ui( ($U > $V) ? $U-$V : $V-$U,  $n );
      if ($f == $n) {
        last if $inloop++;  # We've been here before
      } elsif ($f != 1) {
        push @factors, $f;
        push @factors, int($n/$f);
        croak "internal error in prho" unless ($f * int($n/$f)) == $n;
        return @factors;
      }
    }

  } else {

    for my $i (1 .. $rounds) {
      # U^2+a % n
      $U = _mulmod($U, $U, $n);
      $U = (($n-$U) > $a)  ?  $U+$a  :  $U+$a-$n;
      # V^2+a % n
      $V = _mulmod($V, $V, $n);
      $V = (($n-$V) > $a)  ?  $V+$a  :  $V+$a-$n;
      # V^2+a % n
      $V = _mulmod($V, $V, $n);
      $V = (($n-$V) > $a)  ?  $V+$a  :  $V+$a-$n;
      my $f = _gcd_ui( ($U > $V) ? $U-$V : $V-$U,  $n );
      if ($f == $n) {
        last if $inloop++;  # We've been here before
      } elsif ($f != 1) {
        push @factors, $f;
        push @factors, int($n/$f);
        croak "internal error in prho" unless ($f * int($n/$f)) == $n;
        return @factors;
      }
    }

  }
  push @factors, $n;
  @factors;
}

sub pbrent_factor {
  my($n, $rounds) = @_;
  _validate_positive_integer($n);
  $rounds = 4*1024*1024 unless defined $rounds;

  my @factors = _basic_factor($n);
  return @factors if $n < 4;

  my $a = 1;
  my $Xi = 2;
  my $Xm = 2;

  if ( ref($n) eq 'Math::BigInt' ) {

    $Xi = $n->copy->bzero->badd($Xi);
    $Xm = $n->copy->bzero->badd($Xm);
    for my $i (1 .. $rounds) {
      $Xi->bmul($Xi);  $Xi->badd($a);  $Xi->bmod($n);
      my $f = Math::BigInt::bgcd( ($Xi > $Xm) ? $Xi-$Xm : $Xm-$Xi,  $n);
      if ( ($f != 1) && ($f != $n) ) {
        my $f2 = $n->copy->bdiv($f);
        push @factors, $f;
        push @factors, $f2;
        croak "internal error in pbrent" unless ($f * $f2) == $n;
        return @factors;
      }
      $Xm = $Xi->copy if ($i & ($i-1)) == 0;  # i is a power of 2
    }

  } elsif ($n < $_half_word) {

    for my $i (1 .. $rounds) {
      $Xi = ($Xi * $Xi + $a) % $n;
      my $f = _gcd_ui( ($Xi > $Xm) ? $Xi-$Xm : $Xm-$Xi,  $n );
      if ( ($f != 1) && ($f != $n) ) {
        push @factors, $f;
        push @factors, int($n/$f);
        croak "internal error in pbrent" unless ($f * int($n/$f)) == $n;
        return @factors;
      }
      $Xm = $Xi if ($i & ($i-1)) == 0;  # i is a power of 2
    }

  } else {

    for my $i (1 .. $rounds) {
      # Xi^2+a % n
      $Xi = _mulmod($Xi, $Xi, $n);
      $Xi = (($n-$Xi) > $a)  ?  $Xi+$a  :  $Xi+$a-$n;
      my $f = _gcd_ui( ($Xi > $Xm) ? $Xi-$Xm : $Xm-$Xi,  $n );
      if ( ($f != 1) && ($f != $n) ) {
        push @factors, $f;
        push @factors, int($n/$f);
        croak "internal error in pbrent" unless ($f * int($n/$f)) == $n;
        return @factors;
      }
      $Xm = $Xi if ($i & ($i-1)) == 0;  # i is a power of 2
    }

  }
  push @factors, $n;
  @factors;
}

# This code is bollocks.  See a proper implementation in factor.c
sub pminus1_factor {
  my($n, $rounds) = @_;
  _validate_positive_integer($n);
  $rounds = 4*1024*1024 unless defined $rounds;

  my @factors = _basic_factor($n);
  return @factors if $n < 4;

  if ( ref($n) eq 'Math::BigInt' ) {
    my $kf = $n->copy->bzero->badd(13);
    for my $i (1 .. $rounds) {
      $kf->bmodpow($i,$n);
      $kf = $n if $kf == 0;
      my $f = Math::BigInt::bgcd( $kf-1, $n );
      if ( ($f != 1) && ($f != $n) ) {
        my $f2 = $n->copy->bdiv($f);
        push @factors, $f;
        push @factors, $f2;
        croak "internal error in pminus1" unless ($f * $f2) == $n;
        return @factors;
      }
    }
  } else {
    my $kf = 13;
    for my $i (1 .. $rounds) {
      $kf = _powmod($kf, $i, $n);
      $kf = $n if $kf == 0;
      my $f = _gcd_ui( $kf-1, $n );
      if ( ($f != 1) && ($f != $n) ) {
        push @factors, $f;
        push @factors, int($n/$f);
        croak "internal error in pminus1" unless ($f * int($n/$f)) == $n;
        return @factors;
      }
    }
  }
  push @factors, $n;
  @factors;
}

sub holf_factor {
  my($n, $rounds, $startrounds) = @_;
  _validate_positive_integer($n);
  $rounds = 64*1024*1024 unless defined $rounds;
  $startrounds = 1 unless defined $startrounds;

  my @factors = _basic_factor($n);
  return @factors if $n < 4;

  if ( ref($n) eq 'Math::BigInt' ) {
    for my $i ($startrounds .. $rounds) {
      my $ni = $n->copy->bmul($i);
      my $s = $ni->copy->bsqrt->bfloor;
      $s->binc if ($s * $s) != $ni;
      my $m = $s->copy->bmul($s)->bmod($n);
      # Check for perfect square
      my $mcheck = int(($m & 127)->bstr);
      next if (($mcheck*0x8bc40d7d) & ($mcheck*0xa1e2f5d1) & 0x14020a);
      # ... 82% of non-squares were rejected by the bloom filter
      my $f = $m->copy->bsqrt->bfloor;
      next unless ($f*$f) == $m;
      $f = Math::BigInt::bgcd( ($s > $f) ? $s-$f : $f-$s,  $n);
      last if $f == 1 || $f == $n;   # Should never happen
      my $f2 = $n->copy->bdiv($f);
      push @factors, $f;
      push @factors, $f2;
      croak "internal error in HOLF" unless ($f * $f2) == $n;
      # print "HOLF found factors in $i rounds\n";
      return @factors;
    }
  } else {
    for my $i ($startrounds .. $rounds) {
      my $s = int(sqrt($n * $i));
      $s++ if ($s * $s) != ($n * $i);
      my $m = ($s < $_half_word) ? ($s*$s) % $n : _mulmod($s, $s, $n);
      # Check for perfect square
      my $mcheck = $m & 127;
      next if (($mcheck*0x8bc40d7d) & ($mcheck*0xa1e2f5d1) & 0x14020a);
      # ... 82% of non-squares were rejected by the bloom filter
      my $f = int(sqrt($m));
      next unless $f*$f == $m;
      $f = _gcd_ui( ($s > $f)  ?  $s - $f  :  $f - $s,  $n);
      last if $f == 1 || $f == $n;   # Should never happen
      push @factors, $f;
      push @factors, int($n/$f);
      croak "internal error in HOLF" unless ($f * int($n/$f)) == $n;
      # print "HOLF found factors in $i rounds\n";
      return @factors;
    }
  }
  push @factors, $n;
  @factors;
}




my $_const_euler = 0.57721566490153286060651209008240243104215933593992;
my $_const_li2 = 1.045163780117492784844588889194613136522615578151;

sub ExponentialIntegral {
  my($x, $tol) = @_;
  $tol = 1e-16 unless defined $tol;
  my $sum = 0.0;
  my($y, $t);
  my $c = 0.0;

  croak "Invalid input to ExponentialIntegral:  x must be != 0" if $x == 0;

  $x = new Math::BigFloat "$x"  if defined $bignum::VERSION && ref($x) ne 'Math::BigFloat';

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
# So many terms and digits are used so we can quickly do bignum R.
my @_Riemann_Zeta_Table = (
  '0.64493406684822643647241516664602518921894990',   # zeta(2) - 1
  '0.20205690315959428539973816151144999076498629',
  '0.082323233711138191516003696541167902774750952',
  '0.036927755143369926331365486457034168057080920',
  '0.017343061984449139714517929790920527901817490',
  '0.0083492773819228268397975498497967595998635606',
  '0.0040773561979443393786852385086524652589607906',
  '0.0020083928260822144178527692324120604856058514',
  '0.00099457512781808533714595890031901700601953156',
  '0.00049418860411946455870228252646993646860643576',
  '0.00024608655330804829863799804773967096041608846',
  '0.00012271334757848914675183652635739571427510590',
  '0.000061248135058704829258545105135333747481696169',
  '0.000030588236307020493551728510645062587627948707',
  '0.000015282259408651871732571487636722023237388990',
  '0.0000076371976378997622736002935630292130882490903',
  '0.0000038172932649998398564616446219397304546972190',
  '0.0000019082127165539389256569577951013532585711448',
  '0.00000095396203387279611315203868344934594379418741',
  '0.00000047693298678780646311671960437304596644669478',
  '0.00000023845050272773299000364818675299493504182178',
  '0.00000011921992596531107306778871888232638725499778',
  '0.000000059608189051259479612440207935801227503918837',
  '0.000000029803503514652280186063705069366011844730920',
  '0.000000014901554828365041234658506630698628864788168',
  '0.0000000074507117898354294919810041706041194547190319',
  '0.0000000037253340247884570548192040184024232328930593',
  '0.0000000018626597235130490064039099454169480616653305',
  '0.00000000093132743241966818287176473502121981356795514',
  '0.00000000046566290650337840729892332512200710626918534',
  '0.00000000023283118336765054920014559759404950248298228',
  '0.00000000011641550172700519775929738354563095165224717',
  '0.000000000058207720879027008892436859891063054173122605',
  '0.000000000029103850444970996869294252278840464106981987',
  '0.000000000014551921891041984235929632245318420983808894',
  '0.0000000000072759598350574810145208690123380592648509256',
  '0.0000000000036379795473786511902372363558732735126460284',
  '0.0000000000018189896503070659475848321007300850305893096',
  '0.00000000000090949478402638892825331183869490875386000099',
  '0.00000000000045474737830421540267991120294885703390452991',
  '0.00000000000022737368458246525152268215779786912138298220',
  '0.00000000000011368684076802278493491048380259064374359028',
  '0.000000000000056843419876275856092771829675240685530571589',
  '0.000000000000028421709768893018554550737049426620743688265',
  '0.000000000000014210854828031606769834307141739537678698606',
  '0.0000000000000071054273952108527128773544799568000227420436',
  '0.0000000000000035527136913371136732984695340593429921456555',
  '0.0000000000000017763568435791203274733490144002795701555086',
  '0.00000000000000088817842109308159030960913863913863256088715',
  '0.00000000000000044408921031438133641977709402681213364596031',
  '0.00000000000000022204460507980419839993200942046539642366543',
  '0.00000000000000011102230251410661337205445699213827024832229',
  '0.000000000000000055511151248454812437237365905094302816723551',
  '0.000000000000000027755575621361241725816324538540697689848904',
  '0.000000000000000013877787809725232762839094906500221907718625',
  '0.0000000000000000069388939045441536974460853262498092748358742',
  '0.0000000000000000034694469521659226247442714961093346219504706',
  '0.0000000000000000017347234760475765720489729699375959074780545',
  '0.00000000000000000086736173801199337283420550673429514879071415',
  '0.00000000000000000043368086900206504874970235659062413612547801',
  '0.00000000000000000021684043449972197850139101683209845761574010',
  '0.00000000000000000010842021724942414063012711165461382589364744',
  '0.000000000000000000054210108624566454109187004043886337150634224',
  '0.000000000000000000027105054312234688319546213119497764318887282',
  '0.000000000000000000013552527156101164581485233996826928328981877',
  '0.0000000000000000000067762635780451890979952987415566862059812586',
  '0.0000000000000000000033881317890207968180857031004508368340311585',
  '0.0000000000000000000016940658945097991654064927471248619403036418',
  '0.00000000000000000000084703294725469983482469926091821675222838642',
  '0.00000000000000000000042351647362728333478622704833579344088109717',
  '0.00000000000000000000021175823681361947318442094398180025869417612',
  '0.00000000000000000000010587911840680233852265001539238398470699902',
  '0.000000000000000000000052939559203398703238139123029185055866375629',
  '0.000000000000000000000026469779601698529611341166842038715592556134',
  '0.000000000000000000000013234889800848990803094510250944989684323826',
  '0.0000000000000000000000066174449004244040673552453323082200147137975',
  '0.0000000000000000000000033087224502121715889469563843144048092764894',
  '0.0000000000000000000000016543612251060756462299236771810488297723589',
  '0.00000000000000000000000082718061255303444036711056167440724040096811',
  '0.00000000000000000000000041359030627651609260093824555081412852575873',
  '0.00000000000000000000000020679515313825767043959679193468950443365312',
  '0.00000000000000000000000010339757656912870993284095591745860911079606',
  '0.000000000000000000000000051698788284564313204101332166355512893608164',
  '0.000000000000000000000000025849394142282142681277617708450222269121159',
  '0.000000000000000000000000012924697071141066700381126118331865309299779',
  '0.0000000000000000000000000064623485355705318034380021611221670660356864',
  '0.0000000000000000000000000032311742677852653861348141180266574173608296',
  '0.0000000000000000000000000016155871338926325212060114057052272720509148',
  '0.00000000000000000000000000080779356694631620331587381863408997398684847',
  '0.00000000000000000000000000040389678347315808256222628129858130379479700',
  '0.00000000000000000000000000020194839173657903491587626465673047518903728',
  '0.00000000000000000000000000010097419586828951533619250700091044144538432',
  '0.000000000000000000000000000050487097934144756960847711725486604360898735',
  '0.000000000000000000000000000025243548967072378244674341937966175648398693',
  '0.000000000000000000000000000012621774483536189043753999660777148710632765',
  '0.0000000000000000000000000000063108872417680944956826093943332037500694712',
  '0.0000000000000000000000000000031554436208840472391098412184847972814371270',
  '0.0000000000000000000000000000015777218104420236166444327830159601782237092',
);

# Select n = 55, good for 46ish digits of accuracy.
my $_Borwein_n = 55;
my @_Borwein_dk = (
  '1',
  '6051',
  '6104451',
  '2462539971',
  '531648934851',
  '71301509476803',
  '6504925195108803',
  '429144511928164803',
  '21392068013887742403',
  '832780518854440804803',
  '25977281563850106233283',
  '662753606729324750201283',
  '14062742362385399866745283',
  '251634235316509414702211523',
  '3841603462178827861104812483',
  '50535961819850087101900022211',
  '577730330374203014014104003011',
  '5782012706584553297863989289411',
  '50984922488525881477588707205571',
  '398333597655022403279683908035011',
  '2770992240330783259897072664469955',
  '17238422988353715312442126057365955',
  '96274027751337344115352100618133955',
  '484350301573059857715727453968687555',
  '2201794236784087151947175826243477955',
  '9068765987529892610841571032285864387',
  '33926582279822401059328069515697217987',
  '115535262182820447663793177744255246787',
  '358877507711760077538925500462137369027',
  '1018683886695854101193095537014797787587',
  '2646951832121008166346437186541363159491',
  '6306464665572570713623910486640730071491',
  '13799752848354341643763498672558481367491',
  '27780237373991939435100856211039992177091',
  '51543378762608611361377523633779417047491',
  '88324588911945720951614452340280439890371',
  '140129110249040241501243929391690331218371',
  '206452706984942815385219764876242498642371',
  '283527707823296964404071683165658912154051',
  '364683602811933600833512164561308162744771',
  '441935796522635816776473230396154031661507',
  '508231717051242054487234759342047053767107',
  '559351463001010719709990637083458540691907',
  '594624787018881191308291683229515933311427',
  '616297424973434835299724300924272199623107',
  '628083443816135918099559567176252011864515',
  '633714604276098212796088600263676671320515',
  '636056734158553360761837806887547188568515',
  '636894970116484676875895417679248215794115',
  '637149280289288581322870186196318041432515',
  '637213397278310656625865036925470191411651',
  '637226467136294189739463288384528579584451',
  '637228536449134002301138291602841035366851',
  '637228775173095037281299181461988671775171',
  '637228793021615488494769154535569803469251',
  '637228793670652595811622608101881844621763',
);
# "An Efficient Algorithm for the Riemann Zeta Function", Borwein, 1991.
# About 1.3n terms are needed for n digits of accuracy.
sub _Recompute_Dk {
  my $nterms = shift;
  $_Borwein_n = $nterms;
  @_Borwein_dk = ();
  foreach my $k (0 .. $nterms) {
    my $dsum = Math::BigFloat->bzero;
    $dsum->accuracy(2*$_Borwein_n);
    my $n = Math::BigInt->new($nterms-1)->bfac;
    my $d = Math::BigInt->new($nterms)->bfac;
    foreach my $i (0 .. $k) {
      my $term = Math::BigFloat->bone;
      $term->accuracy(2*$_Borwein_n);
      $term->bmul($n)->bdiv($d);
      $dsum += $term;
      $n->bmul($nterms+$i)->bmul(4);
      $d->bdiv($nterms-$i)->bmul(2*$i+1)->bmul(2*$i+2);
    }
    my $dk = ($nterms * $dsum + 1e-20)->as_int;
    $_Borwein_dk[$k] = $dk;
    #print  "  '$dk',\n";
  }
}

sub RiemannZeta {
  my($x, $tol) = @_;

  if (!defined $bignum::VERSION && ref($x) !~ /^Math::Big/) {
    return 0.0 + $_Riemann_Zeta_Table[int($x)-2]
      if $x == int($x) && defined $_Riemann_Zeta_Table[int($x)-2];
    $tol = 1e-16 unless defined $tol;
    my($y, $t);
    my $sum = 0.0;
    my $c = 0.0;

    for my $k (2 .. 1000000) {
      my $term = (2*$k+1) ** -$x;
      $y = $term-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
      last if $term < abs($tol*$sum);
    }
    my $term = 3 ** -$x;
    $y = $term-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
    $t = 1.0 / (1.0 - (2 ** -$x));
    $sum *= $t;
    $sum += ($t - 1.0);
    return $sum;
  }

  return Math::BigFloat->new($_Riemann_Zeta_Table[int($x)-2])
      if $x == int($x) && defined $_Riemann_Zeta_Table[int($x)-2];

  $tol = 1e-40 unless defined $tol;

  # Trying to work around Math::BigFloat bugs RT 43692 and RT 77105 which make
  # a right mess of things.  Watch this:
  #   my $n = Math::BigFloat->new(11); $n->accuracy(64); say $n**1.1;  # 13.98
  #   my $n = Math::BigFloat->new(11); $n->accuracy(67); say $n**1.1;  # 29.98
  # We can fix some issues with large exponents (e.g. 6^-40.5) by turning it
  # into (6^-(40.5/4))^4  (assuming the base is positive).  Without that hack,
  # none of this would work at all.

  $x = Math::BigFloat->new($x);
  my $superx = 1;
  my $subx = Math::BigFloat->new($x);
  while ($subx > 8) {
    $superx *= 2;
    $subx /= 2;
  }

  # Go with the basic formula for large x, as it best works around the mess,
  # though is unfortunately much slower.
  if ($x > 30) {
    my $sum = 0.0;
    for my $k (4 .. 1000) {
      my $term = ( $k ** -$subx )  ** $superx;
      $sum += $term;
      last if $term < ($sum*$tol);
    }
    for my $k (3, 2) {
      my $term = ( $k ** -$subx )  ** $superx;
      $sum += $term;
    }
    return $sum;
  }
  #if ($x > 25) {
  #  my $sum = 0.0;
  #  my $divisor = 1.0 - ((2 ** -$subx) ** $superx);
  #  for my $k (2 .. 1000) {
  #    my $term = ( (2*$k+1) ** -$subx )  ** $superx;
  #    $sum += $term;
  #    last if $term < ($tol*$divisor);
  #  }
  #  $sum += (3 ** -$subx) ** $superx;
  #  my $t = 1.0 / $divisor;
  #  $sum *= $t;
  #  $sum += ($t - 1.0);
  #  return $sum;
  #}

  # If we wanted to change the Borwein series being used:
  # _Recompute_Dk(55);

  if (ref $_Borwein_dk[0] ne 'Math::BigInt') {
    @_Borwein_dk = map { Math::BigInt->new("$_") } @_Borwein_dk;
  }

  my $n = $_Borwein_n;
  my $intermediate_accuracy = undef;
  my $one = Math::BigFloat->bone;
  $one->accuracy($intermediate_accuracy) if defined $intermediate_accuracy;

  my $d1 = Math::BigFloat->new(2);
  $d1->accuracy($intermediate_accuracy) if defined $intermediate_accuracy;
  # with bignum on, $d1->bpow($one-$x) doesn't change d1 !
  $d1 = $d1 ** ($one - $x);
  my $divisor = $one->copy->bsub($d1)->bmul(-$_Borwein_dk[$n]);

  $tol = $divisor->copy->bmul($tol)->babs();

  my $sum = Math::BigFloat->bzero;
  $sum->accuracy($intermediate_accuracy) if defined $intermediate_accuracy;
  foreach my $k (1 .. $n-1) {
    my $term = Math::BigFloat->new( $_Borwein_dk[$k] - $_Borwein_dk[$n] );
    $term *= -1 if $k % 2;
    $term->accuracy($intermediate_accuracy) if defined $intermediate_accuracy;
    my $den = Math::BigFloat->new($k+1);
    $den->accuracy($intermediate_accuracy) if defined $intermediate_accuracy;
    $den = ($den ** $subx) ** $superx;
    $term /= $den;
    $sum += $term;
    last if $term->copy->babs() < $tol;
  }
  $sum += Math::BigFloat->new( $one - $_Borwein_dk[$n] ); # term k=0
  $sum->bdiv( $divisor );
  $sum->bsub(1);
  return $sum;
}

# Riemann R function
sub RiemannR {
  my($x, $tol) = @_;
  my($y, $t);

  croak "Invalid input to ReimannR:  x must be > 0" if $x <= 0;

  if (!defined $bignum::VERSION && ref($x) !~ /^Math::Big/) {
    $tol = 1e-16 unless defined $tol;
    my $sum = 0.0;
    my $c = 0.0;

    $y = 1.0-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
    my $flogx = log($x);
    my $part_term = 1.0;
    for my $k (1 .. 10000) {
      # Small k from table, larger k from function
      my $zeta = ($k <= $#_Riemann_Zeta_Table)
                 ? $_Riemann_Zeta_Table[$k+1-2]
                 : RiemannZeta($k+1);
      $part_term *= $flogx / $k;
      my $term = $part_term / ($k + $k * $zeta);
      $y = $term-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
      last if abs($term/$sum) < $tol;
    }
    return $sum;
  }

  $x = new Math::BigFloat "$x"  if ref($x) ne 'Math::BigFloat';

  $tol = 1e-35 unless defined $tol;
  my $sum = Math::BigFloat->bzero;
  my $c = Math::BigFloat->bzero;

  $y = 1.0-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
  my $flogx = log($x);
  my $part_term = Math::BigFloat->bone;

  for my $k (1 .. 10000) {
    # Small k from table, larger k from function
    my $zeta = ($k <= $#_Riemann_Zeta_Table)
               ? Math::BigFloat->new($_Riemann_Zeta_Table[$k+1-2])
               : RiemannZeta($k+1);
    $part_term *= $flogx / $k;
    my $term = $part_term / ($k + $k * $zeta);
    $y = $term-$c; $t = $sum+$y; $c = ($t-$sum)-$y; $sum = $t;
    last if abs($term/$sum) < $tol;
  }
  return $sum;
}

1;

__END__


# ABSTRACT: Pure Perl version of Math::Prime::Util

=pod

=encoding utf8


=head1 NAME

Math::Prime::Util::PP - Pure Perl version of Math::Prime::Util


=head1 VERSION

Version 0.13


=head1 SYNOPSIS

The functionality is basically identical to L<Math::Prime::Util>, as this
module is just the Pure Perl implementation.

  # Normally you would just import the functions you are using.
  # Nothing is exported by default.
  use Math::Prime::Util ':all';


  # Get a big array reference of many primes
  my $aref = primes( 100_000_000 );

  # All the primes between 5k and 10k inclusive
  my $aref = primes( 5_000, 10_000 );

  # If you want them in an array instead
  my @primes = @{primes( 500 )};


  # is_prime returns 0 for composite, 2 for prime
  say "$n is prime"  if is_prime($n);

  # is_prob_prime returns 0 for composite, 2 for prime, and 1 for maybe prime
  say "$n is ", qw(composite maybe_prime? prime)[is_prob_prime($n)];


  # step to the next prime (returns 0 if the next one is more than ~0)
  $n = next_prime($n);

  # step back (returns 0 if given input less than 2)
  $n = prev_prime($n);


  # Return Pi(n) -- the number of primes E<lt>= n.
  $primepi = prime_count( 1_000_000 );
  $primepi = prime_count( 10**14, 10**14+1000 );  # also does ranges

  # Quickly return an approximation to Pi(n)
  my $approx_number_of_primes = prime_count_approx( 10**17 );

  # Lower and upper bounds.  lower <= Pi(n) <= upper for all n
  die unless prime_count_lower($n) <= prime_count($n);
  die unless prime_count_upper($n) >= prime_count($n);


  # Return p_n, the nth prime
  say "The ten thousandth prime is ", nth_prime(10_000);

  # Return a quick approximation to the nth prime
  say "The one trillionth prime is ~ ", nth_prime_approx(10**12);

  # Lower and upper bounds.   lower <= nth_prime(n) <= upper for all n
  die unless nth_prime_lower($n) <= nth_prime($n);
  die unless nth_prime_upper($n) >= nth_prime($n);


  # Get the prime factors of a number
  @prime_factors = factor( $n );


  # Precalculate a sieve, possibly speeding up later work.
  prime_precalc( 1_000_000_000 );

  # Free any memory used by the module.
  prime_memfree;

  # Alternate way to free.  When this leaves scope, memory is freed.
  my $mf = Math::Prime::Util::MemFree->new;


  # Random primes
  my $small_prime = random_prime(1000);      # random prime <= limit
  my $rand_prime = random_prime(100, 10000); # random prime within a range
  my $rand_prime = random_ndigit_prime(6);   # random 6-digit prime


=head1 DESCRIPTION

Pure Perl implementations of prime number utilities that are normally
handled with XS or GMP.  Having the Perl implementations (1) provides examples,
(2) allows the functions to run even if XS isn't available, and (3) gives
big number support if L<Math::Prime::Util::GMP> isn't available.  This is a
subset of L<Math::Prime::Util>'s functionality.

All routines should work with native integers or multi-precision numbers.  To
enable big numbers, use bigint or bignum:

    use bigint;
    say prime_count_approx(1000000000000000000000000)'
    # says 18435599767347543283712

This is still experimental, and some functions will be very slow.  The
L<Math::Prime::Util::GMP> module has much faster versions of many of these
functions.  Alternately, L<Math::Pari> has a lot of these types of functions.


=head1 FUNCTIONS

=head2 is_prime

  print "$n is prime" if is_prime($n);

Returns 2 if the number is prime, 0 if not.  For numbers larger than C<2^64>
it will return 0 for composite and 1 for probably prime, using a strong BPSW
test.  Also note there are probabilistic prime testing functions available.


=head2 primes

Returns all the primes between the lower and upper limits (inclusive), with
a lower limit of C<2> if none is given.

An array reference is returned (with large lists this is much faster and uses
less memory than returning an array directly).

  my $aref1 = primes( 1_000_000 );
  my $aref2 = primes( 1_000_000_000_000, 1_000_000_001_000 );

  my @primes = @{ primes( 500 ) };

  print "$_\n" for (@{primes( 20, 100 )});

Sieving will be done if required.  The algorithm used will depend on the range
and whether a sieve result already exists.  Possibilities include trial
division (for ranges with only one expected prime), a Sieve of Eratosthenes
using wheel factorization, or a segmented sieve.


=head2 next_prime

  $n = next_prime($n);

Returns the next prime greater than the input number.  If the input is not a
bigint, then 0 is returned if the next prime is larger than a native integer
type (the last representable primes being C<4,294,967,291> in 32-bit Perl and
C<18,446,744,073,709,551,557> in 64-bit).


=head2 prev_prime

  $n = prev_prime($n);

Returns the prime smaller than the input number.  0 is returned if the
input is C<2> or lower.


=head2 prime_count

  my $primepi = prime_count( 1_000 );
  my $pirange = prime_count( 1_000, 10_000 );

Returns the Prime Count function C<Pi(n)>, also called C<primepi> in some
math packages.  When given two arguments, it returns the inclusive
count of primes between the ranges (e.g. C<(13,17)> returns 2, C<14,17>
and C<13,16> return 1, and C<14,16> returns 0).


=head2 nth_prime

  say "The ten thousandth prime is ", nth_prime(10_000);

Returns the prime that lies in index C<n> in the array of prime numbers.  Put
another way, this returns the smallest C<p> such that C<Pi(p) E<gt>= n>.

This relies on generating primes, so can require a lot of time and space for
large inputs.


=head2 is_strong_pseudoprime

=head2 miller_rabin

  my $maybe_prime = is_strong_pseudoprime($n, 2);
  my $probably_prime = is_strong_pseudoprime($n, 2, 3, 5, 7, 11, 13, 17);

Takes a positive number as input and one or more bases.  The bases must be
greater than 1.  Returns 2 is C<n> is definitely prime, 1 if C<n>
is probably prime, and 0 if C<n> is definitely composite.  Since this is
just the Miller-Rabin test, a value of 2 is only returned for inputs of
2 and 3, which are shortcut.  If 0 is returned, then the number really is a
composite.  If 1 is returned, there is a good chance the number is prime
(depending on the input and the bases), but we cannot be sure.

This is usually used in combination with other tests to make either stronger
tests (e.g. the strong BPSW test) or deterministic results for numbers less
than some verified limit (such as the C<is_prob_prime> function in this module).


=head2 is_strong_lucas_pseudoprime

Takes a positive number as input, and returns 1 if the input is a strong
Lucas pseudoprime using the Selfridge method of choosing D, P, and Q (some
sources call this a strong Lucas-Selfridge pseudoprime).  This is one half
of the BPSW primality test (the Miller-Rabin strong pseudoprime test with
base 2 being the other half).

=head2 is_aks_prime

Takes a positive number as input, and returns 1 if the input can be proven
prime using the AKS primality test.  This code is included for completeness
and as an example, but is impractically slow.


=head1 UTILITY FUNCTIONS

=head2 prime_precalc

  prime_precalc( 1_000_000_000 );

Let the module prepare for fast operation up to a specific number.  It is not
necessary to call this, but it gives you more control over when memory is
allocated and gives faster results for multiple calls in some cases.  In the
current implementation this will calculate a sieve for all numbers up to the
specified number.


=head2 prime_memfree

  prime_memfree;

Frees any extra memory the module may have allocated.  Like with
C<prime_precalc>, it is not necessary to call this, but if you're done
making calls, or want things cleanup up, you can use this.  The object method
might be a better choice for complicated uses.


=head1 FACTORING FUNCTIONS

=head2 factor

  my @factors = factor(3_369_738_766_071_892_021);
  # returns (204518747,16476429743)

Produces the prime factors of a positive number input, in numerical order.
The special cases of C<n = 0> and C<n = 1> will return C<n>, which
guarantees multiplying the factors together will always result in the
input value, though those are the only cases where the returned factors
are not prime.


=head2 trial_factor

  my @factors = trial_factor($n);

Produces the prime factors of a positive number input.  The factors will be
in numerical order.  The special cases of C<n = 0> and C<n = 1> will return
C<n>, while with all other inputs the factors are guaranteed to be prime.
For large inputs this will be very slow.

=head2 fermat_factor

  my @factors = fermat_factor($n);

Currently unimplementated in Perl.

=head2 holf_factor

  my @factors = holf_factor($n);

Produces factors, not necessarily prime, of the positive number input.  An
optional number of rounds can be given as a second parameter.  It is possible
the function will be unable to find a factor, in which case a single element,
the input, is returned.  This uses Hart's One Line Factorization with no
premultiplier.  It is an interesting alternative to Fermat's algorithm,
and there are some inputs it can rapidly factor.  In the long run it has the
same advantages and disadvantages as Fermat's method.

=head2 squfof_factor

  my @factors = squfof_factor($n);

Currently unimplementated in Perl.

=head2 prho_factor

=head2 pbrent_factor

=head2 pminus1_factor

  my @factors = prho_factor($n);

  # Use a very small number of rounds
  my @factors = prho_factor($n, 1000);

Produces factors, not necessarily prime, of the positive number input.  An
optional number of rounds can be given as a second parameter.  These attempt
to find a single factor using one of the probabilistic algorigthms of
Pollard Rho, Brent's modification of Pollard Rho, or Pollard's C<p - 1>.
These are more specialized algorithms usually used for pre-factoring very
large inputs, or checking very large inputs for naive mistakes.  If the
input is prime or they run out of rounds, they will return the single
input value.  On some inputs they will take a very long time, while on
others they succeed in a remarkably short time.



=head1 MATHEMATICAL FUNCTIONS

=head2 ExponentialIntegral

  my $Ei = ExponentialIntegral($x);

Given a non-zero floating point input C<x>, this returns the real-valued
exponential integral of C<x>, defined as the integral of C<e^t/t dt>
from C<-infinity> to C<x>.
Depending on the input, the integral is calculated using
continued fractions (C<x E<lt> -1>),
rational Chebyshev approximation (C< -1 E<lt> x E<lt> 0>),
a convergent series (small positive C<x>),
or an asymptotic divergent series (large positive C<x>).

Accuracy should be at least 14 digits.


=head2 LogarithmicIntegral

  my $li = LogarithmicIntegral($x)

Given a positive floating point input, returns the floating point logarithmic
integral of C<x>, defined as the integral of C<dt/ln t> from C<0> to C<x>.
If given a negative input, the function will croak.  The function returns
0 at C<x = 0>, and C<-infinity> at C<x = 1>.

This is often known as C<li(x)>.  A related function is the offset logarithmic
integral, sometimes known as C<Li(x)> which avoids the singularity at 1.  It
may be defined as C<Li(x) = li(x) - li(2)>.

This function is implemented as C<li(x) = Ei(ln x)> after handling special
values.

Accuracy should be at least 14 digits.

=head2 RiemannZeta

  my $z = RiemannZeta($s);

Given a floating point input C<s> where C<s E<gt>= 0.5>, returns the floating
point value of ζ(s)-1, where ζ(s) is the Riemann zeta function.  One is
subtracted to ensure maximum precision for large values of C<s>.  The zeta
function is the sum from k=1 to infinity of C<1 / k^s>

Accuracy should be at least 14 digits, but currently does not increase
accuracy with big floats.  Small integer values are returned from a table,
values between 0.5 and 5 use rational Chebyshev approximation, and larger
values use a series.


=head2 RiemannR

  my $r = RiemannR($x);

Given a positive non-zero floating point input, returns the floating
point value of Riemann's R function.  Riemann's R function gives a very close
approximation to the prime counting function.

Accuracy should be at least 14 digits.


=head1 LIMITATIONS

The SQUFOF and Fermat factoring algorithms are not implemented yet.

Some of the prime methods use more memory than they should, as the segmented
sieve is not properly used in C<primes> and C<prime_count>.


=head1 PERFORMANCE

Performance compared to the XS/C code is quite poor for many operations.  Some
operations that are relatively close for small and medium-size values:

  next_prime / prev_prime
  is_prime / is_prob_prime
  miller_rabin
  ExponentialIntegral / LogarithmicIntegral / RiemannR
  primearray

Operations that are slower include:

  primes
  random_prime / random_ndigit_prime
  factor / all_factors
  nth_prime
  primecount

Performance improvement in this code is still possible.  The prime sieve is
over 2x faster than anything I was able to find online, but it is still has
room for improvement.

L<Math::Prime::Util::GMP> offers C<C+XS+GMP> support for most of the important
functions, and will be vastly faster for most operations.  If you install that
module, L<Math::Prime::Util> will load it automatically, meaning you should
not have to think about what code is actually being used (C, GMP, or Perl).

Memory use will generally be higher for the PP code, and in some cases B<much>
higher.  Some of this may be addressed in a later release.

For small values (e.g. primes and prime counts under 10M) most of this will
not matter.


=head1 SEE ALSO

L<Math::Prime::Util>

L<Math::Prime::Util::GMP>


=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>


=head1 COPYRIGHT

Copyright 2012 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
