package Math::Prime::Util::ECM;
use strict;
use warnings;
use Carp qw/croak/;

BEGIN {
  $Math::Prime::Util::ECM::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::ECM::VERSION = '0.75';
}

BEGIN {
  use constant OLD_PERL_VERSION=> $] < 5.008;
  use constant MPU_MAXBITS     => (~0 == 4294967295) ? 32 : 64;
  use constant MPU_32BIT       => MPU_MAXBITS == 32;
  use constant INTMAX          => (!OLD_PERL_VERSION || MPU_32BIT) ? ~0 : 562949953421312;
  use constant SINTMAX         => (INTMAX >> 1);
}

*tobigint = \&Math::Prime::Util::_to_bigint;
*getconfig = \&Math::Prime::Util::prime_get_config;

*Mdivint = \&Math::Prime::Util::divint;
*Madd1int = \&Math::Prime::Util::add1int;
*Mmulint = \&Math::Prime::Util::mulint;
*Minvmod = \&Math::Prime::Util::invmod;
*Mgcd = \&Math::Prime::Util::gcd;
*Mis_prime = \&Math::Prime::Util::is_prime;
*Mprimes = \&Math::Prime::Util::primes;
*Msqrtint = \&Math::Prime::Util::sqrtint;

my(%TINY_STAGE1_PLAN, %TINY_STAGE2_PLAN);
Math::Prime::Util::PP::_register_free_sub(sub {
  %TINY_STAGE1_PLAN = ();
  %TINY_STAGE2_PLAN = ();
}) if defined &Math::Prime::Util::PP::_register_free_sub;

my @TINY_ECM_SIGMAS = (
   11,   13,   17,   19,   23,   29,   31,   37,   41,   43,
   47,   53,   59,   61,   67,   71,   73,   79,   83,   89,
  103,  127,  139,  149,  151,  157,  163,  181,  191,  197,
  199,  211,  223,  227,  233,  239,  257,  271,  293,  313,
  331,  337,  347,  359,  379,  389,  397,  401,  409,  421,
  443,  449,  457,  479,  487,  509,  521,  523,  547,  557,
  587,  641,  653,  659,  673,  677,  683,  691,  719,  727,
  739,  751,  769,  797,  809,  853,  919,  929,  941,  997,
 1049, 1051, 1063, 1091, 1093, 1109, 1117, 1129, 1153, 1201,
 1217, 1229, 1283, 1327, 1361, 1381, 1427, 1447, 1459, 1471,
 1481, 1489, 1543, 1549, 1571, 1621, 1709, 1723, 1753, 1759,
 1801, 1811, 1867, 1987, 2039, 2099, 2113, 2131, 2251, 2309,
 2347, 2381, 2399, 2447, 2473, 2551, 2557, 2663, 2677, 2689,
 2713, 2719, 2749, 2857, 2879, 2887, 2939, 3001, 3061, 3067,
 3121, 3137, 3187, 3251, 3259, 3271, 3307, 3359, 3371, 3373,
 3467, 3593, 3607, 3623, 3643, 3709, 3733, 3793, 3851, 3923,
 3989, 4019, 4049, 4129, 4231, 4253, 4283, 4339, 4349, 4441,
 4523, 4649, 4787, 4987, 4999, 5171, 5237, 5273, 5297, 5333,
 5387, 5471, 5479, 5647, 5749, 5791, 6101, 6163, 6257, 6299,
 6337, 6451, 6491, 6659, 6793, 6823, 6967, 7013, 7229, 7253,
 7333, 7369, 7477, 7621, 7793, 7817, 8059, 8167, 8209, 8263,
 8311, 8377, 8573, 8641, 8741, 8837, 8863, 8963, 9001, 9151,
 9203, 9433, 9697, 9743, 9781, 9883,10007,10069,10099,10139,
10163,10193,10267,10429,10457,10487,10691,10837,10949,11087,
11243,11321,11411,11681,11813,11903,12011,12263,12277,12401,
12409,12437,12479,12569,12619,12739,12911,13331,13367,13537,
13721,13789,13841,13873,14051,14149,14221,14419,14431,14827,
14887,15077,15289,15467,15511,15649,15773,15797,15859,15901,
16057,16141,16217,16529,16547,16553,16619,17299,17393,17419,
17449,17737,17921,18049,18223,19073,19183,19477,20021,20323,
);

sub _found_factor {
  my($f, $n, $what) = @_;
  return ($n) if $f == 1 || $f == $n;
  my $f2 = Mdivint($n,$f);
  croak "internal error in $what" unless Mmulint($f,$f2) == $n;
  ($f,$f2) = ($f2,$f) if $f > $f2;
  print "$what found factor $f\n" if getconfig()->{'verbose'} > 0;
  ($f,$f2);
}

sub _check_bounds {
  my($B1, $B2, $ncurves) = @_;
  croak "ecm_factor: B1 must fit in native signed integer" if $B1 > SINTMAX;
  croak "ecm_factor: B2 must fit in native signed integer" if $B2 > SINTMAX;
  croak "ecm_factor: ncurves must fit in native signed integer" if $ncurves > SINTMAX;
}

sub _tiny_double {
  my($x, $z, $a24, $n) = @_;
  my $u = $x - $z;  $u = ($u * $u) % $n;
  my $v = $x + $z;  $v = ($v * $v) % $n;
  my $w = $v - $u;
  my $x2 = ($u * $v) % $n;
  my $z2 = ($w * ($u + (($a24 * $w) % $n))) % $n;
  ($x2, $z2);
}

sub _tiny_dadd {
  my($px, $pz, $qx, $qz, $dx, $dz, $n) = @_;
  my $u = (($px - $pz) * ($qx + $qz)) % $n;
  my $v = (($px + $pz) * ($qx - $qz)) % $n;
  my $s = $u + $v;
  my $d = $u - $v;
  my $rx = ((($s * $s) % $n) * $dz) % $n;
  my $rz = ((($d * $d) % $n) * $dx) % $n;
  ($rx, $rz);
}

sub _tiny_stage1_plan {
  my($B1) = @_;
  return $TINY_STAGE1_PLAN{$B1} if exists $TINY_STAGE1_PLAN{$B1};

  my $sqrt_b1 = int(sqrt($B1));
  my @bprimes = @{ Mprimes(2, $B1) };
  for my $p (@bprimes) {
    last if $p > $sqrt_b1;
    my($k, $pm) = ($p, int($B1 / $p));
    while ($k <= $pm) { $k *= $p; }
    $p = $k;
  }
  $TINY_STAGE1_PLAN{$B1} = \@bprimes;
}

sub _tiny_stage2_plan {
  my($B1, $B2) = @_;
  my $key = "$B1:$B2";
  return $TINY_STAGE2_PLAN{$key} if exists $TINY_STAGE2_PLAN{$key};

  my $D = Msqrtint($B2 >> 1);
  $D = Madd1int($D) if $D % 2;

  my @b2primes = @{ Mprimes($B1+1, $B2) };
  my @windows;
  my $m = 1;
  while ($m < ($B2+$D)) {
    my(@left, @right);
    if ($m+$D > $B1) {
      my($lo, $hi) = ($m-$D, $m+$D);
      for my $p (@b2primes) {
        next if $p < $lo;
        last if $p > $hi;
        if ($p < $m) {
          push @left, $m-$p;
        } elsif ($p > $m) {
          if ($p <= $m+$m) {
            my $mirror = $m+$m-$p;
            next if $mirror > $B1 && $mirror >= $lo && Mis_prime($mirror);
          }
          push @right, $p-$m;
        }
      }
    }
    push @windows, [ $m, \@left, \@right ];
    $m += 2*$D;
  }

  $TINY_STAGE2_PLAN{$key} = { D => $D, windows => \@windows };
}

sub _tiny_mul {
  my($x, $z, $k, $a24, $n) = @_;
  return ($x, $z) if $k == 1;

  my($x0, $z0) = ($x, $z);
  my($x1, $z1) = _tiny_double($x0, $z0, $a24, $n);
  return ($x1, $z1) if $k == 2;

  my $bit = 1;
  $bit <<= 1 while (($bit << 1) <= $k);
  for ($bit >>= 1; $bit; $bit >>= 1) {
    if ($k & $bit) {
      ($x0, $z0) = _tiny_dadd($x0, $z0, $x1, $z1, $x, $z, $n);
      ($x1, $z1) = _tiny_double($x1, $z1, $a24, $n);
    } else {
      ($x1, $z1) = _tiny_dadd($x0, $z0, $x1, $z1, $x, $z, $n);
      ($x0, $z0) = _tiny_double($x0, $z0, $a24, $n);
    }
  }
  ($x0, $z0);
}

sub _tiny_batch_normalize_x {
  my($xs, $zs, $n) = @_;
  my $npoints = scalar @$xs;
  my($acc, @prefix, @xout) = (1);

  return (1, []) if $npoints == 0;

  for my $i (0 .. $npoints-1) {
    $acc = ($acc * $zs->[$i]) % $n;
    $prefix[$i] = $acc;
  }

  my $inv = Minvmod($acc, $n);
  unless (defined $inv) {
    my $g = Mgcd($acc, $n);
    return (0, $g) if $g != 1 && $g != $n;

    # A product can contain different factors whose combined gcd is n.
    # Check each coordinate to recover any proper factor before giving up.
    for my $z (@$zs) {
      $g = Mgcd($z, $n);
      return (0, $g) if $g != 1 && $g != $n;
    }
    return (0, 0);
  }

  for my $i (reverse 0 .. $npoints-1) {
    my $prev = $i ? $prefix[$i-1] : 1;
    my $zinv = ($inv * $prev) % $n;
    $inv = ($inv * $zs->[$i]) % $n;
    $xout[$i] = ($xs->[$i] * $zinv) % $n;
  }

  (1, \@xout);
}

sub _tiny_stage2 {
  my($x, $z, $a24, $n, $plan, $cinfo) = @_;

  my $f = Mgcd($z, $n);
  return _found_factor($f, $n, "ECM S2 $cinfo") if $f != 1 && $f != $n;
  return (0) if $f == $n;

  my $D = $plan->{'D'};

  my $g = 1;
  my(@nqx, @nqz, @Sx, @Sz);

  @nqx = (0, $x);
  @nqz = (0, $z);
  foreach my $i (2 .. 2*$D) {
    if ($i % 2) {
      ($nqx[$i], $nqz[$i]) = _tiny_dadd($nqx[($i-1)/2], $nqz[($i-1)/2],
                                        $nqx[($i+1)/2], $nqz[($i+1)/2],
                                        $x, $z, $n);
    } else {
      ($nqx[$i], $nqz[$i]) = _tiny_double($nqx[$i/2], $nqz[$i/2], $a24, $n);
    }
  }

  @Sx = ($x);
  @Sz = ($z);
  {
    my($xm_x, $xm_z) = ($nqx[2*$D-1], $nqz[2*$D-1]);
    for my $w (1 .. $#{ $plan->{'windows'} }) {
      my($oldx, $oldz) = ($Sx[$w-1], $Sz[$w-1]);
      ($Sx[$w], $Sz[$w]) = _tiny_dadd($nqx[2*$D], $nqz[2*$D],
                                      $oldx, $oldz,
                                      $xm_x, $xm_z, $n);
      ($xm_x, $xm_z) = ($oldx, $oldz);
    }
  }

  {
    # Keep all 1..2D points in the batch.  The upper half is not needed for
    # the final products, but its Z coordinates can expose a factor.  S[0]
    # is the same base point as nQ[1], so normalize it only once.
    my @xs = @nqx[1 .. 2*$D];
    my @zs = @nqz[1 .. 2*$D];
    if (@Sx > 1) {
      push @xs, @Sx[1 .. $#Sx];
      push @zs, @Sz[1 .. $#Sz];
    }
    my($ok, $xout) = _tiny_batch_normalize_x(\@xs, \@zs, $n);
    return (0) if !$ok && $xout == 0;
    return _found_factor($xout, $n, "ECM S2 $cinfo") unless $ok;

    my @allx = @$xout;
    my @nqout = splice(@allx, 0, 2*$D);
    @nqx = (0, @nqout);
    @Sx = ($nqx[1], @allx);
  }

  for my $wi (0 .. $#{ $plan->{'windows'} }) {
    my $window = $plan->{'windows'}[$wi];
    my($m, $left, $right) = @$window;
    if (@$left || @$right) {
      foreach my $i (@$left)  { $g = ($g * (($Sx[$wi] - $nqx[$i]) % $n)) % $n; }
      foreach my $i (@$right) { $g = ($g * (($Sx[$wi] - $nqx[$i]) % $n)) % $n; }
      $f = Mgcd($g, $n);
      last if $f != 1;
    }
  }

  return (0) if $f == $n;
  return _found_factor($f, $n, "ECM S2 $cinfo") if $f != 1;
  ();
}

sub _factor_tiny {
  my($n, $B1, $B2, $ncurves, $sigma_offset) = @_;
  croak "_factor_tiny internal error: parameter not defined"
    unless defined $n && defined $B1 && defined $B2 &&
           defined $ncurves && defined $sigma_offset;

  _check_bounds($B1, $B2, $ncurves);
  $n = tobigint($n);
  my $bprimes = _tiny_stage1_plan($B1);
  my $stage2_plan = ($B2 > $B1) ? _tiny_stage2_plan($B1, $B2) : undef;

  CURVE:
  for my $ci (1 .. $ncurves) {
    my $curve = $sigma_offset + $ci-1;
    my $sigma = $TINY_ECM_SIGMAS[$curve % @TINY_ECM_SIGMAS];
    my $cinfo = "$B1/$B2/$ncurves curve $ci [$curve]";

    my $u = tobigint($sigma*$sigma - 5);
    my $v = tobigint(4 * $sigma);
    my $umv = $u - $v;
    my $t3uv = 3*$u + $v;

    my $u2 = ($u * $u) % $n;
    my $u3 = ($u2 * $u) % $n;
    my $v2 = ($v * $v) % $n;
    my $v3 = ($v2 * $v) % $n;
    my $umv2 = ($umv * $umv) % $n;
    my $umv3 = ($umv2 * $umv) % $n;

    my $num = -(($umv3 * $t3uv) % $n);
    $num %= $n;
    my $den = (($u3 * $v * 16) % $n);
    my $g = Mgcd($den, $n);
    next CURVE if $g == $n;
    return _found_factor($g, $n, "ECM $cinfo") if $g != 1;

    my $den_inv = Minvmod($den, $n);
    next CURVE unless defined $den_inv;
    my $a24 = ($num * $den_inv) % $n;
    my($x, $z) = ($u3, $v3);

    my $j = 0;
    for my $k (@$bprimes) {
      ($x, $z) = _tiny_mul($x, $z, $k, $a24, $n);
      if (($j++ % 64) == 0) {
        $g = Mgcd($z, $n);
        next CURVE if $g == $n;
        return _found_factor($g, $n, "ECM $cinfo") if $g != 1;
      }
    }
    $g = Mgcd($z, $n);
    next CURVE if $g == $n;
    return _found_factor($g, $n, "ECM $cinfo") if $g != 1;

    if (defined $stage2_plan) {
      my @factors = _tiny_stage2($x, $z, $a24, $n, $stage2_plan, $cinfo);
      return @factors if @factors > 1;
      next CURVE if @factors == 1 && $factors[0] == 0;
    }
  }
  ($n);
}

sub ecm_factor_pp {
  my($n, $B1, $B2, $ncurves, $sigma_offset) = @_;
  $ncurves = 10 unless defined $ncurves;
  $sigma_offset = 0 unless defined $sigma_offset;

  if (!defined $B1) {
    my $try_b1 = 500;
    for (1 .. 7) {
      my @factors = _factor_tiny($n, $try_b1, 20*$try_b1, 40, $sigma_offset);
      return @factors if @factors > 1;
      $sigma_offset += 40;
      $try_b1 *= 4;
    }
    return ($n);
  }

  $B2 = 20*$B1 unless defined $B2 && $B2 > 0;
  _factor_tiny($n, $B1, $B2, $ncurves, $sigma_offset);
}

1;

__END__

=pod

=head1 NAME

Math::Prime::Util::ECM - Pure Perl ECM factoring support

=head1 SYNOPSIS

Internal support module for L<Math::Prime::Util>.

=head1 DESCRIPTION

This module contains the pure Perl two-stage ECM implementation used by
L<Math::Prime::Util::PP>.  It is not a public API.

=head1 FUNCTIONS

=head2 ecm_factor_pp

Internal implementation for pure Perl ECM factoring.

=head1 SEE ALSO

L<Math::Prime::Util>

L<Math::Prime::Util::PP>

=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>

=head1 COPYRIGHT

Copyright 2026 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
