package Math::Prime::Util::EllipticCurve;
use strict;
use warnings;
use Carp qw/carp croak confess/;

if (!defined $Math::BigInt::VERSION) {
  eval { require Math::BigInt;   Math::BigInt->import(try=>'GMP,Pari'); 1; }
  or do { croak "Cannot load Math::BigInt"; };
}

BEGIN {
  $Math::Prime::Util::EllipticCurve::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::EllipticCurve::VERSION = '0.26';
}

# Pure perl (with Math::BigInt) manipulation of Elliptic Curves
# in projective coordinates.  Should be split into a point class.

sub new {
  my ($class, $a, $b, $n) = @_;
  $a = Math::BigInt->new("$a") unless ref($a) eq 'Math::BigInt';
  $b = Math::BigInt->new("$b") unless ref($b) eq 'Math::BigInt';
  $n = Math::BigInt->new("$n") unless ref($n) eq 'Math::BigInt';

  my $self = {
    a => $a,
    b => $b,
    n => $n,
  };

  bless $self, $class;
  return $self;
}

sub double_p {
  my ($self, $x, $z) = @_;
  my $n = $self->{'n'};

  my $u = ( ($x+$z) * ($x+$z) ) % $n;
  my $v = ( ($x-$z) * ($x-$z) ) % $n;
  my $w = $u - $v;

  return ( ($u*$v)%$n , ($w*($v+$w*$self->{'b'}))%$n );
}

sub add3_p {
  my ($self, $x1, $z1, $x2, $z2, $xin, $zin) = @_;
  my $n = $self->{'n'};

  my $u = (($x2 - $z2) * ($x1 + $z1) ) % $n;
  my $v = (($x2 + $z2) * ($x1 - $z1) ) % $n;

  my $upv2 = (($u+$v) * ($u+$v)) % $n;
  my $umv2 = (($u-$v) * ($u-$v)) % $n;

  return ( ($upv2*$zin) % $n, ($umv2*$xin) % $n );
}

sub mul_p {
  my ($self, $k, $x, $z) = @_;
  $x = Math::BigInt->new("$x") unless ref($x) eq 'Math::BigInt';
  $z = Math::BigInt->new("$z") unless ref($z) eq 'Math::BigInt';
  my ($x1, $x2, $z1, $z2);

  my $r = --$k;
  my $l = -1;
  while ($r != 1) { $r >>= 1; $l++ }
  if ($k & (1 << $l)) {
    ($x2, $z2) = $self->double_p($x, $z);
    ($x1, $z1) = $self->add3_p($x2, $z2, $x, $z, $x, $z);
    ($x2, $z2) = $self->double_p($x2, $z2);
  } else {
    ($x1, $z1) = $self->double_p($x, $z);
    ($x2, $z2) = $self->add3_p($x, $z, $x1, $z1, $x, $z);
  }
  $l--;
  while ($l >= 1) {
    if ($k & (1 << $l)) {
      ($x1, $z1) = $self->add3_p($x1, $z1, $x2, $z2, $x, $z);
      ($x2, $z2) = $self->double_p($x2, $z2);
    } else {
      ($x2, $z2) = $self->add3_p($x2, $z2, $x1, $z1, $x, $z);
      ($x1, $z1) = $self->double_p($x1, $z1);
    }
    $l--;
  }
  if ($k & 1) {
    ($x, $z) = $self->double_p($x2, $z2);
  } else {
    ($x, $z) = $self->add3_p($x2, $z2, $x1, $z1, $x, $z);
  }
  return ($x, $z);
}

sub add_a {
  my ($self, $P1x, $P1y, $P2x, $P2y) = @_;
  my $n = $self->{'n'};

  if ($P1x == $P2x) {
    my $t = ($P1y + $P2y) % $n;
    return (0, 1) if $t == 0;
  }
  my $deltax = ($P2x - $P1x) % $n;
  $deltax->bmodinv($n);
  return (Math::BigInt->bzero,Math::BigInt->bone) if $deltax eq "NaN";

  my $deltay = ($P2y - $P1y) % $n;
  my $m = ($deltay * $deltax) % $n;   # m = deltay / deltax

  my $x = ($m*$m - $P1x - $P2x) % $n;
  my $y = ($m*($P1x - $x) - $P1y) % $n;
  return ($x,$y);
}

sub double_a {
  my ($self, $P1x, $P1y) = @_;
  my $n = $self->{'n'};

  my $m = 2*$P1y;
  $m->bmodinv($n);
  return (Math::BigInt->bzero,Math::BigInt->bone) if $m eq "NaN";

  $m = ((3*$P1x*$P1x + $self->{'a'}) * $m) % $n;

  my $x = ($m*$m - 2*$P1x) % $n;
  my $y = ($m*($P1x - $x) - $P1y) % $n;
  return ($x,$y);
}

sub mul_a {
  my ($self, $k, $x, $y) = @_;
  $x = Math::BigInt->new("$x") unless ref($x) eq 'Math::BigInt';
  $y = Math::BigInt->new("$y") unless ref($y) eq 'Math::BigInt';
  my $n = $self->{'n'};

  my $Bx = Math::BigInt->bzero;
  my $By = Math::BigInt->bone;
  my $v = 1;

  while ($v && $k > 0) {
    if ( ($k % 2) != 0) {
      $k--;
      my $d = Math::BigInt::bgcd( ($Bx - $x) % $n, $n);
      $v = ($d == 1 || $d == $n);
      if    ($x == 0 && $y == 1)   { }
      elsif ($Bx == 0 && $By == 1) { ($Bx,$By) = ($x,$y); }
      elsif ($v)                   { ($Bx,$By) = $self->add_a($x,$y,$Bx,$By); }
    } else {
      $k >>= 1;
      my $d = Math::BigInt::bgcd( 2*$y % $n, $n);
      $v = ($d == 1 || $d == $n);
      if ($v) { ($x,$y) = $self->double_a($x,$y); }
    }
  }
  return ($Bx, $By);
}

1;

__END__


# ABSTRACT: Elliptic curve operations

=pod

=encoding utf8


=head1 NAME

Math::Prime::Util::EllipticCurve - Elliptic curve operations


=head1 VERSION

Version 0.26


=head1 SYNOPSIS

Todo.

  Todo
  # TO DO
  To do

=head1 DESCRIPTION

This really should just be in Math::EllipticCurve.

To write.


=head1 FUNCTIONS

=head2 new

  $point = Math::Prime::Util::EllipticCurve->new(a, b);

Returns a new curve defined by a and b.


=head1 SEE ALSO

L<Math::EllipticCurve::Prime>


=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>


=head1 COPYRIGHT

Copyright 2012-2013 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
