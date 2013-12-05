package Math::Prime::Util::ECAffinePoint;
use strict;
use warnings;
use Carp qw/carp croak confess/;

if (!defined $Math::BigInt::VERSION) {
  eval { require Math::BigInt;   Math::BigInt->import(try=>'GMP,Pari'); 1; }
  or do { croak "Cannot load Math::BigInt"; };
}

BEGIN {
  $Math::Prime::Util::ECAffinePoint::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::ECAffinePoint::VERSION = '0.35';
}

# Pure perl (with Math::BigInt) manipulation of Elliptic Curves
# in affine coordinates.  Should be split into a point class.

sub new {
  my ($class, $a, $b, $n, $x, $y) = @_;
  $a = Math::BigInt->new("$a") unless ref($a) eq 'Math::BigInt';
  $b = Math::BigInt->new("$b") unless ref($b) eq 'Math::BigInt';
  $n = Math::BigInt->new("$n") unless ref($n) eq 'Math::BigInt';
  $x = Math::BigInt->new("$x") unless ref($x) eq 'Math::BigInt';
  $y = Math::BigInt->new("$y") unless ref($y) eq 'Math::BigInt';

  croak "n must be >= 2" unless $n >= 2;
  $a->bmod($n);
  $b->bmod($n);

  my $self = {
    a => $a,
    b => $b,
    n => $n,
    x => $x,
    y => $y,
    f => $n->copy->bone,
  };

  bless $self, $class;
  return $self;
}

sub _add {
  my ($self, $P1x, $P1y, $P2x, $P2y) = @_;
  my $n = $self->{'n'};

  if ($P1x == $P2x) {
    my $t = ($P1y + $P2y) % $n;
    return (Math::BigInt->bzero,Math::BigInt->bone) if $t == 0;
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

sub _double {
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

sub mul {
  my ($self, $k) = @_;
  my $x = $self->{'x'};
  my $y = $self->{'y'};
  my $n = $self->{'n'};
  my $f = $self->{'f'};

  my $Bx = $n->copy->bzero;
  my $By = $n->copy->bone;
  my $v = 1;

  while ($k > 0) {
    if ( ($k % 2) != 0) {
      $k--;
      $f->bmul($Bx-$x)->bmod($n);
      if    ($x == 0 && $y == 1)   { }
      elsif ($Bx == 0 && $By == 1) { ($Bx,$By) = ($x,$y); }
      else                         { ($Bx,$By) = $self->_add($x,$y,$Bx,$By); }
    } else {
      $k >>= 1;
      $f->bmul(2*$y)->bmod($n);
      ($x,$y) = $self->_double($x,$y);
    }
  }
  $f = Math::BigInt::bgcd( $f, $n);
  $self->{'x'} = $Bx;
  $self->{'y'} = $By;
  $self->{'f'} = $f;
  return $self;
}

sub add {
  my ($self, $other) = @_;
  croak "add takes a EC point"
    unless ref($other) eq 'Math::Prime::Util::ECAffinePoint';
  croak "second point is not on the same curve"
    unless $self->{'a'} == $other->{'a'} &&
           $self->{'b'} == $other->{'b'} &&
           $self->{'n'} == $other->{'n'};

  ($self->{'x'}, $self->{'y'}) = $self->_add($self->{'x'}, $self->{'y'},
                                             $other->{'x'}, $other->{'y'});
  return $self;
}


sub a { return shift->{'a'}; }
sub b { return shift->{'b'}; }
sub n { return shift->{'n'}; }
sub x { return shift->{'x'}; }
sub y { return shift->{'y'}; }
sub f { return shift->{'f'}; }

sub is_infinity {
  my $self = shift;
  return ($self->{'x'}->is_zero() && $self->{'y'}->is_one());
}

1;

__END__


# ABSTRACT: Elliptic curve operations for affine points

=pod

=encoding utf8

=for stopwords mul

=head1 NAME

Math::Prime::Util::ECAffinePoint - Elliptic curve operations for affine points


=head1 VERSION

Version 0.35


=head1 SYNOPSIS

  # Create a point on a curve (a,b,n) with coordinates 0,1
  my $ECP = Math::Prime::Util::ECAffinePoint->new($a, $b, $n, 0, 1);

  # scalar multiplication by k.
  $ECP->mul($k)

  # add two points on the same curve
  $ECP->add($ECP2)

  say "P = O" if $ECP->is_infinity();

=head1 DESCRIPTION

This really should just be in Math::EllipticCurve.

To write.


=head1 FUNCTIONS

=head2 new

  $point = Math::Prime::Util::ECAffinePoint->new(a, b, n, x, y);

Returns a new point at C<(x,y,1)> on the curve C<(a,b,n)>.

=head2 a

=head2 b

=head2 n

Returns the C<a>, C<b>, or C<n> values that describe the curve.

=head2 x

=head2 y

Returns the C<x> or C<y> values that define the point on the curve.

=head2 f

Returns a possible factor found during EC multiplication.

=head2 add

Takes another point on the same curve as an argument and adds it this point.

=head2 mul

Takes an integer and performs scalar multiplication of the point.

=head2 is_infinity

Returns true if the point is (0,1), which is the point at infinity for
the affine coordinates.


=head1 SEE ALSO

L<Math::EllipticCurve::Prime>

This really should just be in a L<Math::EllipticCurve> module.

=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>


=head1 COPYRIGHT

Copyright 2012-2013 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
