package Math::Prime::Util::PrimeArray;
use strict;
use warnings;

BEGIN {
  $Math::Prime::Util::PrimeArray::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::PrimeArray::VERSION = '0.45';
}

# parent is cleaner, and in the Perl 5.10.1 / 5.12.0 core, but not earlier.
# use parent qw( Exporter );
use base qw( Exporter );
our @EXPORT_OK = qw( );
our %EXPORT_TAGS = (all => [ @EXPORT_OK ]);


use Math::Prime::Util qw/nth_prime nth_prime_upper nth_prime_lower primes prime_precalc next_prime prev_prime/;
use Tie::Array;
use Carp qw/carp croak confess/;

use constant SEGMENT_SIZE  =>  10_000;

sub TIEARRAY {
  my $class = shift;
  if (@_) {
    croak "usage: tie ARRAY, '" . __PACKAGE__ . "";
  }
  return bless {
    # used to keep track of shift
    SHIFTINDEX => 0,
    # Remove all extra prime memory when we go out of scope
    MEMFREE    => Math::Prime::Util::MemFree->new,
    # A chunk of primes
    PRIMES     => [2, 3, 5, 7, 11, 13, 17],
    # What's the index of the first one?
    BEG_INDEX  => 0,
    # What's the index of the last one?
    END_INDEX  => 6,
    # positive = forward, negative = backward, 0 = random
    ACCESS_TYPE => 0,
  }, $class;
}
sub STORE     { carp "You cannot write to the prime array"; }
sub DELETE    { carp "You cannot write to the prime array"; }
sub STORESIZE { carp "You cannot write to the prime array"; }
sub EXISTS    { 1 }
#sub EXTEND    { my $self = shift; my $count = shift; prime_precalc($count); }
sub EXTEND    { 1 }
sub FETCHSIZE { 0x7FFF_FFFF }   # Even on 64-bit
# Simple FETCH:
# sub FETCH { return nth_prime($_[1]+1); }

sub FETCH {
  my ($self, $index) = @_;
  # We actually don't get negative indices -- they get turned into big numbers
  croak "Negative index given to prime array" if $index < 0;
  $index += $self->{SHIFTINDEX};  # take into account any shifts
  my $begidx = $self->{BEG_INDEX};
  my $endidx = $self->{END_INDEX};

  if ( $index < $begidx || $index > $endidx ) {

    if ($index == $endidx+1) {       # Forward iteration

      $self->{ACCESS_TYPE}++;
      if ($self->{ACCESS_TYPE} > 2) {
        my $end_prime = nth_prime_upper($index + SEGMENT_SIZE);
        $self->{PRIMES} = primes( $self->{PRIMES}->[-1]+1, $end_prime );
        $begidx = $endidx+1;
      } else {
        push @{$self->{PRIMES}}, next_prime($self->{PRIMES}->[-1]);
      }

    } elsif ($index == $begidx-1) {  # Backward iteration

      $self->{ACCESS_TYPE}--;
      if ($self->{ACCESS_TYPE} < -2) {
        my $beg_prime = $index <= SEGMENT_SIZE
                               ?  2  :  nth_prime_lower($index - SEGMENT_SIZE);
        $self->{PRIMES} = primes($beg_prime, $self->{PRIMES}->[0]-1);
        $begidx -= scalar @{ $self->{PRIMES} };
      } else {
        $begidx--;
        unshift @{$self->{PRIMES}}, prev_prime($self->{PRIMES}->[0]);
      }

    } else {                         # Random access

      $self->{ACCESS_TYPE} = int($self->{ACCESS_TYPE} / 2);
      # Alternately we could get a small window, but that will be quite
      # a bit slower if true random access.
      $begidx = $index;
      $self->{PRIMES} = [nth_prime($begidx+1)];

    }
    $self->{BEG_INDEX} = $begidx;
    $self->{END_INDEX} = $begidx + scalar @{$self->{PRIMES}} - 1;
  }
  return $self->{PRIMES}->[ $index - $begidx ];
}

# Fake out shift and unshift
sub SHIFT {
  my $self = shift;
  my $head = $self->FETCH(0);
  $self->{SHIFTINDEX}++;
  $head;
}
sub UNSHIFT {
  my ($self, $shiftamount) = @_;
  $shiftamount = 1 unless defined $shiftamount;
  $self->{SHIFTINDEX} = ($shiftamount >= $self->{SHIFTINDEX})
                        ? 0
                        : $self->{SHIFTINDEX} - $shiftamount;
  $self->FETCHSIZE;
}
# CLEAR this
# PUSH this, LIST
# POP this
# SPLICE this, offset, len, LIST
# DESTROY this
# UNTIE this

1;

__END__


# ABSTRACT: A tied array for primes

=pod

=head1 NAME

Math::Prime::Util::PrimeArray - A tied array for primes


=head1 VERSION

Version 0.45


=head1 SYNOPSIS

  use Math::Prime::Util::PrimeArray;

  # Create:
  tie my @primes, 'Math::Prime::Util::PrimeArray';

  # Use in a loop by index:
  for my $n (0..9) {
    print "prime $n = $primes[$n]\n";
  }

  # Use in a loop over array:
  for my $p (@primes) {
    last if $p > 1000;   # stop sometime
    print "$p\n";
  }

  # Use via array slice:
  print join(",", @primes[0..49]), "\n";

  # Use via each:
  use 5.012;
  while( my($index,$value) = each @primes ) {
    last if $value > 1000;   # stop sometime
    print "The ${index}th prime is $value\n";
  }

  # Use with shift:
  while ((my $p = shift @primes) < 1000) {
    print "$p\n";
  }


=head1 DESCRIPTION

An array that acts like the infinite set of primes.  This may be more
convenient than using L<Math::Prime::Util> directly, and in some cases it can
be faster than calling C<next_prime> and C<prev_prime>.

If the access pattern is ascending or descending, then a window is sieved and
results returned from the window as needed.  If the access pattern is random,
then C<nth_prime> is used.

Shifting acts like the array is losing elements at the front, so after two
shifts, C<$primes[0] == 5>.  Unshift will move the internal shift index back
one, unless given an argument which is the number to move back.  It will
not shift past the beginning, so C<unshift @primes, ~0> is a useful way to
reset from any shifts.

Example:

  say shift @primes;     # 2
  say shift @primes;     # 3
  say shift @primes;     # 5
  say $primes[0];        # 7
  unshift @primes;       #     back up one
  say $primes[0];        # 5
  unshift @primes, 2;    #     back up two
  say $primes[0];        # 2

If you want sequential primes with low memory, I recommend using
L<Math::Prime::Util/forprimes>.  It is much faster, as the tied array
functionality in Perl is not high performance.  It isn't as flexible as
the prime array, but it is a very common pattern.

If you prefer an iterator pattern, I would recommend using
L<Math::Prime::Util/prime_iterator>.  It will be a bit faster than using this
tied array, but of course you don't get random access.  If you find yourself
using the C<shift> operation, consider the iterator.


=head1 LIMITATIONS

The size of the array will always be shown as 2147483647 (IV32 max), even in
a 64-bit environment where primes through C<2^64> are available.

Some people find the idea of shifting a prime array abhorrent, as after
two shifts, "the second prime is 7?!".  If this bothers you, do not use
C<shift> on the tied array.


=head1 PERFORMANCE

  MPU forprimes:  forprimes { $sum += $_ } nth_prime(100_000);
  MPU iterator:   my $it = prime_iterator; $sum += $it->() for 1..100000;
  MPU array:      $sum = vecsum( @{primes(nth_prime(100_000))} );
  MPUPA:          tie my @primes, ...; $sum += $primes[$_] for 0..99999;
  MNSP:           my $seq = Math::NumSeq::Primes->new;
                  $sum += ($seq->next)[1] for 1..100000;
  MPTA:           tie my @primes, ...; $sum += $primes[$_] for 0..99999;
  List::Gen       $sum = primes->take(100000)->sum

Memory use is comparing the delta between just loading the module and running
the test.  Perl 5.20.0, Math::NumSeq v70, Math::Prime::TiedArray v0.04,
List::Gen 0.974.

Summing the first 0.1M primes via walking the array:

       6ms     56k    Math::Prime::Util      forprimes
       7ms    4 MB    Math::Prime::Util      sum big array
     110ms      0     Math::Prime::Util      prime_iterator
     155ms    644k    Math::Prime::Util::PrimeArray
     120ms   1476k    Math::NumSeq::Primes   sequence iterator
    4656ms   32 MB    List::Gen              sequence
    7540ms   61 MB    Math::Prime::TiedArray (extend 1k)

Summing the first 1M primes via walking the array:

      0.07s   268k    Math::Prime::Util      forprimes
      0.09s  41 MB    Math::Prime::Util      sum big array
      1.4s      0     Math::Prime::Util      prime_iterator
      1.6s    644k    Math::Prime::Util::PrimeArray
      7.1s   2428k    Math::NumSeq::Primes   sequence iterator
     91.3s   93 MB    List::Gen              sequence
    108.4s  760 MB    Math::Prime::TiedArray (extend 1k)

Summing the first 10M primes via walking the array:

      0.7s    432k    Math::Prime::Util      forprimes
      0.9s  394 MB    Math::Prime::Util      sum big array
     16.9s      0     Math::Prime::Util      prime_iterator
     15.4s    772k    Math::Prime::Util::PrimeArray
   3680  s  11.1MB    Math::NumSeq::Primes   sequence iterator
   5555  s  874 MB    List::Gen              sequence
          >5000 MB    Math::Primes::TiedArray (extend 1k)

L<Math::Prime::Util> offers three obvious solutions: a big array, an iterator,
and the C<forprimes> construct.  The big array is fast but uses a B<lot> of
memory, forcing the user to start programming segments.  Using the iterator
avoids all the memory use, but isn't as fast (this may improve in a later
release, as this is a new feature).  The C<forprimes> construct is both fast
and low memory, but it isn't quite as flexible as the iterator (most notably
there is no way to exit early, and it doesn't lend itself to wrapping inside
a filter).

L<Math::NumSeq::Primes> offers an iterator alternative, and works quite well
as long as you don't need lots of primes.  It does not support random access.
It has reasonable performance for the first few hundred thousand, but each
successive value takes much longer to generate, and once past 1 million it
isn't very practical.

L<Math::Primes::TiedArray> is remarkably impractical for anything other
than tiny numbers.

L<List::Gen> includes a built-in prime sequence.  It uses an inefficient
Perl sieve for numbers below 10M and has some performance defects when
getting large numbers of primes.

=head1 SEE ALSO

This module uses L<Math::Prime::Util> to do all the work.  If you're doing
anything but retrieving primes, you should examine that module to see if it
has functionality you can use directly, as it may be a lot faster or easier.

Similar functionality can be had from L<Math::NumSeq>
and L<Math::Prime::TiedArray>.

=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>


=head1 COPYRIGHT

Copyright 2012-2013 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
