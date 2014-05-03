package Math::Prime::Util::PrimeArray;
use strict;
use warnings;

BEGIN {
  $Math::Prime::Util::PrimeArray::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::PrimeArray::VERSION = '0.41';
}

# parent is cleaner, and in the Perl 5.10.1 / 5.12.0 core, but not earlier.
# use parent qw( Exporter );
use base qw( Exporter );
our @EXPORT_OK = qw( );
our %EXPORT_TAGS = (all => [ @EXPORT_OK ]);


use Math::Prime::Util qw/nth_prime nth_prime_upper nth_prime_lower primes prime_precalc next_prime prev_prime/;
use Tie::Array;
use Carp qw/carp croak confess/;

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
  my $self = shift;
  my $index = shift;
  # We actually don't get negative indices -- they get turned into big numbers
  croak "Negative index given to prime array" if $index < 0;
  $index += $self->{SHIFTINDEX};  # take into account any shifts
  my $begidx = $self->{BEG_INDEX};
  my $endidx = $self->{END_INDEX};

  if ( $index < $begidx || $index > $endidx ) {

    if ($index == $endidx+1) {       # Forward iteration

      $self->{ACCESS_TYPE}++;
      if ($self->{ACCESS_TYPE} > 2) {
        my $end_prime = nth_prime_upper($index + 10_000);
        $self->{PRIMES} = primes( $self->{PRIMES}->[-1]+1, $end_prime );
        $begidx = $endidx+1;
      } else {
        push @{$self->{PRIMES}}, next_prime($self->{PRIMES}->[-1]);
      }

    } elsif ($index == $begidx-1) {  # Backward iteration

      $self->{ACCESS_TYPE}--;
      if ($self->{ACCESS_TYPE} < -2) {
        my $num = 10_000;
        my $beg_prime = $index <= $num ? 2 : nth_prime_lower($index - $num );
        $self->{PRIMES} = primes($beg_prime, $self->{PRIMES}->[0]-1);
        $begidx -= scalar @{ $self->{PRIMES} };
      } else {
        $begidx--;
        unshift @{$self->{PRIMES}}, prev_prime($self->{PRIMES}->[0]);
      }

    } else {                         # Random access

      $self->{ACCESS_TYPE} = int($self->{ACCESS_TYPE} / 2);
      # Alternately we could get a small window
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
  my $self = shift;
  my $shiftamount = defined $_[0] ? shift : 1;
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

Version 0.41


=head1 SYNOPSIS

  use Math::Prime::Util::PrimeArray;

  # Create:
  tie my @primes, 'Math::Prime::Util::PrimeArray';

  # Use in a loop by index:
  for my $n (1..10) {
    print "prime $n = $primes[$n]\n";
  }

  # Use in a loop over array:
  for my $p (@primes) {
    print "$p\n";
    last if $p > $limit;   # stop sometime
  }

  # Use via array slice:
  print join(",", @primes[0..50]), "\n";

  # Use via each:
  use 5.012;
  while( my($index,$value) = each @primes ) {
    print "The ${index}th prime is $value\n";
    last if $p > $limit;   # stop sometime
  }

  # Use with shift:
  while ((my $p = shift @primes) < $limit) {
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
one, unless given an argument which is the number to move back (it silently
truncates so it does not shift past the beginning).
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

There are some people that find the idea of shifting a prime array abhorrent,
as after two shifts, "the second prime is 7?!".  If this bothers you, do not
use C<shift> on the tied array.


=head1 PERFORMANCE

  MPU forprimes:  forprimes { $sum += $_ } nth_prime(100_000);
  MPU iterator:   my $it = prime_iterator; $sum += $it->() for 1..100000;
  MPU array:      $sum += $_ for @{primes(nth_prime(100_000))};
  MPUPA:          tie my @primes, ...; $sum += $primes[$_] for 0..99999;
  MNSP:           my $seq = Math::NumSeq::Primes->new;
                  $sum += ($seq->next)[1] for 1..100000;
  MPTA:           tie my @primes, ...; $sum += $primes[$_] for 0..99999;

Memory use is comparing the delta between just loading the module and running
the test.  Perl 5.19.2, Math::NumSeq v61, Math::Prime::TiedArray v0.04.

Summing the first 0.1M primes via walking the array:

       7ms     52k    Math::Prime::Util      forprimes
     140ms      0     Math::Prime::Util      prime_iterator
      12ms   4400k    Math::Prime::Util      sum big array
     220ms    840k    Math::Prime::Util::PrimeArray
     130ms    280k    Math::NumSeq::Primes   sequence iterator
    7560ms   65 MB    Math::Prime::TiedArray (extend 1k)

Summing the first 1M primes via walking the array:

      0.1s    300k    Math::Prime::Util      forprimes
      1.8s      0     Math::Prime::Util      prime_iterator
      0.2s   40 MB    Math::Prime::Util      sum big array
      1.9s   1.1MB    Math::Prime::Util::PrimeArray
      7.5s   1.2MB    Math::NumSeq::Primes   sequence iterator
    110.5s  785 MB    Math::Prime::TiedArray (extend 1k)

Summing the first 10M primes via walking the array:

      0.8s   5.9MB    Math::Prime::Util      forprimes
     22.4s      0     Math::Prime::Util      prime_iterator
      1.5s  368 MB    Math::Prime::Util      sum big array
     19.1s   1.2MB    Math::Prime::Util::PrimeArray
   3680  s  11.1MB    Math::NumSeq::Primes   sequence iterator
          >5000 MB    Math::Primes::TiedArray (extend 1k)

L<Math::Prime::Util> offers three obvious solutions: a big array, an iterator,
and the C<forprimes> construct.  The big array is fast but uses a B<lot> of
memory, forcing the user to start programming segments.  Using the iterator
avoids all the memory use, but isn't as fast (this may improve in a later
release, as this is a new feature).  The C<forprimes> construct is by far
the fastest, but it isn't quite as flexible as the iterator (most notably
there is no way to exit early, and it doesn't lend itself to wrapping inside
a filter).

L<Math::NumSeq::Primes> offers an iterator alternative, and works quite well
for reasonably small numbers.  It does not support random access.  It is
very fast for small values, but is very slow with large counts.

L<Math::Primes::TiedArray> is remarkably impractical for anything other
than very small numbers.


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
