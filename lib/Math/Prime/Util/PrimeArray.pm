package Math::Prime::Util::PrimeArray;
use strict;
use warnings;

BEGIN {
  $Math::Prime::Util::PrimeArray::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::PrimeArray::VERSION = '0.14';
}

# parent is cleaner, and in the Perl 5.10.1 / 5.12.0 core, but not earlier.
# use parent qw( Exporter );
use base qw( Exporter );
our @EXPORT_OK = qw( );
our %EXPORT_TAGS = (all => [ @EXPORT_OK ]);


use Math::Prime::Util qw/nth_prime nth_prime_upper primes prime_precalc/;
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
  if ( ($index < $self->{BEG_INDEX}) || ($index > $self->{END_INDEX}) ) {
    # We're going to get a chunk of primes starting 1000 before the one
    # asked for, and extending at _least_ 1000 past it.  So given a number
    # past index 1000, we would expect to have ~2000 primes with the one
    # requested being about in the middle.  This should be reasonably fast
    # for the siever and amortize the cost over a lot of calls if they're
    # walking us.  Yes, we'll get well over 2000 primes as the numbers get
    # very large, but that's ok (e.g. 80k primes at index 400M).  Calling
    # nth_prime on indices that large is very time consuming.
    my $start_idx = ($index < 1000) ? 0 : $index-1000;
    my $start_prime = nth_prime($start_idx+1);       # This is exact
    my $end_prime   = nth_prime_upper($index+5000);  # This is a bound
    #print "calling primes from $start_idx/$start_prime to $end_prime\n";
    $self->{PRIMES} = primes($start_prime, $end_prime);
    $self->{BEG_INDEX} = $start_idx;
    $self->{END_INDEX} = $start_idx + scalar @{$self->{PRIMES}} - 1;;
  }
  return $self->{PRIMES}->[ $index - $self->{BEG_INDEX} ];
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

Version 0.14


=head1 SYNOPSIS

  use Math::Prime::Util::PrimeArray;

  # Create:
  my @primes;  tie @primes, 'Math::Prime::Util::PrimeArray';

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

Internally when an index is accessed, an area surrounding the index is sieved
if necessary, then the result returned.  This means random access will be a
little slower than using C<nth_prime> directly, but not by very much.
Random access in a small window (1000 or so primes in either direction) will
be very fast, as will sequential access in either direction.

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
  

=head1 LIMITATIONS

The size of the array will always be shown as 2147483647 (IV32 max), even in
a 64-bit environment where primes through C<2^64> are available.


=head1 PERFORMANCE

Summing the first 0.1M primes via walking the array:

      40ms   3.3 MB    $sum += $_ for @{primes(nth_prime(100_000))};
     300ms   0.6 MB    Math::Prime::Util::PrimeArray
     230ms   2.2 MB    Math::NumSeq::Primes
   10300ms  36.2 MB    Math::Primes::TiedArray (extend 1k)

Summing the first 1M primes via walking the array:

     0.3s   35.5 MB    $sum += $_ for @{primes(nth_prime(1_000_000))};
     3.9s    1.3 MB    Math::Prime::Util::PrimeArray
     9.6s    2.2 MB    Math::NumSeq::Primes
   146.7s  444.9 MB    Math::Primes::TiedArray (extend 1k)

Summing the first 10M primes via walking the array:

       4s    365 MB    $sum += $_ for @{primes(nth_prime(10_000_000))};
     2m        2 MB    Math::Prime::Util::PrimeArray
    85m        2 MB    Math::NumSeq::Primes
           >5000 MB    Math::Primes::TiedArray (extend 1k)

Using L<Math::Prime::Util> directly in a naive fashion uses lots of memory,
but is extremely fast.  Sieving segments at a time would control the memory
use, which is one thing the C<PrimeArray> tie is trying to do for you (but
adds more inefficiency than is ideal).

L<Math::NumSeq::Primes> offers an iterator alternative, and works quite well
for reasonably small numbers.  It does not, however, support random access.
There seems to be a 2MB fixed overhead, but after that the memory use is 
is well constrained.  It is very fast for small values, but clearly is
getting slower as we sum to 1 million, and takes well over an hour to count
to 10 million.

L<Math::Primes::TiedArray> is remarkably impractical for anything other
than very small numbers.  I believe the times and memory use in the above
tables illustrate this.


=head1 SEE ALSO

This module uses L<Math::Prime::Util> to do all the work.  If you're doing
anything but retrieving primes, you should examine that module to see if it
has functionality you can use directly, as it may be a lot faster or easier.

Similar functionality can be had from L<Math::NumSeq>
and L<Math::Prime::TiedArray>.

=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>


=head1 COPYRIGHT

Copyright 2012 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
