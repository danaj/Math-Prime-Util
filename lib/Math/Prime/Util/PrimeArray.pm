package Math::Prime::Util::PrimeArray;
use strict;
use warnings;

BEGIN {
  $Math::Prime::Util::PrimeArray::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::PrimeArray::VERSION = '0.73';
}

# parent is cleaner, and in the Perl 5.10.1 / 5.12.0 core, but not earlier.
# use parent qw( Exporter );
use base qw( Exporter );
our @EXPORT_OK = qw(@primes @prime @pr @p $probj);
our %EXPORT_TAGS = (all => [ @EXPORT_OK ]);

# It would be nice to do this dynamically.
our(@primes, @prime, @pr, @p, $probj);
sub import {
  tie @primes, __PACKAGE__ if grep { $_ eq '@primes' } @_;
  tie @prime , __PACKAGE__ if grep { $_ eq '@prime'  } @_;
  tie @pr    , __PACKAGE__ if grep { $_ eq '@pr'     } @_;
  tie @p     , __PACKAGE__ if grep { $_ eq '@p'      } @_;
  $probj = __PACKAGE__->TIEARRAY if grep { $_ eq '$probj' } @_;
  goto &Exporter::import;
}

use Math::Prime::Util qw/nth_prime nth_prime_upper nth_prime_lower primes prime_precalc next_prime prev_prime/;
use Math::Prime::Util::MemFree;
use Tie::Array;
use Carp qw/carp croak confess/;

use constant SEGMENT_SIZE  =>  80_000;
use constant HALFSEG       =>  SEGMENT_SIZE >> 1;
use constant ALLOW_SKIP    =>  3_000;     # Sieve if skipping up to this

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
  $index = 0xFFFFFFFF + $index + 1 if $index < 0;
  $index += $self->{SHIFTINDEX};  # take into account any shifts
  my $begidx = $self->{BEG_INDEX};
  my $endidx = $self->{END_INDEX};

  if ( $index < $begidx || $index > $endidx ) {

    if ($index > $endidx && $index < $endidx + ALLOW_SKIP) { # Forward iteration

      $self->{ACCESS_TYPE}++;
      if ($self->{ACCESS_TYPE} > 2 || $index > $endidx+1) {
        my $prlen = scalar @{$self->{PRIMES}};
        # Keep up to HALFSEG elements from the previous array
        if ($prlen > HALFSEG) {
          @{$self->{PRIMES}} = @{$self->{PRIMES}}[-HALFSEG .. -1];
          $begidx += $prlen - HALFSEG;
        }
        # Add HALFSEG elements to the end
        my $end_prime = nth_prime_upper($index + HALFSEG);
        push @{$self->{PRIMES}}, @{primes($self->{PRIMES}->[-1]+1, $end_prime)};
      } else {
        push @{$self->{PRIMES}}, next_prime($self->{PRIMES}->[-1]);
      }

    } elsif ($index < $begidx && $index > $begidx - ALLOW_SKIP) { # Bk iteration

      $self->{ACCESS_TYPE}--;
      if ($self->{ACCESS_TYPE} < -2 || $index < $begidx-1) {
        my $prlen = scalar @{$self->{PRIMES}};
        my $beg_prime = $index <= HALFSEG
                               ?  2  :  nth_prime_lower($index - HALFSEG);
        unshift @{$self->{PRIMES}}, @{primes($beg_prime, $self->{PRIMES}->[0]-1)};
        my $prnewlen = scalar @{$self->{PRIMES}};
        $begidx -= $prnewlen - $prlen;
        $#{$self->{PRIMES}} = SEGMENT_SIZE-1 if $prnewlen > SEGMENT_SIZE;
      } else {
        $begidx--;
        unshift @{$self->{PRIMES}}, prev_prime($self->{PRIMES}->[0]);
      }

    } else {                         # Random access

      $self->{ACCESS_TYPE} = int($self->{ACCESS_TYPE} / 2);
      # TODO: we are destroying the primes array, just to get $begidx set.
      # We should instead have an additional single last-index-result.
      # return nth_prime($index+1);
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

Version 0.73


=head1 SYNOPSIS

  # Use package and create a tied variable
  use Math::Prime::Util::PrimeArray;
  tie my @primes, 'Math::Prime::Util::PrimeArray';

  # or all in one (allowed: @primes, @prime, @pr, @p):
  use Math::Prime::Util::PrimeArray '@primes';

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

Perl will mask all array arguments to 32-bit, making C<2^32-1> the maximum
prime through the standard array interface.  It will silently wrap after
that.  The only way around this is using the object interface:

    use Math::Prime::Util::PrimeArray;
    my $o = tie my @primes, 'Math::Prime::Util::PrimeArray';
    say $o->FETCH(2**36);

Here we store the object returned by tie, allowing us to call its FETCH
method directly.  This is actually faster than using the array.

Some people find the idea of shifting a prime array abhorrent, as after
two shifts, "the second prime is 7?!".  If this bothers you, do not use
C<shift> on the tied array.


=head1 PERFORMANCE

Performance of tied arrays increased substantially (40% faster) between
Perl v5.18 and 5.24.  It is recommended to use a new-ish Perl.

  sumprimes:      sum_primes(nth_prime(100_000))
  MPU forprimes:  forprimes { $sum += $_ } nth_prime(100_000);
  MPU iterator:   my $it = prime_iterator; $sum += $it->() for 1..100000;
  MPU array:      $sum = vecsum( @{primes(nth_prime(100_000))} );
  MPUPA:          tie my @prime, ...; $sum += $prime[$_] for 0..99999;
  MPUPA-FETCH:    my $o=tie my @pr, ...; $sum += $o->FETCH($_) for 0..99999;
  MNSP:           my $seq = Math::NumSeq::Primes->new;
                  $sum += ($seq->next)[1] for 1..100000;
  MPTA:           tie my @prime, ...; $sum += $prime[$_] for 0..99999;
  List::Gen       $sum = primes->take(100000)->sum

Memory use is comparing the delta between just loading the module and running
the test.  M1 Macbook, Perl 5.42.0, Math::NumSeq v75,
Math::Prime::TiedArray v0.04 with C<extend_step 1000>, List::Gen 0.979.

Summing the first 0.1M primes via walking the array (milliseconds):

       .05      56k    Math::Prime::Util      sumprimes
      1.7       56k    Math::Prime::Util      forprimes
      1.6      4 MB    Math::Prime::Util      sum big array
     12          0     Math::Prime::Util      prime_iterator
     31        3 MB    MPU::PrimeArray        using FETCH
     41        3 MB    MPU::PrimeArray        array
     63        6 MB    List::Gen              sequence
     51        950k    Math::NumSeq::Primes   sequence iterator
   2367ms     78 MB    Math::Prime::TiedArray (extend 1k)

Summing the first 1M primes via walking the array (seconds):

      .0003    268k    Math::Prime::Util      sumprimes
      .018     268k    Math::Prime::Util      forprimes
      .015    41 MB    Math::Prime::Util      sum big array
     0.11        0     Math::Prime::Util      prime_iterator
     0.3       644k    MPU::PrimeArray        using FETCH
     0.4       644k    MPU::PrimeArray        array
     0.8      57 MB    List::Gen              sequence
     4.3      3179k    Math::NumSeq::Primes   sequence iterator
    35.9s    722 MB    Math::Prime::TiedArray (extend 1k)

Summing the first 10M primes via walking the array (seconds):

     0.0015    432k    Math::Prime::Util      sumprimes
     0.19      432k    Math::Prime::Util      forprimes
     0.16    394 MB    Math::Prime::Util      sum big array
     1.2         0     Math::Prime::Util      prime_iterator
     3.0       772k    MPU::PrimeArray        using FETCH
     4.0       772k    MPU::PrimeArray        array
     8.3s    652 MB    List::Gen              sequence
   577       22.8MB    Math::NumSeq::Primes   sequence iterator
           >5000 MB    Math::Prime::TiedArray (extend 1k)

L<Math::Prime::Util> offers four obvious solutions:
the C<sum_primes> function,
summing a big generated array,
an iterator,
and the C<forprimes> construct.
The big array is fast but uses a B<lot> of memory, forcing the user to
start programming segments.  Using the iterator avoids all the memory
use, but isn't as fast.
The C<forprimes> construct is both fast and low memory, but it isn't
quite as flexible as the iterator.

L<Math::NumSeq::Primes> offers an iterator alternative, and works quite well
as long as you don't need lots of primes.  It does not support random access.
It has reasonable performance for the first few hundred thousand, but each
successive value takes much longer to generate, and once past 1 million it
isn't very practical.  Internally it is sieving all primes up to C<n> every
time it makes a new segment which is why it slows down so much.

L<List::Gen> includes a built-in prime sequence.  Version 0.975 will use
this module for primes if it can, which is shown in the above numbers.
It is the odd module out in this comparison, as primes aren't a core feature.
Without this module, it is very slow.

L<Math::Prime::TiedArray> is remarkably impractical for anything other
than tiny numbers.


=head1 SEE ALSO

This module uses L<Math::Prime::Util> to do all the work.  If you're doing
anything but retrieving primes, you should examine that module to see if it
has functionality you can use directly, as it may be a lot faster or easier.

Similar functionality can be had from L<Math::NumSeq>
and L<Math::Prime::TiedArray>.

=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>


=head1 COPYRIGHT

Copyright 2012-2026 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
