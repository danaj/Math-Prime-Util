#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util ":all";
use Benchmark qw/:all/;
my $maxdigits = (~0 <= 4294967295) ? 10 : 20;
my $nnums = 100;

my $count = shift || -5;

srand(29);
my @darray;
push @darray, [gendigits($_)]  for (2 .. 10);
my $sum;

print "Direct sieving:\n";
cmpthese($count,{
    ' 2' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_segment_pi($_) for @{$darray[2-2]} },
    ' 3' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_segment_pi($_) for @{$darray[3-2]} },
    ' 4' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_segment_pi($_) for @{$darray[4-2]} },
    ' 5' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_segment_pi($_) for @{$darray[5-2]} },
    ' 6' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_segment_pi($_) for @{$darray[6-2]} },
    ' 7' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_segment_pi($_) for @{$darray[7-2]} },
    ' 8' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_segment_pi($_) for @{$darray[8-2]} },
    #' 9' => sub { $sum += Math::Prime::Util::_XS_segment_pi($_) for @{$darray[9-2]} },
    #'10' => sub { $sum += Math::Prime::Util::_XS_segment_pi($_) for @{$darray[10-2]} },
});
if (0) {
print "\n";
print "Direct Lehmer:\n";
cmpthese($count,{
    ' 2' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_lehmer_pi($_) for @{$darray[2-2]} },
    ' 3' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_lehmer_pi($_) for @{$darray[3-2]} },
    ' 4' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_lehmer_pi($_) for @{$darray[4-2]} },
    ' 5' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_lehmer_pi($_) for @{$darray[5-2]} },
    ' 6' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_lehmer_pi($_) for @{$darray[6-2]} },
    ' 7' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_lehmer_pi($_) for @{$darray[7-2]} },
    ' 8' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_lehmer_pi($_) for @{$darray[8-2]} },
    ' 9' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_lehmer_pi($_) for @{$darray[9-2]} },
    '10' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_lehmer_pi($_) for @{$darray[10-2]} },
});
}
print "\n";
print "Direct LMO:\n";
cmpthese($count,{
    ' 2' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_LMO_pi($_) for @{$darray[2-2]} },
    ' 3' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_LMO_pi($_) for @{$darray[3-2]} },
    ' 4' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_LMO_pi($_) for @{$darray[4-2]} },
    ' 5' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_LMO_pi($_) for @{$darray[5-2]} },
    ' 6' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_LMO_pi($_) for @{$darray[6-2]} },
    ' 7' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_LMO_pi($_) for @{$darray[7-2]} },
    ' 8' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_LMO_pi($_) for @{$darray[8-2]} },
    ' 9' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_LMO_pi($_) for @{$darray[9-2]} },
    '10' => sub { prime_memfree(); $sum += Math::Prime::Util::_XS_LMO_pi($_) for @{$darray[10-2]} },
});
print "\n";

sub gendigits {
  my $digits = shift;
  die "Digits must be > 0" unless $digits > 0;

  my $base = ($digits == 1) ? 0 : int(10 ** ($digits-1));
  my $max = int(10 ** $digits);
  $max = ~0 if $max > ~0;
  my @nums = map { $base+int(rand($max-$base)) } (1 .. $nnums);
  return @nums;
}
