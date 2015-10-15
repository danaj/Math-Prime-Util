#!/usr/bin/env perl
use warnings;
use strict;
use ntheory ":all";
use bigint;
$|=1;
prime_set_config(verbose=>0);

# Whether to output indices before the values
my $outbstyle = 0;


my $type = shift || die "Must supply type";
my $low = shift || 1;
my $high = shift || 1e9;


my $range = (($high-$low) > 1e15) ? 1e14 : 1+int(($high-$low)/100);
my %patterns = (
   # 2-tuples (twin)
  'A001359' => [2],
   # 3-tuples
  'A022004' => [2,6],
  'A022005' => [4,6],
   # 4-tuples
  'A007530' => [2,6,8],
   # 5-tuples
  'A022007' => [4,6,10,12],
  'A022006' => [2,6,8,12],
   # 6-tuples
  'A022008' => [4,6,10,12,16],
   # 7-tuples
  'A022009' => [2,6,8,12,18,20],
  'A022010' => [2,8,12,14,18,20],
   # 8-tuples
  'A022011' => [2,6,8,12,18,20,26],
  'A022012' => [2,6,12,14,20,24,26],
  'A022013' => [6,8,14,18,20,24,26],
   # 9-tuples
  'A022547' => [4,6,10,16,18,24,28,30],
  'A022548' => [4,10,12,18,22,24,28,30],
  'A022545' => [2,6,8,12,18,20,26,30],
  'A022546' => [2,6,12,14,20,24,26,30],
   # 10-tuples
  'A022569' => [2,6,8,12,18,20,26,30,32],
  'A022570' => [2,6,12,14,20,24,26,30,32],
   # 11-tuples
  'A213646' => [4,6,10,16,18,24,28,30,34,36],
  'A213647' => [2,6,8,12,18,20,26,30,32,36],
   # 12-tuples
  'A213601' => [6,10,12,16,22,24,30,34,36,40,42],
  'A213645' => [2,6,8,12,18,20,26,30,32,36,42],
   # 13-tuples
  'A214947' => [6,12,16,18,22,28,30,36,40,42,46,48],
  'A257137' => [4,6,10,16,18,24,28,30,34,40,46,48],
  'A257138' => [4,6,10,16,18,24,28,30,34,36,46,48],
  'A257139' => [2,6,8,12,18,20,26,30,32,36,42,48],
  'A257140' => [2,8,14,18,20,24,30,32,38,42,44,48],
  'A257141' => [2,12,14,18,20,24,30,32,38,42,44,48],
   # 14-tuples
  'A257167' => [2,6,8,12,18,20,26,30,32,36,42,48,50],
  'A257168' => [2,8,14,18,20,24,30,32,38,42,44,48,50],
   # other
  'A257304' => [2,6,8,12,18,20,26,30,32,36,42,48,50,56],
  'A257375' => [4,6,10,16,18,24,28,30,34,40,46,48,54,58,60,66],
  '5TP39' => [2,6,8,18,20,30,32,36,38],
);

die "Unknown type" unless exists $patterns{$type};

my @cl = @{ $patterns{$type} };

# 30 minutes on Macbook Pro to find first 52 entries of A213601:
# my $high = 25777719656829367;
# my @cl = (6,10,12,16,22,24,30,34,36,40,42);
# which makes it ~3-4x slower than JKA's old hand-tuned code.
#
# 69 seconds on Macbook Pro for the Federighi (5TP39) sequence:
# my $high = 1e14;
# my @cl = (2,6,8,18,20,30,32,36,38);
# which comes out to about 1.5x slower than JKA's 2007 result.


my $i = 0;
my @p;
while ($low < $high) {
  my $chigh = $low + $range - 1;
  $chigh = $high if $chigh > $high;
  # The GMP code will use more residues so favor it with big clusters
  if (scalar(@cl) > 9) {
    @p = Math::Prime::Util::GMP::sieve_prime_cluster($low, $chigh, @cl);
  } else {
    @p = sieve_prime_cluster($low, $chigh, @cl);
  }
  prime_set_config(verbose=>0);
  if ($outbstyle) {
    print ++$i," $_\n" for @p;
  } else {
    print "$_\n" for @p;
  }
  $low += $range;
}
