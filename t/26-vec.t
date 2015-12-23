#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/vecreduce
                         vecextract
                         vecmin vecmax
                         vecsum vecprod factorial
                         vecany vecall vecnotall vecnone vecfirst/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $use64 = Math::Prime::Util::prime_get_config->{'maxbits'} > 32;
$use64 = 0 if $use64 && 18446744073709550592 == ~0;

my @vecmins = (
  [ ],
  [ 1, 1 ],
  [ 0, 0 ],
  [ -1, -1 ],
  [ 1, 1, 2 ],
  [ 1, 2, 1 ],
  [ 1, 2, 1 ],
  [ -6, 0, 4, -5, 6, -6, 0 ],
  [ -6, 0, 4, -5, 7, -6, 0 ],
  [ "27944220269257565027", "81033966278481626507", "27944220269257565027" ],
);
if ($use64) {
  # List::Util::min gets these wrong
  push @vecmins, [ qw/18446744073702958477   18446744073704516093 18446744073706008451 18446744073706436837 18446744073707776433 18446744073702959347 18446744073702958477/ ];
  push @vecmins, [ qw/-9223372036852260731   -9223372036852260673 -9223372036852260731 -9223372036850511139 -9223372036850207017 -9223372036852254557 -9223372036849473359/ ];
  push @vecmins, [ qw/-9223372036853497843   9223372036852278343 -9223372036853497487 -9223372036844936897 -9223372036850971897 -9223372036853497843 9223372036848046999/ ];
}
my @vecmaxs = (
  [ ],
  [ 1, 1 ],
  [ 0, 0 ],
  [ -1, -1 ],
  [ 2, 1, 2 ],
  [ 2, 2, 1 ],
  [ 2, 2, 1 ],
  [  6, 0, 4, -5, 6, -6, 0 ],
  [  7, 0, 4, -5, 7, -8, 0 ],
  [ "81033966278481626507" , "27944220269257565027", "81033966278481626507" ],
);
if ($use64) {
  # List::Util::max gets these wrong
  push @vecmaxs, [ qw/18446744072030630259   18446744070011576186 18446744070972009258 18446744071127815503 18446744072030630259 18446744072030628952 18446744071413452589/ ];
  push @vecmaxs, [ qw/18446744073707508539   18446744073702156661 18446744073707508539 18446744073700111529 18446744073707506771 18446744073707086091 18446744073704381821/ ];
  push @vecmaxs, [ qw/-9223372036847631197   -9223372036853227739 -9223372036847631197 -9223372036851632173 -9223372036847631511 -9223372036852712261 -9223372036851707899/ ];
  push @vecmaxs, [ qw/9223372036846154833   -9223372036846673813 9223372036846154833 -9223372036851103423 9223372036846154461 -9223372036849190963 -9223372036847538803/ ];
}

my @vecsums = (
  [ 0 ],
  [ -1, -1 ],
  [ 0, 1,-1 ],
  [ 0, -1,1 ],
  [ 0, -1,1 ],
  [ 0, -2147483648,2147483648 ],
  [ 0, "-4294967296","4294967296" ],
  [ 0, "-9223372036854775808","9223372036854775808" ],
  [ "18446744073709551615", "18446744073709551615","-18446744073709551615","18446744073709551615" ],
  [ "55340232221128654848", "18446744073709551616","18446744073709551616","18446744073709551616" ],
);
if ($use64) {
  push @vecsums, [ "18446744073709620400", 18446744073709540400, (1000) x 80 ];
}
my @vecprods = (
  [ 1 ],
  [ 1,  1 ],
  [ -1,  -1 ],
  [ 2,  -1, -2 ],
  [ 2,  -1, -2 ],
  [ "-2147385345", 32767, -65535 ],
  [ "-2147385345", 32767, -65535 ],
  [ "-2147450880", 32768, -65535 ],
  [ "-2147483648", 32768, -65536 ],
);

plan tests => 0
            + scalar(@vecmins)
            + scalar(@vecmaxs)
            + scalar(@vecsums)
            + 1 + scalar(@vecprods)
            + 4    # vecreduce
            + 2    # vecextract
            + 3*4  # vec{any,all,notall,none}
            + 5    # vecfirst
            + 0;

###### vecmin
foreach my $r (@vecmins) {
  if (@$r == 0) {
    is(vecmin(), undef, "vecmin() = undef");
  } else {
    my($exp, @vals) = @$r;
    is( vecmin(@vals), $exp, "vecmin(@vals) = $exp" );
  }
}
###### vecmax
foreach my $r (@vecmaxs) {
  if (@$r == 0) {
    is(vecmax(), undef, "vecmax() = undef");
  } else {
    my($exp, @vals) = @$r;
    is( vecmax(@vals), $exp, "vecmax(@vals) = $exp" );
  }
}

###### vecsum
foreach my $r (@vecsums) {
  my($exp, @vals) = @$r;
  is( vecsum(@vals), $exp, "vecsum(@vals) = $exp" );
}
###### vecprod
foreach my $r (@vecprods) {
  my($exp, @vals) = @$r;
  is( vecprod(@vals), $exp, "vecprod(@vals) = $exp" );
}
{
  my(@prod,@fact);
  for my $f (0 .. 50) {
    push @fact, factorial($f);
    push @prod, vecprod(1 .. $f);
  }
  is_deeply(\@prod, \@fact, "vecprod matches factorial for 0 .. 50");
}

##### vecreduce
{
  my $fail = 0;
  is(vecreduce(sub{ $a + $b },()), undef, "vecreduce with empty list is undef");
  is(vecreduce(sub{ $fail = 1; 0; },(15)), 15+$fail, "vecreduce with (a) is a and does not call the sub");
  is(vecreduce(sub{ $a ^ $b },(4,2)), 6, "vecreduce [xor] (4,2) => 6");
  is(vecreduce(sub{ $a * $b**2 },(1, 17, 18, 19)), 17**2 * 18**2 * 19**2, "vecreduce product of squares");
}
###### vecextract
{
  is_deeply([vecextract(['a'..'z'],12345758)], [qw/b c d e h i n o s t u v x/], "vecextract bits");
  is(join("", vecextract(['a'..'z'],[22,14,17,10,18])), "works", "vecextract list");
}

###### vec{any,all,notall,none}
ok(  (vecany { $_ == 1 } 1, 2, 3), 'any true' );
ok( !(vecany { $_ == 1 } 2, 3, 4), 'any false' );
ok( !(vecany { 1 }), 'any empty list' );

ok(  (vecall { $_ == 1 } 1, 1, 1), 'all true' );
ok( !(vecall { $_ == 1 } 1, 2, 3), 'all false' );
ok(  (vecall { 1 }), 'all empty list' );

ok(  (vecnotall { $_ == 1 } 1, 2, 3), 'notall true' );
ok( !(vecnotall { $_ == 1 } 1, 1, 1), 'notall false' );
ok( !(vecnotall { 1 }), 'notall empty list' );

ok(  (vecnone { $_ == 1 } 2, 3, 4), 'none true' );
ok( !(vecnone { $_ == 1 } 1, 2, 3), 'none false' );
ok(  (vecnone { 1 }), 'none empty list' );

###### vecfirst
my $v;
$v = vecfirst { 8 == ($_ - 1) } 9,4,5,6; is($v, 9, "first success");
$v = vecfirst { 0 } 1,2,3,4; is($v, undef, "first failure");
$v = vecfirst { 0 }; is($v, undef, "first empty list");
$v = vecfirst { $_->[1] le "e" and "e" le $_->[2] } [qw(a b c)], [qw(d e f)], [qw(g h i)];
is_deeply($v, [qw(d e f)], 'first with reference args');
$v = vecfirst {while(1) {return ($_>6)} } 2,4,6,12; is($v,12,"first returns in loop");
