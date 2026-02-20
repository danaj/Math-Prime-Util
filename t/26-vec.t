#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/vecreduce
                         vecextract
                         vecequal
                         vecmin vecmax
                         vecsum vecprod factorial
                         vecuniq
                         vecsingleton
                         vecfreq
                         vecslide
                         vecsort vecsorti
                         vecany vecall vecnotall vecnone vecfirst vecfirstidx/;

# vecmex      in t/26-mex.t
# vecpmex     in t/26-mex.t
# vecsample   in t/26-randperm.t

# [related]
# setcontains       return 0 if we are given something NOT in SETA
# setcontainsany    return 1 if we are given anything in SETA

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
  [ "-18446744073709551615", "-9223372036854775807","-9223372036854775807",-1 ],
  [ "-18446744073709551616", "-9223372036854775807","-9223372036854775807",-2 ],
  [ "-18446744073709551617", "-9223372036854775807","-9223372036854775807",-3 ],
  [ "-18446744073709551616", "-9223372036854775808","-9223372036854775808" ],
  [ "-9223372036854775807", "-9223372036854775807",0 ],
  [ "-9223372036854775808", "-9223372036854775807",-1 ],
  [ "-9223372036854775809", "-9223372036854775807",-2 ],
  [ "-9223372036854775808", 0,"-9223372036854775808",0 ],
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

my @vecsorts = (
  [ [], [], "empty input" ],
  [ [0], [0], "single input" ],
  [ [5,0], [0,5], "two positive inputs" ],
  [ [~0,-1], [-1,~0], "-1 and maxuv" ],
  [ ["-9223372036854775808","-9223372036854775809"], ["-9223372036854775809","-9223372036854775808"], "two large negative inputs" ],
  [ [qw/9223372036854775807 18446744073709551615 127 32767 2147483647 1 140737488355327/],
    [qw/1 127 32767 2147483647 140737488355327 9223372036854775807 18446744073709551615/],
    "various 64-bit positive inputs" ],
  [ [qw/18446744073709551615 19342813113834066795298815 4722366482869645213695 65535 4294967295/],
    [qw/65535 4294967295 18446744073709551615 4722366482869645213695 19342813113834066795298815/],
    "large string inputs" ],
  [ [qw/9223372036854775812 9223372036854775809 9223372036854775810 9223372036854775811/],
    [qw/9223372036854775809 9223372036854775810 9223372036854775811 9223372036854775812/],
    "integers over 2^63 broken before 5.26.0" ],
);

plan tests => 1    # vecmin
            + 1    # vecmax
            + 1    # vecsum
            + 1    # vecprod
            + 1    # vecreduce
            + 1    # vecextract
            + 1    # vecequal
            + 1    # vec{any,all,notall,none}
            + 1    # vecfirst
            + 1    # vecfirstidx
            + 1    # vecuniq
            + 1    # vecsingleton
            + 1    # vecfreq
            + 1    # vecsort
            + 1    # vecslide
            + 0;

###### vecmin
subtest 'vecmin', sub {
  foreach my $r (@vecmins) {
    if (@$r == 0) {
      is(vecmin(), undef, "vecmin() = undef");
    } else {
      my($exp, @vals) = @$r;
      is("".vecmin(@vals), $exp, "vecmin(@vals) = $exp");
    }
  }
};
###### vecmax
subtest 'vecmax', sub {
  foreach my $r (@vecmaxs) {
    if (@$r == 0) {
      is(vecmax(), undef, "vecmax() = undef");
    } else {
      my($exp, @vals) = @$r;
      is("".vecmax(@vals), $exp, "vecmax(@vals) = $exp");
    }
  }
};

###### vecsum
subtest 'vecsum', sub {
  foreach my $r (@vecsums) {
    my($exp, @vals) = @$r;
    is( "".vecsum(@vals), $exp, "vecsum(@vals) = $exp" );
  }
};
###### vecprod
subtest 'vecprod', sub {
  foreach my $r (@vecprods) {
    my($exp, @vals) = @$r;
    is( "".vecprod(@vals), $exp, "vecprod(@vals) = $exp" );
  }
  my(@prod,@fact);
  for my $f (0 .. 50) {
    push @fact, "".factorial($f);
    push @prod, "".vecprod(1 .. $f);
  }
  is_deeply(\@prod, \@fact, "vecprod matches factorial for 0 .. 50");
};

##### vecreduce
subtest 'vecreduce', sub {
  my $fail = 0;
  is(vecreduce(sub{ $a + $b },()), undef, "vecreduce with empty list is undef");
  is(vecreduce(sub{ $fail = 1; 0; },(15)), 15+$fail, "vecreduce with (a) is a and does not call the sub");
  is(vecreduce(sub{ $a ^ $b },(4,2)), 6, "vecreduce [xor] (4,2) => 6");
  is(vecreduce(sub{ $a * $b**2 },(1, 17, 18, 19)), 17**2 * 18**2 * 19**2, "vecreduce product of squares");
};
###### vecextract
subtest 'vecextract', sub {
  is_deeply([vecextract(['a'..'z'],12345758)], [qw/b c d e h i n o s t u v x/], "vecextract bits");
  is(join("", vecextract(['a'..'z'],[22,14,17,10,18])), "works", "vecextract list");
};

###### vecequal
subtest 'vecequal', sub {
  is(vecequal([],[]), 1, "vecequal([],[]) = 1");
  is(vecequal([undef],[undef]), 1, "vecequal([undef],[undef]) = 1");
  is(vecequal([0],[0]), 1, "vecequal([0],[0]) = 1");
  is(vecequal([undef],[]), 0, "vecequal([undef],[]) = 0");
  is(vecequal([undef],[0]), 0, "vecequal([undef],[0]) = 0");
  is(vecequal([0],[[]]), 0, "vecequal([0],[[]]) = 0");
  is(vecequal([],[[]]), 0, "vecequal([],[[]]) = 0");
  is(vecequal([0],["a"]), 0, "vecequal([0],[\"a\"]) = 0");

  is(vecequal([1,2,3],[1,2,3]), 1, "vecequal([1,2,3],[1,2,3]) = 1");
  is(vecequal([1,2,3],[3,2,1]), 0, "vecequal([1,2,3],[3,2,1]) = 0");
  is(vecequal([-1,2,3],[-1,2,3]), 1, "vecequal([-1,2,3],[-1,2,3]) = 1");
  is(vecequal([undef,[1,2],"a"],[undef,[1,2],"a"]), 1, "vecequal([undef,[1,2],\"a\"],[undef,[1,2],\"a\"] = 1");

  is(vecequal(\@vecsums, \@vecsums), 1, "vecequal = 1 for vecsums");
  is(vecequal(\@vecsums, \@vecprods), 0, "vecequal = 0 for vecsums");
};

###### vec{any,all,notall,none}
subtest 'vecany vecall vecnotall vecnone', sub {
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
};

###### vecfirst
subtest 'vecfirst', sub {
  my $v;
  $v = vecfirst { 8 == ($_ - 1) } 9,4,5,6; is($v, 9, "first success");
  $v = vecfirst { 0 } 1,2,3,4; is($v, undef, "first failure");
  $v = vecfirst { 0 }; is($v, undef, "first empty list");
  $v = vecfirst { $_->[1] le "e" and "e" le $_->[2] } [qw(a b c)], [qw(d e f)], [qw(g h i)];
  is_deeply($v, [qw(d e f)], 'first with reference args');
  $v = vecfirst {while(1) {return ($_>6)} } 2,4,6,12; is($v,12,"first returns in loop");
};
subtest 'vecfirstidx', sub {
  my $v;
  $v = vecfirstidx { 8 == ($_ - 1) } 9,4,5,6; is($v, 0, "first idx success");
  $v = vecfirstidx { 0 } 1,2,3,4; is($v, -1, "first idx failure");
  $v = vecfirstidx { 0 }; is($v, -1, "first idx empty list");
  $v = vecfirstidx { $_->[1] le "e" and "e" le $_->[2] } [qw(a b c)], [qw(d e f)], [qw(g h i)];  is($v, 1, "first idx with reference args");
  $v = vecfirstidx {while(1) {return ($_>6)} } 2,4,6,12; is($v,3,"first idx returns in loop");
};

###### vecuniq
subtest 'vecuniq', sub {
  my @t = (1..10,1..10);
  my @u = vecuniq @t;
  is_deeply(\@u, [1 .. 10], "vecuniq simple 1..10");
  my $u = vecuniq @t;
  is(10,$u,"vecuniq scalar count correct");

  my @n = map { reverse -5 .. 5 } 0..2;
  my @v = vecuniq @n;
  is_deeply(\@v, [reverse -5 .. 5], "vecuniq simple 5 to -5");

  is_deeply([vecuniq()], [], "vecuniq with empty input returns empty");
  is_deeply([vecuniq(0)], [0], "vecuniq with one input returns it");
  is_deeply([vecuniq(0,"18446744073709551615",0,4294967295,"18446744073709551615",4294967295)], [0,"18446744073709551615",4294967295], "vecuniq with 64-bit inputs");
  is_deeply([vecuniq("-9223372036854775808","9223372036854775807",4294967295,"9223372036854775807",4294967295,"-9223372036854775808")], ["-9223372036854775808","9223372036854775807",4294967295], "vecuniq with signed 64-bit inputs");
};

###### vecsingleton
subtest 'vecsingleton', sub {
  my @t = (15,1..10,4,1..10,-2);
  my @u = vecsingleton @t;
  is_deeply(\@u, [15,-2], "vecsingleton simple");
  my $u = vecsingleton @t;
  is(2,$u,"vecsingleton scalar count correct");

  is_deeply([vecsingleton()], [], "vecsingleton with empty input returns empty");
  is_deeply([vecsingleton(0)], [0], "vecsingleton with one input returns it");
  is_deeply([vecsingleton(0,"18446744073709551615",0,4294967295,-1,4294967295)], ["18446744073709551615",-1], "vecsingleton with 64-bit inputs");
  is_deeply([vecsingleton("-9223372036854775808","9223372036854775807",4294967295,"9223372036854775807",4294967295,"-9223372036854775807")], ["-9223372036854775808","-9223372036854775807"], "vecsingleton with signed 64-bit inputs");

  is_deeply([vecsingleton('a','b','',undef,'b','c','')],['a',undef,'c'],"vecsingleton with strings and one undef");
  is_deeply([vecsingleton('a','b','',undef,'b','c',undef)],['a','','c'],"vecsingleton with strings and two undefs");
};

###### vecfreq
subtest 'vecfreq', sub {
  my @L;
  my %got;
  my %exp;
  is_deeply([vecfreq()], [], "vecfreq on empty list");
  is(0+vecfreq(), 0, "vecfreq on empty list (scalar)");

  is_deeply([vecfreq(1)], [1=>1], "vecfreq one integer");
  is(0+vecfreq(1), 1, "vecfreq one integer (scalar)");

  is_deeply([vecfreq(1,1)], [1=>2], "vecfreq two identical integers");
  is(0+vecfreq(1,1), 1, "vecfreq two identical integers (scalar)");

  @L = (-1,1);
  %got = vecfreq(@L);
  %exp = (-1=>1, 1=>1);
  is_deeply(\%got, \%exp, "vecfreq two integers");
  is(0+vecfreq(@L), 2, "vecfreq two integers (scalar)");

  @L = (-1,14,4,-4,2,2,3,4,3,4,4,1);
  %got = vecfreq(@L);
  %exp = (-1=>1, 14=>1, 4=>4, -4=>1, 2=>2, 3=>2, 1=>1);
  is_deeply(\%got, \%exp, "vecfreq many integers");
  is(0+vecfreq(@L), 7, "vecfreq many integers (scalar)");

  @L = ("hello", 14, "world", "tree", "world");
  %got = vecfreq(@L);
  %exp = (hello=>1, 14=>1, world=>2, tree=>1);
  is_deeply(\%got, \%exp, "vecfreq strings");
  is(0+vecfreq(@L), 4, "vecfreq strings (scalar)");

  { # from List::MoreUtils::frequency test
    @L = ('a', 'b', '', undef, 'b', 'c', '', undef);
    my %e = (a=>1, b=>2, ''=>2, c=>1);
    my @f = vecfreq(@L);
    my $seen_undef;
    # This works because we always put undef at the end
    ref $f[-2] and ref $f[-2] eq "SCALAR" and not defined ${$f[-2]} and (undef, $seen_undef) = splice @f, -2, 2, ();
    my %f = @f;
    is_deeply(\%f, \%e, "vecfreq mixed with undef");
    is($seen_undef, 2, "vecfreq counts two undefs");
  }

  is(scalar(vecfreq(-1,~0)),2,"vecfreq doesn't confuse -1 and ~0");
};

###### vecsort
subtest 'vecsort', sub {
  foreach my $r (@vecsorts) {
    my($in, $out, $str) = @$r;
    my @got1 = vecsort(@$in);
    my @got2 = vecsort($in);
    vecsorti($in);
    is_deeply( [ \@got1, \@got2, $in ],
               [ $out, $out, $out ],
               "vecsort list, ref, in-place [$str]" );
  }

  my @s = ("5",2,1,3,4);

  my $in0_beg = length( do { no if $] >= 5.022, "feature", "bitwise"; no warnings "numeric"; $s[0] & "" }) ? "number" : "string";
  my @t = vecsort(\@s);
  my $in0_end = length( do { no if $] >= 5.022, "feature", "bitwise"; no warnings "numeric"; $s[0] & "" }) ? "number" : "string";

  is_deeply([[@s],[@t]], [[5,2,1,3,4],[1,2,3,4,5]], "vecsort sorts without modifying input");

  my @ivd = (qw/-3937 4322 -3619 -390 2039 2123 -1614 -879 -4372 1793 4404 4229 286 -3613 2707 -4166 4025 2450 -2003 3390 4498 -3094 -4854 3441 3501 -2871 -1206 315 71 -2101 4881 -3141 10 -2545 -2825 -519 3534 -4904 -3523 -1170 -3 3 -2 2 -1 1 0/);
  my @sivd = (qw/-4904 -4854 -4372 -4166 -3937 -3619 -3613 -3523 -3141 -3094 -2871 -2825 -2545 -2101 -2003 -1614 -1206 -1170 -879 -519 -390 -3 -2 -1 0 1 2 3 10 71 286 315 1793 2039 2123 2450 2707 3390 3441 3501 3534 4025 4229 4322 4404 4498 4881/);
  is_deeply([vecsort(@ivd)], \@sivd, "vecsort list of negative integers");

  # Both of these should be "string".  XS doesn't copy for validation.
  if ($extra) {
    diag "vecsort input:   $in0_beg => $in0_end" if $in0_beg ne $in0_end;
  }

  my @actx = return_sort(12,13,14,11);
  my $sctx = return_sort(12,13,14,11);
  is($sctx, scalar(@actx), "returning vecsort(\@L) gives the number of items");
};
sub return_sort { return vecsort(@_); }

###### vecslide
subtest 'vecslide', sub {
  is_deeply([vecslide {$a+$b} ()],[],"vecslide with empty array returns empty");
  is_deeply([vecslide {$a+$b} 1..1],[],"vecslide with 1 element returns empty");

  is_deeply([vecslide {$a+$b} 1..5],[3,5,7,9],"vecslide {\$a+\$b} 1..5");
  is_deeply([vecslide { "$a->[0] $b->[1]" } ["hello","world"], ["goodbye","friends"], ["love","hate"]], ["hello friends","goodbye hate"], "vecslide with array refs");
  is(join(", ", vecslide { "$a and $b" } 0..3), "0 and 1, 1 and 2, 2 and 3", "vecslide example from LMU");
};
