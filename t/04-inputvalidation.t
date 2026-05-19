#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/toint next_prime/;
use Math::BigInt try=>"GMP,Pari";
use Math::BigFloat;
use Config;
use Carp;

my @incorrect = (
  -4,
  '-',
  '+',
  '++4',
  '+-4',
  '0-0',
  '-0004',
  'a',
  '5.6',
  '4e',
  '1.1e12',
  '1e8',
  'NaN',
  Math::BigInt->new("-4"),
  Math::BigFloat->new("15.6"),
);
push @incorrect, Math::BigInt->bnan() if $Config{d_isnan};

my %correct = (
  4       => 5,
  '+4'    => 5,
  '0004'  => 5,
  '+0004' => 5,
  5.0     => 7,
  1e8     => 100000007,
  Math::BigInt->new("10000000000000000000000012") => Math::BigInt->new("10000000000000000000000013"),
  Math::BigFloat->new("9") => 11,
);

plan tests => 2                      # undefined and empty string
            + scalar(@incorrect)     # values that should be rejected
            + scalar(keys(%correct)) # values that should be accepted
            + 2                      # infinity and nan
            + 2                      # internal validators normalize scalar ints
            + 1                      # internal validators reject coderefs
            + 7                      # internal canonicalization helper
            + 1;                     # long invalid string

my $qrnn = qr/ must be a (non-negative|positive) integer/;

eval { next_prime(undef); };
like($@, qr/^Parameter must be defined/, "Gives Error:  next_prime(undef)");

eval { next_prime(""); };
like($@, $qrnn, "Gives Error:  next_prime('')");

foreach my $v (@incorrect) {
  $v = "$v" if $] < 5.008 && ref($v) eq 'Math::BigFloat';
  eval { next_prime($v); };
  like($@, $qrnn, "Gives Error:  next_prime($v)");
}

while (my($v, $expect) = each (%correct)) {
  is("".next_prime($v), $expect, "Correct:      next_prime($v)");
}

# The next two tests really are not critical, but are nice to check.
SKIP: {
  skip "Your machine does not have infinity", 1 unless $Config{d_isinf};
  my $infinity = ($^O ne 'MSWin32') ? 0+'inf' : '1.#INF';
  $infinity = Math::BigInt->binf()->numify() if 65535 > $infinity;
  $infinity = +(20**20**20) if 65535 > $infinity;
  skip "Your machine seems to not have infinity", 1 if 65535 > $infinity;
  eval { next_prime($infinity); };
  like($@, $qrnn, "Gives Error:  next_prime( infinity )");
}

SKIP: {
  skip "Your machine does not have NaN", 1 unless $Config{d_isnan};
  no warnings 'numeric';
  my $nan = ($^O ne 'MSWin32') ? 0+'nan' : '1.#IND';
  $nan      = Math::BigInt->bnan()->numify() if $nan >= 0;
  $nan      = -sin('inf') if $nan >= 0;
  skip "Your machine seems to not have NaN", 1 if $nan >= 0 || $nan =~ /^\d*$/;
  eval { next_prime($nan); };
  like($@, $qrnn, "Gives Error:  next_prime( nan ) [nan = '$nan']");
}

{
  my @vals = ("42", "+42", "007", "-0");
  my @got;
  for my $v (@vals) {
    my $x = $v;
    Math::Prime::Util::_validate_integer($x);
    push @got, [ "$x", ref($x) ? ref($x) : "" ];
  }
  is_deeply(\@got,
            [ ["42",""], ["42",""], ["7",""], ["0",""] ],
            "_validate_integer normalizes in-range scalar integers");
}

{
  my @vals = ("42", "+42", "007");
  my @got;
  for my $v (@vals) {
    my $x = $v;
    Math::Prime::Util::_validate_integer_nonneg($x);
    push @got, [ "$x", ref($x) ? ref($x) : "" ];
  }
  is_deeply(\@got,
            [ ["42",""], ["42",""], ["7",""] ],
            "_validate_integer_nonneg normalizes in-range scalar integers");
}

{
  my $ok = 0;
  $ok++ if !defined eval { Math::Prime::Util::_validate_integer(sub { 123 }); 1 } &&
            $@ =~ /must be an integer/;
  $ok++ if !defined eval { Math::Prime::Util::_validate_integer_nonneg(sub { 123 }); 1 } &&
            $@ =~ /must be a non-negative integer/;
  is($ok, 2, "_validate_integer* reject coderefs");
}

################################################################################

# How to find the default bigint class?
# 1. assume Math::BigInt.  That is Wrong.
# 2. ref(powint(2,80)).  Works but many things being tested here.
# 3. ref(toint("1208925819614629174706176")).  Works, almost there.
# 4. _load_bigint().  Very low level private function.  Works.
#
# In later tests we would want to use (3).  Here we'll use (4).
#
my $bigint_class = Math::Prime::Util::_load_bigint();
#diag "Default bigint type: $bigint_class";

{
  my $x = Math::BigInt->new(123);
  Math::Prime::Util::_canonicalize_integers(\$x);
  is(ref($x), "", "_canonicalize_integers scalar ref down-converts small bigint");
}

{
  my $x = Math::BigInt->new("123456789012345678901234567890");
  Math::Prime::Util::_canonicalize_integers(\$x);
  is("".ref($x).":$x", "$bigint_class:123456789012345678901234567890",
     "_canonicalize_integers scalar ref keeps large bigint");
}

{
  my $x = Math::BigInt->new(123);
  my $r = Math::Prime::Util::_canonicalized_integer($x);
  is_deeply([ref($r), "$r", ref($x), "$x"],
            ["", "123", "Math::BigInt", "123"],
            "_canonicalized_integer returns native without mutating source");
}

{
  my $x = Math::BigInt->new("123456789012345678901234567890");
  my $r = Math::Prime::Util::_canonicalized_integer($x);
  is("".ref($r).":$r", "$bigint_class:123456789012345678901234567890",
     "_canonicalized_integer returns canonical large bigint");
}

{
  my @v = ("42", Math::BigInt->new(7),
           Math::BigInt->new("123456789012345678901234567890"),
           undef);
  Math::Prime::Util::_canonicalize_integers(\@v);
  is_deeply([map { defined($_) ? ((ref($_)||"").":$_") : "undef" } @v],
            [":42", ":7", "$bigint_class:123456789012345678901234567890", "undef"],
            "_canonicalize_integers array ref canonicalizes elements");
}

{
  my @v;
  $v[2] = Math::BigInt->new(5);
  Math::Prime::Util::_canonicalize_integers(\@v);
  ok(!exists($v[0]) && !exists($v[1]) && exists($v[2]) && !ref($v[2]) && $v[2] == 5,
     "_canonicalize_integers array ref preserves holes");
}

{
  my $ok = 0;
  $ok++ if !eval { Math::Prime::Util::_canonicalize_integers(5); 1 };
  my $x = "abc";
  $ok++ if !eval { Math::Prime::Util::_canonicalize_integers(\$x); 1 };
  is($ok, 2, "_canonicalize_integers rejects bad inputs");
}


SKIP: {
  skip "Perl $], Carp $Carp::VERSION.  We need a minimum of 5.8 or Carp 1.17 to avoid segfaults.", 1 if $] < 5.008 && $Carp::VERSION < 1.17;
  eval { next_prime("11111111111111111111111111111111111111111x"); };
  like($@, $qrnn, "Gives Error:  next_prime('111...111x')");
}
