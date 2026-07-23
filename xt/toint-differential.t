#!/usr/bin/env perl
use strict;
use warnings;

use Math::BigInt lib => 'Calc';
use Test::More;
use Math::Prime::Util qw/toint/;
use Math::Prime::Util::PP;

sub outcome {
  my ($func, $input) = @_;
  my ($value, $error);
  {
    local $@;
    $value = eval { $func->($input) };
    $error = $@;
  }
  return ['error'] if length($error);
  return ['value', defined($value) ? "$value" : '<undef>'];
}

sub display {
  my ($s) = @_;
  return '<undef>' unless defined $s;
  $s = "$s";
  $s =~ s/\\/\\\\/g;
  $s =~ s/\n/\\n/g;
  $s =~ s/\r/\\r/g;
  $s =~ s/\t/\\t/g;
  $s =~ s/\f/\\f/g;
  $s =~ s/\x0b/\\v/g;
  $s =~ s/([\x00-\x08\x0e-\x1f\x7f])/sprintf("\\x{%02x}", ord($1))/eg;
  return "'$s'";
}

sub exact_integer {
  my ($input) = @_;
  my $s = "$input";
  $s =~ s/\A[ \t\r\n\f\x0b]+//;
  $s =~ s/[ \t\r\n\f\x0b]+\z//;
  $s =~ s/_//g;

  my $negative = ($s =~ s/\A-//);
  $s =~ s/\A\+//;

  my $base = 10;
  if    ($s =~ s/\A0[xX]//) { $base = 16; }
  elsif ($s =~ s/\A0[oO]//) { $base = 8;  }
  elsif ($s =~ s/\A0[bB]//) { $base = 2;  }

  my $n = Math::BigInt->new(0);
  for my $c (split //, lc($s)) {
    my $digit = index('0123456789abcdef', $c);
    $n->bmul($base)->badd($digit);
  }
  $n->bneg if $negative && !$n->is_zero;
  return $n->bstr;
}

sub compare_valid {
  my ($errors, $input, $expected, $label) = @_;
  my $xs = outcome(\&toint, $input);
  my $pp = outcome(\&Math::Prime::Util::PP::toint, $input);
  if ($xs->[0] ne 'value' || $pp->[0] ne 'value' ||
      $xs->[1] ne $expected || $pp->[1] ne $expected) {
    push @$errors, "$label input=" . display($input) .
                   " expected=$expected XS=@$xs PP=@$pp";
  }
}

sub compare_invalid {
  my ($errors, $input) = @_;
  my $xs = outcome(\&toint, $input);
  my $pp = outcome(\&Math::Prime::Util::PP::toint, $input);
  if ($xs->[0] ne 'error' || $pp->[0] ne 'error') {
    push @$errors, "input=" . display($input) . " XS=@$xs PP=@$pp";
  }
}

sub separated {
  my ($s) = @_;
  return join('_', split //, $s);
}

subtest 'generated integer grammar' => sub {
  my @errors;
  my %seen;
  my @space = (
    ['', ''],
    [' ', ' '],
    ["\t", "\r\n"],
  );
  my @spec = (
    [10, [''], [
      qw/0 00 7 00042
         2147483647 2147483648 4294967295 4294967296
         9223372036854775807 9223372036854775808
         18446744073709551615 18446744073709551616/,
      ('1234567890' x 12),
    ]],
    [2, [qw/0b 0B/], [
      qw/0 1 000101/,
      ('1' x 31), ('1' x 32), ('1' x 63), ('1' x 64),
      ('1' . ('0' x 128)), ('1011001' x 30),
    ]],
    [8, [qw/0o 0O/], [
      qw/0 7 000123 17777777777 20000000000/,
      ('7' x 22), ('7' x 44), ('12345670' x 20),
    ]],
    [16, [qw/0x 0X/], [
      qw/0 f 000abc 7fffffff 80000000 ffffffff 100000000
         7fffffffffffffff 8000000000000000 ffffffffffffffff
         10000000000000000/,
      ('abcdef0123456789' x 10),
    ]],
  );

  for my $spec (@spec) {
    my ($base, $prefixes, $bodies) = @$spec;
    for my $prefix (@$prefixes) {
      for my $body (@$bodies) {
        my @body_forms = ($body);
        push @body_forms, separated($body) if length($body) > 1;
        for my $body_form (@body_forms) {
          for my $sign ('', '+', '-') {
            for my $space (@space) {
              my $input = $space->[0] . $sign . $prefix . $body_form .
                          $space->[1];
              next if $seen{$input}++;
              compare_valid(\@errors, $input, exact_integer($input),
                            "base $base");
            }
          }
        }
      }
    }
  }

  is(scalar(@errors), 0, 'XS and PP accept generated integer forms');
  if (@errors) {
    my $last = @errors < 10 ? $#errors : 9;
    diag($_) for @errors[0 .. $last];
  }
};

subtest 'floating-point strings' => sub {
  my @errors;
  my @cases = (
    ['.5', '0'],
    ['-.5', '0'],
    ['+1.', '1'],
    ['1.9', '1'],
    ['-1.9', '-1'],
    ['1e3', '1000'],
    ['1E+3', '1000'],
    ['1e-3', '0'],
    ['-1e-3', '0'],
    ['1_2.3_4e+2', '1234'],
    ['-1_2.3_4e+2', '-1234'],
    ['  1_000_000.9_9  ', '1000000'],
    ['123456789012345678901234567890.99',
     '123456789012345678901234567890'],
    ['-123456789012345678901234567890.99',
     '-123456789012345678901234567890'],
    ['1e30', '1000000000000000000000000000000'],
    ['-1e30', '-1000000000000000000000000000000'],
  );
  compare_valid(\@errors, $_->[0], $_->[1], 'floating point') for @cases;
  is(scalar(@errors), 0, 'XS and PP agree on floating-point strings');
  if (@errors) {
    diag($_) for @errors;
  }
};

subtest 'invalid grammar' => sub {
  my @invalid = (
    ' ', "\t\r\n", '+', '-', '.', '+.', '-.',
    '_', '_1', '1_', '1__0', '1_.0', '1._0',
    '1e', '1e+', '1e-', '1e_2', '1_e2', '1e+_2',
    '1.2.3', '1 2', '--1', '++1', '+-1', '-+1',
    '0x', '0x_', '0x_FF', '0xFF_', '0xF__F', '0x1g', '0x1.0',
    '0b', '0b_', '0b_1', '0b1_', '0b1__0', '0b102',
    '0o', '0o_', '0o_7', '0o7_', '0o7__0', '0o128',
    '0d10', 'NaN', 'Inf', '-Infinity', "12\0" . '34',
  );
  my @errors;
  compare_invalid(\@errors, $_) for @invalid;
  is(scalar(@errors), 0, 'XS and PP reject malformed numeric strings');
  if (@errors) {
    diag($_) for @errors;
  }
};

subtest 'special scalar forms' => sub {
  my @errors;
  compare_valid(\@errors, undef, '0', 'undef');
  compare_valid(\@errors, '', '0', 'empty string');
  compare_valid(\@errors, 0, '0', 'native zero');
  compare_valid(\@errors, 1.75, '1', 'positive NV');
  compare_valid(\@errors, -1.75, '-1', 'negative NV');

  {
    my $source = '59049';
    $source =~ /(\d+)/;
    compare_valid(\@errors, $1, '59049', 'magical capture');
  }

  for my $case (
    ['0xFF',  '255'],
    ['0b101', '5'],
    ['0o777', '511'],
    ['123456789012345678901234567890',
     '123456789012345678901234567890'],
  ) {
    my $input = $case->[0];
    {
      no warnings 'numeric';
      my $unused = 0 + $input;  # Add cached numeric flags to the scalar.
    }
    compare_valid(\@errors, $input, $case->[1], 'previously numified string');
  }

  is(scalar(@errors), 0, 'XS and PP agree on special scalar forms');
  if (@errors) {
    diag($_) for @errors;
  }
};

{
  package Local::IntegerString;
  use overload '""' => sub { $_[0]->[0] }, fallback => 1;
  sub new { bless [$_[1]], $_[0] }
}

subtest 'integer objects and input preservation' => sub {
  my @errors;
  for my $value (qw/0 -42 18446744073709551616
                    -123456789012345678901234567890/) {
    compare_valid(\@errors, Local::IntegerString->new($value), $value,
                  'overloaded integer object');
  }

  my $input = " \t-0xF_F\r\n";
  my $copy = $input;
  compare_valid(\@errors, $input, '-255', 'normalization copy');
  is($input, $copy, 'normalization does not modify its input');
  is(scalar(@errors), 0, 'XS and PP agree on integer objects');
  if (@errors) {
    diag($_) for @errors;
  }
};

done_testing();
