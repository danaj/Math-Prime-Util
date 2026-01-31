#!/usr/bin/perl

use strict;
use warnings;

use ntheory qw/vecequal irand irand64/;
use Math::Prime::Util::PP;
use Array::Compare;
use Benchmark qw( cmpthese );
use List::AllUtils qw( each_arrayref );
use Data::Cmp qw/cmp_data/;
use Algorithm::Diff qw/LCS_length/;
use List::Compare::Functional qw/is_LequivalentR/;
use FreezeThaw qw/cmpStr/;
use Storable qw/freeze/;
use Sereal qw/encode_sereal/;
use match::smart;

#my @x = 1 .. 1_000;
#my @y = map { "$_" } 1 .. 1_000;
my @x = 1 .. 5000, map { irand64 } 1 .. 100;
my @y = 1 .. 5000, map { irand64 } 1 .. 100;

my $comp = Array::Compare->new;

cmpthese -2, {
    iterator => sub { my $r = elementwise_eq(\(@x, @y)) },
    array_comp => sub { my $r = $comp->compare(\(@x, @y)) },
    my_comp => sub { my $r = my_comp(\(@x, @y)) },
    vecequal => sub { my $r = vecequal(\@x, \@y) },
    vecequal_pp => sub { my $r = Math::Prime::Util::PP::vecequal(\@x, \@y) },
    msmart => sub { my $r = \@x |M| \@y; },
    data_cmp => sub { my $r = cmp_data(\@x, \@y) },
    alg_diff => sub { my $r = LCS_length(\@x, \@y) },
    list_compare => sub { my $r = is_LequivalentR([\@x, \@y]) },
    freezethaw => sub { my $r = 0==cmpStr(\@x, \@y); },
    storable => sub { my $r = freeze(\@x) eq freeze(\@y); },
    sereal => sub { my $r = encode_sereal(\@x) eq encode_sereal(\@y); },
};



sub elementwise_eq {
    my ($xref, $yref) = @_;
    return unless  @$xref == @$yref;

    my $it = each_arrayref($xref, $yref);
    while ( my ($x, $y) = $it->() ) {
        return unless $x eq $y;
    }
    return 1;
}

sub my_comp {
    my ($xref, $yref) = @_;
    return unless  @$xref == @$yref;

    my $i;
    for my $e (@$xref) {
        return unless $e eq $yref->[$i++];
    }
    return 1;
}

__END__
                Rate 
msmart         204/s
list_compare   417/s
data_cmp       606/s
freezethaw     718/s
iterator      1058/s
array_comp    1099/s
vecequal_pp   1321/s
alg_diff      1729/s
my_comp       3973/s
storable      4773/s
sereal       12669/s
vecequal     15689/s
