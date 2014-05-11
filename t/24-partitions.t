#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/partitions forpart/;
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @parts = qw/
1
1
2
3
5
7
11
15
22
30
42
56
77
101
135
176
231
297
385
490
627
792
1002
1255
1575
1958
2436
3010
3718
4565
5604
6842
8349
10143
12310
14883
17977
21637
26015
31185
37338
44583
53174
63261
75175
89134
105558
124754
147273
173525
204226
/;

my %bparts = (
    101 => "214481126",
    256 => "365749566870782",
    501 => "2431070104309287327876",
   1001 => "25032297938763929621013218349796",
   2347 => "56751384003004060684283391440819878903446789803099",
   4128 => "13036233928924552978434294180871407270098426394166677221003078079504",
   #9988 => "31043825285346179203111322344702502691204288916782299617140664920755263693739998376431336412511604846065386",
  #13337 => "4841449229081281114351180373774137636239639013054790559544724995314398354517477085116206336008004971541987422037760634642695",
  #37373 => "885240148270777711759915557428752066370785294706979437063536090533501018735098279767013023483349639513395622225840616033227700794918506274833787569446519667398089943122156454986205555766363295867812094833219935",
);
if (!$extra) {
  my @ns = grep { $_ > 300 } keys %bparts;
  foreach my $n (@ns) { delete $bparts{$n} }
}

plan tests => scalar(@parts) + scalar(keys(%bparts)) + 14;


foreach my $n (0..$#parts) {
  is( partitions($n), $parts[$n], "partitions($n)" );
}

while (my($n, $epart) = each (%bparts)) {
  is( partitions($n), $epart, "partitions($n)" );
}

################### forpart

{ my @p=(); forpart { push @p, [@_] } 0;
  is_deeply( [@p], [[]], "forpart 0" ); }

{ my @p=(); forpart { push @p, [@_] } 1;
  is_deeply( [@p], [[1]], "forpart 1" ); }

{ my @p=(); forpart { push @p, [@_] } 2;
  is_deeply( [@p], [[2],[1,1]], "forpart 2" ); }

{ my @p=(); forpart { push @p, [@_] } 3;
  is_deeply( [@p], [[3],[2,1],[1,1,1]], "forpart 3" ); }

{ my @p=(); forpart { push @p, [@_] } 3;
  is_deeply( [@p], [[3],[2,1],[1,1,1]], "forpart 3" ); }

{ my @p=(); forpart { push @p, [@_] } 6;
  is_deeply( [@p], [[6],[5,1],[4,2],[4,1,1],[3,3],[3,2,1],[3,1,1,1],[2,2,2],[2,2,1,1],[2,1,1,1,1],[1,1,1,1,1,1]], "forpart 6" ); }

{ my @p=(); forpart { push @p, [@_] } 17,{n=>2};
  is_deeply( [@p], [[16,1],[15,2],[14,3],[13,4],[12,5],[11,6],[10,7],[9,8]], "forpart 17 restricted n=[2,2]" ); }

{ my @p1 = (); my @p2 = ();
  forpart { push @p1, [@_] if @_ <= 5 } 27;
  forpart { push @p2, [@_] } 27, {nmax=>5};
  is_deeply( [@p1], [@p2], "forpart 27 restricted nmax 5" ); }

{ my @p1 = (); my @p2 = ();
  forpart { push @p1, [@_] if @_ >= 20 } 27;
  forpart { push @p2, [@_] } 27, {nmin=>20};
  is_deeply( [@p1], [@p2], "forpart 27 restricted nmin 20" ); }

{ my @p1 = (); my @p2 = ();
  forpart { push @p1, [@_] if @_ >= 10 && @_ <= 13 } 19;
  forpart { push @p2, [@_] } 19, {nmin=>10,nmax=>13};
  is_deeply( [@p1], [@p2], "forpart 19 restricted n=[10..13]" ); }

{ my @p1 = (); my @p2 = ();
  forpart { push @p1, [@_] unless scalar grep { $_ > 4 } @_ } 20;
  forpart { push @p2, [@_] } 20, {amax=>4};
  is_deeply( [@p1], [@p2], "forpart 20 restricted amax 4" ); }

{ my @p1 = (); my @p2 = ();
  forpart { push @p1, [@_] unless scalar grep { $_ < 4 } @_ } 15;
  forpart { push @p2, [@_] } 15, {amin=>4};
  is_deeply( [@p1], [@p2], "forpart 15 restricted amin 4" ); }

{ my @p1 = (); my @p2 = ();
  forpart { push @p1, [@_] unless scalar grep { $_ < 3 || $_ > 6 } @_ } 21;
  forpart { push @p2, [@_] } 21, {amin=>3,amax=>6};
  is_deeply( [@p1], [@p2], "forpart 21 restricted a=[3..6]" ); }

#{ my @p1 = (); my @p2 = ();
#  forpart { push @p1, [@_] unless @_ != 4 || scalar grep { $_ < 2 || $_ > 8 } @_ } 22;
#  forpart { push @p2, [@_] } 22, {amin=>2,amax=>8,n=>4};
#  is_deeply( [@p1], [@p2], "forpart 22 restricted n=4 and a=[3..6]" ); }
{ my @p=(); forpart { push @p, [@_] } 22, {amin=>2,amax=>8,n=>4};
  is_deeply( [@p], [[8,8,4,2],[8,8,3,3],[8,7,5,2],[8,7,4,3],[8,6,6,2], [8,6,5,3], [8,6,4,4], [8,5,5,4], [7,7,6,2], [7,7,5,3], [7,7,4,4], [7,6,6,3], [7,6,5,4], [7,5,5,5], [6,6,6,4], [6,6,5,5]], "forpart 22 restricted n=4 and a=[3..6]" ); }
