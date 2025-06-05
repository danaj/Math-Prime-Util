#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util qw/is_congruent_number/;
my $usexs = Math::Prime::Util::prime_get_config->{'xs'};
#my $usegmp = Math::Prime::Util::prime_get_config->{'gmp'};
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @cn200 = (5,6,7,13,14,15,20,21,22,23,24,28,29,30,31,34,37,38,39,41,45,46,47,52,53,54,55,56,60,61,62,63,65,69,70,71,77,78,79,80,84,85,86,87,88,92,93,94,95,96,101,102,103,109,110,111,112,116,117,118,119,120,124,125,126,127,133,134,135,136,137,138,141,142,143,145,148,149,150,151,152,154,156,157,158,159,161,164,165,166,167,173,174,175,180,181,182,183,184,188,189,190,191,194,197,198,199);
my @cn1e6 = (1,5,6,7,9);

my @ncfilters = (
  ["p3",          [19,43,59]],
  ["(1+sqrt2,p)", [17,73,89]],
  ["2p Bastien",  [26,58,74,82]],
  ["p3q3",        [33,57,129]],
  ["p1q3",        [51,123,187]],
  ["p5q7",        [35,91,115]],
  ["2p5q5",       [130,290,370]],
  ["2p3q3",       [66,114,258]],
  ["2p1q5",       [170,730,970]],
  ["2p3q7",       [42,186,266]],
  ["2p7q7",       [1442,2338,2786,4354,4738,6146,7042]],
  ["2p1q1",       [1394,2482,7922,7954]],
  ["p1q3r3",      [969,2193,2409]],
  ["p3q5r7",      [105,273,345]],
  ["p3q7r7",      [483,1491,2387,3619]],
  ["p1q3r1",      [2091,4539,5763]],
  ["p1q3r1 (2)",  [3723,4947,9843,11931,12291,13243]],
  ["p1p5p7",      [595,1435,1547,1955]],
  ["p3p5p5",      [195,435,555,715,795]],
  ["p3p3p3",      [1947,2211,2451,2739,3363]],
  ["p1p1p1",      [78761,95489,120377]],
  ["2p1q3r3",     [1938,4386,4818]],
  ["2p1q5r5",     [2210,5330,9010]],
  ["2p3q5r7",     [690,770,930,1218]],
  ["2p1q1r1",     [135218,157522,240754,269042]],
  ["2p5q7r7",     [4186,4970,5530,7130]],
  ["2p1q1r5",     [6970,12410,16490,18122,19210]],
  ["2p3q3r5",     [570,858,1290,1914,2010,2090]],
  ["2p5q5r5",     [3770,7930,13130,15370,17690]],
  ["p5q5r7s7",    [10465,14105,29785,31465,32305,40145]],
);

my @ncfilters2 = (
  # Once we added the extra Lagrange 1974 filters, many fewer are found here
  [1,"Iskra 1996",[qw/340689 681378 693633 1291449 1387266 1517169 1575993/]],
  [1,"Reinholz 2013",[qw/203433 360393 615201 717369 924801 1085793 1368609/]],
  [1,"Cheng 2019",[qw/161601 255057 323202 358809 380721 441969 490809/]],
  [1,"Cheng 2018 T1.1",[qw/2679170 4102930 7181330 7681570 8118370/]],
  [1,"Cheng 2018 T1.2",[qw/106530 110922 218139 225834 276963 395970 449970/]],
  [1,"Cheng 2018 T1.3",[qw/714 1722 3162 3738 8058 8602 10794 10906 11594/]],
  [1,"Cheng 2018 T1.4",[qw/505155 1147755 1404795 1529745 2129505 2172345/]],
  [1,"Das 2020",[qw/30849 51865 63609 80185 89585 100113 103385 110929/]],

  # These only work if the factor ordering is allowed to permute.
  # The 0 indicates we're going to skip them.
  [0,"Iskra 1996",[qw/161601 255057 293073 317361 323202 359841 360393/]],
  [0,"Reinholz 2013",[qw/26961 36993 83721 87153 107844 112233 132297/]],
  [0,"Cheng 2019",[qw/36993 42009 52041 73986 82137 83721 84018 104082/]],
  #[0,"Cheng 2018 T1.1",[qw//]],
  [0,"Cheng 2018 T1.2",[qw/47082 60027 73491 90651 91443 127347 131970/]],
  [0,"Cheng 2018 T1.3",[qw/8715 23835 29667 32538 34860 49035 54570 59115/]],
  [0,"Cheng 2018 T1.4",[qw/405195 504339 1620780 2017356 2239755 3403995/]],
  [0,"Das 2020",[qw/31465 125860 214401 263465 275065 283185 315665 316217/]],
);


plan tests => 0
            + 3
            + scalar(@ncfilters)
            + scalar(@ncfilters2)
            ;

# This covers all simple special cases in the code
is_deeply([grep { is_congruent_number($_) } 1..200], \@cn200, "congruent numbers to 200");

is_deeply([grep { is_congruent_number(1000000+$_) } 1..10], \@cn1e6, "congruent numbers 10^6 + (1..10)");

# Skip the filters and directly test the Tunnell loop
SKIP: {
  skip "PP doesn't have test interfaces", 1 unless $usexs;
  is_deeply([grep { Math::Prime::Util::_is_congruent_number_tunnell($_) } 1..200], \@cn200, "congruent numbers to 200 (no filtering)");
}


# Test the various filters finding non-congruent families
SKIP: {
  skip "PP doesn't have all the NC families",scalar(@ncfilters) unless $usexs;
  for my $td (@ncfilters) {
    my($name,$narr) = @$td;
    my @got = map {Math::Prime::Util::_is_congruent_number_filter($_)} @$narr;
    my @exp = map {0} @$narr;
    # Any return value of -1 means the filter didn't work.
    is_deeply(\@got, \@exp, "Non-congruent family of $name");
  }
}
SKIP: {
  skip "PP doesn't have all the NC families",scalar(@ncfilters2) unless $usexs;
  for my $td (@ncfilters2) {
    my($order,$name,$narr) = @$td;
    SKIP: {
      skip "NC filters don't use arbitrary factor order", 1 unless $order;
      my @got = map {Math::Prime::Util::_is_congruent_number_filter($_)} @$narr;
      my @exp = map {0} @$narr;
      # Any return value of -1 means the filter didn't work.
      is_deeply(\@got, \@exp, "Non-congruent family of $name");
    }
  }
}

  #my @das = (17*3*409*19,17*3*859*3697,19*409*3697*859,17*3*409*19*3697*859, 5*7*29*79,5*7*821*151,29*79*821*151,5*7*29*79*821*151);
  #@das = grep { $_ <= ~0 } @das;  # Only test native integers
  #is_deeply([map { is_congruent_number($_) } @das],
  #          [map { 0 } @das],
  #          "Non-congruent examples from Das and Saikia 2020 section 5");
  # 30849, 52865, 63609, .... 4 factors
  # 432502235                 6 factors
  # no 8/10/12 factor results below 5e11
  # 311199575628433           8 factors  (probably not the smallest)
