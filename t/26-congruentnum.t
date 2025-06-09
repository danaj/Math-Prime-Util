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

  ["pq 33",       [33,57,129]],
  ["pq 13",       [51,123,187]],
  ["pq 57",       [35,91,115]],

  ["2pq 55",      [130,290,370]],
  ["2pq 33",      [66,114,258]],
  ["2pq 15",      [170,730,970]],
  ["2pq 37",      [42,186,266]],
  ["2pq 77",      [1442,2338,2786,4354,4738,6146,7042]],
  ["2pq 11",      [1394,2482,7922,7954]],

  ["pqr 133",     [969,2193,2409]],
  ["pqr 357",     [105,273,345]],
  ["pqr 377",     [483,1491,2387,3619]],
  ["pqr 131",     [2091,4539,5763]],
  ["pqr 131 (2)", [3723,4947,9843,11931,12291,13243]],
  ["pqr 157",     [595,1435,1547,1955]],
  ["pqr 355",     [195,435,555,715,795]],
  ["pqr 333",     [1947,2211,2451,2739,3363]],
  ["pqr 111",     [78761,95489,120377]],

  ["2pqr 133",    [1938,4386,4818]],
  ["2pqr 155",    [2210,5330,9010]],
  ["2pqr 357",    [690,770,930,1218]],
  ["2pqr 111",    [135218,157522,240754,269042]],
  ["2pqr 577",    [4186,4970,5530,7130]],
  ["2pqr 115",    [6970,12410,16490,18122,19210]],
  ["2pqr 335",    [570,858,1290,1914,2010,2090]],
  ["2pqr 555",    [3770,7930,13130,15370,17690]],

  ["pqrs 5577",   [10465,14105,29785,31465,32305,40145]],
  ["pqrs 3357",   [8715,23835,29667,29715,32235,49035,54915]],
  ["pqrs 1333",   [60027,73491,90651,91443,127347,156579,194667]],
  ["pqrs 3333",   [161601,255057,293073,317361,340689,348513]],
  ["pqrs 1357",   [30849,31161,51865,63609,89585,100113,104673]],
  ["pqrs 1133",   [214401,316217,330033,386097,396321,419577]],

  ["2pqrs 1133",  [79458,81906,126786,141474,179826,187986]],
  ["2pqrs 3557",  [36330,61530,75530,142170,162330,163770,167370]],
  ["2pqrs 1335",  [32538,54570,78474,83130,115770,135546,144330]],
  ["2pqrs 3337",  [47082,110922,169818,181146,205674,225834]],
  ["2pqrs 1555",  [466570,928330,1697930,1742330,1869530,2198170]],
  ["2pqrs 3333",  [323202,510114,586146,634722,681378,697026]],
  ["2pqrs 5555",  [2679170,4102930,7181330,7681570,8118370,10167170]],
  ["2pqrs 3355",  [106530,131970,193314,395970,449970,555330]],

);

my @ncfilters2 = (
  # Most of the 1/2/3/4-factor results are in the first filter set.
  [1,"Iskra 1996",[qw/161601 255057 293073 317361 323202 340689 348513/]],
  #[1,"Reinholz 2013",[qw/26961 36993 42009 52041 67089 82137 83721 87153/]],
  [1,"Reinholz 2013",[qw/203433 615201 717369 924801 1085793 1988673/]],
  [1,"Cheng 2019",[qw/761442 981618 1155858 1201794 1490082 1550274/]],
  [1,"Cheng 2018 T1.1",[qw/2679170 4102930 7181330 7681570 8118370/]],
  [1,"Cheng 2018 T1.2",[qw/15073995 24459690 24649707 26534010 28219395/]],
  [1,"Cheng 2018 T1.3",[qw/4451145 5322345 5488305 5676555 9122505 12063594/]],
  [1,"Cheng 2018 T1.4",[qw/505155 1147755 1404795 1529745 2129505 2172345/]],
  [1,"Das 2020",[qw/432502235 1085236971 1332878635 1445059707 1641579115/]],

  # These only work if the factor ordering is allowed to permute.
  # The 0 indicates we're going to skip them.
  [0,"Reinholz 2013",[qw/26961 36993 42009 52041 82137 83721 87153 112233/]],
  [0,"Cheng 2019 (ord)",[qw/53922 73986 84018 104082 164274 167442 174306/]],
  [0,"Cheng 2018 T1.1",[qw/71885379 83674803 106112739 112768059 159654099/]],
  [0,"Cheng 2018 T1.2 (ord)",[qw/5221227 7981347 11477499 14789643 15709683/]],
  [0,"Cheng 2018 T1.3 (ord)",[qw/1334058 1978305 4069905 4072593 4190298/]],
  [0,"Cheng 2018 T1.4 (ord)",[qw/405195 504339 2239755 3403995/]],
  [0,"Das 2020 (ord)",[qw/332868235 682507315 881126547 932968715 1180617347/]],
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
