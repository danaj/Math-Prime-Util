use ExtUtils::MakeMaker;
use lib 'inc'; # load our bundled version of Devel::CheckLib
use Devel::CheckLib;

my %require_mpugmp;
my $have_gmp = check_lib(lib => 'gmp', header => 'gmp.h');
if ($have_gmp) {
  warn "\n   It looks like you have the GMP C library.\n";
  warn "   Adding Math::Prime::Util::GMP to dep list.\n\n";
  $require_mpugmp{'Math::Prime::Util::GMP'} = '0.53';
} else {
  warn "\n  It looks like you don't have the GMP library.  Sad face.\n";
}

my $broken64 = (18446744073709550592 == ~0);
if ($broken64) {
  warn <<EOW;

Your Perl, version $], has a broken 64-bit implementation.
Perl's arithmetic operations on numbers larger than 2^53 are wrong.

Much core functionality will work, and we try to work around this in
our functions, but there are still exceptions.

I highly recommend upgrading to a newer version of Perl.

EOW
}

WriteMakefile1(
    NAME         => 'Math::Prime::Util',
    ABSTRACT     => 'Utilities related to prime numbers, including fast sieves and factoring',
    VERSION_FROM => 'lib/Math/Prime/Util.pm',
    LICENSE      => 'perl',
    AUTHOR       => 'Dana A Jacobsen <dana@acm.org>',

    OBJECT       => 'XS.o '           .
                    'util.o '         .
                    'factor128.o '    .
                    'factor.o '       .
                    'rootmod.o '      .
                    'strops.o '       .
                    'primality.o '    .
                    'sieve_cluster.o '.
                    'congruent_numbers.o ' .
                    'totients.o '     .
                    'sieve.o '        .
                    'lucky_numbers.o '.
                    'factmod.o '      .
                    'lmo.o '          .
                    'almost_primes.o '.
                    'prime_counts.o ' .
                    'legendre_phi.o ' .
                    'lucas_seq.o '    .
                    'aks.o '          .
                    'ds_iset.o '      .
                    'moebius.o '      .
                    'semi_primes.o '  .
                    'twin_primes.o '  .
                    'ramanujan_primes.o ' .
                    'powerful.o '     .
                    'powerfree.o '    .
                    'real.o '         .
                    'sort.o '         .
                    'chacha.o '       .
                    'xs_internal.o '  .
                    'xs_set.o '       .
                    'rational.o '     .
                    'prime_sums.o '   .
                    'prime_powers.o ' .
                    'perfect_powers.o ' .
                    'prime_count_cache.o ' .
                    'lehmer.o '       .
                    'omega_primes.o ' .
                    'inverse_interpolate.o ' .
                    'intcomplexity.o '.
                    'goldbach.o '     .
                    'cache.o '        .
                    'csprng.o '       .
                    'random_prime.o ' .
                    'entropy.o',
    LIBS         => ['-lm'],

    EXE_FILES    => ['bin/primes.pl', 'bin/factor.pl'],

    TEST_REQUIRES=> {
                      'Test::More'       => '0.96',
                      # See issue #61 for reasons why we want 0.65 vs 0.22.
                      'bignum'           => '0.65',
                    },
    PREREQ_PM    => {
                      'Exporter'         => '5.57',
                      'XSLoader'         => '0.01',
                      'Carp'             => ($] < 5.008) ? '1.17' : 0,
                      'Tie::Array'       => 0,
                      'base'             => 0,
                      'constant'         => 0,
                      'Config'           => 0,
                      # Has from_base
                      'Math::BigInt'     => '1.999814',
                      'Math::BigFloat'   => '1.59',
                      # Add in MPU::GMP if we can
                      %require_mpugmp,
                    },
    META_MERGE   => {
               'meta-spec'  => {
                  version => '2',
                  url     => 'http://search.cpan.org/perldoc?CPAN::Meta::Spec',
               },
               release_status => 'testing',
               dynamic_config => 1,  # Check for GMP on install
               resources    => {
                  license    => [ 'http://dev.perl.org/licenses/' ],
                  homepage   => 'https://github.com/danaj/Math-Prime-Util',
                  repository => {
                    url  => 'https://github.com/danaj/Math-Prime-Util',
                  },
               },
               provides => {
                 'ntheory' => {
                    version => '0.75', file => 'lib/ntheory.pm',
                 },
                 'Math::Prime::Util' => {
                    version => '0.75', file => 'lib/Math/Prime/Util.pm',
                 },
                 'Math::Prime::Util::MemFree' => {
                    version => '0.75', file => 'lib/Math/Prime/Util/MemFree.pm',
                 },
                 'Math::Prime::Util::PP' => {
                    version => '0.75', file => 'lib/Math/Prime/Util/PP.pm',
                 },
                 'Math::Prime::Util::PrimeArray' => {
                    version => '0.75', file => 'lib/Math/Prime/Util/PrimeArray.pm',
                 },
                 'Math::Prime::Util::PrimeIterator' => {
                    version => '0.75', file => 'lib/Math/Prime/Util/PrimeIterator.pm',
                 },
                 'Math::Prime::Util::Entropy' => {
                    version => '0.75', file => 'lib/Math/Prime/Util/Entropy.pm',
                 },
                 'Math::Prime::Util::ChaCha' => {
                    version => '0.75', file => 'lib/Math/Prime/Util/ChaCha.pm',
                 },
                 # Skip: PPFE, PrimalityProving, RandomPrimes, ZetaBigFloat,
                 #       ECAffinePoint, ECProjectivePoint
               },
               prereqs => {
                 runtime => {
                   recommends => {
                     'Math::Prime::Util::GMP' => 0.53,
                     'Math::BigInt::GMP'      => 0,
                     'Math::GMPz'             => 0.68,
                     'Digest::SHA'            => 5.87,
                   },
                 },
                 test  => {
                   suggests => {
                     'Test::Warn' => 0,
                   },
                 },
               },
             },

    MIN_PERL_VERSION => 5.006002,
);

sub WriteMakefile1 {   # Cribbed from eumm-upgrade by Alexandr Ciornii
  my %params = @_;
  my $eumm_version = $ExtUtils::MakeMaker::VERSION;
  $eumm_version = eval $eumm_version;

  if ($params{TEST_REQUIRES} and $eumm_version < 6.6303) {
      $params{BUILD_REQUIRES}={ %{$params{BUILD_REQUIRES} || {}} , %{$params{TEST_REQUIRES}} };
      delete $params{TEST_REQUIRES};
  }
  if ($params{BUILD_REQUIRES} and $eumm_version < 6.5503) {
      #EUMM 6.5502 has problems with BUILD_REQUIRES
      $params{PREREQ_PM}={ %{$params{PREREQ_PM} || {}} , %{$params{BUILD_REQUIRES}} };
      delete $params{BUILD_REQUIRES};
  }
  delete $params{CONFIGURE_REQUIRES} if $eumm_version < 6.52;
  delete $params{MIN_PERL_VERSION} if $eumm_version < 6.48;
  delete $params{META_MERGE} if $eumm_version < 6.46;
  delete $params{META_ADD} if $eumm_version < 6.46;
  delete $params{LICENSE} if $eumm_version < 6.31;
  delete $params{AUTHOR} if $] < 5.005;
  delete $params{ABSTRACT_FROM} if $] < 5.005;
  delete $params{BINARY_LOCATION} if $] < 5.005;

  WriteMakefile(%params);
}
