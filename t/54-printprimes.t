#!/usr/bin/env perl
use strict;
use warnings;

use File::Temp qw/tempfile/;
use IO::Handle ();
use POSIX ();
use Test::More;
use Math::Prime::Util qw/print_primes/;

sub capture_fds (&) {
  my ($code) = @_;
  my ($outfh) = tempfile();
  my ($errfh) = tempfile();

  STDOUT->flush;
  STDERR->flush;
  my $oldout = POSIX::dup(fileno(STDOUT));
  my $olderr = POSIX::dup(fileno(STDERR));
  die "Cannot dup STDOUT: $!" unless defined $oldout;
  die "Cannot dup STDERR: $!" unless defined $olderr;

  my $ok = eval {
    POSIX::dup2(fileno($outfh), fileno(STDOUT)) or die "Cannot redirect STDOUT: $!";
    POSIX::dup2(fileno($errfh), fileno(STDERR)) or die "Cannot redirect STDERR: $!";
    STDOUT->autoflush(1);
    STDERR->autoflush(1);
    $code->();
    STDOUT->flush;
    STDERR->flush;
    1;
  };
  my $err = $@;

  POSIX::dup2($oldout, fileno(STDOUT)) or die "Cannot restore STDOUT: $!";
  POSIX::dup2($olderr, fileno(STDERR)) or die "Cannot restore STDERR: $!";
  POSIX::close($oldout);
  POSIX::close($olderr);

  die $err unless $ok;

  seek($outfh, 0, 0) or die "Cannot seek STDOUT capture: $!";
  seek($errfh, 0, 0) or die "Cannot seek STDERR capture: $!";

  local $/;
  my $stdout = <$outfh>;
  my $stderr = <$errfh>;
  return ($stdout || '', $stderr || '');
}

{
  my ($stdout, $stderr) = capture_fds { print_primes(11); };
  is($stdout, "2\n3\n5\n7\n11\n", "print_primes(hi) writes primes to STDOUT");
  is($stderr, "", "print_primes(hi) leaves STDERR alone");
}

{
  my ($stdout, $stderr) = capture_fds { print_primes(10, 19); };
  is($stdout, "11\n13\n17\n19\n", "print_primes(lo, hi) writes primes to STDOUT");
  is($stderr, "", "print_primes(lo, hi) leaves STDERR alone");
}

{
  my ($stdout, $stderr) = capture_fds { print_primes(20, 10); };
  is($stdout, "", "print_primes empty range writes nothing to STDOUT");
  is($stderr, "", "print_primes empty range writes nothing to STDERR");
}

{
  my ($stdout, $stderr) = capture_fds {
    print_primes("99999999999999999999999999900",
                 "100000000000000000000000000378");
  };
  is($stdout, join("", map { "$_\n" }
                   "99999999999999999999999999947",
                   "99999999999999999999999999973",
                   "100000000000000000000000000319"),
     "print_primes bigint range writes primes to STDOUT");
  is($stderr, "", "print_primes bigint range leaves STDERR alone");
}

{
  my ($stdout, $stderr) = capture_fds {
    print_primes(5, 11, fileno(STDOUT));
  };
  is($stdout, "5\n7\n11\n", "print_primes(..., fileno(STDOUT)) writes to STDOUT");
  is($stderr, "", "print_primes(..., fileno(STDOUT)) leaves STDERR alone");
}

{
  my ($stdout, $stderr) = capture_fds {
    print_primes(5, 11, fileno(STDERR));
  };
  is($stdout, "", "print_primes(..., fileno(STDERR)) leaves STDOUT alone");
  is($stderr, "5\n7\n11\n", "print_primes(..., fileno(STDERR)) writes to STDERR");
}

{
  for my $case (
    [ undef, "explicit undef fd", qr/defined/ ],
    [    -1, "negative fd",       qr/non-negative integer/ ],
    [   1.5, "non-integer fd",    qr/non-negative integer/ ],
    [ "2foo", "invalid fd",       qr/non-negative integer/ ],
  ) {
    my($fd, $name, $error) = @$case;
    my $ok = eval { print_primes(2, 10, $fd); 1 };
    ok(!$ok, "print_primes rejects $name");
    like($@, $error, "print_primes $name error");
  }
}

{
  my $fd = POSIX::dup(fileno(STDOUT));
  die "Cannot dup STDOUT for closed fd test: $!" unless defined $fd;
  POSIX::close($fd) or die "Cannot close duplicated fd: $!";

  my $ok = eval { print_primes(2, 10, $fd); 1 };
  ok(!$ok, "print_primes rejects closed fd");
  like($@, qr/print_primes: open fd \Q$fd\E failed:/, "print_primes closed fd error");
}

{
  ok(!eval { &print_primes(2, 10, 1, 0); 1 }, "print_primes rejects too many arguments");
  like($@, qr/expected|Usage/, "print_primes too many arguments error");
}

done_testing();
