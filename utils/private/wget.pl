#Retrieve contents of file from URI
#
# SYNOPSIS:
#   wget.pl $file $url [$k1 $v1 [$k2 $v2 [ ... $kn $vn]]]
#
# PARAMETERS:
#   $file   - Local filename into which contents of URL will be stored.
#
#   $url    - Uniform resource locator for remote file resource
#
#   $ki/$vi - Key/value pairs that will be supplied in the GET request.
#
# NOTE:
#   This is *NOT* efficient or idiomatic.  We use it as a fall-back for 
#   function WEBSAVE in older releases of MATLAB only.

# Copyright 2009-2016 SINTEF ICT, Applied Mathematics.
#
# This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
#
# MRST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MRST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MRST.  If not, see <http://www.gnu.org/licenses/>.

use LWP::UserAgent;

my $file = shift @ARGV;
my $url  = shift @ARGV;
my %args = @ARGV;

my $args = '';
while (($k, $v) = each %args) {
    $args .= "$k=$v{$k} ";
}

$args =~ s# $##;
$args =~ s# +#&#g;

my $u = $url . '?' . $args;

my $ua  = LWP::UserAgent->new;
my $req = HTTP::Request->new(GET => $u);
my $res = $ua->request($req, $file);

die $res->status_line unless $res->is_success;
