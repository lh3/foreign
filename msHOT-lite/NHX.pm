package NHX;

=head1 NAME

LHPM::NHX - Parser for NHX/TFF format

=head1 SYNOPSIS

  use LHPM::NHX qw(get_leaves);

  # read NHX/TFF string from STDIN and parse the first tree.
  my $nhx = LHPM::NHX->new;
  $nhx->parse(join('', <STDIN>));

  # visit each node and do some modifications
  foreach my $p ($nhx->node_array) {
      print "length: $p->{_L}; " if ($p->{_L});
      if ($p->{_C}) {
          print "internal node with ", scalar(@{$p->{_C}}), " children ";
      } else {
          print "external node '$p->{N}' ";
          $p->{N} .= "-modified";
      }
      print "attributes: ";
      foreach my $key (keys %$p) {
          print ":$key=$p->{$key}";
      }
      $p->{S} = 'unknown' unless ($defined{$p->{S}});
      print "\n";
  }

  # export the tree as an NHX string
  my $nhx_str = $nhx->string;

  # get a list of external nodes (or leaves)
  my @leaves;
  get_leaves($nhx_str, \@leaves);

=head1 DESCRIPTION

This module parses trees in NHX format and stores the result in a tree-like
recursive structure. Various nodes attributes, including topological
connections, are stored in a hash which can be visited by L<node_array> method
or by recursively traversing from the B<root> node. Predefined keys and
attributes are:

    N        node name (display name)
    _P       reference to parent node
    _C       list of references to children (empty if the node is a leaf)
    _L       branch length
    _index   DFS finish time of the node
    _lindex  DFS finish time of the left-most leaf
    B        bootstrap value
    S        taxon name
    D        'Y' for duplication, 'N' for speciation, or not defined
    O        sequence ID
    G        gene ID
    E        gene loss stored in a string like "$-Eutheria-DROME"
    Com      'Y' if the branch exists in compared tree
    Col      color of the descended subtree

=head2 Methods

=cut

use strict;
use Exporter;

use vars qw(@ISA @EXPORT_OK @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw/tff2nhx nhx2tff/;
@EXPORT_OK = qw/get_leaves/;

=head3 new

  Arg [0]     : NONE
  ReturnType  : LHPM::NHX
  Example     : $nhx = LHPM::NHX->new;
  Description : Get a LHPM::NHX object.

=cut

sub new
{
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = { @_ };
	bless($self, $class);
	return $self;
}

=head3 parse

  Arg [1]     : string $nhx_str
  ReturnType  : NONE
  Example     : $nhx->parse($nhx_str);
  Description : Parse the first the tree stored in the string $str.

=cut

sub parse
{
	my ($self, $str) = @_;
	return $self->parse_nhx($str);
}

=head3 parse_nhx

  Arg [1]     : string $str
  ReturnType  : NONE
  Example     : $nhx->parse_nhx($nhx_str);
  Description : Parse the first NHX-formatted tree stored in the string $str.

=cut

sub parse_nhx {
  my ($self, $str) = @_;
  my ($array, @stack, $pos);
  $self->DESTROY if ($self->{_root});
  $self->{_error} = 0;
  @{$self->{_node}} = ();
  $array = \@{$self->{_node}};
  $_ = (($pos = index($str, ';')) >= 0)? substr($str, 0, $pos) : $str;
  s/\s//g;
  while (/(\(|((\)?[^,;:\[\]\(\)]+|\))(:[\d.eE\-]+)?(\[[^\[\]]*\])?))/g) {
	&parse_aux($self, $array, \@stack, $1, $3, $4, $5);
  }
  if (@stack != 1) {
	my $count = @stack;
	warn(qq/[LHPM::NHX::parse] unmatched left bracket ($count)/);
	$self->{_error} = 1;
	@stack = ();
  }
  if ($self->{_error} == 0) {
	$self->{_root} = shift(@stack);
  } else {
	@{$self->{_node}} = ();
	delete($self->{_root});
	return;
  }
  $self->set_index;
  return $self->{_root};
}

sub parse_aux {
  my ($self, $array, $stack, $str, $name, $dist, $nhx) = @_;
  if ($str eq '(') {
	push(@$stack, $str);
  } elsif ($name) {
	my %hash;
	if ($name =~ /^\)/) {
	  my (@s, $t);
	  while (($t = pop(@$stack))) {
		last if (ref($t) ne 'HASH');
		push(@s, $t);
	  }
	  unless (defined($t)) {
		warn('[LHPM::NHX::parse_aux] unmatched right bracket');
		$self->{_error} = 1;
		return;
	  }
	  foreach (reverse @s) {
		push(@{$hash{_C}}, $_);
		$_->{_P} = \%hash;
	  }
	  $hash{N} = substr($name, 1) if (length($name) > 1);
	} else {
	  $hash{N} = $name;
	}
	$hash{_L} = substr($dist, 1) if ($dist);
	$nhx =~ s/:([^=:]+)=([^:=\[\]]+)/$hash{$1}=$2,''/eg if ($nhx && $nhx =~ /^\[&&NHX/);
	push(@$stack, \%hash);
	push(@$array, \%hash);
  }
}

=head3 parse_tff

  Arg [1]     : string $str
  ReturnType  : NONE
  Example     : $nhx->parse_nhx($nhx_str);
  Description : Parse the first TFF-formatted tree stored in the string $str.

=cut

sub parse_tff {
  my ($self, $str) = @_;
  my (@data, %conv, @node, $root_index);
  my $i = 0;
  $self->DESTROY if ($self->{_root});
  @{$data[$i++]} = split("\t") foreach (split("\n", $str));
  $i = 0;
  foreach my $p (sort {$a->[0]<=>$b->[0]} @data) { # initialization
	my %hash;
	push(@node, \%hash);
	$conv{$p->[0]} = $i++;
	$root_index = $p->[0] if (!$root_index && $p->[0] == $p->[1]);
  }
  if (!$root_index) {
	warn('[LHPM::NHX::parse_tff] no root');
	$self->{_error} = 1;
	return;
  }
  $root_index = $conv{$root_index};
  foreach my $p (@data) {
	my ($chi, $par) = ($conv{$p->[0]}, $conv{$p->[1]});
	my $q = \%{$node[$chi]};
	$q->{_P} = $node[$par];
	push(@{$node[$par]->{_C}}, $q) if ($chi != $par);
	$q->{_L} = $p->[2];
	foreach (my $i = 3; $i < @$p; ++$i) {
	  $q->{$1} = $2 if ($p->[$i] =~ /^\s*([A-Za-z][A-Za-z0-9]*)\s*=\s*([^\n\t]+)\s*$/);
	}
  }
  $self->{_root} = $node[$root_index];
  $self->update_node_array;
  $i = 0;
  for my $p ($self->node_array) {
	if (!$p->{N} && !$p->{_C}) {
	  $p->{N} = $i;
	  warn("[LHPM::NHX::parse_tff] missing name for external leaf!\n");
	}
	++$i;
  }
  return $self->{_root};
}

=head3 node_array

  Arg [0]     : NONE
  ReturnType  : array
  Example     : @array = @{$nhx->node_array};
  Description : Return the list of nodes in suffix order.

=cut

sub node_array
{
	my ($self) = @_;
	if (defined($_[0])) {
		return $self->{_node}->[$_[0]];
	}
	return @{$self->{_node}};
}

=head3 ext_name

  Arg [0]     : NONE
  ReturnType  : array
  Example     : @leaves = @{$nhx->ext_name};
  Description : Return the list of names of external nodes.

=cut

sub ext_name
{
	my $self = shift;
	my @array;
	foreach my $x ($self->node_array) {
		push(@array, $x->{N}) unless($x->{_C});
	}
	return @array;
}

=head3 n_leaf

  Arg [0]     : NONE
  ReturnType  : int
  Example     : $n_leaf = $nhx->n_leaf;
  Description : Return the number of external nodes.

=cut

sub n_leaf
{
	my $self = shift;
	return $self->{_n_leaf};
}

=head3 root

  Arg [0]     : NONE
  ReturnType  : node reference
  Example     : $root = $nhx->root;
  Description : Return the reference to the root node.

=cut

sub root
{
	my $self = shift;
	return $self->{_root};
}
sub DESTROY
{
	my $self = shift;
	foreach my $p (@{$self->{_node}}) {
		delete($p->{$_}) foreach (keys %$p);
	}
	delete($self->{_leafhash}) if (defined($self->{_leafhash}));
	delete($self->{_node});
	delete($self->{_root});
	delete($self->{_n_leaf});
}

sub string
{
	my $self = shift;
	return $self->string_nhx;
}

=head3 string_nhx

  Arg [0]     : NONE
  ReturnType  : string
  Example     : $nhx_str = $nhx->string_nhx;
  Description : Export the tree to a string in NHX format.

=cut

sub string_nhx # No recursion is used.
{
	my $self = shift;
	# calculate depth
	$self->{_root}->{_depth} = 0;
	foreach my $p (reverse $self->node_array) {
		next if ($p == $self->{_root});
		$p->{_depth} = $p->{_P}->{_depth} + 1;
	}
	# generate string
	my ($cur_depth, $is_first, $str) = (0, 1, '');
	foreach my $p ($self->node_array) {
		my $n = $p->{_depth} - $cur_depth;
		if ($n > 0) {
			if ($is_first) {
				$is_first = 0;
			} else {
				$str .= ",\n";
			}
			$str .= '('x$n;
		} elsif ($n < 0) { # $n == -1
			$str .= "\n)";
		} else { # $n == 0
			$str .= ",\n";
		}
		$str .= $p->{N} if ($p->{N}); # node name
		$str .= ":" . $p->{_L} if (defined($p->{_L}) && $p->{_L} >= 0.0); # length
		my $s = '';
		foreach my $q (keys %$p) { $s .= ":$q=".$p->{$q} if ($q !~ /^_/ && $q ne 'N'); }
		$str .= "[&&NHX$s]" if ($s);
		$cur_depth = $p->{_depth};
		delete($p->{_depth}); # useless now
	}
	$str .= ";\n";
}

=head3 string_tff

  Arg [0]     : NONE
  ReturnType  : string
  Example     : $tab_str = $nhx->string_tft;
  Description : Export the tree to a string in TFF format.

=cut

sub string_tff
{
	my $self = shift;
	my ($str, $s) = ('', '');
	# print out
	foreach my $p ($self->node_array) {
		$s = '';
		foreach my $q (sort keys %$p) { $s .= "\t$q=".$p->{$q} if ($q !~ /^_/); }
		$str .= "$p->{_index}\t" . ($p == $self->root ? $p->{_index} : $p->{_P}->{_index});
		$str .= "\t" . (defined($p->{_L}) ? $p->{_L} : 0) . "$s\n";
	}
	return $str;
}

=head3 get_leaf_by_name

  Arg [1]     : string $name
  ReturnType  : HASH reference
  Example     : my $node_ref = $nhx->get_leaf_by_name('foo');
  Description : Get the reference to a leaf node.

=cut

sub get_leaf_by_name
{
	my ($self, $leafname) = @_;
	if (!defined($self->{_leafhash})) { # build a hash to accelerate node fetching
		my $hash = \%{$self->{_leafhash}};
		%$hash = ();
		foreach my $q (@{$self->{_node}}) {
			unless ($q->{_C}) {
				if ($hash->{$q->{N}}) {
					warn("[LHPM::NHX::get_leaf_by_name] duplicated leaf names: $q->{N}\n");
				}
				$hash->{$q->{N}} = $q;
			}
		}
	}
	return $self->{_leafhash}->{$leafname};
}

=head3 get_lca

  Arg [..]    : array of node references or leaf names
  ReturnType  : HASH reference
  Example     : my $node_ref = $nhx->get_lca('foo', $nhx->get_leaf_by_name('bar'));
  Description : Get the reference to the last common ancestror (LCA) of all the
                input nodes. Three or more nodes can be applied at the same
                time.

=cut

# This is not the most efficient way to find LCA, but it is the simplest way.
sub get_lca
{
	my $self = shift;
	my $p = shift;
	$p = $self->get_leaf_by_name($p) if (ref($p) ne 'HASH');
	foreach my $r (@_) {
		my $q = (ref($r) ne 'HASH')? $self->get_leaf_by_name($r) : $r;
		next unless ($q);
		while ($p != $q) {
			if ($p->{_index} < $q->{_index}) { $p = $p->{_P}; }
			else { $q = $q->{_P}; }
		}
	}
	return $p;
}

=head3 get_tips

  Arg [1]     : HASH reference $q
  ReturnType  : array of leaf names
  Example     : my @leaves = $nhx->get_tips($nhx->get_lca('foo', 'bar'));
  Description : Return an array of names of leaves that descand from the node $q.

=cut

sub get_tips
{
	my ($self, $q) = @_;
	my @array;
	my $nodes = \@{$self->{_node}};
	for (my $i = $q->{_lindex}; $i <= $q->{_index}; ++$i) {
		my $p = $nodes->[$i];
		push(@array, $p->{N}) unless ($p->{_C});
	}
	return @array;
}

sub set_index
{
	my $self = shift;
	my ($j, $k) = (0, 0);
	foreach my $p (@{$self->{_node}}) {
		++$j unless ($p->{_C});
		$p->{_lindex} = ($p->{_C})? $p->{_C}[0]->{_lindex} : $k;
		$p->{_index} = $k++;
	}
	$self->{_n_leaf} = $j;
}

sub get_DFS_array
{
	my $self = shift;
	my $q = shift;
	my $is_node = (@_)? shift : 0;
	my (@stack, @array);
	my $k = 0;
	my $p = \%{$stack[0]};
	$p->{_P} = $q; $p->{i} = 0;
	for (;;) {
		while ($p->{_P}{_C} && $p->{i} != @{$p->{_P}{_C}}) {
			$stack[++$k]{i} = 0;
			$stack[$k]{_P} = $p->{_P}{_C}[$p->{i}];
			$p = \%{$stack[$k]};
		}
		push(@array, $p->{_P});
		$p = \%{$stack[--$k]};
		if ($k >= 0) { ++$p->{i}; }
		else { last; }
	}
	return @array;
}

=head3 update_node_array

  Arg [0]     : NONE
  ReturnType  : NONE
  Example     : $nhx->update_node_array;
  Description : Update @{$self->{_node}}. Nodes appear in suffix order in array
                @{$self->{_node}}. When branches are swapped or nodes are added
                or deleted, @{$self->{_node}} must be updated accordingly. This
                array is the key to various algorthms.

=cut

sub update_node_array
{
	my $self = shift;
	@{$self->{_node}} = $self->get_DFS_array($self->{_root});
	$self->set_index;
}

=head3 get_leaves

  Arg [1|2]   : string $string, [ref $leaves]
  ReturnType  : int
  Example     : $n_leaf = get_leaves($nhx_str, \@leaves);
  Description : Fetch the leaves of the tree $string and store the name in
                %$leaves or @$leaves. If $leaves is not specified, only the
                number of leaves will be returned.

=cut

sub get_leaves
{
	my ($string, $leaves) = @_;
	my $nhx = LHPM::NHX->new;
	$nhx->parse($string);
	return $nhx->n_leaf unless(defined($leaves));
	if (ref($leaves) eq 'HASH') {
		%$leaves = ();
		foreach my $p ($nhx->ext_name) { $leaves->{$p} = 1; }
	} elsif (ref($leaves) eq 'ARRAY') {
		@$leaves = ();
		foreach my $p ($nhx->ext_name) { push(@$leaves, $p); }
	}
	return $nhx->n_leaf;
}

=head3 tff2nhx

=cut

sub tff2nhx
{
	my ($tff) = @_;
	my $nhx = LHPM::NHX->new;
	$nhx->parse_tff($tff);
	return $nhx->string_nhx;
}

=head3 nhx2tff

=cut

sub nhx2tff
{
	my ($nhx_str) = @_;
	my $nhx = LHPM::NHX->new;
	$nhx->parse_nhx($nhx_str);
	return $nhx->string_tff;
}

1;

=head1 LIMITATIONS

This module only supports very simple manipulations on trees. Although it is
not hard to implement advanced functions, I prefer to use my own C library
to achieve these. Perl is not good at algorithmic things after all. Please
check out my NJTREE program, which is distributed under GPL, if you want
more.

=head1 AUTHOR

Heng Li <lh3@sanger.ac.uk>

=cut
