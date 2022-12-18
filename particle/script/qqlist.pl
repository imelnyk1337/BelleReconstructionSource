#!/usr/bin/perl

### struct pdt

package pdt;

sub new
{
	my($type) = @_;

	my($self) = {};

	$self->{m_gen     } = '';
	$self->{m_name    } = '';
	$self->{m_pdgid   } = 0;
	$self->{m_stable  } = 0;
	$self->{m_mass    } = 0.0;
	$self->{m_charge  } = 0.0;
	$self->{m_spin    } = 0.0;
	$self->{m_ctau    } = 0.0;
	$self->{m_width   } = 0.0;
	$self->{m_mass_min} = 0.0;
	$self->{m_mass_max} = 0.0;

	return bless($self, $type);
}


sub dump
{
  my($self) = @_;

	print(
		$self->{m_name    }," ",
		$self->{m_pdgid   }," ",
		$self->{m_stable  }," ",
		$self->{m_mass    }," ",
		$self->{m_charge  }," ",
		$self->{m_spin    }," ",
		$self->{m_ctau    }," ",
		$self->{m_width   }," ",
		$self->{m_mass_min}," ",
		$self->{m_mass_max}," ",
		$self->{m_gen},     "\n"
	);
}



### crt0.o

package main;

my($argc, @argv) = ($#ARGV, @ARGV);
unshift(@argv, $0);
$argc += 2;

exit(&main($argc, @argv));



### user.c


my($Prog);


sub main
{
	my($argc, @argv) = @_;

	$Prog = shift(@argv);

	if( $argc<3 ){
		&usage($Prog);
	}

	my($sortby) = shift(@argv);
	if( $sortby ne '-sortbyname' && $sortby ne '-sortbyid' ){
		&usage($Prog);
	}


	my(@pdt_list);

	my($arg);
	foreach $arg (@argv)
	{
		if     ( $arg =~ m/^-qq=(.*)/     ){
			push(@pdt_list, &process_qq($1));
		} elsif( $arg =~ m/^-evtgen=(.*)/ ){
			push(@pdt_list, &process_evtgen($1));
		} else {
			&usage($Prog);
		}
	}


	my($pdt);
	my(%evtgen_name_hash);
	my(%qq_name_hash);

	# EvtGen info. will override QQ info.
	for $pdt (@pdt_list){
		my($pdgid) = $pdt->{m_pdgid};
		my($name) = $pdt->{m_name};
		$evtgen_name_hash{$pdgid} .= $name."\0" if( $pdt->{m_gen} eq 'evtgen' );
		$qq_name_hash{$pdgid} .= $name."\0" if( $pdt->{m_gen} eq 'qq' );
	}

	my(@pdt_dump_list);
	for $pdt (@pdt_list){
		my($gen) = $pdt->{m_gen};
		my($pdgid) = $pdt->{m_pdgid};

		# to be replaced with particle info in evtgen
		if( $gen eq 'qq' && exists($evtgen_name_hash{$pdgid}) ){
			next;
		}

		if( $gen eq 'evtgen' && exists($qq_name_hash{$pdgid}) ){
			my($name);
			my(@namelist) = split(/\0/, $qq_name_hash{$pdgid});
			foreach $name (@namelist)
			{
				my(%pdttmp) = %$pdt;
				my($pdtobj) = \%pdttmp;
				bless($pdtobj,pdt);
				$pdtobj->{m_name} = $name;
				push(@pdt_dump_list, $pdtobj);
			}
		} else {
			push(@pdt_dump_list, $pdt);
		}
	}


	if( $sortby eq '-sortbyname' ){ @pdt_dump_list = sort({ $a->{m_name} cmp $b->{m_name};} @pdt_dump_list);}
	if( $sortby eq '-sortbyid'   ){ @pdt_dump_list = sort({ $a->{m_pdgid} <=> $b->{m_pdgid};} @pdt_dump_list);}

	for $pdt (@pdt_dump_list){ $pdt->dump()}
	return 0;
}


sub usage
{
	my($prog) = @_;
	print(STDERR
		"usage: $prog -sortbyname| -sortbyid -qq=file1 -evtgen=file2 ...\n");
	exit(1);
}


sub process_qq
{
	my($file) = @_;
	unless( open(F,$file) ){
		die("process_qq: open(\"$file\"): $!\n");
	}

	my(@pdt_list);
	my(%pdt_ref_hash);

	while(<F>){
		m/^\s*(QQBAR|PARTICLE)\s+/i || next;

		my(@tokens) = split;
		my($pname) = $tokens[1];

		my($pdt) = new pdt;
		$pdt_ref_hash{$pname} = $pdt;

		push(@pdt_list, $pdt);

		if     ( $tokens[0] =~ /QQBAR/i ){
			$pdt->{m_gen     } = 'qq';
			$pdt->{m_name    } = $tokens[1];
		} elsif( $tokens[0] =~ /PARTICLE/i ){
			$pdt->{m_gen     } = 'qq';
			$pdt->{m_name    } = $tokens[1];
			$pdt->{m_stable  } = ($tokens[3] != -1) + 0;
			$pdt->{m_mass    } = $tokens[4]+0;
			$pdt->{m_charge  } = $tokens[5]+0;
			$pdt->{m_spin    } = $tokens[6]+0;
			$pdt->{m_ctau    } = $tokens[7]+0;
			$pdt->{m_width   } = defined($tokens[8]) ? $tokens[ 8]+0 : -1;
			$pdt->{m_mass_min} = defined($tokens[8]) ? $tokens[ 9]+0 :  0;
			$pdt->{m_mass_max} = defined($tokens[8]) ? $tokens[10]+0 :  0;
		}
	}

	seek(F, 0, 0);

	while(<F>){
		m/^\s*PDG\s+/i || next;

		my(@tokens) = split;
		my($pname) = $tokens[1];

		my($pdt) = $pdt_ref_hash{$pname};

		$pdt->{m_pdgid   } = $tokens[2];
	}

	close(F);

	return @pdt_list;
}


sub process_evtgen
{
	my($file) = @_;
	unless( open(F,$file) ){
		die("process_evtgen: open(\"$file\"): $!\n");
	}

	my(@pdt_list);
	my(%pdt_ref_hash);

	while(<F>){
		m/^\s*add\s+/ || next;

		my(@tokens) = split;
		my($pname) = $tokens[3];

		my($pdt) = new pdt;
		$pdt_ref_hash{$pname} = $pdt;

		push(@pdt_list, $pdt);

		$pdt->{m_gen     } = 'evtgen';
		$pdt->{m_name    } = $tokens[3];
		$pdt->{m_pdgid   } = $tokens[4];
		$pdt->{m_mass    } = $tokens[5]+0;

		$pdt->{m_charge  } = $tokens[8]/3;
		$pdt->{m_spin    } = $tokens[9]/2;
		$pdt->{m_ctau    } = $tokens[10]/1e3;

		$pdt->{m_width   } = $tokens[6]+0;
		$pdt->{m_mass_min} = $tokens[5]-3*$tokens[6];
		$pdt->{m_mass_max} = $tokens[5]+3*$tokens[6];
	}

	seek(F, 0, 0);

	while(<F>){
		m/^\s*sets\s+p\s+(\S+)\s+isStable\s+(\d+)/ || next;
		if( exists($pdt_ref_hash{$1}) ){
			$pdt_ref_hash{$1}->{m_stable} = $2;
		}
	}

	return @pdt_list;
}

