// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:f;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/apps/pilot/phil/test1.cc
/// @brief  Some simple examples of how to use basic functionality + some DNA functionality
/// @author Phil Bradley (pbradley@fhcrc.org)

// libRosetta headers
#ifndef INCLUDED_apps_pilot_phil_topology_HH
#define INCLUDED_apps_pilot_phil_topology_HH

#include <apps/pilot/phil/phil.hh>
#include <apps/pilot/phil/phil_options.hh>
// #include <apps/pilot/phil/phil_io.hh>
#include <apps/pilot/phil/sym_basic.hh>
#include <apps/pilot/phil/types.hh>


static basic::Tracer TR_TOPOLOGY_HH( "apps.pilot.phil.topology_hh" );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
get_backbone_hbond_interactions(
	bools const & subset,
	Pose const & pose,
	SizePairs & helices,
	SizePairs & strands,
	SizePairs & antiparallel_bridges, // (i,j) with i<j
	SizePairs & parallel_bridges      // (i,j) i hbonded to j-1 and j+1, possibly i>j
)
{
	Real const bb_hbond_maxdis2( 3.4 * 3.4 );
	Size const min_helix_length( 4 ), min_strand_length( 3 ); // arbitrary


	Size const nres( min( pose.total_residue(), subset.size() ) );


	vector1< Sizes > bbn_hbond_partners( pose.total_residue() );

	{ /// figure out all backbone N-->O hbond interactions in the subset
		for ( Size i=1; i<= nres; ++i ) {
			if ( !subset[i] ) continue;
			Residue const irsd( pose.residue(i));
			if ( !irsd.is_protein() ) continue;
			if ( irsd.aa() == aa_pro ) continue; // skip prolines here !!!!!!!!!!!!!!!!!!!!
			for ( Size j=1; j<= nres; ++j ) {
				if ( !subset[j] ) continue;
				Residue const jrsd( pose.residue(j));
				if ( j+2>=i && i+2>= j ) continue; // rule out j=i, j=i+-1, j=i+-2
				if ( !jrsd.is_protein() ) continue;

				if ( irsd.xyz("N").distance_squared( jrsd.xyz("O") ) <= bb_hbond_maxdis2 ) {
					bbn_hbond_partners[i].push_back( j );
				}
			}
		}
	}


	/// look for helical bridges:
	/// bb == 'A' for 4 rsds, i<--i+4 hbond

	string allbb( string("X")+torsion2big_bin_string( 1, nres, pose, true ) );

	bools is_hbonded_helix( pose.total_residue(), false );
	for ( Size i=1; i+4<= nres; ++i ) {
		if ( subset[i] && subset[i+1] && subset[i+2] && subset[i+3] && subset[i+4] &&
				allbb.substr(i,4) == "AAAA" &&
				has_element( bbn_hbond_partners[i+4], i ) ) {
			for ( Size j=i; j<= i+4; ++j ) is_hbonded_helix[j] = true;
		}
		// pi helix
		if ( subset[i] && subset[i+1] && subset[i+2] && subset[i+3] &&
				allbb.substr(i,3) == "AAA" &&
				has_element( bbn_hbond_partners[i+3], i ) ) {
			for ( Size j=i; j<= i+3; ++j ) is_hbonded_helix[j] = true;
		}
	}

	/// look for beta bridges
	string const beta_bridge_bbs( "BE" );
	antiparallel_bridges.clear();
	parallel_bridges.clear();
	bools is_beta_paired( nres,false );
	for ( Size i=1; i<= nres; ++i ) {
		if ( !subset[i] ) continue;
		if ( beta_bridge_bbs.find( allbb[i] ) == string::npos ) continue;
		Sizes const & partners( bbn_hbond_partners[i] );
		for ( Sizes::const_iterator p= partners.begin(); p!= partners.end(); ++p ) {
			Size const j( *p );
			if ( !subset[j] ) continue;
			// antiparallel (i:j)
			if ( beta_bridge_bbs.find( allbb[j] ) != string::npos &&
					has_element( bbn_hbond_partners[j], i ) && i < j ) {
				is_beta_paired[i] = is_beta_paired[j] = true;
				antiparallel_bridges.push_back( make_pair( i, j ) );
			}
			// parallel (i:j+1)  i-->j, j+2-->i
			if ( j+2<= nres && subset[j+2] && has_element( bbn_hbond_partners[j+2], i ) &&
					beta_bridge_bbs.find( allbb[j+1] ) != string::npos ) {
				is_beta_paired[i] = is_beta_paired[j+1] = true;
				parallel_bridges.push_back( make_pair( i, j+1 ) );
			}
		}
	}

	// locate helices, strands
	helices.clear();
	{
		for ( Size hbegin=1; hbegin< nres; ++hbegin ) {
			if ( allbb[hbegin] == 'A' && ( hbegin==1 || allbb[hbegin-1] != 'A' ) ) {
				Size hend( hbegin );
				while ( hend <nres && allbb[hend+1] == 'A' ) ++hend; // go to end of 'A' stretch
				Size const hlen( hend-hbegin+1 );
				if ( hlen<min_helix_length ) continue;
				bool has_helical_hbonds( false );
				for ( Size i=hbegin; i<= hend; ++i ) has_helical_hbonds = ( has_helical_hbonds || is_hbonded_helix[i] );
				if ( has_helical_hbonds ) helices.push_back( make_pair( hbegin, hend ) );
			}
		}
	}


	strands.clear();
	{
		for ( Size sbegin=1; sbegin< nres; ++sbegin ) {
			if ( allbb[sbegin] == 'B' && ( sbegin==1 || allbb[sbegin-1] != 'B' ) ) {
				Size send( sbegin );
				while ( send <nres && ( ( allbb[send+1] == 'B' ) ||
						( allbb[send+1] == 'E' && send+1<nres && allbb[send+2]=='B' && is_beta_paired[send+1] &&
						( is_beta_paired[send+2] || ( send+2<nres && allbb[send+3]=='B' &&
						is_beta_paired[ send+3 ] ) ) ) ) ) ++send;
				Size const slen( send-sbegin+1 );
				if ( slen<min_strand_length ) continue;
				bool has_beta_hbonds( false );
				for ( Size i=sbegin; i<= send; ++i ) has_beta_hbonds = ( has_beta_hbonds || is_beta_paired[i] );
				if ( has_beta_hbonds ) strands.push_back( make_pair( sbegin, send ) );
			}
		}
	}

}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct SSE {

	int begin, end;
	SegmentType type;

	bool
	contains( SSE const & other ) const;

	bool
	contains( int const pos ) const;

};

bool
SSE::contains( SSE const & other ) const {
	return begin <= other.begin && end >= other.end;
}

bool
SSE::contains( int const pos ) const {
	return begin <= pos && end >= pos;
}

bool
operator<( SSE const & a, SSE const & b ) {
	// int const
	//  ab( a.begin < a.end ? a.begin : a.end ),
	//  ae( a.begin > a.end ? a.begin : a.end ),
	//  bb( b.begin < b.end ? b.begin : b.end ),
	//  be( b.begin > b.end ? b.begin : b.end );
	// return ( ab<bb || ab==bb && ae<be );
	return ( a.begin < b.begin || ( a.begin == b.begin && a.end < b.end ) );
}

bool
operator==( SSE const & a, SSE const & b ) {
	return ( a.begin == b.begin && a.end == b.end && a.type == b.type );
}

typedef vector1< SSE > SSEs;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// this assumes the pose is a repeat pose, ie calls deduce_repeatlen...
//
bool // return TRUE on success
parse_topology_parameters_from_pose(
	Pose const & pose,
	string & topology,
	Size & nrepeat,
	Size & repeatlen,
	string & repeatseq,
	string & repeatbb,
	string & repeatss,
	Sizes & sse_lens,
	SegmentTypes & sse_types,
	strings & turns,
	vector1< SizePairs > & potential_jump_points,
	vector1< Sizes > & potential_cutpoints,
	BetaPairingTypes & jump_types_single_repeat,
	string const tag = string()
)
{
	Size nres( pose.total_residue() );
	while ( pose.residue(nres).name() == "VRT" ) --nres;
	runtime_assert( pose.residue( nres ).is_protein() );

	for ( Size i=1; i< nres; ++i ) {
		if ( ( pose.residue( i  ).has_variant_type( CUTPOINT_LOWER ) &&
				pose.residue( i+1).has_variant_type( CUTPOINT_UPPER ) ) ) runtime_assert( pose.fold_tree().is_cutpoint(i));
	}

	string const fullseq( pose.sequence().substr(0,nres));

	repeatlen = deduce_repeatlen( fullseq ); // returns 0 on failure

	if ( !repeatlen ) return false;
	runtime_assert( nres%repeatlen == 0 );

	nrepeat = nres / repeatlen;
	//runtime_assert( pose.total_residue()%nrepeat == 0 );
	//int const repeatlen( pose.total_residue()/nrepeat );

	string const fake_repeatbb( torsion2big_bin_string( repeatlen+1, 2*repeatlen, pose ) ); // 2nd repeat
	string const fullbb( torsion2big_bin_string( 1, nres, pose ) ); // 0-indexed!

	bools subset( pose.total_residue(), true );
	if ( nres<pose.total_residue() ) for ( Size i=nres+1; i<= pose.total_residue(); ++i ) subset[ i ] = false;

	SizePairs helices, strands, antiparallel_bridges, parallel_bridges;
	get_backbone_hbond_interactions( subset, pose, helices, strands, antiparallel_bridges, parallel_bridges );


	SSEs sses;
	for ( SizePairs::const_iterator h= helices.begin(); h!= helices.end(); ++h ) {
		if ( h->first <= 2 || h->second >= nres-1 ) continue; // exclude potentially truncated sses
		SSE sse;
		sse.begin = ( int(h->first )-1 )%repeatlen+1;
		sse.end   = ( int(h->second)-1 )%repeatlen+1;
		sse.type = alpha_type;
		if ( !has_element( sses, sse ) ) sses.push_back( sse );
	}
	for ( SizePairs::const_iterator s= strands.begin(); s!= strands.end(); ++s ) {
		if ( s->first <= 2 || s->second >= nres-1 ) continue; // exclude potentially truncated sses
		SSE sse;
		sse.begin = ( int(s->first )-1 )%repeatlen+1;
		sse.end   = ( int(s->second)-1 )%repeatlen+1;
		sse.type = beta_type;
		if ( !has_element( sses, sse ) ) sses.push_back( sse );
	}
	for ( SSEs::iterator sse=sses.begin(); sse!= sses.end(); ++sse ) {
		runtime_assert( sse->begin >= 1 && sse->end <= (int)repeatlen );
		if ( sse->begin > sse->end ) {
			// should it be at the beginning or at the end?
			if ( sse->end - 1 > (int)repeatlen - sse->begin ) { // more at the Nterminus
				sse->begin -= int(repeatlen);
			} else {
				sse->end += int(repeatlen);
			}
			runtime_assert( sse->begin < sse->end );
		}
	}


	std::sort( sses.begin(), sses.end() );

	{ // get nr set
		SSEs nr_sses;
		for ( Size i=1; i<= sses.size(); ++i ) {
			// check for redundancy
			bool addme( true );
			for ( Size j=1; j<= nr_sses.size(); ++j ) {
				if ( sses[i].contains( nr_sses[j] ) ) {
					nr_sses[j] = sses[i];
					addme = false;
				} else if ( nr_sses[j].contains( sses[i] ) ) addme = false;
			}
			if ( addme ) nr_sses.push_back( sses[i] );
		}
		sses.swap( nr_sses );
	}

	sse_lens.clear();
	sse_types.clear();
	turns.clear();

	string bbtag;

	for ( Size i=1; i<= sses.size(); ++i ) {
		SSE const & sse( sses[i] );
		runtime_assert( sse.begin <= sse.end );
		TR_TOPOLOGY_HH.Trace << "SSE: " << i << ' ' << sse.type << ' ' << sse.begin << ' ' << sse.end << ' ' <<
			fake_repeatbb << ' ' << tag << endl;
	}
	for ( Size i=1; i<= sses.size(); ++i ) {
		SSE const & sse( sses[i] );
		runtime_assert( sse.begin <= sse.end );
		TR_TOPOLOGY_HH.Trace << "SSE: " << i << ' ' << sse.type << ' ' << sse.begin << ' ' << sse.end << ' ' <<
			fake_repeatbb << ' ' << tag << endl;
		sse_lens.push_back( sse.end - sse.begin+ 1 );
		sse_types.push_back( sse.type );
		int next_begin(0);
		if ( i<sses.size() ) next_begin = sses[i+1].begin;
		else next_begin = sses[1].begin + int(repeatlen);
		TR_TOPOLOGY_HH.Trace << sses.size() << ' ' << next_begin << ' ' << sse.end << endl;
		runtime_assert( next_begin > sse.end );
		runtime_assert( next_begin < (int)fullbb.size() );
		runtime_assert( sse.end > 1 ); // fullbb at pos 1 is bogus
		int const turnlen( next_begin - sse.end - 1 );
		if ( turnlen>0 ) turns.push_back( fullbb.substr( sse.end, turnlen ) );
		else turns.push_back( "-" );
		bbtag += string_of( sse_lens.back() )+"."+turns.back()+".";
	}
	bbtag.erase( bbtag.size()-1 );

	if ( tag.size() && tag.find( bbtag ) == string::npos ) {
		TR_TOPOLOGY_HH.Trace << "parse_topology_parameters_from_pose bbtagchange: " << bbtag << ' ' << tag << endl;
		TR_TOPOLOGY_HH.Trace << "parse_topology_parameters_from_pose fullbb: " << fullbb.substr(0,repeatlen) << ' ' <<
			fullbb.substr(repeatlen,2*repeatlen) << ' ' << fullbb.substr(2*repeatlen) << ' ' << tag << endl;
	}

	int const shift( 1 - sses.front().begin ); // we add this to residue numbers in the pose to get resampled rsd #s

	///  figure out where the potential jump_points, cutpoints could be
	///
	/// look for intra- and inter- repeat jumps

	topology.clear();
	Size const num_sses( sses.size() );
	Size num_strands(0);
	Sizes sse2strand;
	for ( Size i=1; i<= num_sses; ++i ) {
		if ( sses[i].type == beta_type ) {
			++num_strands;
			sse2strand.push_back( num_strands );
			topology.push_back('B');
		} else {
			sse2strand.push_back( 0 );
			topology.push_back('A');
		}
	}

	potential_jump_points.clear();
	potential_cutpoints.clear();

	for ( Size rr=1; rr<= 2; ++rr ) { // 1st time= intra-, 2nd=inter
		int const sse2_offset( rr == 1 ? 0 : repeatlen );

		for ( Size ii=1; ii<= num_sses; ++ii ) {
			SSE const & sse1( sses[ii] );
			if ( sse1.type != beta_type ) continue;
			Size const strand1( sse2strand[ii] );
			for ( Size jj=1; jj<= num_sses; ++jj ) {
				SSE const & sse2( sses[jj] );
				if ( sse2.type != beta_type ) continue;
				if ( sse2_offset==0 && jj<=ii ) continue;
				Size strand2( sse2strand[jj] );
				if ( sse2_offset ) strand2 += num_strands;

				SizePairs jumps_this_pairing;

				BetaPairingType jump_type( none_type );

				// look for beta bridges
				for ( SizePairs::const_iterator b= antiparallel_bridges.begin(); b!= antiparallel_bridges.end(); ++b ) {
					if ( sse1.contains( b->first ) && sse2.contains( int(b->second)-sse2_offset ) ) {
						if ( ii == jj ) { // only include symmetric jump points if same-strand pairing
							if ( b->first != b->second - repeatlen ) continue;
							jump_type = symmetrized_antiparallel_type;
						} else {
							jump_type = antiparallel_type;
						}
						jumps_this_pairing.push_back( make_pair( int( b->first )+shift, int( b->second )+shift ) );
						TR_TOPOLOGY_HH.Trace << "newpairing: " << jump_type << " A " << ii << ' ' << jj << ' ' <<
							strand1 << ' ' << strand2 << ' ' <<
							I(4,b->first) << I(4,b->second) <<
							I(4,jumps_this_pairing.back().first) << I(4,jumps_this_pairing.back().second) << I(4,repeatlen) <<
							I(2,shift) << ' ' << tag << endl;
					}
				}

				for ( SizePairs::const_iterator b= parallel_bridges.begin(); b!= parallel_bridges.end(); ++b ) {
					if ( sse1.contains( b->first ) && sse2.contains( int(b->second)-sse2_offset ) ) {
						jumps_this_pairing.push_back( make_pair( int( b->first )+shift, int( b->second )+shift ) );
						TR_TOPOLOGY_HH.Trace << "newpairing: " << jump_type << " P " << ii << ' ' << jj << ' ' <<
							strand1 << ' ' << strand2 << ' ' <<
							I(4,b->first) << I(4,b->second) <<
							I(4,jumps_this_pairing.back().first) << I(4,jumps_this_pairing.back().second) << I(4,repeatlen) <<
							I(2,shift) << ' ' << tag << endl;
						runtime_assert( jump_type != antiparallel_type );
						jump_type = parallel_type;
					}
				}

				if ( jumps_this_pairing.size() ) {
					potential_jump_points.push_back( jumps_this_pairing );
					runtime_assert( jump_type != none_type );
					jump_types_single_repeat.push_back( jump_type );

					Sizes cutpoints;
					if ( sse2_offset == 0 ) { // intra-repeat
						runtime_assert( jj == ii+1 ); // just hairpins at the moment
						runtime_assert( jump_type == antiparallel_type );
						for ( int cut= sse1.end; cut<sse2.begin; ++cut ) {
							cutpoints.push_back( cut+shift );
						}
					} else {
						cutpoints.push_back( repeatlen );
					}
					potential_cutpoints.push_back( cutpoints );
					if ( jump_type == antiparallel_type || jump_type == symmetrized_antiparallel_type ) {
						topology += string("a") + string_of(strand1)+string_of(strand2);
					} else topology += string("p")+string_of(strand1)+string_of(strand2);
				}
			}
		}
	}


	/// repeat strings
	repeatseq.clear();
	repeatbb.clear();
	repeatss.clear();
	for ( int i=1; i<= (int)repeatlen; ++i ) {
		int posepos( i-shift );
		if ( posepos <= 1 ) posepos += (int) repeatlen;
		runtime_assert( posepos<(int)nres );
		repeatseq.push_back( pose.residue( posepos ).name1() );
		repeatbb.push_back( fullbb[ posepos-1 ] );
		repeatss.push_back( pose.secstruct( posepos ) );
	}

	return true;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////
//
// this also sets up the symmetrized_antiparallel_type transforms
//
//
//
void
load_beta_transforms( vector1< RTs> & beta_transforms )
{

	// load beta transforms
	beta_transforms.clear();
	beta_transforms.resize( none_type );
	ifstream data( option[ my_options::beta_transforms_file ]().c_str() );

	runtime_assert( data.good() );

	strings ab_turns, ba_turns;
	if ( option[ my_options::ab_turns ].user() ) ab_turns = option[ my_options::ab_turns ]();
	if ( option[ my_options::ba_turns ].user() ) ba_turns = option[ my_options::ba_turns ]();

	string line, tag;
	while ( getline( data, line ) ) {
		istringstream l( line );
		l >> tag;
		if ( l.fail() ) continue;
		if ( tag == "parallel_beta_transform" ||
				tag == "antiparallel_beta_transform" ||
				tag == "antiparallel2_beta_transform" ) {
			RT rt;
			l >> rt;
			if ( !l.fail() ) {
				if (      tag ==     "parallel_beta_transform" ) beta_transforms[      parallel_type ].push_back( rt );
				else if ( tag == "antiparallel_beta_transform" ) beta_transforms[  antiparallel_type ].push_back( rt );
				else                                             beta_transforms[ antiparallel2_type ].push_back( rt );
			}
		} else if ( tag == "turn_transform:" ) {
			RT rt;
			Size turnlen;
			string turnbb, turnseq;
			l >> tag >> turnlen >> turnbb >> turnseq >> rt;
			if ( l.fail() ) {
				TR_TOPOLOGY_HH.Trace << "bad line: " << line << endl;
			} else if ( tag == "alpha_beta" && ( ab_turns.empty() || has_element( ab_turns, turnbb ) ) ) {
				beta_transforms[ ab_turn_type ].push_back( rt );
			} else if ( tag == "beta_alpha" && ( ba_turns.empty() || has_element( ba_turns, turnbb ) ) ) {
				beta_transforms[ ba_turn_type ].push_back( rt );
			}
		}
	}
	TR_TOPOLOGY_HH.Trace << "read " << beta_transforms[     parallel_type ].size() <<     " parallel_beta_transforms" << endl;
	TR_TOPOLOGY_HH.Trace << "read " << beta_transforms[ antiparallel_type ].size() << " antiparallel_beta_transforms" << endl;
	TR_TOPOLOGY_HH.Trace << "read " << beta_transforms[ antiparallel2_type ].size() << " antiparallel2_beta_transforms" << endl;
	TR_TOPOLOGY_HH.Trace << "read " << beta_transforms[ ab_turn_type ].size() << " ab_turn_transforms" << endl;
	TR_TOPOLOGY_HH.Trace << "read " << beta_transforms[ ba_turn_type ].size() << " ba_turn_transforms" << endl;
	// symmetrize the antiparallel_beta_transforms

	foreach_ ( BetaPairingType bptype, make_vector1( antiparallel_type, antiparallel2_type ) ) {
		BetaPairingType const symmtype( bptype == antiparallel_type ? symmetrized_antiparallel_type :
			symmetrized_antiparallel2_type );
		for ( Size i=1; i<= beta_transforms[ bptype ].size(); ++i ) {
			RT const & rt( beta_transforms[ bptype ][i] );
			Real theta;
			Vector n = numeric::rotation_axis( rt.get_rotation(), theta );
			RT rtnew;
			rtnew.set_rotation( numeric::rotation_matrix( n, numeric::constants::d::pi ) ); // 180 degree rotation
			rtnew.set_translation( rt.get_translation() - n.dot( rt.get_translation() ) * n );
			runtime_assert( fabs( rtnew.get_translation().dot( n ) )<1e-3 );
			//TR.Trace << "symtransformdis: " << F(9,3,sqrt( rtnew.distance_squared( rt ) ) ) << endl;
			beta_transforms[ symmtype ].push_back( rtnew );
		}
	}
} //scope for loading transforms

////////////////////////////////////////////
void
get_beta_pairing_geometry(
	Vectors const & coords,
	Size const pos1,
	Size const pos2,
	Real & orientation,
	Reals & distances
)
{
	distances.clear();
	orientation = ( ( coords[pos1+1] - coords[pos1] ).cross( coords[pos1+2] - coords[pos1] ).dot( coords[pos2+1] -
		coords[pos1+1] ) );
	for ( Size i=pos1; i<= pos1+2; ++i ) {
		for ( Size j=pos2; j<= pos2+2; ++j ) {
			distances.push_back( coords[i].distance( coords[j] ) );
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// scores the triplets {pos1,pos1+1,pos1+2} aligned with {pos2,pos2+1,pos2+2}
//
Real
compute_parallel_beta_pairing_score(
	Vectors const & coords,
	Size const pos1,
	Size const pos2
)
{
	///
	static vector1< Reals > all_pos_distances, all_neg_distances;
	if ( all_pos_distances.empty() ) {
		string const filename( option[ my_options::beta_pairings_file ] );
		ifstream data( filename.c_str());
		runtime_assert( data.good() );
		string line;
		while ( getline( data, line ) ) {
			strings const l( split_to_vector1( line ) );
			if ( l.size() == 11 && l[1] == "beta_pairing" ) {
				Reals geoms;
				for ( Size i=2; i<= 11; ++i ) {
					runtime_assert( is_float( l[i] ) );
					geoms.push_back( float_of( l[i] ) );
				}
				if ( geoms[1]<0 ) all_neg_distances.push_back( geoms );
				else all_pos_distances.push_back( geoms );
			}
		}
		data.close();

		{ // debugging
			for ( Size r=1; r<= 2; ++r ) {
				vector1< Reals > const & ref_distances( r == 1 ? all_pos_distances : all_neg_distances );

				for ( Size k=1; k<= ref_distances.size(); ++k ) {
					Real const orientation( ref_distances[k][1] );
					Reals distances;
					for ( Size j=1; j<= 9; ++j ) distances.push_back( ref_distances[k][j+1] );

					{ // compute meandev for this guy
						Real meandev(0);

						for ( Size i=1; i<= ref_distances.size(); ++i ) {
							Real dev( numeric::square( ( ref_distances[i][1] - orientation )/10 ) );
							for ( Size j=2; j<= 10; ++j ) {
								dev += numeric::square( ref_distances[i][j] - distances[j-1] );
							}
							dev = sqrt( dev/10 );
							meandev += dev;
						}
						meandev /= ref_distances.size();
						//TR.Trace << "meandev_parallel " << F(9,3,meandev) << I(3,r) << I(6,k) << I(6,ref_distances.size()) << endl;
					}
				}
			}
		}
	}


	Real orientation;
	Reals distances;
	get_beta_pairing_geometry( coords, pos1, pos2, orientation, distances );

	vector1< Reals > const & ref_distances( orientation < 0 ? all_neg_distances : all_pos_distances );

	Real meandev(0);

	for ( Size i=1; i<= ref_distances.size(); ++i ) {
		Real dev( numeric::square( ( ref_distances[i][1] - orientation )/10 ) );
		for ( Size j=2; j<= 10; ++j ) {
			dev += numeric::square( ref_distances[i][j] - distances[j-1] );
		}
		dev = sqrt( dev/10 );
		meandev += dev;
	}
	meandev /= ref_distances.size();

	return meandev;

}

/// scores pos1,pos1+1,pos1+2 vs pos2,pos2+1,pos2+2
///
/// so pos1+1 and pos2+1 should be in alignment (and pos1,pos2+2 and pos1+2,pos2)
///

Real
compute_antiparallel_beta_pairing_score(
	Vectors const & coords,
	Size const pos1,
	Size const pos2
)
{
	///
	static vector1< Reals > all_neg_distances;
	if ( all_neg_distances.empty() ) {
		string const filename( option[ my_options::beta_pairings_file ] );
		ifstream data( filename.c_str());
		runtime_assert( data.good() );
		string line;
		while ( getline( data, line ) ) {
			strings const l( split_to_vector1( line ) );
			if ( l.size() == 11 && l[1] == "antiparallel_beta_pairing" ) {
				Reals geoms;
				for ( Size i=2; i<= 11; ++i ) {
					runtime_assert( is_float( l[i] ) );
					geoms.push_back( float_of( l[i] ) );
				}
				if ( geoms[1]<0 ) all_neg_distances.push_back( geoms );
				else {
					//TR.Trace << "ignoring wonky ap beta pairing w/ positive orientation" << endl;
				}
			}
		}
		data.close();

		{ // debugging
			vector1< Reals > const & ref_distances( all_neg_distances );

			for ( Size k=1; k<= ref_distances.size(); ++k ) {
				Real const orientation( ref_distances[k][1] );
				Reals distances;
				for ( Size j=1; j<= 9; ++j ) distances.push_back( ref_distances[k][j+1] );

				{ // compute meandev for this guy
					Real meandev(0);

					for ( Size i=1; i<= ref_distances.size(); ++i ) {
						Real dev( numeric::square( ( ref_distances[i][1] - orientation )/10 ) );
						for ( Size j=2; j<= 10; ++j ) {
							dev += numeric::square( ref_distances[i][j] - distances[j-1] );
						}
						dev = sqrt( dev/10 );
						meandev += dev;
					}
					meandev /= ref_distances.size();
					//TR.Trace << "meandev_antiparallel " << F(9,3,meandev) << I(6,k) << I(6,ref_distances.size()) << endl;
				}
			}
		}
	}


	Real orientation;
	Reals distances;
	get_beta_pairing_geometry( coords, pos1, pos2, orientation, distances );

	if ( orientation > 0 ) return 100.0; // only score pairings that could be a true antiparallel beta bridge

	vector1< Reals > const & ref_distances( all_neg_distances );

	Real meandev(0);

	for ( Size i=1; i<= ref_distances.size(); ++i ) {
		Real dev( numeric::square( ( ref_distances[i][1] - orientation )/10 ) );
		for ( Size j=2; j<= 10; ++j ) {
			dev += numeric::square( ref_distances[i][j] - distances[j-1] );
		}
		dev = sqrt( dev/10 );
		meandev += dev;
	}
	meandev /= ref_distances.size();

	return meandev;

}

inline
int
get_strand1_pos_antiparallel(
	int const strand2_pos,
	int const offset,
	int const strand1_begin,
	int const strand2_end
)
{
	//                     |---------------- i -----------------|
	return strand1_begin + ( strand2_end - strand2_pos - offset );
}

inline
int
get_strand2_pos_antiparallel(
	int const strand1_pos,
	int const offset,
	int const strand1_begin,
	int const strand2_end
)
{
	//                     |---------- i ------------|
	return strand2_end - ( strand1_pos - strand1_begin + offset ); // = strand2_pos
}

void
get_antiparallel_strand_pairing_scores(
	Vectors const & coords,
	int const strand1_begin,
	int const strand1_end,
	int const strand2_begin,
	int const strand2_end,
	Reals & strand_pairing_scores,
	int & strand_pairing_offset,
	std::pair< Size, Size > & best_beta_bridge // would be forming double hbond pair
)
{
	// initialize
	strand_pairing_scores.clear();
	strand_pairing_offset = 0;
	best_beta_bridge = make_pair(0,0);


	int const strand1_len( strand1_end-strand1_begin+1 ), strand2_len( strand2_end-strand2_begin+1 );
	/// choose these so we have a minimum overlap of three
	int const min_overlap( 3 );
	runtime_assert( strand1_len >= min_overlap && strand2_len >= min_overlap );
	int const min_offset( -1 * ( strand1_len - min_overlap ) ), max_offset( strand2_len - min_overlap );
	/// look for three consecutive coords that look "paired"
	Real best_bscore( 1e6 );
	//Reals best_bscores;
	Size const max_dis2( 5.7 * 5.7 );
	int pos1,pos2;
	// Size best_count(0);
	// int best_offset(0);
	Reals bscores; bscores.reserve( 10 );
	for ( int offset = min_offset; offset <= max_offset; offset+= 1 ) {
		Size count(0);
		bool this_offset_is_best_so_far( false );
		for ( int i=max(0,-1*offset),ie=min(strand1_len,strand2_len-offset); i<ie; ++i ) {
			pos1 = strand1_begin + i;
			pos2 = strand2_end - ( i + offset ); // i+offset is guaranteed to be positive
			runtime_assert( pos2 >= strand2_begin && pos2 <= strand2_end ); // take this out later, sanity check
			if ( coords[pos1].distance_squared( coords[pos2] ) < max_dis2 ) {
				++count;
				if ( count >= 3 ) { // 3 in a row, look at geometry
					Real const bscore( compute_antiparallel_beta_pairing_score( coords, pos1-2, pos2 ) );
					bscores.push_back( bscore );
					if ( bscore < best_bscore ) {
						best_bscore = bscore;
						this_offset_is_best_so_far = true;
						best_beta_bridge = make_pair( Size(pos1-1), Size(pos2+1) );
						runtime_assert( get_strand2_pos_antiparallel( pos1-1, offset, strand1_begin, strand2_end ) == pos2+1 );
						runtime_assert( get_strand1_pos_antiparallel( pos2+1, offset, strand1_begin, strand2_end ) == pos1-1 );
					}
				}
			} else {
				count =0;
			}
		}
		if ( this_offset_is_best_so_far ) {
			strand_pairing_scores.swap( bscores );
			strand_pairing_offset = offset;
		}
		bscores.clear(); // get ready for next round
	} // offset

}



void
get_parallel_strand_pairing_scores(
	Vectors const & coords,
	int const strand1_begin,
	int const strand1_end,
	int const strand2_begin,
	int const strand2_end,
	Reals & strand_pairing_scores,
	int & strand_pairing_offset
)
{
	int const strand1_len( strand1_end-strand1_begin+1 ), strand2_len( strand2_end-strand2_begin+1 );
	int const min_overlap( 3 );
	runtime_assert( strand1_len >= min_overlap && strand2_len >= min_overlap );
	int const min_offset( -1*( strand1_len - min_overlap ) ), max_offset( strand2_len - min_overlap );
	/// look for three consecutive coords that look "paired"
	Real best_bscore( 1e6 );
	//Reals best_bscores;
	Size const max_dis2( 5.7 * 5.7 );
	int pos1,pos2;
	// Size best_count(0);
	// int best_offset(0);
	Reals bscores;
	for ( int offset = min_offset; offset <= max_offset; offset+= 1 ) {
		Size count(0);
		for ( int i=max(0,-1*offset),ie=min(strand1_len,strand2_len-offset); i<ie; ++i ) {
			pos1 = strand1_begin + i;
			pos2 = strand2_begin + i + offset;
			if ( coords[pos1].distance_squared( coords[pos2] ) < max_dis2 ) {
				++count;
				if ( count >= 3 ) { // 3 in a row, look at geometry
					bscores.push_back( compute_parallel_beta_pairing_score( coords, pos1-2, pos2-2 ) );
				}
			} else {
				count =0;
			}
		}
		if ( !bscores.empty() ) {
			Real const minbscore( min( bscores ) );
			if ( minbscore < best_bscore ) {
				best_bscore = minbscore;
				strand_pairing_scores.swap( bscores );
				strand_pairing_offset = offset;
			}
			bscores.clear(); // get ready for next round
		}
	} // offset

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string
get_hbond_ss( Pose const & pose )
{
	string hbond_ss( pose.total_residue(), 'L' );

	SizePairs helices, strands, antiparallel_bridges, parallel_bridges;

	get_backbone_hbond_interactions( bools( pose.total_residue(), true ), pose, helices, strands, antiparallel_bridges,
		parallel_bridges );

	foreach_( SizePair p, helices ) {
		for ( Size i=p.first; i<= p.second; ++i ) hbond_ss[i-1] = 'H';
	}

	foreach_( SizePair p, strands ) {
		for ( Size i=p.first; i<= p.second; ++i ) hbond_ss[i-1] = 'E';
	}

	return hbond_ss;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif
