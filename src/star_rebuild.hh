// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief
#ifndef INCLUDED_apps_pilot_phil_star_rebuild_HH
#define INCLUDED_apps_pilot_phil_star_rebuild_HH

#include <apps/pilot/phil/phil.hh>
#include <apps/pilot/phil/types.hh>
#include <apps/pilot/phil/phil_io.hh>
#include <apps/pilot/phil/rms.hh>
#include <apps/pilot/phil/phil_options.hh>
#include <apps/pilot/phil/fragments.hh>
#include <apps/pilot/phil/rotamer_correctness.hh>

#include <apps/pilot/phil/constraints.hh> // MinMaxFunc


#include <core/pack/pack_rotamers.hh>
// #include <core/pack/pack_rotamers_envdep.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.hh>
#include <core/conformation/util.hh>
#include <core/scoring/Energies.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/dna/water.hh>
#include <protocols/loops/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_KIC.hh>

#include <devel/dna/relax_util.hh>
#include <devel/dna/DoubleFragment.hh>
#include <devel/blab/move_map_util.hh>
// #include <devel/blab/loops/util.hh>
// #include <devel/blab/loops/protocols.hh>
#include <devel/blab/classic_frags/legacy_loop_util.hh>

// #include <devel/blab/opte/sidechain_relax.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

using namespace std; // try it out
using id::AtomID;

basic::Tracer TR_REBUILD( "apps.pilot.phil.star_rebuild" );



///////////////////////////////////////////////////////////////////////////////////////
//
// uses the (dssp-based) secondary structure to cut the pose into segments
//
//
void
setup_protein_segments(
	Pose const & pose,
	vector1< Sizes > &protein_segments,
	bool const use_randomness = true,
	bools const & is_frozen = empty_bools
)
{
	if ( !is_frozen.empty() ) { runtime_assert( is_frozen.size() == pose.total_residue() ); }
	protein_segments.clear();


	//  if ( has_motif_data( pose ) && get_motif_data( pose ).has_segment( "PROTEIN_SEGMENTS") ) {
	//   TR_REBUILD.Trace << "setup_protein_segments: using MotifData" << endl;
	//   Sizes const & segs( get_motif_data( pose ).segment("PROTEIN_SEGMENTS") );
	//   pbassert( segs.size()%2==0 );
	//   for ( Size ii=0; ii< segs.size()/2; ++ii ) {
	//    protein_segments.push_back( make_vector1( segs[2*ii+1], segs[2*ii+2] ) );
	//   }
	//   return;
	//  }


	Size const min_loop_segment_size( 6 ), max_loop_segment_size( 13 ), target_loop_segment_size( 8 ); // split loops>=14 into 8+x
	Size const min_segment_size( is_frozen.empty() ? 0 : 3 );


	for ( Size c=1; c<= num_chains( pose ); ++c ) {
		Size const cb( chain_begin(pose,c) ), ce( chain_end(pose,c ) ), chainlen( ce-cb+1 );
		if ( pose.residue( ce ).is_protein() ) {
			string ss( string("X") + pose.secstruct().substr( cb-1, ce-cb+1 ) ); // 1-indexed!
			if ( is_frozen.size() ) {
				for ( Size i=1; i<= chainlen; ++i ) {
					if ( is_frozen[cb+i-1] ) ss[i] = 'F';
				}
			}

			// find contiguous stretches...
			vector1< Sizes > segments;
			for ( Size i=1; i<= chainlen; ) {
				Size j(i);
				while ( j<chainlen && ss[j+1] == ss[i] ) ++j;
				segments.push_back( make_vector1( i,j) );
				i = j+1;
			}

			// look for short loop segments
			for ( Size ii=1; ii<= segments.size(); ++ii ) {
				Size const i( segments[ii].front() ), j( segments[ii].back() ), seglen( j-i+1 );
				//pbassert( ss[ i ] == ss[j] );
				string const segss( ss.substr(i,seglen) );
				bool const isloop( segss == string( seglen, 'L' ) ), isflexible( segss != string( seglen, 'F' ) );
				if ( ( isloop || ( isflexible && seglen<min_segment_size ) ) &&
						seglen < min_loop_segment_size && segments.size()>1 ) {
					TR_REBUILD.Trace << "shortseg: " << segss << ' ' << i << ' ' << j << endl;
					/// merge onto ss elements
					bool delete_me( false ),
						next_segment_is_flexible( ii<segments.size() && ss[ segments[ii+1].front() ] != 'F' ),
						previous_segment_is_flexible( ii>1 && ss[ segments[ii-1].front() ] != 'F' );
					if ( ii == 1 ) {
						// n-terminal loop
						if ( next_segment_is_flexible ) {
							segments[ii+1].front() = i;
							delete_me = true;
						}

					} else if ( ii == segments.size() ) {
						// c-terminal loop
						if ( previous_segment_is_flexible ) {
							segments[ii-1].back() = j;
							delete_me = true;
						}

					} else {
						// internal loop
						Size cutpoint(0);
						if ( next_segment_is_flexible && previous_segment_is_flexible ) {
							cutpoint = ( use_randomness ? (i-1 + int( numeric::random::uniform()* ( seglen+1 ) ) ) : ( i-1+j)/2 );
						} else if ( next_segment_is_flexible ) {
							cutpoint = i-1;
						} else if ( previous_segment_is_flexible ) {
							cutpoint = j;
						}
						if ( cutpoint ) {
							delete_me = true;
							pbassert( cutpoint>=i-1 && cutpoint <= j );
							segments[ii-1].back() = cutpoint;
							segments[ii+1].front() = cutpoint+1;
						}
					}
					if ( delete_me ) {
						segments.erase( segments.begin()+(ii-1) );
						--ii;
					}
				}

				if ( isloop && seglen > max_loop_segment_size ) {
					// split to two loops; make the first one target_loop_segment_size
					TR_REBUILD.Trace   << "Long loop: " << seglen << endl;
					Sizes const new_seg( make_vector1( i, i+target_loop_segment_size-1 ) );
					segments.insert( segments.begin()+(ii-1), new_seg );
					pbassert( segments[ii].front() == i  && segments[ii+1].front() == i );
					segments[ii+1].front() = i+target_loop_segment_size;
				}
			} //ii=1,segments.size()

			// status output
			for ( Size ii=1; ii<= segments.size(); ++ii ) {
				Size const i( segments[ii].front() ), j( segments[ii].back() ), seglen( j-i+1 );
				TR_REBUILD.Trace << "SEG " << I(4,cb) << I(4,i) << I(4,j) << ss.substr(i,seglen) << endl;
				protein_segments.push_back( make_vector1( cb + i-1, cb + j-1 ) );
			}
		}
	}
	if ( !is_frozen.empty() ) {
		// split segments that contain frozen parts
		vector1< Sizes > const old_protein_segments( protein_segments );
		protein_segments.clear();
		foreach_ ( Sizes const & seg, old_protein_segments ) {
			Size const segbegin( seg[1] ), segend( seg[2] ), seglen( segend-segbegin+1 );
			string const ss( pose.secstruct().substr(segbegin-1,seglen));
			string fs;
			for ( Size i=segbegin; i<= segend; ++i ) {
				if ( is_frozen[i] ) fs.push_back('1');
				else fs.push_back('0');
			}
			for ( Size i=segbegin; i<= segend; ++i ) {
				// i is the start of a segment
				Size j(i);
				while ( j<segend && is_frozen[i] == is_frozen[j+1] ) ++j;
				protein_segments.push_back( make_vector1( i,j ) );
				if ( i>segbegin ||j<segend ) {
					TR_REBUILD.Trace << "SPLITSEG " << ss << ' ' << fs <<
						" old: " << I(4,segbegin) << I(4,segend) << " new: " << I(4,i) << I(4,j) <<
						" is_frozen: " << is_frozen[i] << endl;
				}
				i = j;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
setup_dna_segments(
	Pose const & pose,
	SizePairs & dna_segments
	//bool const use_randomness = true
)
{
	scoring::dna::BasePartner const & partner( scoring::dna::retrieve_base_partner_from_pose( pose ) );
	Size const max_dna_seglen( 4 );
	// first identify continuous segments that are either perfect duplex or not

	bools in_segment( pose.total_residue(), false);

	//bool inseg( false );
	dna_segments.clear();
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( !pose.residue(i).is_DNA() || in_segment[i] ) continue;

		// starting a new segment with i
		Size j(i), ce( chain_end( pose.chain(i), pose  ) );
		if ( partner[i] ) {
			// look ahead to end of segment
			while ( j< ce && !pose.residue(j).is_upper_terminus() && partner[j+1] &&
					( partner[i] - partner[j+1] == j+1-i ) ) ++j;
		} else { // unpaired segment
			while ( j< ce && !pose.residue(j).is_upper_terminus() && !partner[j+1] ) ++j;
		}
		dna_segments.push_back( make_pair(i,j) );
		for ( Size ii=i; ii<=j; ++ii ) {
			in_segment[ii] = true;
			if ( partner[ii] ) in_segment[partner[ii]] = true;
		}
		i = j+1;
	}

	// check that every dna residue is in a segment
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		TR_REBUILD.Trace << "inseg: " << in_segment[i] << ' ' <<  pose.residue(i).is_DNA() << endl;
		runtime_assert( in_segment[i] == pose.residue(i).is_DNA() );
	}

	/// now go through and break the segments up into pieces
	for ( Size i=1, i_end = dna_segments.size(); i<= i_end; ++i ) {
		SizePair & p( dna_segments[i] );
		Size const plen( p.second - p.first + 1 );
		Size ncuts(0);
		while ( Real(plen) / (ncuts+1) > max_dna_seglen ) ++ncuts;
		if ( ncuts ) {
			TR_REBUILD.Trace << "break segment " << p.first << ' ' << p.second << ' ' << plen << endl;
			Real const cutsize( Real(plen)/(ncuts+1) );
			SizePair const oldp(p);
			for ( Size i=1; i<= ncuts+1; ++i ) {
				SizePair newp;
				if ( i==1 ) newp.first = oldp.first;
				else newp.first  = Size( floor( 0.5 + ( oldp.first-1 ) + (i-1)*cutsize ) ) + 1; // +1
				if ( i==ncuts+1 ) newp.second = oldp.second;
				else newp.second = Size( floor( 0.5 + ( oldp.first-1 ) +     i*cutsize ) );
				TR_REBUILD.Trace << "newp " << newp.first << ' ' << newp.second << endl;
				if ( i== 1 ) p = newp;
				else dna_segments.push_back( newp );
			}
		}
	}


	foreach_ ( SizePair p, dna_segments ) {
		// sanity
		Size const pb( p.first), pe( p.second );
		for ( Size i=pb; i<= pe; ++i ) {
			runtime_assert( pose.residue(i).is_DNA() );
			if ( partner[pb] ) {
				runtime_assert( partner[i] && partner[pb] - partner[i] == i-pb );
			} else {
				runtime_assert( partner[i] == 0 );
			}
		}
		TR_REBUILD.Trace << "dnaseg: " << pb << ":" << partner[pb] << " --> " << pe << ":" << partner[pe] << endl;
	}

}


//using devel::blab::append_virtual_residue; // in pose_tools.cc
// void
// append_virtual_residue( Pose & pose )
// {
//  if ( pose.residue( pose.total_residue() ).name() == "VRT" ) return; // already done

//  Vector xyz_orig;
//  {
//   /// get the dnachains
//   Sizes dnachains;
//   for ( Size i=1; i<= pose.num_chains(); ++i ) {
//    if ( pose.residue(pose.chain_begin(i)).is_DNA()) dnachains.push_back(i );
//   }
//   if ( dnachains.empty() ) {
//    // no DNA
//    xyz_orig = pose.residue( pose.total_residue()/2 ).nbr_atom_xyz();
//   } else {
//    //pbassert( dnachains.size() == 2 && dnachains[1]+1 == dnachains[2] ); // could relax?
//    Size const dna_root( ( pose.chain_begin( dnachains.front() ) + pose.chain_end( dnachains.front() ) )/2 );
//    xyz_orig = pose.residue( dna_root ).xyz("P");
//   }
//  }

//  ResidueOP vrtrsd
//   ( conformation::ResidueFactory::create_residue( pose.residue(1).residue_type_set().name_map( "VRT" ) ) );
//  translate_residue( *vrtrsd, xyz_orig );
//  pose.append_residue_by_jump( *vrtrsd, pose.total_residue()/2 );
//  pose.conformation().insert_chain_ending( pose.total_residue()-1 );

// }


/// this assumes that the segment is base-paired
///
void
add_dna_segment_base_jumps_to_star_fold_tree(
	Size segbegin,
	Size segend,
	Pose const & pose,
	kinematics::FoldTree & f,
	bool const stochastic_flip = true
)
{
	scoring::dna::BasePartner const & partner( scoring::dna::retrieve_base_partner_from_pose( pose ) );

	runtime_assert( pose.chain( segbegin ) == pose.chain( segend ) );
	for ( Size i=segbegin; i<= segend; ++i ) {
		runtime_assert( partner[ i ] && ( partner[ segbegin ] - partner[i] == i - segbegin ) );
	}

	if ( partner[segbegin] < segbegin ) {
		Size const tmp( segbegin );
		segbegin = partner[segend];
		segend = partner[tmp];
	}

	runtime_assert( partner[segbegin] > segbegin );

	Size const vrtpos( pose.total_residue() );
	runtime_assert( pose.residue(vrtpos).name() == "VRT" );


	Size const npairs( segend - segbegin + 1 );

	bool const flip( stochastic_flip ? ( numeric::random::uniform() < 0.5 ) : false );

	for ( Size i=segbegin; i<= segend; ++i ) {
		Size const ii( i - segbegin + 1 ); // basepair index
		Size const p_i( partner[i] );

		// should the jump between base-steps be on strand I (i-strand) or strand II
		bool const base_step_jump_on_strand_I( flip ? ( i%2 == 0 ) : ( i%2 == 1 ) );

		// the basepair jump
		int basepair_cut( base_step_jump_on_strand_I ? p_i - 1 : i );

		if ( ii==npairs ) basepair_cut = segend;

		f.new_jump( i, p_i, basepair_cut );

		// the basestep jump
		if ( ii < npairs ) {
			if ( base_step_jump_on_strand_I ) {
				f.new_jump( i, i+1, i );
			} else {
				f.new_jump( p_i - 1, p_i, p_i - 1 );
			}
		}
	}

	// connection from the segment to the vrtpos
	Size anchor( (segbegin+segend)/2 );
	if ( stochastic_flip && ( numeric::random::uniform() < 0.5 ) ) anchor = partner[ anchor ];
	// Sizes vrtcuts( make_vector1( segbegin-1, segend, partner[ segbegin], partner[segend]-1 ) );
	// std::sort( vrtcuts.begin(), vrtcuts.end() );
	// std::reverse( vrtcuts.begin(), vrtcuts.end() ); // prefer later cuts
	// Size vrtcut(0);
	// foreach_( Size cut, vrtcuts ) {
	//  if ( cut >=1 && cut < pose.total_residue() && !f.is_cutpoint(cut) ) {
	//   vrtcut = cut;
	//   break;
	//  }
	// }
	// runtime_assert( vrtcut );
	Size const vrtcut( partner[segbegin ] );
	TR_REBUILD.Trace << "vrtcut " << segbegin << ' ' << segend << ' ' << vrtcut << endl;
	runtime_assert( !f.is_cutpoint(vrtcut) );
	f.new_jump( anchor, vrtpos, vrtcut );

}







//// uses MotifData::DNA_MOTIF to set dna root
void
setup_star_fold_tree_cutpoint_variants_and_virtual_residue(
	bool const anchor_dna_jumps_in_backbone,
	vector1< Sizes > const & protein_segments,
	Pose & pose,
	bool const freeze_dna = false
)
{
	/// get the dnachains
	Sizes dnachains;
	for ( Size i=1; i<= num_chains( pose ); ++i ) {
		if ( pose.residue(chain_begin(pose,i)).is_DNA() ) dnachains.push_back( i );
	}
	bool const pose_has_dna( dnachains.size() );
	Size dna_root(0);
	if ( pose_has_dna && !freeze_dna ) { // figure out where to root the dna
		pbassert( dnachains.size() == 2 && dnachains[1]+1 == dnachains[2] ); // could relax?

		if ( option[ my_options::dna_root_shift_window ].user() ) {
			/// shift dna_root randomly within a window of size "window" and on either strand
			Size const chain1_begin( chain_begin(pose, dnachains.front() ) ), chain1_end( chain_end(pose, dnachains.front() )),
				chain1_midpoint( (chain1_end + chain1_begin)/2 ), chain1_len( chain1_end - chain1_begin + 1 ),
				window( min( chain1_len, Size(option[ my_options::dna_root_shift_window ] ) ) ),
				window_midpoint( (1 + window)/2 );
			/// align chain1_midpoint with window_midpoint
			Size const window_start( chain1_midpoint - window_midpoint+1 );
			for ( Size i=window_start; i< window_start+window; ++i ) {pbassert( retrieve_base_partner_from_pose( pose )[i] );}
			dna_root = window_start + int( numeric::random::uniform() * window );
			if ( numeric::random::uniform()<0.5 ) dna_root = retrieve_base_partner_from_pose( pose )[ dna_root ];

		} else {
			dna_root = ( chain_begin(pose, dnachains.front() ) + chain_end(pose, dnachains.front() ) )/2;
			dna_root = int( dna_root ) + option[ my_options::dna_root_shift_hack ];
			if ( option[ my_options::dna_root_strand_flip_hack ] ) { /// THIS IS A SILLY HACK !!!
				dna_root = retrieve_base_partner_from_pose( pose )[ dna_root ];
			}
		}
		pbassert( has_element( dnachains, Size(pose.chain( dna_root ) ) ) );
	}

	// find the virtual residue, if it exists
	Size const vrtpos( get_vrt_pos( pose ) );
	pbassert( vrtpos == pose.total_residue() ); /// assuming we've got no waters, just protein and dna
	pbassert( chain_begin(pose, num_chains( pose ) ) == chain_end(pose, num_chains( pose ) ) ); // last chain: vrtpos


	kinematics::FoldTree f( pose.total_residue() );

	Size const root( vrtpos );

	// jump from root to dna
	if ( pose_has_dna ) {
		if ( freeze_dna ) { // each dna chain individually
			for ( Sizes::const_iterator ch = dnachains.begin(); ch != dnachains.end(); ++ch ) {
				Size const ce( chain_end(pose, *ch ) );
				TR_REBUILD.Trace << "setup_star_fold_tree_cutpoint_variants_and_virtual_residue: frozen_dna_chain_jump: ce= " << ce << " root= " << root << endl;
				f.new_jump( ce, root, ce );
			}
		} else {
			f.new_jump( dna_root, root, chain_end(pose, dnachains.back() ) );
		}
	}


	// jump from root to protein
	for ( Size ii=1; ii<= protein_segments.size(); ++ii ) {
		Size const i( protein_segments[ii].front() ), j( protein_segments[ii].back() );
		TR_REBUILD.Trace << "setup_star_fold_tree_cutpoint_variants_and_virtual_residue: protein_segment_jump: " << i << ' ' << j << ' ' << (i+j)/2 << endl;
		f.new_jump( root, (i+j)/2, j );
	}

	// now the intra-base jumps
	bool const basetreeflip( numeric::random::uniform() < 0.5 );
	if ( pose_has_dna && !freeze_dna ) devel::dna::add_dna_base_jumps_to_fold_tree( pose, f, basetreeflip ); // STOCHASTIC


	f.reorder( root );

	// now set the dna atoms
	if ( pose_has_dna && !freeze_dna ) devel::dna::set_dna_jump_atoms_in_fold_tree( pose, anchor_dna_jumps_in_backbone, f);

	// set pose fold_tree equal to f
	pose.fold_tree( f );

	// now setup protein cutpoint variants
	for ( Size ii=1; ii<= protein_segments.size(); ++ii ) {
		Size const j( protein_segments[ii].back() );
		if ( !( pose.residue(j).is_upper_terminus() ||
				j == pose.total_residue() ||
				pose.residue(j+1).is_lower_terminus() ) ) {
			if ( !pose.residue(j).has_variant_type( CUTPOINT_LOWER ) ) {
				pbassert( !pose.residue(j+1).has_variant_type( CUTPOINT_UPPER ) );
				add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, j   );
				add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, j+1 );
			}
		}
	}

	// sanity check:
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).is_protein() && pose.residue(i).has_variant_type( CUTPOINT_LOWER ) ) {
			pbassert( f.is_cutpoint(i) );
		}
	}

	// now set cutpoint variants between all bonded DNA residues
	if ( pose_has_dna && !freeze_dna ) devel::dna::setup_dna_cutpoint_variants( pose );


	/// confirm correctness of stub atoms
	/// will this work?
	if ( pose_has_dna && !freeze_dna && !anchor_dna_jumps_in_backbone ) {
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const & rsd( pose.residue(i) );
			if ( rsd.is_DNA() ) {
				// just check that the correct tree has been set:
				kinematics::tree::Atom const & jatom( pose.atom_tree().atom( AtomID( rsd.chi_atoms(1)[4], i ) ) );
				if ( jatom.stub_atom2_id() != AtomID( rsd.chi_atoms(1)[3], i ) ||
						jatom.stub_atom3_id() != AtomID( rsd.chi_atoms(1)[2], i ) ||
						!jatom.is_jump() ) {
					utility_exit_with_message( "bad atom_tree setup in set_dna_jump_atoms!" );
				}
			}
		}
	}


}

/// this is a little hacky right now
//// uses MotifData::DNA_MOTIF to set dna root
void
setup_star_fold_tree_cutpoint_variants_and_virtual_residue(
	bool const anchor_dna_jumps_in_backbone,
	vector1< Sizes > const & protein_segments,
	SizePairs const & dna_segments,
	Pose & pose
)
{


	// find the virtual residue, if it exists
	Size const vrtpos( get_vrt_pos( pose ) );
	pbassert( vrtpos == pose.total_residue() ); /// assuming we've got no waters, just protein and dna
	pbassert( chain_begin(pose, num_chains( pose ) ) == chain_end(pose, num_chains( pose ) ) ); // last chain: vrtpos


	kinematics::FoldTree f( pose.total_residue() );

	Size const root( vrtpos );


	// jump from root to protein
	for ( Size ii=1; ii<= protein_segments.size(); ++ii ) {
		Size const i( protein_segments[ii].front() ), j( protein_segments[ii].back() );
		f.new_jump( root, (i+j)/2, j );
	}

	// intra-dna zig-zag jumps and jump to vrtpos
	foreach_ ( SizePair const & seg, dna_segments ) {
		add_dna_segment_base_jumps_to_star_fold_tree( seg.first, seg.second, pose, f );
	}

	f.reorder( root );

	// now set the dna atoms
	if ( dna_segments.size() ) devel::dna::set_dna_jump_atoms_in_fold_tree( pose, anchor_dna_jumps_in_backbone, f);

	// set pose fold_tree equal to f
	pose.fold_tree( f );

	// now setup protein cutpoint variants
	for ( Size ii=1; ii<= protein_segments.size(); ++ii ) {
		Size const j( protein_segments[ii].back() );
		if ( !pose.residue(j).is_upper_terminus() ) {
			if ( !pose.residue(j).has_variant_type( CUTPOINT_LOWER ) ) {
				pbassert( !pose.residue(j+1).has_variant_type( CUTPOINT_UPPER ) );
				add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, j   );
				add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, j+1 );
			}
		}
	}

	// sanity check:
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).is_protein() && pose.residue(i).has_variant_type( CUTPOINT_LOWER ) ) {
			pbassert( f.is_cutpoint(i) );
		}
	}

	// now set cutpoint variants between all bonded DNA residues
	if ( dna_segments.size() ) devel::dna::setup_dna_cutpoint_variants( pose );


	/// confirm correctness of stub atoms
	/// will this work?
	if ( dna_segments.size() && !anchor_dna_jumps_in_backbone ) {
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const & rsd( pose.residue(i) );
			if ( rsd.is_DNA() ) {
				// just check that the correct tree has been set:
				kinematics::tree::Atom const & jatom( pose.atom_tree().atom( AtomID( rsd.chi_atoms(1)[4], i ) ) );
				if ( jatom.stub_atom2_id() != AtomID( rsd.chi_atoms(1)[3], i ) ||
						jatom.stub_atom3_id() != AtomID( rsd.chi_atoms(1)[2], i ) ||
						!jatom.is_jump() ) {
					utility_exit_with_message( "bad atom_tree setup in set_dna_jump_atoms!" );
				}
			}
		}
	}


}

/// this is a little hacky right now
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
superimpose_segments(
	Vectors const & target_coords, // CA coords of target pose (only works for protein right now)
	SizePairs const & segments,
	Sizes const & segment_jumps,
	Pose & pose
)
{
	conformation::symmetry::SymmetryInfoCOP symminfo(0);
	if ( pose::symmetry::is_symmetric( pose ) ) symminfo = pose::symmetry::symmetry_info( pose );
	for ( Size seg=1; seg<= segments.size(); ++seg ) {
		Size const segbegin( segments[seg].first ), segend( segments[seg].second );
		if ( symminfo && !symminfo->jump_is_independent( segment_jumps[seg] ) ) continue;
		Vectors movcoords, fixcoords;
		for ( Size i=segbegin; i<= segend; ++i ) {
			movcoords.push_back( pose.residue(i).xyz("CA"));
			fixcoords.push_back( target_coords[i] );
		}
		if ( movcoords.size()<3 ) {
			TR_REBUILD.Trace << "superimpose_segments: short segment? " << movcoords.size() << ' ' << segbegin << " --- " <<
				segend << endl;
			continue;
		}
		runtime_assert( movcoords.size()>=3 );
		Stub const oldstub( movcoords[1], movcoords[2], movcoords[3] );
		superimpose_coords( fixcoords, movcoords );
		Stub const newstub( movcoords[1], movcoords[2], movcoords[3] );

		Stub const   upstub( pose.conformation().  upstream_jump_stub( segment_jumps[seg] ) );
		Stub const downstub( pose.conformation().downstream_jump_stub( segment_jumps[seg] ) );
		// we want to set a new jump so that downstub got moved by the transform that mapped oldstub to newstub
		Matrix const R( newstub.M * oldstub.M.transposed() );
		Vector const v( newstub.v - R*oldstub.v );
		Stub const newdownstub( R * downstub.M, R * downstub.v + v );
		pose.set_jump( segment_jumps[seg], kinematics::Jump( upstub, newdownstub ) );

		if ( TR_REBUILD.Trace.visible() ) {
			Real dev(0);
			for ( Size i=segbegin; i<= segend; ++i ) {
				dev += movcoords[i-segbegin+1].distance_squared( pose.residue(i).xyz("CA") );
			}
			dev = sqrt( dev / movcoords.size() );
			TR_REBUILD.Trace << "superimpose_segments: " << I(3,seg) << I(4,segbegin) << I(4,segend) <<
				I(4,segment_jumps[seg]) << F(9,3,dev) << endl;
		}
	}

}


void
parse_target_coords_and_segments(
	Pose const & pose,
	bools const & is_flexible_protein,
	Vectors & target_coords,
	SizePairs & segments,
	Sizes & segment_jumps
)

{ // analysis of constraints set
	Vector const dummy_xyz( 1, 2, 3 );
	target_coords.clear(); target_coords.resize( pose.total_residue(), dummy_xyz );
	segment_jumps.clear();
	segments.clear();

	using namespace scoring::constraints;
	//ConstraintSetCOP cst_set( pose.constraint_set() );

	ConstraintCOPs all_csts( pose.constraint_set()->get_all_constraints() );
	for ( ConstraintCOPs::const_iterator cst= all_csts.begin(); cst != all_csts.end(); ++cst ) {
		ostringstream out;
		(*cst)->show_def( out, pose );
		string const outstring( out.str() );
		TR_REBUILD.Trace << "cstdef: " << outstring << endl;
		strings const l( split_to_vector1( out.str() ) );
		if ( l[1] == "CoordinateConstraint" ) {
			// runtime_assert( l[2] == "CA" ); // not necessary, actually
			runtime_assert( is_int( l[3] ) );
			Size const pos( int_of( l[3] ) );
			runtime_assert( l[4] == "ORIG" );
			runtime_assert( is_int( l[5] ) );
			Size const vrtpos( int_of( l[5] ) );
			if ( pose.residue(pos).is_protein() ) {
				runtime_assert( !pose.residue(vrtpos).is_protein() );
				runtime_assert( is_float( l[6] ) );
				runtime_assert( is_float( l[7] ) );
				runtime_assert( is_float( l[8] ) );
				target_coords[ pos ] = Vector( float_of( l[6] ), float_of( l[7] ), float_of( l[8] ) );
			}
		}
	}
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).is_protein() ) { runtime_assert( target_coords[i].distance( dummy_xyz ) > 1e-3 );}
	}

	Size seg_begin(0);
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( !is_flexible_protein[i] ) continue;
		TR_REBUILD.Trace << i << ' ' << pose.residue(i).name() << endl;
		if ( !pose.residue(i).is_protein() ) {
			runtime_assert( !seg_begin );
			continue;
		}
		if ( pose.residue(i).has_variant_type( CUTPOINT_UPPER ) || pose.residue(i).is_lower_terminus() ) {
			runtime_assert( i==1 || pose.fold_tree().is_cutpoint(i-1) );
			seg_begin = i;
		}
		if ( pose.residue(i).has_variant_type( CUTPOINT_LOWER ) || pose.residue(i).is_upper_terminus() ) {
			runtime_assert( i==pose.total_residue() || pose.fold_tree().is_cutpoint(i) );
			runtime_assert( seg_begin );
			segments.push_back( make_pair( seg_begin, i ) );
			Size root(0);
			for ( Size j=seg_begin; j<=i; ++j ) {
				if ( pose.fold_tree().is_jump_point(j) ) {
					runtime_assert( !root );
					root = j;
				}
			}
			runtime_assert( root );
			segment_jumps.push_back( pose.fold_tree().get_jump_that_builds_residue( root ) );
			seg_begin = 0;
		} else {
			runtime_assert( seg_begin );
		}
	}
}

void
star_fragment_rebuild(
	ScoreFunction const & scorefxn_in,
	bools is_flexible_dna, // ensure that it includes partner positions
	bools const & is_flexible_protein,
	devel::blab::classic_frags::FragLibCOP protein_frag_lib,
	Pose & pose,
	Size const nouter = 10,
	Size const nouter_9mers_only = 3, // 5-bp frags <= this
	Size const nouter_3mers_only = 8, // 2-bp frags >= this (3 bp in middle)
	bool const use_superimpose_segments = false,
	Size const big_frag_size = 9,
	Sizes dna_frag_sizes = make_vector1( 2,3,5 ),
	string mintype = string()
)
{
	runtime_assert( nouter_9mers_only <= nouter_3mers_only );
	runtime_assert( nouter_9mers_only <= nouter );
	runtime_assert( nouter_3mers_only <= nouter );

	std::sort( dna_frag_sizes.begin(), dna_frag_sizes.end() );
	runtime_assert( dna_frag_sizes.size() == 3 );

	if ( mintype.empty() ) {
		if ( option[ my_options::star_fragment_rebuild_mintype ].user() ) {
			mintype = option[ my_options::star_fragment_rebuild_mintype ]();
		} else {
			mintype = "dfpmin_armijo_nonmonotone"; // new default, old was linmin
		}
	}

	protocols::moves::PyMOLMoverOP pymol(0);
	if ( option[ my_options::show_pymol ] ) {
		static Size counter(0);
		++counter;
		pymol = protocols::moves::PyMOLMoverOP( new protocols::moves::PyMOLMover());
		pymol->pymol_name("star_fragment_rebuild_"+lead_zero_string_of(counter,2));
		pymol->keep_history(true);
	}

	//PROF_START( util::CENTROID_SIMULATION );
	using namespace  protocols::moves;
	using namespace protocols::minimization_packing;
	using namespace devel::blab::classic_frags;
	using namespace devel::dna;
	using numeric::random::uniform;
	using ObjexxFCL::format::I;
	using ObjexxFCL::format::F;
	//bool const debug( false );

	ScoreFunctionOP scorefxn( scorefxn_in.clone() );

	bool pose_has_flexible_dna( false ), pose_has_flexible_protein( false );
	Size const vrtpos( get_vrt_pos( pose ) ); // will be 0 if there's not VRT residue
	runtime_assert( vrtpos == pose.total_residue() ); // should be, right?
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).is_DNA() && is_flexible_dna[i] ) {
			pose_has_flexible_dna = true;
			Size const ip( scoring::dna::retrieve_base_partner_from_pose( pose )[i] );
			if ( ip ) is_flexible_dna[ip] = true;
		}
		if ( pose.residue(i).is_protein() && is_flexible_protein[i] ) pose_has_flexible_protein = true;
	}

	// setup the DNA fragments
	devel::dna::DoubleFragmentLibraryOP dna_frag_lib(0);
	if ( pose_has_flexible_dna ) {
		Size const frag_nn( 30 );
		//Sizes frag_sizes( make_vector1( 2, 3, 5 ) );
		dna_frag_lib = devel::dna::setup_dna_frags( pose, is_flexible_dna, frag_nn, dna_frag_sizes );
	}

	Vectors target_coords;
	SizePairs segments;
	Sizes segment_jumps;
	if ( use_superimpose_segments ) {
		parse_target_coords_and_segments( pose, is_flexible_protein, target_coords, segments, segment_jumps );
	}

	/// monte carlo setup
	Real const temp( 2.0 );
	MonteCarloOP mc( new MonteCarlo( pose, *scorefxn, temp ) );

	/// setup the move_map
	kinematics::MoveMapOP mm_star_jumps( new kinematics::MoveMap ), mm_full( new kinematics::MoveMap );
	{

		bools is_chi_flexible( pose.total_residue(), false ), is_bb_flexible( pose.total_residue(), false ),
			is_jump_flexible( pose.total_residue(), false );
		for ( Size i=1; i< pose.total_residue(); ++i ) {
			if ( pose.residue(i).is_DNA() ) {
				bool tmpbool;
				is_jump_flexible[i] = ( is_flexible_dna[i] && pose.fold_tree().get_parent_residue( i, tmpbool ) == vrtpos );
				if ( is_jump_flexible[i] ) {
					runtime_assert( tmpbool );
					TR_REBUILD.Trace << "star_fragment_rebuild:: mm_star_jumps-dna-jump-flex: " << i << endl;
				}
			} else if ( pose.residue(i).is_protein() ) {
				is_jump_flexible[i] = is_flexible_protein[i];
			}
		}

		bool const vary_omega( false ), flex_dna_sugar( false ), set_jump_flex_by_downstream_only( true );
		devel::blab::setup_move_map( is_chi_flexible, is_bb_flexible, is_jump_flexible, vary_omega, flex_dna_sugar, pose,
			*mm_star_jumps, set_jump_flex_by_downstream_only );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( pose.residue(i).is_DNA() && is_flexible_dna[i] ) {
				is_bb_flexible[i] = is_chi_flexible[i] = is_jump_flexible[i] = true;
			} else if ( pose.residue(i).is_protein() && is_flexible_protein[i] ) {
				is_bb_flexible[i] = is_chi_flexible[i] = true;
			}
		}
		devel::blab::setup_move_map( is_chi_flexible, is_bb_flexible, is_jump_flexible, vary_omega, flex_dna_sugar, pose,
			*mm_full, set_jump_flex_by_downstream_only );
	}




	{ /// randomize the starting position
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( pose.residue(i).is_protein() && is_flexible_protein[i] ) {
				conformation::idealize_position( i, pose.conformation() );
			}
		}
		if ( pose_has_flexible_dna ) {
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				if ( pose.residue(i).is_DNA() && is_flexible_dna[i] ) {
					// this is new:
					// it's not so bad if we assume that there are zig-zag jumps anchored in the backbone
					// so the bases will move around a bit but other stuff should be OK...
					// if we were folding through the bb then there might be trouble...
					conformation::idealize_position( i, pose.conformation() );
				}
			}
			insert_random_fragments_in_flexible_dna_regions( sizes_from_bools( is_flexible_dna ), *dna_frag_lib, pose );
		}
		if ( pose_has_flexible_protein ) {
			insert_random_fragments_in_flexible_protein_regions( sizes_from_bools( is_flexible_protein ), *protein_frag_lib,
				pose );
		}
	}
	mc->reset( pose );

	/// setup the Movers
	TorsionFragmentMoverOP torsion_frag_mover( 0 );
	if ( pose_has_flexible_protein ) {
		torsion_frag_mover =
			TorsionFragmentMoverOP( new TorsionFragmentMover( protein_frag_lib, mm_full ) );
	}


	// min mover
	Real const mintol( 0.1 );
	MinMoverOP min_mover( new MinMover( mm_star_jumps, scorefxn, mintype, mintol, true ) );

	/// setup trials
	// MoverOP torsion_frag_trial(0), torsion_frag_min_trial(0);
	// if ( pose_has_flexible_protein ) {
	//  torsion_frag_trial = new TrialMover( torsion_frag_mover, mc );
	//  torsion_frag_min_trial = new TrialMover( new SequenceMover( torsion_frag_mover, min_mover ), mc );
	// }
	TrialMoverOP min_trial( new TrialMover( min_mover, mc ) );

	MoverOP double_frag_trial(0), double_frag_min_trial(0);
	DoubleFragmentMoverOP double_frag_mover( 0 );
	if ( pose_has_flexible_dna ) {
		double_frag_mover = DoubleFragmentMoverOP( new DoubleFragmentMover( dna_frag_lib, mm_full, true ) );
		double_frag_trial = MoverOP( new TrialMover(  double_frag_mover, mc ) );
		double_frag_min_trial = MoverOP( new TrialMover( MoverOP( new SequenceMover(  double_frag_mover, min_mover ) ), mc ));
	}

	// { // hacking
	//  double_frag_mover->frag_size( 3 );
	//  pose.dump_pdb("pre_idl.pdb");

	//  for ( Size i=1; i<= pose.total_residue(); ++i ) {
	//   if ( pose.residue(i).is_DNA() && is_flexible_dna[i] ) {
	//    conformation::idealize_position( i, pose.conformation() );
	//   }
	//  }
	//  pose.dump_pdb("post_idl.pdb");

	//  Pose const start_pose( pose );

	//  for ( Size n=1; n<= 10; ++n ) {
	//   pose = start_pose;
	//   cout << "test " << n << endl;
	//   double_frag_mover->apply( pose );
	//   pose.dump_pdb("test"+lead_zero_string_of(n,4)+".pdb");
	//  }
	//  exit(0);
	// }


	bool const min_moves( true ); // try something different
	Size ninner( dry_run() ? 2 : 150 ); // zf uses 250/finger, but no minmoves until the end...
	// Size const nouter( 10 );
	// Size const nouter_9mers_only( 3 );
	// Size const nouter_3mers_only( 8 ); // should give: 9 9 9 3 9 3 9 3 3 3

	Real const //min_chainbreak_dna_weight( 0.5 ), max_chainbreak_dna_weight( 1.0 ),
		min_chainbreak_protein_weight( 0.1 ), max_chainbreak_protein_weight( 1.0 ); // changing min from 0 to 0.1, now w/ dna

	//scorefxn->set_weight( chainbreak, 0.0 ); // set them independently

	for ( Size n=1; n<= nouter; ++n ) {
		Real const scale( Real(n-1) / (nouter-1 ) ); // goes from 0 to 1
		//   Real const chainbreak_dna_weight( min_chainbreak_dna_weight +
		//                    scale * ( max_chainbreak_dna_weight - min_chainbreak_protein_weight ) );
		Real const chainbreak_protein_weight( min_chainbreak_protein_weight +
			scale * ( max_chainbreak_protein_weight - min_chainbreak_protein_weight ) );

		//scorefxn->set_weight( chainbreak_dna, chainbreak_dna_weight );
		scorefxn->set_weight( chainbreak, chainbreak_protein_weight );

		mc->score_function( *scorefxn );
		mc->reset( pose );

		// now set everything flexible
		if ( n == nouter ) min_mover->movemap( mm_full );



		// choose fragment size
		Size frag_size, dna_frag_size;
		//bool min_moves( n == nouter );
		if ( n <= nouter_9mers_only ) { // early
			frag_size = big_frag_size;
			dna_frag_size = dna_frag_sizes[3];
		} else if ( n >= nouter_3mers_only ) { // late
			frag_size = 3;
			dna_frag_size = dna_frag_sizes[1];
		} else { // middle
			frag_size = ( ( n - nouter_9mers_only )%2 == 0 ? big_frag_size : 3 );
			dna_frag_size = dna_frag_sizes[2];
		}
		if ( pose_has_flexible_protein ) torsion_frag_mover->frag_size( frag_size );
		if ( pose_has_flexible_dna ) double_frag_mover->frag_size( dna_frag_size );

		if ( min_moves ) min_trial->apply( pose );

		// now inner cycles
		for ( Size m=1; m<= ninner; ++m ) {
			if ( !min_moves ) {
				if ( pose_has_flexible_protein ) { // cant use trial movers since we want to use superimpose_segments
					torsion_frag_mover->apply( pose );
					if ( use_superimpose_segments ) superimpose_segments( target_coords, segments, segment_jumps, pose );
					mc->boltzmann( pose );
					if ( pymol && mc->mc_accepted() == protocols::moves::MCA_accepted_score_beat_low ) pymol->apply( pose );
				}
				if ( pose_has_flexible_dna ) double_frag_trial->apply( pose );
			} else {
				if ( pose_has_flexible_protein ) { // cant use trial movers since we want to use superimpose_segments
					torsion_frag_mover->apply( pose );
					if ( use_superimpose_segments ) superimpose_segments( target_coords, segments, segment_jumps, pose );
					min_mover->apply( pose );
					mc->boltzmann( pose );
					if ( pymol && mc->mc_accepted() == protocols::moves::MCA_accepted_score_beat_low ) pymol->apply( pose );
				}
				if ( pose_has_flexible_dna ) double_frag_min_trial->apply( pose );
			}
			TR_REBUILD.Trace << "outer " << I(3,n) << " inner " << I(4,m) << " min_moves: " << min_moves <<
				" chnbrk: " << F(9,3,mc->last_accepted_pose().energies().total_energies()[ chainbreak] ) <<
				" crdcst: " << F(9,3,mc->last_accepted_pose().energies().total_energies()[ coordinate_constraint ] ) <<
				" last: " << F(9,3,mc->last_accepted_score()) << " low: " << F(9,3,mc->lowest_score()) << endl;
		}
		mc->show_counters();
		mc->recover_low( pose );
	}
	if ( pymol ) pymol->apply( pose ); // final pose

	//basic::Tracer status( "blab_status.motif_refine_loops" ); mc->show_counters( status );
	//PROF_STOP( util::CENTROID_SIMULATION );
}



/// Right now:
/// include_current = TRUE
/// bump_check = FALSE

void
envdep_pack_rotamers_wrapper(
	kinematics::MoveMap const & mm,
	ScoreFunction const & env_indep_scorefxn,
	ScoreFunction const &, // env_dep_scorefxn,
	bool const use_envdep_packing,
	Pose & pose,
	bool const allow_new_waters = true
)
{
	runtime_assert( !use_envdep_packing );
	// runtime_assert( !pack::score_function_has_local_context_dependent_methods( env_indep_scorefxn ) );

	Size const pack_rotamers_nloop( 25 );
	/// if linmem or other on-the-fly igs/annealers are used, loop takes forever!
	runtime_assert( ! basic::options::option[ basic::options::OptionKeys::packing::linmem_ig ].user() );

	/// note -- using bump check?
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ) );
	{ // setup the task
		using namespace pack::rotamer_set;
		//bool const pose_has_water( has_water_packing_info( pose ));
		task->initialize_from_command_line(); //.restrict_to_repacking().restrict_to_residues(allow_repack);
		task->or_include_current( true );
		task->set_bump_check( false );
		protocols::dna::SetupFlexibilityIncludingWaters( pose.total_residue(), mm, allow_new_waters ).apply( pose, *task );
		// bools is_flexible_water( pose.total_residue(), false );
		// for ( Size i=1; i<= pose.total_residue(); ++i ) {
		//  if ( !mm.get_chi( i ) ) {
		//   task->nonconst_residue_task(i).prevent_repacking();
		//  } else {
		//   if ( pose_has_water && get_water_packing_info( pose ).has(i) ) { // flexible water
		//    is_flexible_water[ i ] = true;
		//   } else {
		//    task->nonconst_residue_task(i).restrict_to_repacking();
		//   }
		//  }
		// }
		// if ( pose_has_water ) {
		//  protocols::dna::add_water_design_flexibility_to_packer_task( pose, base_is_flexible_water_anchor, *task,
		//                                 allow_new_waters );
		// }
	}

	// clock_t starttime( clock() );
	//cerr << "hacking" << endl;
	pack::pack_rotamers_loop( pose, env_indep_scorefxn, task, pack_rotamers_nloop );

	// if ( use_envdep_packing ) {
	//  Real const pack_scorefxn_score1( env_indep_scorefxn( pose ) );
	//  Real const full_scorefxn_score1(   env_dep_scorefxn( pose ) );
	//  Real const pack1time( ( (double) clock() - starttime )/( CLOCKS_PER_SEC*60 ) );
	//  starttime = clock();
	//  pack::task::PackerTaskOP task2( pack::task::TaskFactory::create_packer_task( pose ) );
	//  { // setup the task
	//   bools allow_repack( pose.total_residue(), false );
	//   for ( Size i=1; i<= pose.total_residue(); ++i ) if ( mm.get_chi( i ) ) allow_repack[i] = true;
	//   task2->initialize_from_command_line().restrict_to_repacking().restrict_to_residues(allow_repack);
	//   task2->or_include_current( true );
	//   task2->set_bump_check( false );
	//  }
	//  pack::pack_rotamers_envdep_lowtemp( pose, env_dep_scorefxn, task2 ); // dont re-use task
	//  Real const pack_scorefxn_score2( env_indep_scorefxn( pose ) );
	//  Real const full_scorefxn_score2(   env_dep_scorefxn( pose ) );
	//  Real const pack2time( ( (double) clock() - starttime )/( CLOCKS_PER_SEC*60 ) );
	//  TR_REBUILD.Trace << "envdep_refinement: " <<
	//   " scores_after_pack1 " << F(9,3,pack_scorefxn_score1) << F(9,3,full_scorefxn_score1) << F(9,3,pack1time) <<
	//   " scores_after_pack2 " << F(9,3,pack_scorefxn_score2) << F(9,3,full_scorefxn_score2) << F(9,3,pack2time) <<
	//   endl;
	// }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Sizes const empty_poslist;

void
mike_fast_relax(
	ScoreFunction const & scorefxn_in,
	MoveMap const & mm_in,
	Pose & pose,
	Size const n_outer = 5, // default in FastRelax is 5
	Sizes const mutable_dna_positions = empty_poslist
)
{
	Reals const fa_rep_scales( make_vector1( 0.02, 0.25, 0.55, 1.0 ) );
	Reals const mintols( make_vector1( 0.01, 0.01, 0.01, 0.00001 ) );

	Real best_score( scorefxn_in( pose ) );
	Pose best_pose( pose );

	kinematics::MoveMapOP mm( new MoveMap( mm_in ) );// no mm clone in blab mini ?!

	//Size const n_outer( 5 );
	Size const n_inner( mintols.size() );

	for ( Size n=1; n<= n_outer; ++n ) {

		for ( Size m=1; m<= n_inner; ++m ) {
			Real const fa_rep_scale( fa_rep_scales[m] );
			Real const mintol( mintols[m] );

			/// setup scorefxn
			ScoreFunctionOP scorefxn( scorefxn_in.clone() );

			scorefxn->set_weight( fa_rep, fa_rep_scale * scorefxn_in.get_weight( fa_rep ) );

			ScoreFunctionOP pack_scorefxn( scorefxn );
			Real startscore( (*scorefxn)(pose) ), repackscore, minscore, repacktime, mintime;

			{ // repack
				/// setup flexible positions for packing
				pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ) );
				task->initialize_from_command_line();
				task->or_include_current( true );

				if ( !mutable_dna_positions.empty() ) { // setup residue couplings, allow design of DNA
					using namespace pack::rotamer_set;
					RotamerCouplingsOP couplings( new RotamerCouplings() );
					couplings->resize( pose.total_residue() );
					for ( Size i=1; i<= pose.total_residue(); ++i ) {
						if ( retrieve_base_partner_from_pose(pose)[i] ) {
							(*couplings)[i].first = retrieve_base_partner_from_pose(pose)[i];
							(*couplings)[i].second = ResidueMatcherCOP( new conformation::WatsonCrickResidueMatcher() );
						}
					}

					task->rotamer_couplings( couplings );
					for ( Size i=1; i<= pose.total_residue(); ++i ) {
						bool const is_mm_flexible( mm->get( id::TorsionID( i, id::CHI, 1) ) );
						bool const is_mutable( has_element( mutable_dna_positions, i ) ||
							has_element( mutable_dna_positions, retrieve_base_partner_from_pose( pose )[i] ) );
						if ( is_mutable && is_mm_flexible ) {
							pbassert( pose.residue(i).is_DNA() );
							task->nonconst_residue_task( i ).allow_aa( na_ade );
							task->nonconst_residue_task( i ).allow_aa( na_thy );
							task->nonconst_residue_task( i ).allow_aa( na_gua );
							task->nonconst_residue_task( i ).allow_aa( na_cyt );
							pbassert( task->design_residue(i) );
						} else if ( is_mm_flexible ) {
							task->nonconst_residue_task( i ).restrict_to_repacking();
						} else {
							task->nonconst_residue_task( i ).prevent_repacking();
						}
					}
				} else { // not designing the dna
					bools allow_repack( pose.total_residue(), false );
					for ( Size i=1; i<= pose.total_residue(); ++i ) { // use chi1 as proxy
						if ( mm->get( id::TorsionID( i, id::CHI, 1)) ) allow_repack[i] = true;
					}
					task->restrict_to_repacking().restrict_to_residues(allow_repack);
				}

				clock_t starttime = clock();
				// if ( !dry_run() ) pack::pack_rotamers_loop( pose, *pack_scorefxn, task, 25 ); // hacking
				if ( !dry_run() ) pack::pack_rotamers_loop( pose, *pack_scorefxn, task, 25 ); // hacking
				//if ( !dry_run() ) pack::pack_rotamers( pose, *pack_scorefxn, task );

				repackscore = (*scorefxn)(pose);
				repacktime = ( (double) clock() - starttime )/( CLOCKS_PER_SEC*60 ); // in minutes
			}


			{ // minimize
				//bools is_interface;
				string const mintype( "dfpmin_armijo" );

				protocols::minimization_packing::MinMoverOP min_mover
					( new protocols::minimization_packing::MinMover( mm, scorefxn, mintype, mintol, true ) );

				clock_t starttime = clock();
				if ( !dry_run() ) min_mover->apply( pose );
				minscore = (*scorefxn)(pose);
				mintime = ( (double) clock() - starttime )/( CLOCKS_PER_SEC*60 ); // in minutes
			}

			cout << "FAST_RELAX_CYCLE " << I(3,n) << I(3,m) <<
				F(9,3,scorefxn->get_weight( fa_rep ) ) <<
				F(12,9,mintol) << ' ' << F(9,3,startscore) <<
				F(9,3,repacktime ) << ' ' << F(9,3,repackscore) <<
				F(9,3,mintime ) << ' ' << F(9,3,minscore) << ' ' << F(9,3,best_score ) << endl;

			//basic::prof_show_oneliner();

		} // m=1,n_inner

		///
		Real const score( scorefxn_in( pose ) );
		if ( score < best_score ) {
			best_score = score;
			best_pose = pose;
		}
	}

	/// recover best pose
	pose = best_pose;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
mike_slow_relax(
	ScoreFunction const & scorefxn_in,
	MoveMap const & mm_in,
	Pose & pose,
	Size const n_outer = 5, // default in FastRelax is 5
	Sizes const mutable_dna_positions = empty_poslist,
	Size const packer_nloop = 15
)
{
	protocols::moves::PyMOLMoverOP pymol(0);
	if ( option[ my_options::show_pymol ] ) {
		static Size counter(0);
		++counter;
		pymol = protocols::moves::PyMOLMoverOP( new protocols::moves::PyMOLMover());
		pymol->pymol_name("mike_slow_relax_"+lead_zero_string_of(counter,2));
		pymol->keep_history(true);
	}

	Reals const fa_rep_scales( make_vector1( 0.02, 0.25, 0.55, 1.0 ) );
	Reals const mintols( make_vector1( 0.01, 0.01, 0.01, 0.00001 ) );

	Real best_score( scorefxn_in( pose ) );
	Pose best_pose( pose );

	kinematics::MoveMapOP mm( new MoveMap( mm_in ) );// no mm clone in blab mini ?!

	//Size const n_outer( 5 );
	Size const n_inner( mintols.size() );

	for ( Size n=1; n<= n_outer; ++n ) {

		for ( Size m=1; m<= n_inner; ++m ) {
			Real const fa_rep_scale( fa_rep_scales[m] );
			Real const mintol( mintols[m] );

			/// setup scorefxn
			ScoreFunctionOP scorefxn( scorefxn_in.clone() );

			scorefxn->set_weight( fa_rep, fa_rep_scale * scorefxn_in.get_weight( fa_rep ) );

			ScoreFunctionOP pack_scorefxn( scorefxn );
			Real startscore( (*scorefxn)(pose) ), repackscore, minscore, repacktime, mintime;

			PoseOPs repacked_poses;
			{ // repack
				/// setup flexible positions for packing
				pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ) );
				task->initialize_from_command_line();
				task->or_include_current( true );

				if ( !mutable_dna_positions.empty() ) { // setup residue couplings, allow design of DNA
					using namespace pack::rotamer_set;
					RotamerCouplingsOP couplings( new RotamerCouplings() );
					couplings->resize( pose.total_residue() );
					for ( Size i=1; i<= pose.total_residue(); ++i ) {
						if ( retrieve_base_partner_from_pose(pose)[i] ) {
							(*couplings)[i].first = retrieve_base_partner_from_pose(pose)[i];
							(*couplings)[i].second = ResidueMatcherCOP( new conformation::WatsonCrickResidueMatcher() );
						}
					}

					task->rotamer_couplings( couplings );
					for ( Size i=1; i<= pose.total_residue(); ++i ) {
						bool const is_mm_flexible( mm->get( id::TorsionID( i, id::CHI, 1) ) );
						bool const is_mutable( has_element( mutable_dna_positions, i ) ||
							has_element( mutable_dna_positions, retrieve_base_partner_from_pose( pose )[i] ) );
						if ( is_mutable && is_mm_flexible ) {
							pbassert( pose.residue(i).is_DNA() );
							task->nonconst_residue_task( i ).allow_aa( na_ade );
							task->nonconst_residue_task( i ).allow_aa( na_thy );
							task->nonconst_residue_task( i ).allow_aa( na_gua );
							task->nonconst_residue_task( i ).allow_aa( na_cyt );
							pbassert( task->design_residue(i) );
						} else if ( is_mm_flexible ) {
							task->nonconst_residue_task( i ).restrict_to_repacking();
						} else {
							task->nonconst_residue_task( i ).prevent_repacking();
						}
					}
				} else { // not designing the dna
					bools allow_repack( pose.total_residue(), false );
					for ( Size i=1; i<= pose.total_residue(); ++i ) { // use chi1 as proxy
						if ( mm->get( id::TorsionID( i, id::CHI, 1)) ) allow_repack[i] = true;
					}
					task->restrict_to_repacking().restrict_to_residues(allow_repack);
				}

				clock_t starttime = clock();
				// if ( !dry_run() ) pack::pack_rotamers_loop( pose, *pack_scorefxn, task, 25 ); // hacking

				vector1< std::pair< Real, string > > results;
				if ( !dry_run() ) pack::pack_rotamers_loop( pose, *pack_scorefxn, task, packer_nloop, results, repacked_poses );

				// look for duplicates
				while ( true ) {
					Real const epsilon( 1e-3 );
					bool found_a_duplicate( false );
					for ( Size i=1; i<= results.size() && !found_a_duplicate; ++i ) {
						for ( Size j=1; j<i && !found_a_duplicate; ++j ) {
							if ( results[i].second == results[j].second ) {
								Real const scoredev( fabs( results[i].first - results[j].first ) );
								if ( scoredev<epsilon ) {
									// delete i
									found_a_duplicate = true;
									results.erase( results.begin()+i-1 );
									repacked_poses.erase( repacked_poses.begin()+i-1 );
									TR_REBUILD.Trace << "found_a_duplicate: " << i << endl;
									break;
								}
							}
						}
					}
					if ( !found_a_duplicate ) break;
				}

				repackscore = (*scorefxn)(pose); // the best by total score after nloop
				repacktime = ( (double) clock() - starttime )/( CLOCKS_PER_SEC*60 ); // in minutes
			}


			{ // minimize
				//bools is_interface;
				string const mintype( "dfpmin_armijo" );

				protocols::minimization_packing::MinMoverOP min_mover
					( new protocols::minimization_packing::MinMover( mm, scorefxn, mintype, mintol, true ) );

				clock_t starttime = clock();

				/// try minimizing each of the poses generated by the simannealer...
				Real bestenergy(0);
				runtime_assert( repacked_poses.size() <= packer_nloop );
				for ( Size i=1; i<= repacked_poses.size(); ++i ) {
					Pose & repacked_pose( *repacked_poses[i] );
					if ( !dry_run() ) min_mover->apply( repacked_pose );
					Real const score( (*scorefxn)( repacked_pose ) );
					TR_REBUILD.Trace << "minloop " << I(4,i) << ' ' << F(9,3,score) << std::endl;
					if ( i==1 || score < bestenergy ) {
						bestenergy = score;
						pose = repacked_pose; // save the repacked_pose
						if ( pymol ) pymol->apply( pose );
					}
				}
				minscore = (*scorefxn)(pose);
				mintime = ( (double) clock() - starttime )/( CLOCKS_PER_SEC*60 ); // in minutes
			}

			cout << "FAST_RELAX_CYCLE " << I(3,n) << I(3,m) <<
				F(9,3,scorefxn->get_weight( fa_rep ) ) <<
				F(12,9,mintol) << ' ' << F(9,3,startscore) <<
				F(9,3,repacktime ) << ' ' << F(9,3,repackscore) <<
				F(9,3,mintime ) << ' ' << F(9,3,minscore) << ' ' << F(9,3,best_score ) << ' ' << repacked_poses.size() << endl;

			basic::prof_show_oneliner();

		} // m=1,n_inner

		///
		Real const score( scorefxn_in( pose ) );
		if ( score < best_score ) {
			best_score = score;
			best_pose = pose;
		}
	}

	/// recover best pose
	pose = best_pose;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
//
// ramp up the repulsive and the chainbreak weights
//
//
void
fast_star_relax(
	ScoreFunction const & scorefxn_in,
	bools const & is_frozen,
	bool const dry_run,
	bool const use_coordinate_constraint_weight,
	Pose & pose,
	bool const keep_coordinate_constraints_constant = false,
	Size nouter = 3,
	bool const allow_new_waters_in_early_rounds = false
)
{
	using namespace protocols::moves;
	using namespace protocols::minimization_packing;
	using devel::blab::get_ramped_weight;

	Size const vrtpos( get_vrt_pos( pose ) ); // will be 0 if there's not VRT residue

	// params
	Size ninner( 4 );
	//Size nouter( 3 );
	bool const debug( false );
	if ( debug ) nouter = ninner = 1;

	Real const mc_temp( 0.0 ); // only accept poses that get better...
	//bool const flex_water_jumps( true );
	bool const vary_omega( true );
	bool const flex_dna_sugar( true );
	//Size const pack_mover_nloop( 25 );

	// make a local copy
	ScoreFunctionOP scorefxn( scorefxn_in.clone() );

	/// move map -- includes both protein and dna
	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	{
		Size const nres( pose.total_residue() );
		bools is_chi_flexible( nres, true ), is_bb_flexible( nres, true ), is_jump_flexible( nres, true );
		if ( vrtpos ) is_chi_flexible[ vrtpos ] = is_bb_flexible[ vrtpos ] = is_jump_flexible[ vrtpos ] = false;
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( is_frozen[i] ) is_chi_flexible[ i ] = is_bb_flexible[ i ] = is_jump_flexible[ i ]= false;
		}
		devel::blab::setup_move_map( is_chi_flexible, is_bb_flexible, is_jump_flexible, vary_omega, flex_dna_sugar, pose,
			*mm );
	}


	/// MC
	MonteCarloOP mc( new MonteCarlo( pose, *scorefxn, mc_temp ) );

	/// min mover
	MinMoverOP min_mover( new MinMover( mm, scorefxn, "dfpmin_armijo_atol", 0.05, true ) ); // tol-unused,use_nblist

	/// task factories
	// pack::task::TaskFactoryOP pack_task_factory
	//  ( devel::blab::setup_task_factory_from_move_map( pose, *mm, false, true ) ); // use_bump_check, inc_current
	// PackRotamersMoverOP repack_mover( new PackRotamersMover( scorefxn, 0 /* task */, pack_mover_nloop ) );
	// repack_mover->task_factory( pack_task_factory );


	Real const rep_weight_full( scorefxn->get_weight( fa_rep ) );
	Real const dun_weight_full( scorefxn->get_weight( fa_dun ) );
	Real const chainbreak_protein_full( scorefxn->get_weight( chainbreak ) );//_protein ) );
	//Real const chainbreak_dna_full( scorefxn->get_weight( chainbreak_dna ) );
	Real coordinate_constraint_full( use_coordinate_constraint_weight ? 1.0 : 0.0 );


	/// now do the ramping thing
	Real const atol_b( 0.5 ), atol_e( 0.02 );
	Real const rep_b( 0.02 ), rep_e( 1.0 );
	Real const cbp_b( 0.25 ), cbp_e( 1.0 ); // chainbreak protein, well actually all chainbreaks
	//Real const cbd_b( 0.75 ), cbd_e( 1.0 );
	Real coord_b( 0.2 ), coord_e( 1.0 ); // ramp them up, not down
	Real dun_b( 1.0 ), dun_e( 1.0 );
	if ( option[ my_options::ramp_fa_dun ] ) dun_b = 0.1;
	if ( keep_coordinate_constraints_constant ) {
		coordinate_constraint_full = scorefxn->get_weight( coordinate_constraint );
		coord_b = coord_e = 1.0;
		pbassert( coordinate_constraint_full > 0.01 ); // assert that we have set this in scorefxn
	}

	scorefxn->set_weight( fa_rep, rep_weight_full * rep_e );
	scorefxn->set_weight( chainbreak/*_protein*/, chainbreak_protein_full * cbp_e );
	//scorefxn->set_weight( chainbreak_dna, chainbreak_dna_full * cbd_e );
	scorefxn->set_weight( coordinate_constraint, coordinate_constraint_full * coord_e );
	mc->score_function( *scorefxn );

	static Size const debug_N( int( numeric::random::uniform()*1000.0 ) );
	static Size debug_M( 0 );
	++debug_M;

	basic::Tracer status( "blab_status" );
	for ( Size ii=1; ii<= nouter; ++ii ) {


		for ( Size jj=1; jj<= ninner; ++jj ) {
			string const debug_tag( "fast_star_relax_"+
				lead_zero_string_of(debug_N,4) +"_" +lead_zero_string_of(debug_M,4) + "_"+
				lead_zero_string_of(ii,4) +"_" +lead_zero_string_of(jj,4) );
			if ( debug && ii==1 && jj == 1 ) pose.dump_pdb( debug_tag + "_begin.pdb" );

			/// set the repulsive, chainbreaks
			scorefxn->set_weight( fa_rep, rep_weight_full * get_ramped_weight( jj, ninner, rep_b, rep_e, false ) );
			scorefxn->set_weight( fa_dun, dun_weight_full * get_ramped_weight( jj, ninner, dun_b, dun_e, false ) );
			scorefxn->set_weight( chainbreak/*_protein*/, chainbreak_protein_full *
				get_ramped_weight( jj, ninner, cbp_b, cbp_e, false ) );
			//    scorefxn->set_weight( chainbreak_dna, chainbreak_dna_full *
			//               get_ramped_weight( jj, ninner, cbd_b, cbd_e, false ) );
			scorefxn->set_weight( coordinate_constraint, coordinate_constraint_full *
				get_ramped_weight( jj, ninner, coord_b, coord_e, false ) );

			Real const min_atol( get_ramped_weight( jj, ninner, atol_b, atol_e, true ) ); // geometric ramping
			min_mover->tolerance( min_atol );

			//
			ScoreFunctionOP pack_scorefxn( scorefxn );
			bool const use_envdep_packing( false );//pack::score_function_has_local_context_dependent_methods( *pack_scorefxn ) );
			//if ( use_envdep_packing ) pack_scorefxn = create_env_indep_score_function( *scorefxn );

			Real const start_score( (*scorefxn)( pose ) );

			clock_t starttime = clock();
			//if ( !dry_run ) repack_mover->apply( pose );
			if ( !dry_run ) {
				bool const use_envdep_packing_this_round( use_envdep_packing && jj >= ninner-1 );
				bool const allow_new_waters_this_round( allow_new_waters_in_early_rounds || jj == ninner );
				envdep_pack_rotamers_wrapper( *mm, *pack_scorefxn, *scorefxn, use_envdep_packing_this_round, pose,
					allow_new_waters_this_round );
			}
			clock_t stoptime = clock();

			Real const repack_time(  ((double) stoptime - starttime )/( CLOCKS_PER_SEC*60 ) );
			Real const repack_score( (*scorefxn)( pose ) );

			starttime = clock();
			if ( !dry_run ) min_mover->apply( pose );
			stoptime = clock();

			Real const min_time(  ((double) stoptime - starttime )/( CLOCKS_PER_SEC*60 ) );
			Real const min_score( (*scorefxn)( pose ) );

			Real const mc_score( mc->score_function()( pose ) );

			status << "FAST_RELAX_CYCLE " << I(3,ii) << I(3,jj) <<
				F(9,3,scorefxn->get_weight( fa_rep ) ) <<
				F(9,3,scorefxn->get_weight( chainbreak/*_protein*/ ) ) <<
				//F(9,3,scorefxn->get_weight( chainbreak_dna ) ) <<
				F(9,3,scorefxn->get_weight( coordinate_constraint ) ) <<
				F(9,3,min_atol) << ' ' << F(9,3,start_score) <<
				F(9,3,repack_time ) << ' ' << F(9,3,repack_score) <<
				F(9,3,min_time ) << ' ' << F(9,3,min_score) << ' ' << F(9,3,mc_score ) << endl;

			if ( debug ) pose.dump_pdb( debug_tag + "_end.pdb" );
		}
		mc->boltzmann( pose ); // check energy at this point...
		mc->recover_low( pose );

	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
my_fast_relax(
	ScoreFunction const & scorefxn_in,
	MoveMap const & mm_in,
	bool const cartesian,
	Pose & pose,
	Size const n_outer = 5 // default in FastRelax is 5
)
{
	protocols::moves::PyMOLMoverOP pymol(0);
	if ( option[ my_options::show_pymol ] ) {
		static Size counter(0);
		++counter;
		pymol = protocols::moves::PyMOLMoverOP( new protocols::moves::PyMOLMover());
		pymol->pymol_name("my_fast_relax_"+lead_zero_string_of(counter,2));
		pymol->keep_history(true);
	}

	if ( cartesian ) {
		runtime_assert( scorefxn_in.get_weight( cart_bonded )>1e-3 || scorefxn_in.get_weight( cart_bonded_length ) >1e-3 );
	}
	Reals const fa_rep_scales( make_vector1( 0.02, 0.25, 0.55, 1.0 ) );
	Reals const mintols( make_vector1( 0.01, 0.01, 0.01, 0.00001 ) );

	Real best_score( scorefxn_in( pose ) );
	Pose best_pose( pose );

	MoveMapOP mm( mm_in.clone() );

	//Size const n_outer( 5 );
	Size const n_inner( mintols.size() );

	/// we will only use envdep packing in the final inner iterations...
	bool const use_envdep_packing( false );//pack::score_function_has_local_context_dependent_methods( scorefxn_in ) );

	for ( Size n=1; n<= n_outer; ++n ) {

		for ( Size m=1; m<= n_inner; ++m ) {
			Real const fa_rep_scale( fa_rep_scales[m] );
			Real const mintol( mintols[m] );
			//bool const cartesian( false ); // for the moment

			/// setup scorefxn
			ScoreFunctionOP scorefxn( scorefxn_in.clone() );

			scorefxn->set_weight( fa_rep, fa_rep_scale * scorefxn_in.get_weight( fa_rep ) );

			ScoreFunctionOP pack_scorefxn( scorefxn );

			runtime_assert( !use_envdep_packing );
			//if ( use_envdep_packing ) pack_scorefxn = create_env_indep_score_function( *scorefxn );

			Real startscore( (*scorefxn)(pose) ), repackscore, minscore, repacktime, mintime;

			{ // repack
				/// setup flexible positions for packing
				// pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ) );

				// bools allow_repack( pose.total_residue(), false );
				// for ( Size i=1; i<= pose.total_residue(); ++i ) if ( mm->get_chi( i ) ) allow_repack[i] = true;
				// task->initialize_from_command_line().restrict_to_repacking().restrict_to_residues(allow_repack);
				// task->or_include_current( true );

				clock_t starttime = clock();
				bool const use_envdep_packing_this_round( use_envdep_packing && m>= n_inner-1 );
				if ( !dry_run() ) {
					envdep_pack_rotamers_wrapper( *mm, *pack_scorefxn, *scorefxn, use_envdep_packing_this_round,
						pose );
				}
				//  pack::pack_rotamers( pose, *pack_scorefxn, task );
				//  if ( use_envdep_packing && m >= n_inner-1 ) { // refine with the true scorefxn...
				//   Real const pack_scorefxn_score1( (*pack_scorefxn)( pose ) );
				//   Real const full_scorefxn_score1( (*scorefxn)( pose ) );
				//   Real const pack1time( ( (double) clock() - starttime )/( CLOCKS_PER_SEC*60 ) );
				//   clock_t starttime2 = clock();
				//   pack::pack_rotamers_envdep_lowtemp( pose, *scorefxn, task ); // OK to re-use the task? since not designing...
				//   Real const pack_scorefxn_score2( (*pack_scorefxn)( pose ) );
				//   Real const full_scorefxn_score2( (*scorefxn)( pose ) );
				//   Real const pack2time( ( (double) clock() - starttime2 )/( CLOCKS_PER_SEC*60 ) );
				//   cout << "envdep_refinement: " <<
				//    " scores_after_pack1 " << F(9,3,pack_scorefxn_score1) << F(9,3,full_scorefxn_score1) << F(9,3,pack1time) <<
				//    " scores_after_pack2 " << F(9,3,pack_scorefxn_score2) << F(9,3,full_scorefxn_score2) << F(9,3,pack2time) <<
				//    endl;
				//  }
				// }

				repackscore = (*scorefxn)(pose);
				repacktime = ( (double) clock() - starttime )/( CLOCKS_PER_SEC*60 ); // in minutes
				basic::prof_show_oneliner();
				if ( pymol ) pymol->apply( pose );
			}


			{ // minimize
				//bools is_interface;
				string const mintype( cartesian ? "lbfgs_armijo_nonmonotone" : "dfpmin_armijo_nonmonotone" );

				protocols::minimization_packing::MinMoverOP min_mover
					( new protocols::minimization_packing::MinMover( mm, scorefxn, mintype, mintol, true ) );

				min_mover->cartesian( cartesian );
				clock_t starttime = clock();
				if ( !dry_run() ) min_mover->apply( pose );
				minscore = (*scorefxn)(pose);
				mintime = ( (double) clock() - starttime )/( CLOCKS_PER_SEC*60 ); // in minutes
				if ( pymol ) pymol->apply( pose );
			}

			cout << "FAST_RELAX_CYCLE " << I(3,n) << I(3,m) <<
				F(9,3,scorefxn->get_weight( fa_rep ) ) <<
				F(12,9,mintol) << ' ' << F(9,3,startscore) <<
				F(9,3,repacktime ) << ' ' << F(9,3,repackscore) <<
				F(9,3,mintime ) << ' ' << F(9,3,minscore) << ' ' << F(9,3,best_score ) << endl;

			basic::prof_show_oneliner();

		} // m=1,n_inner

		///
		Real const score( scorefxn_in( pose ) );
		if ( score < best_score ) {
			best_score = score;
			best_pose = pose;
		}
	}

	/// recover best pose
	pose = best_pose;

	if ( pymol ) pymol->apply( pose );
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size // actual min flexible seglen
setup_protein_segments_from_flexible_positions(
	bools was_flexible, // before expanding
	bools const & is_frozen,
	Pose const & pose,
	bools & is_flexible,
	vector1< Sizes > & protein_segments,
	bool const dont_extend_flex_at_start = false
)
{
	/// mod was_flexible to reflect current length and to freeze out non-protein
	was_flexible.resize( pose.total_residue(), false ); // in case vrt rsd at end added after was_flexible
	is_flexible = was_flexible;
	bools is_protein( pose.total_residue(), false );
	for ( Size i=1; i<= pose.total_residue(); ++i ) is_protein[i]= pose.residue(i).is_protein();


	Size const min_flex_seglen( 6 ), min_fix_seglen(4), max_flex_seglen( 10 );

	Size nres_flexible_old(0);

	// add residues on either side of flexible sidechains
	if ( !dont_extend_flex_at_start ) {
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( was_flexible[i] ) ++nres_flexible_old;
			if ( !is_protein[i] ) continue;
			Size const ichain( pose.chain(i) );
			if ( i>1 && pose.chain(i-1) == ichain && was_flexible[i-1] && !is_frozen[i] ) is_flexible[i]= true;
			if ( i<pose.total_residue() && pose.chain(i+1) == ichain && was_flexible[i+1] && !is_frozen[i] ) {
				is_flexible[i]= true;
			}
			if ( was_flexible[i] ) ++nres_flexible_old;
		}
	}


	SizePairs fix_segments, flex_segments;
	bool made_a_change( true );
	while ( made_a_change ) {
		made_a_change = false;
		fix_segments.clear();
		flex_segments.clear();
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			// at the start of a segment
			Size j(i), j_end( chain_end( pose.chain(i), pose ) );
			while ( j<j_end && is_flexible[i] == is_flexible[j+1] ) ++j;
			// segment running from i to j
			if ( is_flexible[i] && is_protein[i] ) {
				flex_segments.push_back( make_pair( i,j ) );
				TR_REBUILD.Trace << "FLEXPROT " << I(6,i) << I(6,j) << " len: " << I(4,j-i+1) << " ch: " << pose.chain(i) << endl;
			} else {
				fix_segments.push_back( make_pair(i,j) );
				TR_REBUILD.Trace << "FIXORDNA " << I(6,i) << I(6,j) << " len: " << I(4,j-i+1) << " ch: " << pose.chain(i) << endl;
			}
			i = j;
		}

		// look for short fixed segments surrounded on both sides by flexible segments in the same chain
		foreach_ ( SizePair p, fix_segments ) {
			runtime_assert( pose.chain(p.first) == pose.chain(p.second) );
			if ( p.second - p.first + 1 < min_fix_seglen &&
					p.first > 1 && pose.chain( p.first ) == pose.chain( p.first-1 ) && is_flexible[ p.first-1] &&
					p.second < pose.total_residue() && pose.chain( p.second ) == pose.chain( p.second+1 ) &&
					is_flexible[ p.second+1 ] ) {
				bool frozen( false );
				for ( Size i=p.first; i<= p.second; ++i ) frozen = ( frozen || is_frozen[i] );
				if ( !frozen ) {
					for ( Size i=p.first; i<= p.second; ++i ) is_flexible[i] = true;
					made_a_change = true;
					break;
				}
			}
		}
		if ( made_a_change ) continue;

		// look for short flexible segments
		foreach_ ( SizePair p, flex_segments ) {
			Size const plen( p.second - p.first + 1 ), pcb( chain_begin( pose.chain(p.first), pose ) ),
				pce( chain_end( pose.chain(p.first), pose ) );
			runtime_assert( pose.chain(p.first) == pose.chain(p.second) );
			if ( plen < min_flex_seglen ) { // extend, deterministically
				bool const
					cant_go_backward( p.first == pcb || is_frozen[p.first-1] ),
					cant_go_forward( p.second == pce || is_frozen[p.second+1] ),
					go_backward( !cant_go_backward && ( plen%2==0 || cant_go_forward ) ),
					go_forward ( !cant_go_forward  && ( plen%2==1 || cant_go_backward ) );
				TR_REBUILD.Trace << "go? " << go_forward << ' ' << go_backward << ' ' << p.first << ' ' << p.second << ' ' <<
					pose.chain(p.first) << ' ' << pcb << ' ' << pce << endl;
				if      ( go_forward  ) { is_flexible[p.second+1] = true; made_a_change = true; }
				else if ( go_backward ) { is_flexible[ p.first-1] = true; made_a_change = true; }
				else {
					TR_REBUILD.Trace << "unable to extend FLEX segment " <<
						go_forward << ' ' << go_backward << ' ' << p.first << ' ' << p.second << ' ' <<
						pose.chain(p.first) << ' ' << pcb << ' ' << pce << endl;
				}
				if ( made_a_change ) break;
			}
		}
		if ( made_a_change ) continue;
	} // while ( made_a_change )


	// now we've got a set of fixed and flexible segments that cover the protein
	// we may want to break up the flexible segments that are too long

	for ( Size i=1, i_end = flex_segments.size(); i<= i_end; ++i ) { // skip any new ones we create
		SizePair & p( flex_segments[i] );
		Size const plen( p.second - p.first + 1 );
		Size ncuts(0);
		while ( Real(plen) / (ncuts+1) > max_flex_seglen ) ++ncuts;
		if ( ncuts ) {
			TR_REBUILD.Trace << "break segment " << p.first << ' ' << p.second << ' ' << plen << endl;
			Real const cutsize( Real(plen)/(ncuts+1) );
			SizePair const oldp(p);
			for ( Size i=1; i<= ncuts+1; ++i ) {
				SizePair newp;
				if ( i==1 ) newp.first = oldp.first;
				else newp.first  = Size( floor( 0.5 + ( oldp.first-1 ) + (i-1)*cutsize ) ) + 1; // +1
				if ( i==ncuts+1 ) newp.second = oldp.second;
				else newp.second = Size( floor( 0.5 + ( oldp.first-1 ) +     i*cutsize ) );
				TR_REBUILD.Trace << "newp " << newp.first << ' ' << newp.second << endl;
				if ( i== 1 ) p = newp;
				else flex_segments.push_back( newp );
			}
		}
	}

	{ // debug
		bools is_flex( pose.total_residue(), false ), is_fix( pose.total_residue(), false );
		foreach_ ( SizePair p, flex_segments ) {
			for ( Size i=p.first; i<= p.second; ++i ) {
				runtime_assert( !is_flex[i] );
				is_flex[i] = true;
			}
		}
		foreach_ ( SizePair p, fix_segments ) {
			for ( Size i=p.first; i<= p.second; ++i ) {
				runtime_assert( !is_fix[i] );
				is_fix[i] = true;
			}
		}
		Size nres_flexible_new(0);
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			// TR.Trace << "fixorflex " << i << ' ' << is_fix[i] << ' ' << is_flex[i] << endl;
			runtime_assert( is_flex[i] || is_fix[i] );
			runtime_assert( !( is_flex[i] && is_fix[i] ) );
			nres_flexible_new += is_flex[i];
		}
		// TR.Trace << "flex_status: nres_flexible_old:  "<< nres_flexible_old << " nres_flexible_new: " <<
		//  nres_flexible_new << " n_flex_segs: " << flex_segments.size() << " n_fix_segs: " << fix_segments.size() <<
		//  ' ' << filename << endl;
	}
	// {
	//  ostringstream out;
	//  for ( Size i=1; i<= pose.total_residue(); ++i ) {
	//   if ( is_flexible[i] ) {
	//    if ( out.str().size()) out << '+';
	//    out << i;
	//   }
	//  }
	//  TR.Trace << "select is_flex, " << out.str() << "/" << endl;
	// }

	SizePairs all_segs;
	Size actual_min_flexible_seglen(1000);
	foreach_ ( SizePair p, flex_segments ) {
		all_segs.push_back(p);
		actual_min_flexible_seglen = min( actual_min_flexible_seglen, (p.second-p.first+1) );
	}
	foreach_ ( SizePair p, fix_segments ) all_segs.push_back(p);
	std::sort( all_segs.begin(), all_segs.end() );

	//vector1< Sizes > protein_segments;
	foreach_ ( SizePair p, all_segs ) {
		if ( !pose.residue(p.first).is_protein() ) continue; // eg VRT residue... they are "protein_segments" after all...
		runtime_assert( pose.residue(p.second).is_protein() );
		protein_segments.push_back( make_vector1( p.first, p.second ) );
	}

	return actual_min_flexible_seglen;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// legacy version
Size // actual min flexible seglen
setup_protein_segments_from_flexible_positions(
	bools const & was_flexible, // before expanding
	Pose const & pose,
	bools & is_flexible,
	vector1< Sizes > & protein_segments,
	bool const dont_extend_flex_at_start = false
)
{
	bools const is_frozen( pose.total_residue(), false );
	return setup_protein_segments_from_flexible_positions( was_flexible, is_frozen, pose, is_flexible, protein_segments,
		dont_extend_flex_at_start );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
centroid_star_rebuild_around_isolated_positions(
	bools is_core_flexible, // const except for resizing
	bools is_frozen, // ditto
	bool const restore_original_nres_fold_tree_and_cutpoint_variants,
	bools & is_flexible, // communicate additional bb flex to outside world
	Pose & pose,
	Real const coordcst_distol = 0.25
)
{
	Pose start_pose( pose );

	bool pose_has_dna( false ), flexdna( false );
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).is_DNA() ) {
			pose_has_dna = true;
			if ( is_core_flexible[i] ) {
				flexdna = true;
				runtime_assert( retrieve_base_partner_from_pose( pose )[i] ); // assuming flex dna posns are paired
			}
		}
	}

	append_virtual_residue( pose ); // array bounds problems?
	if ( !restore_original_nres_fold_tree_and_cutpoint_variants ) append_virtual_residue( start_pose );

	is_core_flexible.resize( pose.total_residue(), false );
	is_frozen.resize( pose.total_residue(), false );

	// pick segments
	vector1< Sizes > protein_segments;
	Size const min_flex_seglen
		( setup_protein_segments_from_flexible_positions( is_core_flexible, is_frozen, pose, is_flexible, protein_segments ) );

	//
	SizePairs dna_segments;
	Size min_flex_dna_seglen(1000);
	if ( flexdna ) {
		scoring::dna::BasePartner const & partner
			( scoring::dna::retrieve_base_partner_from_pose( pose ));
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( pose.residue(i).is_DNA() ) {
				Size const ip( partner[i] );
				if ( is_core_flexible[i] || ( ip && is_core_flexible[ip] ) ) {
					runtime_assert( ip );
					TR_REBUILD.Trace << "centroid_star_rebuild_around_isolated_positions: flexdna: " << i << ' ' << ip << endl;
					is_flexible[i] = is_flexible[ip] = true;
				} else {
					is_flexible[i] = false;
					if ( ip ) is_flexible[ip] = false;
				}
			}
		}
		bools in_segment( pose.total_residue(), false );
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( !pose.residue(i).is_DNA() || in_segment[i] ) continue;

			// starting a new segment with i
			Size j(i), ce( chain_end( pose.chain(i), pose  ) );
			if ( partner[i] ) {
				// look ahead to end of segment
				while ( j< ce && !pose.residue(j).is_upper_terminus() && partner[j+1] && is_core_flexible[i] == is_core_flexible[j+1] &&
						( partner[i] - partner[j+1] == j+1-i ) ) ++j;
			} else { // unpaired segment
				while ( j< ce && !pose.residue(j).is_upper_terminus() && !partner[j+1] ) ++j;
			}
			dna_segments.push_back( make_pair(i,j) );
			TR_REBUILD.Trace << "dna_segment: " << dna_segments.size() << ' ' << i << ' ' << j << ' ' <<
				partner[i] << ' ' << partner[j] << ' '<< pose.residue(i).name1() << ' '<<
				pose.residue(j).name1() << ' ' << ( is_core_flexible[i] ? "FLEX" : "FIX" ) << endl;
			if ( is_core_flexible[i] ) {
				min_flex_dna_seglen = min( min_flex_dna_seglen, j-i+1 );
			}
			for ( Size ii=i; ii<=j; ++ii ) {
				in_segment[ii] = true;
				if ( partner[ii] ) in_segment[partner[ii]] = true;
			}
			i = j;
		}
	}


	if ( min_flex_seglen < 3  ) {
		foreach_ ( Sizes const & seg, protein_segments ) {
			if ( is_flexible[seg[1]] && seg[2]-seg[1]+1<3 ) {
				TR_REBUILD.Trace << "switching short segment to fixed " << seg[1] << ' ' << seg[2] << endl;
				for ( Size i=seg[1]; i<= seg[2]; ++i ) {
					runtime_assert( is_flexible[i] );
					is_flexible[i] = false;
				}
			}
		}
	}
	Size nres_flexible(0), nres_core_flexible(0);
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( is_flexible[i] ) ++nres_flexible;
		if ( is_core_flexible[i] ) ++nres_core_flexible;
	}

	TR_REBUILD.Trace << "centroid_star_rebuild_around_isolated_positions:: min_flex_seglen: " << min_flex_seglen <<
		" nres_flexible: " << nres_flexible << " nres_core_flexible: " << nres_core_flexible << endl;

	// pick fragments, taken from sctest ////////////////////////////////////
	devel::blab::classic_frags::FragLibOP protein_fraglib(0);
	{ // NOW setup frags, when we know what backbone regions are flexible; taken from setup_cheating_fragments
		Sizes vall_homs;
		Real const min_torsion_dev( option[ my_options::cheating_frags_min_torsion_dev ] ), // dflt 30
			max_torsion_dev( option[ my_options::cheating_frags_max_torsion_dev ] ); // dflt 90
		string const ss( pose.secstruct() );
		{ // confirm that ss has been initialized
			bool all_L( true );
			for ( Size i=0; i<ss.size(); ++i ) if ( ss[i] != 'L' ) all_L = false;
			if ( all_L ) utility_exit_with_message("Need to initialize pose's secondary structure.");
		}
		kinematics::MoveMap mm;
		for ( Size i=1; i<= pose.total_residue(); ++i ) mm.set_bb( i, is_flexible[i] );
		vector1< Size > frag_sizes;
		frag_sizes.push_back( 6 );
		frag_sizes.push_back( 3 );
		Size const nfrags( 200 );
		Real const seq_weight( 1.0 ), ss_weight( 3.0 ), torsion_weight( 5 );

		TR_REBUILD.Trace << "setup_cheating_fragments: ss= " << ss << endl;

		if ( !dry_run() ) {
			protein_fraglib =
				devel::blab::classic_frags::setup_vall_cheating_fragments( frag_sizes, nfrags, pose, mm, ss,
				seq_weight, ss_weight, torsion_weight,
				min_torsion_dev, max_torsion_dev,
				vall_homs );
		}
	} // scope ///////////////////


	// switch to centroid
	ScoreFunctionOP cen_scorefxn( get_centroid_score_function_from_command_line() );
	cen_scorefxn->set_weight( coordinate_constraint, 1.0 );
	if ( pose_has_dna ) {
		scoring::methods::EnergyMethodOptions options( cen_scorefxn->energy_method_options() );
		options.atom_vdw_atom_type_set_name( CENTROID_DNA );
		cen_scorefxn->set_energy_method_options( options );
	}
	// cen_scorefxn->set_weight( atom_pair_constraint, 1.0 ); // used for disulfides
	//runtime_assert( std::abs( fa_scorefxn->get_weight( chainbreak ) ) > 1e-3 );

	if ( false ) { // hack testing -- can we do these things in fullatom?
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			cout << "old name1: " << pose.residue(i).name() << endl;
			remove_variant_type_from_pose_residue( pose, VIRTUAL_DNA_PHOSPHATE, i );
			cout << "new name1: " << pose.residue(i).name() << endl;
		}
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( pose.residue(i).is_DNA() && !pose.residue(i).is_upper_terminus() ) {
				cout << "old name2: " << i << ' ' << pose.residue(i).name() << endl;
				add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, i );
				cout << "new name2: " << i << ' ' << pose.residue(i).name() << endl;
			}
		}
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( pose.residue(i).is_DNA() && !pose.residue(i).is_lower_terminus() ) {
				cout << "old name2: " << i << ' ' << pose.residue(i).name() << endl;
				add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, i );
				cout << "new name2: " << i << ' ' << pose.residue(i).name() << endl;
			}
		}
		cout << "done testing" << endl;
	}


	if ( pose_has_dna ) {
		devel::blab::classic_frags::switch_pose_to_residue_type_set( pose, CENTROID_DNA );
	} else {
		devel::blab::classic_frags::switch_pose_to_residue_type_set( pose, CENTROID );
	}

	// simulation
	// setup foldtree, cutpoints

	if ( flexdna ) {
		bool const anchor_dna_jumps_in_backbone( true );
		setup_star_fold_tree_cutpoint_variants_and_virtual_residue( anchor_dna_jumps_in_backbone, protein_segments, dna_segments, pose );

		if ( !restore_original_nres_fold_tree_and_cutpoint_variants ) {
			setup_star_fold_tree_cutpoint_variants_and_virtual_residue( anchor_dna_jumps_in_backbone, protein_segments, dna_segments, start_pose );
		}
	} else {
		bool const anchor_dna_jumps_in_backbone( true ), freeze_dna( true );
		setup_star_fold_tree_cutpoint_variants_and_virtual_residue( anchor_dna_jumps_in_backbone, protein_segments,
			pose, freeze_dna );
		if ( !restore_original_nres_fold_tree_and_cutpoint_variants ) {
			setup_star_fold_tree_cutpoint_variants_and_virtual_residue( anchor_dna_jumps_in_backbone, protein_segments,
				start_pose, freeze_dna ); // otherwise recover fixed sidechains fails
		}
	}

	// set constraints
	pose.constraint_set(0);
	bool const tether_phosphates( true ), tether_c1stars( false ); //option[ my_options::tether_centroid_c1stars ] );
	add_simple_cartesian_constraints( start_pose, pose, coordcst_distol, tether_phosphates, tether_c1stars );

	// add_centroid_disulfide_constraints( start_pose, pose );

	if ( !dry_run() ) {
		// do some centroid rebuilding
		bool const use_superimpose_segments( true );
		bools is_flexible_dna( pose.total_residue(), false ), is_flexible_protein( is_flexible );
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			is_flexible_protein[i] = is_flexible[i] && pose.residue(i).is_protein();
			is_flexible_dna[i] = flexdna && is_flexible[i] && pose.residue(i).is_DNA();
		}
		Size const big_frag_size(6);
		if ( !flexdna ) {
			star_fragment_rebuild( *cen_scorefxn, is_flexible_dna, is_flexible_protein, protein_fraglib, pose,
				10, 3, 8, use_superimpose_segments, big_frag_size );
		} else {
			TR_REBUILD.Trace << "centroid_star_rebuild_around_isolated_positions:: min_flex_dna_seglen= " << min_flex_dna_seglen << endl;
			runtime_assert( min_flex_dna_seglen >= 2 );
			Sizes dna_frag_sizes( make_vector1(2,3,5) );
			foreach_ ( Size & fragsize, dna_frag_sizes ) fragsize = min( fragsize, min_flex_dna_seglen );
			star_fragment_rebuild( *cen_scorefxn, is_flexible_dna, is_flexible_protein, protein_fraglib, pose,
				10, 3, 8, use_superimpose_segments, big_frag_size, dna_frag_sizes );
		}
	}

	// switch back to fullatom
	devel::blab::classic_frags::switch_pose_to_residue_type_set( pose, FA_STANDARD );
	//pose.conformation().fix_disulfides( ssbonds );

	// reset the constraints
	pose.constraint_set(0);


	// restore original fold tree
	if ( restore_original_nres_fold_tree_and_cutpoint_variants ) {
		if ( pose.total_residue() == start_pose.total_residue()+1 ) {
			FoldTree f( pose.fold_tree() );
			for ( Size i=1; i<= f.num_jump(); ++i ) {
				f.set_jump_atoms( i, "", "" );
			}
			pose.fold_tree(f);
			pose.conformation().delete_residue_slow( pose.total_residue() ); // remove virtual rsd
		}
		pose.fold_tree( start_pose.fold_tree() );
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( pose.residue(i).has_variant_type( CUTPOINT_LOWER ) &&
					!start_pose.residue(i).has_variant_type( CUTPOINT_LOWER ) ) {
				remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, i );
			}
			if ( pose.residue(i).has_variant_type( CUTPOINT_UPPER ) &&
					!start_pose.residue(i).has_variant_type( CUTPOINT_UPPER ) ) {
				remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, i );
			}
		}
	}

	// recover fixed sidechains:
	recover_frozen_protein_sidechains( start_pose, is_flexible, pose );

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
centroid_rebuild_loop(
	Size const old_begin,
	Size const old_end,
	string const & newseq,
	bool add_cheating_frags, // only if size is the same
	bool const add_sequence_frags,
	Pose & pose,
	Sizes frag_sizes_in = Sizes()
)
{
	using devel::blab::motif::has_motif_data;
	using devel::blab::motif::get_nonconst_motif_data;

	TR_REBUILD.Trace << "centroid_rebuild_loop: " << old_begin << ' ' << old_end << ' ' << newseq <<
		" old/new-lens: " << old_end-old_begin+1 << ' ' << newseq.size();
	for ( Size i=old_begin; i<= old_end; ++i ) {
		TR_REBUILD.Trace <<' ' << pose.residue(i).name();
	}
	TR_REBUILD.Trace << endl;
	Pose const pose_in( pose );

	string pymol_prefix;
	protocols::moves::PyMOLMoverOP pymol(0);

	if ( option[ my_options::show_pymol ] ) {
		static Size counter(0);
		++counter;
		pymol = protocols::moves::PyMOLMoverOP( new protocols::moves::PyMOLMover());
		pymol_prefix = "centroid_rebuild_loop_"+lead_zero_string_of(counter,2);
		pymol->keep_history(true);
	}

	Size loop_begin( old_begin ), loop_end( old_end );

	/// add loop residues
	while ( newseq.size() > loop_end-loop_begin+1 ) {
		ResidueOP rsd( get_vanilla_protein_residue( 'A' ) );
		Size const insertpos( (loop_begin+loop_end)/2 );
		pose.conformation().append_polymer_residue_after_seqpos( *rsd, insertpos, true );
		pose.set_secstruct( insertpos+1, 'L' );
		++loop_end;
		add_cheating_frags = false;
		if ( has_motif_data( pose ) ) get_nonconst_motif_data( pose ).insert_residue( insertpos );
	}

	/// delete loop residues
	while ( newseq.size() < loop_end-loop_begin+1 ) {
		Size const deletepos( (loop_begin+loop_end)/2 );
		pose.conformation().delete_residue_slow( deletepos );
		if ( has_motif_data( pose ) ) get_nonconst_motif_data( pose ).delete_residue( deletepos );
		--loop_end;
		add_cheating_frags = false;
	}

	runtime_assert( add_cheating_frags || add_sequence_frags ); // have to have some frags!

	if ( pymol ) {
		pymol->pymol_name(pymol_prefix+"_set_looplen");
		pymol->apply( pose );
	}


	// choose a cutpoint
	Size const cutpoint( random_range( loop_begin, loop_end-1 ) );

	{
		FoldTree f( pose.fold_tree() );
		f.new_jump( loop_begin-1, loop_end+1, cutpoint );
		f.reorder( f.root() );
		f.put_jump_stubs_intra_residue(); // could just do this for the new jump
		pose.fold_tree(f);
		add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, cutpoint );
		add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, cutpoint+1 );
	}
	if ( pymol ) {
		pymol->pymol_name(pymol_prefix+"_set_fold_tree");
		pymol->apply( pose );
	}


	// idealize the loop
	for ( Size i=loop_begin; i<= loop_end; ++i ) {
		conformation::idealize_position( i, pose.conformation() );
	}
	if ( pymol ) {
		pymol->pymol_name(pymol_prefix+"_idealize");
		pymol->apply( pose );
	}

	/// mutate to the correct sequence
	for ( Size i=loop_begin; i<= loop_end; ++i ) {
		make_sequence_change( i, aa_from_oneletter_code( newseq[i-loop_begin] ), pose );
	}

	/// initialize the torsions
	for ( Size i=loop_begin; i<= loop_end; ++i ) {
		pose.set_phi  (i, -150.0 );
		pose.set_psi  (i,  150.0 );
		pose.set_omega(i,  180.0 );
	}

	if ( pymol ) {
		pymol->pymol_name(pymol_prefix+"_extend");
		pymol->apply( pose );
	}

	runtime_assert( pose.sequence().substr( loop_begin-1, loop_end-loop_begin+1 ) == newseq );

	// loops classes
	protocols::loops::Loop loop( loop_begin, loop_end, cutpoint, 0.0 /*skip_rate*/, true /*extended*/ );
	protocols::loops::LoopsOP loops( new protocols::loops::Loops() );
	loops->add_loop( loop );

	// pick fragments
	devel::blab::classic_frags::FragLibOP fraglib( new devel::blab::classic_frags::FragLib() );
	if ( !dry_run() ) {
		kinematics::MoveMap mm;
		for ( Size i=loop_begin; i<= loop_end; ++i ) mm.set_bb( i, true );

		Sizes fragsizes;
		Size make_1mers_from(0);
		if ( frag_sizes_in.empty() ) {
			fragsizes = make_vector1( 3, 6 );
			make_1mers_from = 3;
		} else {
			runtime_assert( has_element( frag_sizes_in, Size(3) ) ); // because we use 3mers when inserting random frags
			fragsizes = frag_sizes_in;
			std::sort( fragsizes.begin(), fragsizes.end() );
			if ( fragsizes.front() == 1 ) {
				fragsizes.erase( fragsizes.begin() );
				runtime_assert( !fragsizes.empty() );
				make_1mers_from = fragsizes.front();
			}
		}

		Real const cheating_frags_ratio( option[ my_options::cheating_frags_ratio ] ); // default is 3.0
		Size const n_cheating_frags( add_sequence_frags? 100 : 100 ),
			n_sequence_frags( add_cheating_frags ? int( Real(n_cheating_frags)/cheating_frags_ratio) : 100 );
		Real const seq_weight( 1.0 ), ss_weight( 3.0 ), torsion_weight( 100.0 );

		Real const
			min_torsion_dev( option[ my_options::cheating_frags_min_torsion_dev ] ),
			max_torsion_dev( option[ my_options::cheating_frags_max_torsion_dev ] );


		if ( add_cheating_frags ) {
			fraglib =
				( devel::blab::classic_frags::setup_vall_cheating_fragments( fragsizes, n_cheating_frags, pose_in, mm,
				pose_in.secstruct(),
				seq_weight, ss_weight, torsion_weight,
				min_torsion_dev, max_torsion_dev ) );
		}

		if ( add_sequence_frags ) {
			TR_REBUILD.Trace << "add_sequence_frags: " << pose.secstruct().substr( loop_begin-1, loop_end-loop_begin+1 ) <<
				" n_frags: " << n_sequence_frags << ' ' << seq_weight << ' ' << ss_weight << " frag_sizes: ";
			for ( Size sz : fragsizes ) TR_REBUILD.Trace << ' '<< sz;
			TR_REBUILD.Trace << endl;

			devel::blab::classic_frags::add_vall_fragments( fragsizes, n_sequence_frags, pose, mm, pose.secstruct(),
				seq_weight, ss_weight, *fraglib );
		}

		if ( make_1mers_from ) {
			fraglib->library( 1 ).derive_from_src_lib( 1, make_1mers_from, fraglib->library( make_1mers_from ) );
		}
	}

	Pose const pose_before_loop_building( pose ); // still fullatom

	// switch to centroid
	devel::blab::classic_frags::switch_pose_to_residue_type_set( pose, CENTROID );
	ScoreFunctionOP cen_scorefxn( protocols::loops::get_cen_scorefxn() );

	// loop model
	if ( !dry_run() ) {
		bool const use_all_frags( true );
		devel::blab::classic_frags::chu_perturb_one_loop( pose, loop, fraglib, *cen_scorefxn, use_all_frags );
		if ( pymol ) {
			pymol->pymol_name(pymol_prefix+"_after_loopbuild");
			pymol->apply( pose );
		}
	}

	// switch to fullatom
	devel::blab::classic_frags::switch_pose_to_residue_type_set( pose, FA_STANDARD );

	// recover sidechains except in loop region
	devel::blab::classic_frags::recover_template_sidechains( pose_before_loop_building, *loops, pose );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
centroid_rebuild_loop_using_KIC(
	Size const old_begin,
	Size const old_end,
	string const & newseq,
	Pose & pose
)
{
	using devel::blab::motif::has_motif_data;
	using devel::blab::motif::get_nonconst_motif_data;

	TR_REBUILD.Trace << "centroid_rebuild_loop: " << old_begin << ' ' << old_end << ' ' << newseq <<
		" old/new-lens: " << old_end-old_begin+1 << ' ' << newseq.size();
	for ( Size i=old_begin; i<= old_end; ++i ) {
		TR_REBUILD.Trace <<' ' << pose.residue(i).name();
	}
	TR_REBUILD.Trace << endl;
	Pose const pose_in( pose );

	string pymol_prefix;
	protocols::moves::PyMOLMoverOP pymol(0);

	if ( option[ my_options::show_pymol ] ) {
		static Size counter(0);
		++counter;
		pymol = protocols::moves::PyMOLMoverOP( new protocols::moves::PyMOLMover());
		pymol_prefix = "centroid_rebuild_loop_"+lead_zero_string_of(counter,2);
		pymol->keep_history(true);
	}

	Size loop_begin( old_begin ), loop_end( old_end );

	/// add loop residues
	while ( newseq.size() > loop_end-loop_begin+1 ) {
		ResidueOP rsd( get_vanilla_protein_residue( 'A' ) );
		Size const insertpos( (loop_begin+loop_end)/2 );
		pose.conformation().append_polymer_residue_after_seqpos( *rsd, insertpos, true );
		pose.set_secstruct( insertpos+1, 'L' );
		++loop_end;
		if ( has_motif_data( pose ) ) get_nonconst_motif_data( pose ).insert_residue( insertpos );
	}

	/// delete loop residues
	while ( newseq.size() < loop_end-loop_begin+1 ) {
		Size const deletepos( (loop_begin+loop_end)/2 );
		pose.conformation().delete_residue_slow( deletepos );
		if ( has_motif_data( pose ) ) get_nonconst_motif_data( pose ).delete_residue( deletepos );
		--loop_end;
	}

	if ( pymol ) {
		pymol->pymol_name(pymol_prefix+"_set_looplen");
		pymol->apply( pose );
	}


	// choose a cutpoint
	Size const cutpoint( random_range( loop_begin, loop_end-1 ) );

	{
		FoldTree f( pose.fold_tree() );
		f.new_jump( loop_begin-1, loop_end+1, cutpoint );
		f.reorder( f.root() );
		f.put_jump_stubs_intra_residue(); // could just do this for the new jump
		pose.fold_tree(f);
		add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, cutpoint );
		add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, cutpoint+1 );
	}
	if ( pymol ) {
		pymol->pymol_name(pymol_prefix+"_set_fold_tree");
		pymol->apply( pose );
	}


	// idealize the loop
	for ( Size i=loop_begin; i<= loop_end; ++i ) {
		conformation::idealize_position( i, pose.conformation() );
	}
	if ( pymol ) {
		pymol->pymol_name(pymol_prefix+"_idealize");
		pymol->apply( pose );
	}

	/// mutate to the correct sequence
	for ( Size i=loop_begin; i<= loop_end; ++i ) {
		make_sequence_change( i, aa_from_oneletter_code( newseq[i-loop_begin] ), pose );
	}

	/// initialize the torsions
	for ( Size i=loop_begin; i<= loop_end; ++i ) {
		pose.set_phi  (i, -150.0 );
		pose.set_psi  (i,  150.0 );
		pose.set_omega(i,  180.0 );
	}

	if ( pymol ) {
		pymol->pymol_name(pymol_prefix+"_extend");
		pymol->apply( pose );
	}

	runtime_assert( pose.sequence().substr( loop_begin-1, loop_end-loop_begin+1 ) == newseq );

	// loops classes
	protocols::loops::Loop loop( loop_begin, loop_end, cutpoint, 0.0 /*skip_rate*/, true /*extended*/ );
	protocols::loops::LoopsOP loops( new protocols::loops::Loops() );
	loops->add_loop( loop );

	Pose const pose_before_loop_building( pose ); // still fullatom

	// switch to centroid
	devel::blab::classic_frags::switch_pose_to_residue_type_set( pose, CENTROID );
	ScoreFunctionOP cen_scorefxn( protocols::loops::get_cen_scorefxn() );

	// loop model
	if ( !dry_run() ) {
		bool const oldval( option[OptionKeys::loops::fast] );
		option[OptionKeys::loops::fast].value( true ); // hack!
		runtime_assert( option[OptionKeys::loops::fast] );
		protocols::loops::loop_mover::perturb::LoopMover_Perturb_KIC loopmover( loops, cen_scorefxn );
		loopmover.apply( pose );
		if ( pymol ) {
			pymol->pymol_name(pymol_prefix+"_after_loopbuild");
			pymol->apply( pose );
		}
		option[OptionKeys::loops::fast].value( oldval ); // unhack!
	}

	// switch to fullatom
	devel::blab::classic_frags::switch_pose_to_residue_type_set( pose, FA_STANDARD );

	// recover sidechains except in loop region
	devel::blab::classic_frags::recover_template_sidechains( pose_before_loop_building, *loops, pose );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





// ResidueOP
// get_vanilla_protein_residue( char const name1 )
// {
//  ResidueTypeSetCAP rsd_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
//  ResidueTypeCAP new_rsd_type( ResidueSelector().set_name1( name1 ).exclude_variants().select( *rsd_set )[1] );
//  pbassert( new_rsd_type->is_protein() );
//  return ResidueFactory::create_residue( *new_rsd_type );
// }

// ResidueOP
// get_vanilla_dna_residue( char const name1 )
// {
//  ResidueTypeSetCAP rsd_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
//  ResidueTypeCAP new_rsd_type( ResidueSelector().set_name1( name1 ).exclude_variants().select( *rsd_set )[1] );
//  pbassert( new_rsd_type->is_DNA() );
//  return ResidueFactory::create_residue( *new_rsd_type );
// }

// ResidueOP
// get_vanilla_dna_residue_centroid( char const name1 )
// {
//  ResidueTypeSetCAP rsd_set( ChemicalManager::get_instance()->residue_type_set( CENTROID ) );
//  ResidueTypeCAP new_rsd_type( ResidueSelector().set_name1( name1 ).exclude_variants().select( *rsd_set )[1] );
//  pbassert( new_rsd_type->is_DNA() );
//  return ResidueFactory::create_residue( *new_rsd_type );
// }

// ResidueOP
// get_vanilla_protein_residue_centroid( char const name1 )
// {
//  ResidueTypeSetCAP rsd_set( ChemicalManager::get_instance()->residue_type_set( CENTROID ) );
//  ResidueTypeCAP new_rsd_type( ResidueSelector().set_name1( name1 ).exclude_variants().select( *rsd_set )[1] );
//  return ResidueFactory::create_residue( *new_rsd_type );
// }
// static Real const init_phi( -150.0 ), init_psi( 150.0 ), init_omega( 180.0 );

// ///////////////////////////////////////////////////////////////////////////////////////
// bool
// init_torsions_still_present(
//               bools const & is_flexible,
//               Pose const & pose
//               )
// {
//  using util::subtract_degree_angles;

//  bool done( true );
//  for ( Size i=1; i<= pose.total_residue(); ++i ) {
//   if ( is_flexible[i] && pose.residue(i).is_protein() ) {
//    Residue const & rsd( pose.residue(i) );
//    if ( ( abs( subtract_degree_angles( pose.phi  (i), init_phi   ) ) < 1e-3 || rsd.is_lower_terminus() ) &&
//       ( abs( subtract_degree_angles( pose.psi  (i), init_psi   ) ) < 1e-3 || rsd.is_upper_terminus() ) &&
//       ( abs( subtract_degree_angles( pose.omega(i), init_omega ) ) < 1e-3 || rsd.is_upper_terminus() ) ) {
//     done = false;
//     break;
//    }
//   }
//  }
//  return !done;
// }


///////////////////////////////////////////////////////////////////////////////////////

// using namespace devel::dna;
// using namespace protocols::moves;
// using namespace id;

// static vector1< DoubleFragmentCOP > double_frags;

// Reals
// optimize_newly_added_base_pair_geometry(ScoreFunction & scorefxn, Pose & pose, Size const cycles,
//                     Real const gamma, Real & temperature,
//                     Size const pos1, Size const pos2)
// {
//  if ( double_frags.empty() ) {
//   Size const fragsize( 3 );
//   DoubleFragmentVall const vall
//    ( core::options::option[ core::options::OptionKeys::dna::specificity::double_frag_vall ] );

//   vall.pick_fragments( fragsize, "xxx", vall.nfrags( fragsize ), double_frags );

//   TR_REBUILD.Trace << "Got " << double_frags.size() << " double_frags" << endl;
//  }

//  MonteCarloOP mc( new MonteCarlo( pose, scorefxn, temperature ) );
//  Real beginning_score;
//  for ( Size n=1; n<= cycles; ++n ) {
//   temperature *= gamma;
//   mc->set_temperature( temperature );

//   /// insert randomly chosen fragment torsion angles
//   DoubleFragment const & frag( *numeric::random::random_element( double_frags ) );

//   //Real random_unif_num = numeric::random::uniform();
//   //TR_REBUILD.Trace << "Randome uniform number: " << random_unif_num << endl;
//   bool move_strand1( false ), move_strand2( false );
//   bool const move_both_strands( numeric::random::uniform() < 0.5 );
//   if ( move_both_strands ) move_strand1 = move_strand2 = true;
//   else {
//    if ( numeric::random::uniform() < 0.5 ) move_strand1 = true;
//    else move_strand2 = true;
//   }
//   if ( move_strand1 ) {
//    // this does pos1-1 torsion 5 and 6 as well
//    set_frag_pos_torsions( pos1, pose.conformation(), frag.torsions()[1][2] );
//   }
//   if ( move_strand2 ) {
//    set_frag_pos_torsions( pos2, pose.conformation(), frag.torsions()[1][2] );
//    // also want to do torsions 5 and 6, and pos2+1 torsion 1
//    pose.conformation().set_torsion( TorsionID( pos2  , BB, 5 ), frag.torsions()[1][3][1] ); // epsilon
//    pose.conformation().set_torsion( TorsionID( pos2  , BB, 6 ), frag.torsions()[1][3][2] ); // epsilon
//    pose.conformation().set_torsion( TorsionID( pos2+1, BB, 1 ), frag.torsions()[1][3][3] ); // alpha (pos2+1)
//   }

//   Real const score( scorefxn( pose ) ); // superfluous score call for debugging
//   bool const pass_mc( mc->boltzmann( pose ) );

//   if ( n == 1 )
//    beginning_score = mc->last_accepted_score();

//   TR_REBUILD.Trace << "MC " << I(6,n) << " mc_temp: " << F(9,3,mc->temperature() ) <<
//    " pass: " << pass_mc << ' ' << F(9,3,score) <<
//    F(9,3,mc->last_accepted_score()) << F(9,3,mc->lowest_score() ) << endl;
//  }
//  mc->show_counters( TR_REBUILD.Trace, false );
//  mc->recover_low( pose );
//  Reals scores;
//  scores.add_back(beginning_score);
//  scores.add_back(mc->lowest_score());
//  return scores;
// }


// /**
// Extend DNA duplex to be longer, based on DNA fragments
// **/

// void
// extend_pose_dna_one_bp_using_fragment_torsions(
//                         Size const strand1_chain,
//                         AA const dna_aa,
//                         Pose & pose
//                         )
// {
//  using namespace devel::dna;
//  using namespace protocols::moves;
//  using namespace id;

//  if ( double_frags.empty() ) {
//   Size const fragsize( 3 );
//   DoubleFragmentVall const vall
//    ( core::options::option[ core::options::OptionKeys::dna::specificity::double_frag_vall ] );

//   vall.pick_fragments( fragsize, "xxx", vall.nfrags( fragsize ), double_frags );

//   TR_REBUILD.Trace << "Got " << double_frags.size() << " double_frags" << endl;
//  }

//  bool const dry_run( option[ OK::out::dry_run ] || option[ OK::run::dry_run ] );
//  Size const cycles( dry_run ? 25 : 2500 );

//  Size pairing_base_pos(retrieve_base_partner_from_pose( pose )[ pose.chain_begin(strand1_chain ) ] );
//  if (pairing_base_pos == 0) { // move this check up or we get segfault/assert-fail in vector lookup
//   utility_exit_with_message("Base pairing information in error");
//  }
//  Size const strand2_chain( pose.chain( pairing_base_pos ));
//  TR_REBUILD.Trace << "Strand 1: " << strand1_chain << " Strand 2: " << strand2_chain << endl;

//  pbassert(strand1_chain != strand2_chain); // should be two chains

//  ScoreFunctionOP scorefxn( new ScoreFunction() );
//  scorefxn->set_energy_method_options
//   ( scoring::methods::EnergyMethodOptions().exclude_DNA_DNA( false ) );

//  scorefxn->set_weight( dna_bp, 0.01 );
//  scorefxn->set_weight( dna_bs, 0.01 );

//  TR_REBUILD.Trace << "Pose size is " << pose.total_residue() << endl;


//  // add the new bases
//  remove_upper_terminus_type_from_pose_residue( pose, pose.chain_end  (strand1_chain) );
//  remove_lower_terminus_type_from_pose_residue( pose, pose.chain_begin(strand2_chain) );
//  ResidueOP rsd1( get_vanilla_dna_residue( 'a' ) ), rsd2( get_vanilla_dna_residue( 't' ) );
//  retrieve_nonconst_base_partner_from_pose( pose ).insert_residue_after_seqpos( pose.chain_end( strand1_chain ) );
//  pose.append_polymer_residue_after_seqpos( *rsd1, pose.chain_end(strand1_chain), true );
//  retrieve_nonconst_base_partner_from_pose( pose ).insert_residue_after_seqpos( pose.chain_begin( strand2_chain )-1 );
//  pose.prepend_polymer_residue_before_seqpos( *rsd2, pose.chain_begin(strand2_chain), true );
//  add_upper_terminus_type_to_pose_residue( pose, pose.chain_end  ( strand1_chain ) );
//  add_lower_terminus_type_to_pose_residue( pose, pose.chain_begin( strand2_chain ) );

//  Size const pos1( pose.chain_end( strand1_chain ) ), pos2( pose.chain_begin( strand2_chain ) );

//  // manually set A and T to be base-paired
//  BasePartner & partner( retrieve_nonconst_base_partner_from_pose( pose ) );
//  partner[pos1] = pos2;
//  partner[pos2] = pos1;
//  //pose.data().set( util::BASE_PARTNER, new BasePartner( partner ) );

//  make_base_pair_mutation( pose, pos1, dna_aa );

//  TR_REBUILD.Trace << "Extended pose size is " << pose.total_residue() << endl;

//  // update base partner to declare that these are paired
//   //TR_REBUILD.Trace << "Pos1 " << pos1 << " pairs with " << retrieve_base_partner_from_pose( pose )[pos1] << endl
//  //  << "Pos2 " << pos2 << " pairs with " << retrieve_base_partner_from_pose( pose )[pos2] << endl;
//  /// debugging
//  BasePartner const & partner_full( retrieve_base_partner_from_pose( pose ) );
//  for ( Size i=1; i<= pose.total_residue(); ++i ) {
//   Size base_partner( partner_full[i] );
//   if ( pose.residue(i).is_DNA() && base_partner > i ) {
//    pbassert( partner_full[ base_partner ] == i );
//    TR_REBUILD.Trace << "DNA " << pose.chain(i) << " "
//      << i << " " << pose.residue(i).name3() << " "
//      << pose.chain( base_partner ) << " "
//      << base_partner << " "
//      << pose.residue( base_partner ).name3() << endl;
//   }
//  }

//  pbassert( retrieve_base_partner_from_pose( pose )[pos1] == pos2 );
//  pbassert( retrieve_base_partner_from_pose( pose )[pos2] == pos1 );

//  Real const init_temp( 2.0 ), final_temp( 0.5 );
//  Real const gamma = std::pow( (final_temp/init_temp), Real(1.0/( cycles-1 ) ) );
//  Real temperature( init_temp / gamma );

//  Reals scores;
//  scorefxn->set_weight( dna_bp, 0.00 );
//  scorefxn->set_weight( dna_bs, 0.01 );
//  scores = optimize_newly_added_base_pair_geometry(*scorefxn, pose, cycles,
//                           gamma, temperature,
//                           pos1, pos2);

//  temperature = init_temp / gamma; // reset temperature
//  scorefxn->set_weight( dna_bp, 0.01 );
//  scorefxn->set_weight( dna_bs, 0.01 );
//  scores = optimize_newly_added_base_pair_geometry(*scorefxn, pose, cycles,
//                           gamma, temperature,
//                           pos1, pos2);

//  if ( scores.size() == 2 && (scores[2] / scores[1]) > 0.8 ) {
//   TR_REBUILD.Trace << "WARNING: Potentially bad take-off angle. Check structure to verify. " << endl;

//   // try a new take-off angle?
//  }


//  // debugging
//  string outputname( "temp_extended_" + string_of(pos1) + "_" + string_of(pos2) + ".pdb" );
//  dump_complete_motif_pdb( pose, outputname);

// }



#endif

