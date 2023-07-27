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

#ifndef INCLUDED_apps_pilot_phil_interface
#define INCLUDED_apps_pilot_phil_interface


#include <apps/pilot/phil/phil.hh>
// #include <apps/pilot/phil/phil_core.hh>
// #include <apps/pilot/phil/phil_protocols.hh>
#include <apps/pilot/phil/sasa.hh>

#include <apps/pilot/phil/phil_io.hh>
#include <apps/pilot/phil/phil_options.hh>

#include <utility/graph/Graph.hh>


#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.hh>

#include <devel/dna/util.hh>

#include <sstream>

static basic::Tracer TRI( "apps.pilot.phil.interface" );


bool // returnval: did we find any chains to *include* ?
setup_exclude_bools_from_chains_list(
	Pose const & pose, // with pdb info, hopefully
	string const & filename,
	string const & listfile,
	bools & exclude
)
{
	using namespace std;
	string const my_pdb_id( ObjexxFCL::uppercased( filebase( filename ).substr( 0, 4 ) ) );

	vector1< char > chains;

	ifstream data( listfile.c_str() );
	string line,tag;
	while ( getline( data, line ) ) {
		istringstream l( line );
		l >> tag;
		if ( tag.size() == 5 && ObjexxFCL::uppercased( tag.substr(0,4) ) == my_pdb_id ) {
			TRI.Trace << "Found chain for my pdb id: " << tag << endl;
			chains.push_back( tag[4] );
		}
	}

	if ( chains.empty() ) return false;

	exclude.clear(); exclude.resize( pose.total_residue(), false );

	Size included(0), excluded(0);
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).is_protein() ) {
			if ( !has_element( chains, pose.pdb_info()->chain( i ) ) ) {
				TRI.Trace << "Excluding position: " << i << ' ' << pose.pdb_info()->chain(i) << endl;
				exclude[i]= true;
				++excluded;
			} else ++included;
		}
	}

	TRI.Trace << "included: " << included << " excluded: " << excluded << ' ' << filename << endl;


	return ( included > 0 );

}

bool
protein_position_could_interact_with_dna_base_pair(
	Size const seqpos,
	Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	Size const dnapos,
	bool const dna_sidechain_only = false,
	bool const fix_dna_sequence = false
)
{
	using namespace pack;
	using namespace pack::task;
	using namespace pack::rotamer_set;

	Real const heavyatom_distance_threshold( 6.5 ), // allow a little slop in rotamers... do we need more?
		heavyatom_dis2_threshold( heavyatom_distance_threshold * heavyatom_distance_threshold );

	Real const dna_nbr_radius( 7.2 );

	Size const dnapos_partner( retrieve_base_partner_from_pose( pose )[ dnapos ] );

	if ( pose.residue( seqpos ).nbr_atom_xyz().distance_squared( pose.residue( dnapos ).nbr_atom_xyz() ) >
			numeric::square( pose.residue( seqpos ).nbr_radius() + dna_nbr_radius + heavyatom_distance_threshold ) &&
			( !dnapos_partner ||
			( pose.residue( seqpos ).nbr_atom_xyz().distance_squared( pose.residue( dnapos_partner ).nbr_atom_xyz() ) >
			numeric::square( pose.residue( seqpos ).nbr_radius() + dna_nbr_radius + heavyatom_distance_threshold ) ) ) ) {
		return false;
	}


	PackerTaskOP ptask( TaskFactory::create_packer_task( pose ) );
	ptask->set_bump_check( false );
	//ptask->temporarily_set_pack_residue( seqpos, true );
	ptask->initialize_from_command_line();
	ptask->or_include_current( true );

	// just repacking rotamers at the protein position
	ptask->nonconst_residue_task( seqpos ).restrict_to_repacking();

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( pose.residue(ii).is_DNA() ) {
			if ( fix_dna_sequence ) {
				ptask->nonconst_residue_task( ii ).restrict_to_repacking();
			} else {
				ptask->nonconst_residue_task( ii ).allow_aa( na_ade );
				ptask->nonconst_residue_task( ii ).allow_aa( na_thy );
				ptask->nonconst_residue_task( ii ).allow_aa( na_gua );
				ptask->nonconst_residue_task( ii ).allow_aa( na_cyt );
				assert( ptask->design_residue(ii) );
			}
		}
	}

	// unnecessary here, yet also required
	utility::graph::GraphOP dummygraph( new utility::graph::Graph( pose.total_residue() ) );

	RotamerSetFactory rsf;
	RotamerSetOP rotset( rsf.create_rotamer_set( pose ) ),
		rotset_strand1( rsf.create_rotamer_set( pose ) ),
		rotset_strand2;
	rotset->set_resid( seqpos );
	rotset->build_rotamers( pose, scorefxn, *ptask, dummygraph, false );

	rotset_strand1->set_resid( dnapos );
	rotset_strand1->build_rotamers( pose, scorefxn, *ptask, dummygraph, false );

	if ( dnapos_partner ) {
		rotset_strand2 = rsf.create_rotamer_set( pose );
		rotset_strand2->set_resid( dnapos_partner );
		rotset_strand2->build_rotamers( pose, scorefxn, *ptask, dummygraph, false );
	}

	// TRI.Trace << "rotamer_contact_screen: n_protein_rots= " << rotset->num_rotamers() <<
	//  " n_strand1_rots: " << rotset_strand1->num_rotamers() <<
	//  " n_strand2_rots: " << ( dnapos_partner ? rotset_strand2->num_rotamers() : 0 ) << std::endl;

	for ( Rotamers::const_iterator p_rotamer( rotset->begin() ); p_rotamer != rotset->end(); ++p_rotamer ) {
		Residue const & rsd1( **p_rotamer );
		assert( rsd1.aa() == pose.residue( seqpos ).aa() );

		for ( Size strand=1; strand<=2; ++strand ) {
			if ( strand == 2 && !dnapos_partner ) continue;
			RotamerSetOP const & rotset_dna( strand==1 ? rotset_strand1 : rotset_strand2 );
			for ( Rotamers::const_iterator d_rotamer= rotset_dna->begin(); d_rotamer != rotset_dna->end(); ++d_rotamer ) {
				Residue const & rsd2( **d_rotamer );
				if ( rsd1.nbr_atom_xyz().distance_squared( rsd2.nbr_atom_xyz() ) >
						numeric::square( rsd1.nbr_radius() + rsd2.nbr_radius() + heavyatom_distance_threshold ) ) continue;

				Size jj_begin( 1 );
				if ( dna_sidechain_only ) jj_begin = rsd2.first_sidechain_atom();

				for ( Size ii=1; ii<= rsd1.nheavyatoms(); ++ii ) {
					for ( Size jj=jj_begin; jj<= rsd2.nheavyatoms(); ++jj ) {
						if ( rsd1.xyz( ii ).distance_squared( rsd2.xyz( jj ) ) < heavyatom_dis2_threshold ) {
							// TRI.Trace << "found_contact: protein_pos= " << seqpos << " strand= " << strand << " dnapos= " <<
							//  rsd2.seqpos() << " distance= " << rsd1.xyz(ii).distance( rsd2.xyz(jj) ) << endl;
							return true; // IN CONTACT !!!!!!!!
						}
					}
				}
			}
		}
	}

	return false;
}


/////////////////// for unanchored waters
bool
unanchored_water_position_could_interact_with_dna_base_pair(
	Vector const & wat_center,
	Real const wat_radius,
	Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	Size const dnapos,
	bool const dna_sidechain_only = false,
	bool const fix_dna_sequence = false,
	Real const heavyatom_distance_threshold = 6.5
)
{
	using namespace pack;
	using namespace pack::task;
	using namespace pack::rotamer_set;

	Real const heavyatom_dis2_threshold( numeric::square( heavyatom_distance_threshold + wat_radius ) );
	Real const dna_nbr_radius( 7.2 );
	Size const dnapos_partner( retrieve_base_partner_from_pose( pose )[ dnapos ] );

	if ( wat_center.distance_squared( pose.residue( dnapos ).nbr_atom_xyz() ) >
			numeric::square( wat_radius + dna_nbr_radius + heavyatom_distance_threshold ) &&
			( !dnapos_partner ||
			( wat_center.distance_squared( pose.residue( dnapos_partner ).nbr_atom_xyz() ) >
			numeric::square( wat_radius + dna_nbr_radius + heavyatom_distance_threshold ) ) ) ) {
		return false;
	}


	PackerTaskOP ptask( TaskFactory::create_packer_task( pose ) );
	ptask->set_bump_check( false );
	ptask->initialize_from_command_line();
	ptask->or_include_current( true );

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( pose.residue(ii).is_DNA() ) {
			if ( fix_dna_sequence ) {
				ptask->nonconst_residue_task( ii ).restrict_to_repacking();
			} else {
				ptask->nonconst_residue_task( ii ).allow_aa( na_ade );
				ptask->nonconst_residue_task( ii ).allow_aa( na_thy );
				ptask->nonconst_residue_task( ii ).allow_aa( na_gua );
				ptask->nonconst_residue_task( ii ).allow_aa( na_cyt );
				assert( ptask->design_residue(ii) );
			}
		}
	}

	// unnecessary here, yet also required
	utility::graph::GraphOP dummygraph( new utility::graph::Graph( pose.total_residue() ) );

	RotamerSetFactory rsf;
	RotamerSetOP rotset_strand1( rsf.create_rotamer_set( pose ) ),
		rotset_strand2;

	rotset_strand1->set_resid( dnapos );
	rotset_strand1->build_rotamers( pose, scorefxn, *ptask, dummygraph, false );

	if ( dnapos_partner ) {
		rotset_strand2 = rsf.create_rotamer_set( pose );
		rotset_strand2->set_resid( dnapos_partner );
		rotset_strand2->build_rotamers( pose, scorefxn, *ptask, dummygraph, false );
	}

	for ( Size strand=1; strand<=2; ++strand ) {
		if ( strand == 2 && !dnapos_partner ) continue;
		RotamerSetOP const & rotset_dna( strand==1 ? rotset_strand1 : rotset_strand2 );
		for ( Rotamers::const_iterator d_rotamer= rotset_dna->begin(); d_rotamer != rotset_dna->end(); ++d_rotamer ) {
			Residue const & rsd2( **d_rotamer );
			if ( wat_center.distance_squared( rsd2.nbr_atom_xyz() ) >
					numeric::square( wat_radius + rsd2.nbr_radius() + heavyatom_distance_threshold ) ) continue;

			Size jj_begin( 1 );
			if ( dna_sidechain_only ) jj_begin = rsd2.first_sidechain_atom();

			for ( Size jj=jj_begin; jj<= rsd2.nheavyatoms(); ++jj ) {
				if ( wat_center.distance_squared( rsd2.xyz( jj ) ) < heavyatom_dis2_threshold ) return true; // contact
			}
		}
	}

	return false;
}
bool
unanchored_water_position_could_interact_with_dna_base_pair(
	Size const wpos,
	Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	Size const dnapos,
	bool const dna_sidechain_only = false,
	bool const fix_dna_sequence = false,
	Real const heavyatom_distance_threshold = 6.5
)
{
	using namespace pack::rotamer_set;
	runtime_assert( has_water_packing_info( pose ) && get_water_packing_info( pose ).has( wpos ) );
	WaterAnchorInfo const & watinfo( get_water_packing_info( pose )[wpos] );
	runtime_assert( watinfo.anchor_residue() == 0 ); // only unanchored posns right now...

	return unanchored_water_position_could_interact_with_dna_base_pair(
		watinfo.anchor_sphere_center(), watinfo.anchor_sphere_radius(), pose, scorefxn, dnapos, dna_sidechain_only,
		fix_dna_sequence, heavyatom_distance_threshold );
}
bool
protein_positions_could_interact(
	Size const pos1,
	bool const design1,
	bool const sc_only1,
	Size const pos2,
	bool const design2,
	bool const sc_only2,
	Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	Real const heavyatom_distance_threshold
)
{
	using namespace pack;
	using namespace pack::task;
	using namespace pack::rotamer_set;

	Real const heavyatom_dis2_threshold( heavyatom_distance_threshold * heavyatom_distance_threshold );
	Real const max_protein_rsd_radius( 6.1209 ); // ARG
	Real const radius1( design1 ? max_protein_rsd_radius : pose.residue( pos1 ).nbr_radius() );
	Real const radius2( design2 ? max_protein_rsd_radius : pose.residue( pos2 ).nbr_radius() );

	if ( pose.residue( pos1 ).nbr_atom_xyz().distance_squared( pose.residue( pos2 ).nbr_atom_xyz() ) >
			numeric::square( radius1 + radius2 + heavyatom_distance_threshold ) ) {
		return false;
	}

	PackerTaskOP ptask( TaskFactory::create_packer_task( pose ) );
	ptask->set_bump_check( false );
	//ptask->temporarily_set_pack_residue( seqpos, true );
	// ptask->initialize_from_command_line();
	ptask->or_include_current( true );

	// just repacking rotamers at the protein position
	if ( !design1 ) ptask->nonconst_residue_task( pos1 ).restrict_to_repacking();
	if ( !design2 ) ptask->nonconst_residue_task( pos2 ).restrict_to_repacking();


	// unnecessary here, yet also required
	utility::graph::GraphOP dummygraph( new utility::graph::Graph( pose.total_residue() ) );

	RotamerSetFactory rsf;
	RotamerSetOP
		rotset1( rsf.create_rotamer_set( pose ) ),
		rotset2( rsf.create_rotamer_set( pose ) );
	rotset1->set_resid( pos1 );
	rotset2->set_resid( pos2 );
	rotset1->build_rotamers( pose, scorefxn, *ptask, dummygraph, false );
	rotset2->build_rotamers( pose, scorefxn, *ptask, dummygraph, false );

	// TRI.Trace << "rotamer_contact_screen: n_pos1_rots= " << rotset1->num_rotamers() <<
	//  " n_pos2_rots: " << rotset2->num_rotamers() << endl;

	for ( Rotamers::const_iterator rot1( rotset1->begin() ); rot1 != rotset1->end(); ++rot1 ) {
		Residue const & rsd1( **rot1 );
		runtime_assert( design1 || rsd1.aa() == pose.residue( pos1 ).aa() );

		for ( Rotamers::const_iterator rot2( rotset2->begin() ); rot2 != rotset2->end(); ++rot2 ) {
			Residue const & rsd2( **rot2 );
			runtime_assert( design2 || rsd2.aa() == pose.residue( pos2 ).aa() );

			if ( rsd1.nbr_atom_xyz().distance_squared( rsd2.nbr_atom_xyz() ) >
					numeric::square( rsd1.nbr_radius() + rsd2.nbr_radius() + heavyatom_distance_threshold ) ) continue;

			for ( Size ii=1; ii<= rsd1.nheavyatoms(); ++ii ) {
				if ( sc_only1 && rsd1.atom_is_backbone( ii ) ) continue;
				for ( Size jj=1; jj<= rsd2.nheavyatoms(); ++jj ) {
					if ( sc_only2 && rsd2.atom_is_backbone( jj ) ) continue;
					if ( rsd1.xyz( ii ).distance_squared( rsd2.xyz( jj ) ) < heavyatom_dis2_threshold ) {
						// TRI.Trace << "found_contact: pos1= " << pos1 << " pos2= " << pos2 <<
						//  " distance= " << rsd1.xyz(ii).distance( rsd2.xyz(jj) ) << endl;
						return true; // IN CONTACT !!!!!!!!
					}
				}
			}
		}
	}

	return false;
}






///////////////////////////////////////////////////////////////////////////////////////
///
/// this assumes protein backbone and water's anchor_sphere_center stay fixed
///
///
bool
protein_position_could_interact_with_water_position(
	Size const ppos,
	Vector const & wat_center,
	Real const wat_radius,
	Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	Real const heavyatom_distance_threshold
)
{
	using namespace pack;
	using namespace pack::task;
	using namespace pack::rotamer_set;

	Real const dis2_threshold( numeric::square( heavyatom_distance_threshold + wat_radius ) );
	Real const max_protein_rsd_radius( 6.1209 ); // ARG
	bool const designing( false );
	if ( designing ) { runtime_assert( !pose.residue( ppos ).is_DNA() ); } // otherwise radius is messed up
	Real const pradius( designing ? max_protein_rsd_radius : pose.residue( ppos ).nbr_radius() );

	if ( pose.residue( ppos ).nbr_atom_xyz().distance_squared( wat_center ) >
			numeric::square( pradius + wat_radius + heavyatom_distance_threshold ) ) {
		return false;
	}

	PackerTaskOP ptask( TaskFactory::create_packer_task( pose ) );
	ptask->set_bump_check( false );
	//ptask->temporarily_set_pack_residue( seqpos, true );
	// ptask->initialize_from_command_line();
	ptask->or_include_current( true );

	// just repacking rotamers at the protein position
	if ( !designing ) ptask->nonconst_residue_task( ppos ).restrict_to_repacking();


	// unnecessary here, yet also required
	utility::graph::GraphOP dummygraph( new utility::graph::Graph( pose.total_residue() ) );

	RotamerSetFactory rsf;
	RotamerSetOP rotset( rsf.create_rotamer_set( pose ) );
	rotset->set_resid( ppos );
	rotset->build_rotamers( pose, scorefxn, *ptask, dummygraph, false );

	// TRI.Trace << "rotamer_contact_screen: n_ppos_rots= " << rotset->num_rotamers() << endl;

	for ( Rotamers::const_iterator rot( rotset->begin() ); rot != rotset->end(); ++rot ) {
		Residue const & prsd( **rot );
		runtime_assert( designing || prsd.aa() == pose.residue( ppos ).aa() );

		for ( Size ii=1; ii<= prsd.nheavyatoms(); ++ii ) {
			if ( prsd.xyz( ii ).distance_squared( wat_center ) <= dis2_threshold ) {
				return true; // IN CONTACT !!!!!!!!
			}
		}
	}

	return false;
}


bool
protein_position_could_interact_with_water_position(
	Size const ppos,
	Size const wpos,
	Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	Real const heavyatom_distance_threshold
)
{
	using namespace pack;
	using namespace pack::task;
	using namespace pack::rotamer_set;

	runtime_assert( has_water_packing_info( pose ) && get_water_packing_info( pose ).has( wpos ) );
	WaterAnchorInfo const & watinfo( get_water_packing_info( pose )[wpos] );
	runtime_assert( watinfo.anchor_residue() == 0 ); // only unanchored posns right now...

	return protein_position_could_interact_with_water_position( ppos,
		watinfo.anchor_sphere_center(), watinfo.anchor_sphere_radius(), pose, scorefxn, heavyatom_distance_threshold );
}





///////////////////////////////////////////////////////////////////////////////////////
bools
get_protein_positions_potentially_interacting_with_base_pair(
	Size const dnapos,
	Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	bool const dna_sidechain_only = false
)
{
	bools is_interacting( pose.total_residue(), false );
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).is_protein() ) {
			is_interacting[ i ] = protein_position_could_interact_with_dna_base_pair( i, pose, scorefxn, dnapos,
				dna_sidechain_only );
		}
	}

	return is_interacting;
}

///////////////////////////////////////////////////////////////////////////////////////
/// both protein and dna positions are included in the interface
///

bools
get_protein_dna_interface_by_rotamer_contacts(
	Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	Sizes & contact_counts_for_dna_by_basepair
)
{
	bools is_interface( pose.total_residue(), false );
	contact_counts_for_dna_by_basepair.clear(); contact_counts_for_dna_by_basepair.resize( pose.total_residue(), 0 );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( !pose.residue(i).is_DNA() ) continue;

		Size const ppos( retrieve_base_partner_from_pose( pose )[i] );
		if ( ppos && ppos < i ) continue; // only do once for each basepair

		for ( Size j=1; j<= pose.total_residue(); ++j ) {
			if ( !pose.residue(j).is_protein() ) continue;

			if ( protein_position_could_interact_with_dna_base_pair( j, pose, scorefxn, i ) ) {
				is_interface[j] = true;
				is_interface[i] = true;
				++contact_counts_for_dna_by_basepair[i];
				if ( ppos ) {
					is_interface[ ppos ] = true;
					++contact_counts_for_dna_by_basepair[ ppos ];
				}
			}
		}
	}

	return is_interface;
}


bools
get_protein_dna_interface_by_rotamer_contacts(
	Pose const & pose,
	scoring::ScoreFunction const & scorefxn
)
{
	Sizes contact_counts_for_dna_by_basepair;
	return get_protein_dna_interface_by_rotamer_contacts( pose, scorefxn, contact_counts_for_dna_by_basepair );
}


///////////////////////////////////////////////////////////////////////////////////////
/// uses nbr_atom distances and residue nbr_radii
///
///
void
find_neighbors(
	bools const & is_core,
	Pose const & pose,
	bools & is_nbr,
	Real const heavyatom_distance_threshold = 6.0
)
{
	//Real const heavyatom_distance_threshold( 6.0 );

	Size const nres( pose.total_residue() );
	is_nbr = is_core;

	for ( Size i=1; i<= nres; ++i ) {
		Residue const & rsd1( pose.residue(i) );
		if ( rsd1.is_virtual_residue() ) continue;
		for ( Size j=1; j<= nres; ++j ) {
			if ( !is_core[j] ) continue;
			if ( is_nbr[i] ) break;
			Residue const & rsd2( pose.residue(j) );
			if ( rsd2.is_virtual_residue() ) continue;
			if ( rsd1.nbr_atom_xyz().distance_squared( rsd2.nbr_atom_xyz() ) <=
					numeric::square( rsd1.nbr_radius() + rsd2.nbr_radius() + heavyatom_distance_threshold ) ) {
				is_nbr[i] = true;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// void
// create_subset_pose(
//           Pose const & pose_full,
//           bools const & subset,
//           sequence::SequenceMapping & mapping,
//           Pose & pose
//           )
// {
//  pose = pose_full;

//  mapping = sequence::SequenceMapping::identity( pose_full.total_residue() );

//  /// delete the extra data, setup the mapping
//  for ( Size i=pose_full.total_residue(); i>= 1; --i ) {
//   if ( !subset[i] ) {
//    // need to delete this position
//    //pose.conformation().delete_residue_slow( i );
//    mapping.delete_target_residue( i );
//    if ( has_motif_data( pose ) ) get_nonconst_motif_data( pose ).delete_residue( i );
//    if ( has_base_partner( pose ) ) retrieve_nonconst_base_partner_from_pose( pose ).delete_residue( i );
//   }
//  }

//  /// now delete the pose residues
//  for ( Size i=pose_full.total_residue(); i>= 1; --i ) {
//   if ( !subset[i] && ( i == 1 || subset[i-1] ) ) {
//    // first residue of an excluded segment
//    Size j(i);
//    while ( j< pose.total_residue() && !subset[j+1] ) ++j;
//    assert( j <= pose.total_residue() && !subset[j] );
//    //TRI.Trace << "delete: " << i << ' ' << j << endl;
//    pose.conformation().delete_residue_range_slow( i, j );
//   }
//  }

//  //set_base_partner( pose ); // set new base partner info


//  for ( Size i=1; i<= pose.num_chains(); ++i ) {
//   add_lower_terminus_type_to_pose_residue( pose, pose.chain_begin( i ) );
//   add_upper_terminus_type_to_pose_residue( pose, pose.chain_end( i ) );
//  }
// }


string
get_symmetric_name( Size const ii, Residue const & rsd )
{
	string name( stripped( rsd.atom_name(ii) ) ); // added "stripped" 3/18/14 (doh!)

	if (      rsd.aa() == aa_asp && ( name == "OD1" || name == "OD2" ) ) name = "OD";
	else if ( rsd.aa() == aa_glu && ( name == "OE1" || name == "OE2" ) ) name = "OE";
	else if ( rsd.aa() == aa_arg && ( name == "NH1" || name == "NH2" ) ) name = "NH";
	return name;
}

// Size
// get_example_atom_index_from_symmetric_name( string const & name, Residue const & rsd )
// {
//  if (      rsd.aa() == aa_asp && name == "OD" ) return rsd.atom_index("OD1" );
//  else if ( rsd.aa() == aa_glu && name == "OE" ) return rsd.atom_index("OE1" );
//  else if ( rsd.aa() == aa_arg && name == "NH" ) return rsd.atom_index("NH1" );
//  return rsd.atom_index( name );
// }


struct PolarContact {
	PolarContact( Size const ii, Residue const & rsd1, Size const jj, Residue const & rsd2 ):
		atom1( get_symmetric_name( ii, rsd1 ) ),
		atype1( rsd1.atom_type(ii).name() ),
		resname1( rsd1.name1() ),
		pos1( rsd1.seqpos() ),
		atom2( get_symmetric_name( jj, rsd2 ) ),
		atype2( rsd2.atom_type(jj).name() ),
		resname2( rsd2.name1() ),
		pos2( rsd2.seqpos() )
	{
		pbassert( rsd1.seqpos() < rsd2.seqpos() );
	}

	string
	string_of() const
	{
		std::ostringstream os;
		os << "PC " <<
			pos1 << ' ' << resname1 << ' ' << atom1 << ' ' << atype1 << ' ' <<
			pos2 << ' ' << resname2 << ' ' << atom2 << ' ' << atype2;
		return os.str();
	}

	string
	comma_separated_string_of() const
	{
		std::ostringstream os;
		os <<
			pos1 << ',' << resname1 << ',' << atom1 << ',' << atype1 << ',' <<
			pos2 << ',' << resname2 << ',' << atom2 << ',' << atype2;
		return os.str();
	}


	string atom1, atype1;
	char resname1;
	Size pos1;
	string atom2, atype2;
	char resname2;
	Size pos2;
};

bool operator==( PolarContact const & a, PolarContact const & b ) {
	return ( a.atom1 == b.atom1 &&
		a.atype1 == b.atype1 &&
		a.resname1 == b.resname1 &&
		a.pos1 == b.pos1 &&
		a.atom2 == b.atom2 &&
		a.atype2 == b.atype2 &&
		a.resname2 == b.resname2 &&
		a.pos2 == b.pos2 );
}

typedef utility::vector1< PolarContact > PolarContacts;



struct ProteinDNA_Contact {
	ProteinDNA_Contact( Size const dp, Size const da, Size const pp, Size const pa, Real const d, bool const p ):
		dpos( dp ),
		datm( da ),
		ppos( pp ),
		patm( pa ),
		distance( d ),
		polar( p )
	{}


	string
	string_of( Pose const & pose ) const
	{
		std::ostringstream out;
		out << " CONTACT: " << dpos << ' ' << ppos << ' ' << distance << ' ' << polar << ' ' <<
			pose.residue(dpos).name1() << ' ' <<
			pose.residue(dpos).atom_name(datm) << ' ' <<
			pose.residue(ppos).name1() << ' ' <<
			pose.residue(ppos).atom_name(patm);
		return out.str();
	}

	string
	compact_string_of( Pose const & pose ) const
	{
		std::ostringstream out;
		out << polar << ' ' << F(4,2,distance) << ' ' <<
			pose.residue(dpos).name1() << ' ' << pose.residue(dpos).atom_name(datm) << ' ' <<
			pose.residue(ppos).name1() << ' ' << pose.residue(ppos).atom_name(patm);
		return out.str();
	}

	Size dpos;
	Size datm;
	Size ppos;
	Size patm;
	Real distance;
	bool polar;
};

/// NOTE -- ignoring distance!!

bool operator==( ProteinDNA_Contact const & a, ProteinDNA_Contact const & b ) {
	return ( a.dpos == b.dpos &&
		a.datm == b.datm &&
		a.ppos == b.ppos &&
		a.patm == b.patm &&
		a.polar == b.polar );
}

typedef utility::vector1< ProteinDNA_Contact > PDContacts;

void
get_protein_dna_base_contacts(
	Size const dpos,
	Size const ppos,
	Pose const & pose,
	PDContacts & polar_contacts,
	PDContacts & nonpolar_contacts
)
{
	using numeric::square;

	Real polar_distance_threshold( 3.3 ), nonpolar_distance_threshold( 4.2 );

	polar_contacts.clear(); nonpolar_contacts.clear();

	Residue const & prsd( pose.residue( ppos ) );
	pbassert( prsd.is_protein() );
	Residue const & drsd( pose.residue( dpos ) );
	pbassert( drsd.is_DNA() );

	if ( prsd.nbr_atom_xyz().distance_squared( drsd.nbr_atom_xyz() ) <=
			square( prsd.nbr_radius() + drsd.nbr_radius() + nonpolar_distance_threshold ) ) {
		for ( Size ii=1; ii<= prsd.nheavyatoms(); ++ii ) {
			AtomType const & itype( prsd.atom_type(ii) );
			for ( Size jj=drsd.first_sidechain_atom(); jj<= drsd.nheavyatoms(); ++jj ) {
				AtomType const & jtype( drsd.atom_type(jj) );
				Real const distance( prsd.xyz(ii).distance( drsd.xyz(jj) ) );
				if ( ( ( itype.is_donor() && jtype.is_acceptor() ) ||
						( jtype.is_donor() && itype.is_acceptor() ) ) &&
						distance <= polar_distance_threshold ) {
					polar_contacts.push_back( ProteinDNA_Contact( dpos, jj, ppos, ii, distance, true ) );
					TRI.Trace << "polar_contact: " << prsd.name() << ' ' << prsd.atom_name( ii ) << ' ' <<
						drsd.name() << ' ' << drsd.atom_name( jj ) << F(9,3,distance) << endl;
				}
				if ( !itype.is_donor() && !itype.is_acceptor() &&
						!jtype.is_donor() && !jtype.is_acceptor() &&
						distance <= nonpolar_distance_threshold ) {
					nonpolar_contacts.push_back( ProteinDNA_Contact( dpos, jj, ppos, ii, distance, false ) );
					TRI.Trace << "nonpolar_contact: " << prsd.name() << ' ' << prsd.atom_name( ii ) << ' ' <<
						drsd.name() << ' ' << drsd.atom_name( jj ) << F(9,3,distance) << endl;
				}
			}
		}
	}
}

void
find_interface_by_potential_protein_dna_contacts(
	Pose const & pose,
	ScoreFunction const & fa_scorefxn,
	bools & is_interface_protein,
	bools & is_interface_dna,
	bool const restrict_to_dna_bases,
	Real const distance_threshold,
	Sizes restrict_to_dna_positions = Sizes()
)
{
	Real const dis2_threshold( distance_threshold * distance_threshold );
	is_interface_protein.clear(); is_interface_protein.resize( pose.total_residue(), false );
	is_interface_dna.clear(); is_interface_dna.resize( pose.total_residue(), false );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & irsd( pose.residue(i) );
		if ( !irsd.is_protein() ) continue;
		/// build rotamers for this position
		pack::rotamer_set::RotamerSetOP rotset(0);

		for ( Size j=1; j<= pose.total_residue(); ++j ) {
			if ( is_interface_protein[i] && is_interface_dna[j] ) continue; // nothing to learn
			if ( !restrict_to_dna_positions.empty() && !has_element( restrict_to_dna_positions, j ) ) continue;

			Residue const & jrsd( pose.residue(j) );
			if ( !jrsd.is_DNA() ) continue;

			if ( numeric::square( irsd.nbr_radius() + jrsd.nbr_radius() + distance_threshold ) <
					irsd.nbr_atom_xyz().distance_squared( jrsd.nbr_atom_xyz() ) ) continue;

			if ( !rotset ) {
				using namespace pack::task;
				using namespace pack::rotamer_set;
				PackerTaskOP ptask( TaskFactory::create_packer_task( pose ) );
				ptask->set_bump_check( false );
				//ptask->initialize_from_command_line(); /// no cmdline dependence
				ptask->or_include_current( true );
				ptask->nonconst_residue_task( i ).or_ex1( true );
				ptask->nonconst_residue_task( i ).or_ex1aro_sample_level( pack::task::ExtraRotSample(6) );
				ptask->nonconst_residue_task( i ).restrict_to_repacking();

				// unnecessary here, yet also required
				utility::graph::GraphOP dummygraph = utility::graph::GraphOP( new utility::graph::Graph( pose.total_residue() ) );

				RotamerSetFactory rsf;
				rotset = RotamerSetFactory().create_rotamer_set( pose );
				rotset->set_resid( i );
				// use_neighbor_context = false => auto-buried
				rotset->build_rotamers( pose, fa_scorefxn, *ptask, dummygraph, false );
			}


			{
				using pack::rotamer_set::Rotamers;
				Real mindis2( 1e6 );
				for ( Rotamers::const_iterator rot= rotset->begin(); rot != rotset->end(); ++rot ) {
					Residue const & rotrsd( **rot );
					for ( Size ii=1; ii<= rotrsd.nheavyatoms(); ++ii ) {
						Size const jj_begin( restrict_to_dna_bases ? jrsd.first_sidechain_atom() : 1 );
						for ( Size jj=jj_begin; jj<= jrsd.nheavyatoms(); ++jj ) {
							mindis2 = min( mindis2, rotrsd.xyz(ii).distance_squared( jrsd.xyz(jj) ) );
							if ( mindis2 <= dis2_threshold ) break;
						}
						if ( mindis2 <= dis2_threshold ) break;
					}
					if ( mindis2 <= dis2_threshold ) break;
				}
				if ( mindis2 <= dis2_threshold ) {
					if ( !is_interface_protein[i] ) TRI.Trace << "is_interface_protein: " << i << F(9,3,sqrt(mindis2) ) << endl;
					if ( !is_interface_dna[j] ) TRI.Trace << "is_interface_dna: " << j << F(9,3,sqrt(mindis2) ) << endl;
					is_interface_protein[i] = is_interface_dna[j] = true;
				}
			} // scope
		} // j = dna rsd
	} // i = protein rsd
}


void
find_interface_by_protein_dna_contacts(
	Pose const & pose,
	bools & is_interface_protein,
	bools & is_interface_DNA
)
{
	is_interface_protein.clear(); is_interface_protein.resize( pose.total_residue(), false );
	is_interface_DNA    .clear(); is_interface_DNA    .resize( pose.total_residue(), false );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).is_protein() ) {
			for ( Size j=1; j<= pose.total_residue(); ++j ) {
				if ( pose.residue(j).is_DNA() ) {
					PDContacts polar_contacts, nonpolar_contacts;
					get_protein_dna_base_contacts( j, i, pose, polar_contacts, nonpolar_contacts );
					if ( ( !polar_contacts.empty() ) || ( !nonpolar_contacts.empty() ) ) {
						if ( !is_interface_protein[i] ) TRI.Trace << "is_interface_protein: " << i << std::endl;
						if ( !is_interface_DNA[j] ) TRI.Trace << "is_interface_DNA: " << j << std::endl;
						is_interface_protein[i] = is_interface_DNA[j] = true;
					}
				}
			}
		}
	}
}

void
get_flexible_polar_contacts_old(
	bools const & is_flexible,
	Pose const & pose,
	PDContacts & polar_contacts,
	bool const only_protein_or_dna = true // to exclude bogus water ones
)
{
	using numeric::square;

	Real polar_distance_threshold( 3.3 ); //, nonpolar_distance_threshold( 4.2 );

	polar_contacts.clear(); //nonpolar_contacts.clear();

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( !is_flexible[i] ) continue;
		Residue const & rsd1( pose.residue( i ) );

		if ( only_protein_or_dna && (!rsd1.is_protein()) && (!rsd1.is_NA() ) ) continue;

		for ( Size j=1; j<= pose.total_residue(); ++j ) {
			if ( is_flexible[j] && j <= i ) continue; // we'll hit this one in the other order
			Residue const & rsd2( pose.residue( j ) );

			if ( only_protein_or_dna && (!rsd2.is_protein()) && (!rsd2.is_NA() ) ) continue;

			if ( rsd1.nbr_atom_xyz().distance_squared( rsd2.nbr_atom_xyz() ) <=
					square( rsd1.nbr_radius() + rsd2.nbr_radius() + polar_distance_threshold ) ) {
				for ( Size ii=1; ii<= rsd1.nheavyatoms(); ++ii ) {
					AtomType const & itype( rsd1.atom_type(ii) );
					for ( Size jj=1; jj<= rsd2.nheavyatoms(); ++jj ) {
						if ( rsd1.atom_is_backbone( ii ) && ( ( !is_flexible[j] ) || rsd2.atom_is_backbone(jj) ) ) continue;
						AtomType const & jtype( rsd2.atom_type(jj) );
						Real const distance( rsd1.xyz(ii).distance( rsd2.xyz(jj) ) );
						if ( ( ( itype.is_donor() && jtype.is_acceptor() ) ||
								( jtype.is_donor() && itype.is_acceptor() ) ) &&
								distance <= polar_distance_threshold ) {
							polar_contacts.push_back( ProteinDNA_Contact( i, ii, j, jj, distance, true ) );
							// Would outputting residue number (in .cpdb files) be helpful?
							// added by AL
							// TRI.Trace << "polar_contact: " << rsd1.name() << ' ' << rsd1.atom_name( ii ) << ' ' <<
							//  rsd2.name() << ' ' << rsd2.atom_name( jj ) << endl;
							TRI.Trace << "polar_contact: " << i << " " << rsd1.name() << ' ' << rsd1.atom_name( ii )
								<< ' ' << j << " " << rsd2.name() << ' ' << rsd2.atom_name( jj ) << endl;
						}
						//       if ( !itype.is_donor() && !itype.is_acceptor() &&
						//          !jtype.is_donor() && !jtype.is_acceptor() &&
						//          distance <= nonpolar_distance_threshold ) {
						//        nonpolar_contacts.push_back( ProteinDNA_Contact( dpos, jj, ppos, ii, distance, false ) );
						//        TRI.Trace << "nonpolar_contact: " << rsd1.name() << ' ' << rsd1.atom_name( ii ) << ' ' <<
						//         rsd2.name() << ' ' << rsd2.atom_name( jj ) << endl;
						//       }
					} // jj
				} // ii
			} // nbrs
		} // j
	} // i
}

void
get_flexible_polar_contacts(
	bools const & is_flexible,
	bool const backbone_is_flexible,
	Pose const & pose,
	PolarContacts & polar_contacts,
	bool const force_protdna = false
)
{
	using numeric::square;

	Real const polar_distance_threshold( 3.4 ), polar_dis2_threshold( numeric::square( polar_distance_threshold ) );
	bool const exclude_intrachain_backbone_backbone( false );

	polar_contacts.clear(); //nonpolar_contacts.clear();

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & irsd( pose.residue( i ) );
		if ( irsd.name() == "VRT1" || irsd.name() == "TP3" ) continue;

		for ( Size j=i+1; j<= pose.total_residue(); ++j ) {
			Residue const & jrsd( pose.residue( j ) );
			if ( force_protdna && !( ( irsd.is_protein() && jrsd.is_DNA() ) || ( irsd.is_DNA() && jrsd.is_protein() ) ) ) continue;
			if ( jrsd.name() == "VRT1" || jrsd.name() == "TP3" ) continue;

			if ( !is_flexible[i] && !is_flexible[j] ) continue; // neither is flexible

			if ( irsd.nbr_atom_xyz().distance_squared( jrsd.nbr_atom_xyz() ) <=
					square( irsd.nbr_radius() + jrsd.nbr_radius() + polar_distance_threshold ) ) {
				for ( Size ii=1; ii<= irsd.nheavyatoms(); ++ii ) {
					AtomType const & itype( irsd.atom_type(ii) );
					for ( Size jj=1; jj<= jrsd.nheavyatoms(); ++jj ) {
						if ( exclude_intrachain_backbone_backbone && irsd.atom_is_backbone( ii ) && jrsd.atom_is_backbone(jj) &&
								pose.chain(i) == pose.chain(j) ) continue;
						if ( !backbone_is_flexible &&
								( ( irsd.atom_is_backbone( ii ) && jrsd.atom_is_backbone(jj) ) ||
								( irsd.atom_is_backbone( ii ) && !is_flexible[j] ) ||
								( jrsd.atom_is_backbone( jj ) && !is_flexible[i] ) ) ) continue;
						AtomType const & jtype( jrsd.atom_type(jj) );

						if ( ( itype.is_donor() && jtype.is_acceptor() ) ||
								( jtype.is_donor() && itype.is_acceptor() ) ) {
							Real const dis2( irsd.xyz(ii).distance_squared( jrsd.xyz(jj) ) );
							if ( dis2 <= polar_dis2_threshold ) {
								polar_contacts.push_back( PolarContact( ii, irsd, jj, jrsd ) );
								TRI.Trace << "polar_contact: " << irsd.name() << ' ' << irsd.atom_name( ii ) << ' ' <<
									jrsd.name() << ' ' << jrsd.atom_name( jj ) << endl;
							}
						}
					} // jj
				} // ii
			} // nbrs
		} // j
	} // i
}

Size
get_symmetric_atom(
	Size const atom,
	Residue const & rsd
)
{

	if ( rsd.aa() == aa_phe || rsd.aa() == aa_tyr ) {
		if ( atom == rsd.atom_index("CD1") ) return rsd.atom_index("CD2");
		if ( atom == rsd.atom_index("CD2") ) return rsd.atom_index("CD1");
		if ( atom == rsd.atom_index("CE1") ) return rsd.atom_index("CE2");
		if ( atom == rsd.atom_index("CE2") ) return rsd.atom_index("CE1");
	}
	return 0;
}

void
count_nonpolar_contacts_recovered(
	bools const & is_flexible,
	bool const backbone_is_flexible,
	Pose const & refpose,
	Pose const & modpose,
	Size & total_contacts,
	Size & recovered_contacts,
	bool const force_protdna = false
)
{
	using numeric::square;
	// runtime_assert( refpose.total_residue() == modpose.total_residue() );

	total_contacts = recovered_contacts = 0;

	Real const distance_threshold( 4.4 ), dis2_threshold( numeric::square( distance_threshold ) ),
		distance_threshold_wslop( distance_threshold + 0.6 ),
		dis2_threshold_wslop( numeric::square( distance_threshold_wslop ) );

	bool const exclude_intrachain_backbone_backbone( true );

	Size const min_nres( min( refpose.total_residue(), modpose.total_residue() ) );

	for ( Size i=1; i<= min_nres; ++i ) {
		Residue const & irsd( refpose.residue( i ) ), &irsdmod( modpose.residue(i) );

		for ( Size j=i+1; j<= min_nres; ++j ) {
			Residue const & jrsd( refpose.residue( j ) ), &jrsdmod( modpose.residue(j) );

			if ( !is_flexible[i] && !is_flexible[j] ) continue; // neither is flexible

			if ( force_protdna && !( ( irsd.is_protein() && jrsd.is_DNA() ) || ( irsd.is_DNA() && jrsd.is_protein() ) ) ) continue;

			if ( irsd.nbr_atom_xyz().distance_squared( jrsd.nbr_atom_xyz() ) <=
					square( irsd.nbr_radius() + jrsd.nbr_radius() + distance_threshold ) ) {
				for ( Size ii=1; ii<= irsd.nheavyatoms(); ++ii ) {
					AtomType const & itype( irsd.atom_type(ii) );
					string const & iiname( irsd.atom_name(ii) );
					for ( Size jj=1; jj<= jrsd.nheavyatoms(); ++jj ) {
						if ( exclude_intrachain_backbone_backbone && irsd.atom_is_backbone( ii ) && jrsd.atom_is_backbone(jj) &&
								refpose.chain(i) == refpose.chain(j) ) continue;
						if ( !backbone_is_flexible &&
								( ( irsd.atom_is_backbone( ii ) && jrsd.atom_is_backbone(jj) ) ||
								( irsd.atom_is_backbone( ii ) && !is_flexible[j] ) ||
								( jrsd.atom_is_backbone( jj ) && !is_flexible[i] ) ) ) continue;
						AtomType const & jtype( jrsd.atom_type(jj) );
						string const & jjname( jrsd.atom_name(jj) );

						if ( !itype.is_donor() && !itype.is_acceptor() && !jtype.is_donor() && !jtype.is_acceptor() ) { /// nonpolar
							Real const dis2( irsd.xyz(ii).distance_squared( jrsd.xyz(jj) ) );
							if ( dis2 <= dis2_threshold ) {
								++total_contacts;

								if ( irsdmod.aa() == irsd.aa() && jrsdmod.aa() == jrsd.aa() &&
										irsdmod.has( iiname ) && jrsdmod.has( jjname ) ) {
									Size const iimod( irsdmod.atom_index( iiname ) ), jjmod( jrsdmod.atom_index( jjname ) );
									Real mindis2( 1e6 );
									for ( Size ir=1; ir<= 2; ++ir ) {
										Size const iim( ir == 1 ? iimod : get_symmetric_atom( iimod, irsdmod ) );
										for ( Size jr=1; jr<= 2; ++jr ) {
											Size const jjm( jr == 1 ? jjmod : get_symmetric_atom( jjmod, jrsdmod ) );
											if ( iim && jjm ) mindis2 = min( mindis2, irsdmod.xyz(iim).distance_squared( jrsdmod.xyz(jjm) ) );
										}
									}
									if ( mindis2 <= dis2_threshold_wslop ) ++recovered_contacts;
								}
							}
						}
					} // jj
				} // ii
			} // nbrs
		} // j
	} // i
}

string
trim_off_first_character_or_return_dash_if_empty( string const & s )
{
	if ( s.size() ) return s.substr(1);
	else return string("-");
}


void
get_flexible_polar_contacts_status(
	Pose const & mod_pose,
	Pose const & ref_pose,
	bools const & is_flexible,
	bool const backbone_is_flexible,
	std::ostream & out,
	bool const terse = false,
	bool const force_protdna = false
)
{
	PolarContacts mod_contacts, ref_contacts;
	get_flexible_polar_contacts( is_flexible, backbone_is_flexible, mod_pose, mod_contacts, force_protdna );
	get_flexible_polar_contacts( is_flexible, backbone_is_flexible, ref_pose, ref_contacts, force_protdna );

	Size recovered_polar_contacts(0);
	std::ostringstream correct_contacts_sstream, missed_contacts_sstream;
	for ( PolarContacts::const_iterator c = ref_contacts.begin(); c != ref_contacts.end(); ++c ) {
		if ( std::find( mod_contacts.begin(), mod_contacts.end(), *c ) != mod_contacts.end() ) {
			++recovered_polar_contacts;
			correct_contacts_sstream << ';' << c->comma_separated_string_of();
		} else {
			missed_contacts_sstream << ';' << c->comma_separated_string_of();
		}
	}
	std::ostringstream incorrect_contacts_sstream;
	for ( PolarContacts::const_iterator c = mod_contacts.begin(); c != mod_contacts.end(); ++c ) {
		if ( std::find( ref_contacts.begin(), ref_contacts.end(), *c ) == ref_contacts.end() ) {
			incorrect_contacts_sstream << ';' << c->comma_separated_string_of();
		}
	}

	string const flexible_sasa_status( get_flexible_sasa_status_slow( is_flexible, mod_pose ) );

	out << "recovered_polar_contacts: " << recovered_polar_contacts <<
		" total_ref_polar_contacts: " << ref_contacts.size() <<
		" total_decoy_polar_contacts: " << mod_contacts.size();
	if ( !terse ) {
		out << " correct_contacts: " << trim_off_first_character_or_return_dash_if_empty( correct_contacts_sstream.str() ) <<
			" incorrect_contacts: " << trim_off_first_character_or_return_dash_if_empty( incorrect_contacts_sstream.str() ) <<
			" missed_contacts: " << trim_off_first_character_or_return_dash_if_empty( missed_contacts_sstream.str() );
		if ( !force_protdna ) out << " sasa_status: " << flexible_sasa_status;
	}

}



////
void
get_flexible_nonpolar_contacts_status(
	Pose const & mod_pose,
	Pose const & ref_pose,
	bools const & is_flexible,
	bool const backbone_is_flexible,
	std::ostream & out,
	bool const force_protdna = false
)
{
	Size total_ref_contacts, recovered_ref_contacts, total_mod_contacts, correct_mod_contacts;

	count_nonpolar_contacts_recovered( is_flexible, backbone_is_flexible, ref_pose, mod_pose,
		total_ref_contacts, recovered_ref_contacts, force_protdna );

	count_nonpolar_contacts_recovered( is_flexible, backbone_is_flexible, mod_pose, ref_pose,
		total_mod_contacts, correct_mod_contacts, force_protdna );


	out <<
		" recovered_nonpolar_contacts: " << recovered_ref_contacts <<
		" total_ref_nonpolar_contacts: " << total_ref_contacts <<
		" correct_nonpolar_contacts: " << correct_mod_contacts <<
		" total_decoy_nonpolar_contacts: " << total_mod_contacts;

}



////
// returns is_interface, which includes positions in partner1 and partner2
//
bools
find_interface_by_heavyatom_contacts(
	Real const polar_distance_threshold,
	Real const nonpolar_distance_threshold,
	bools const & is_partner1,
	bools const & is_partner2,
	Pose const & pose
)
{
	bools is_interface( pose.total_residue(), false );
	Real const polar_dis2_threshold( numeric::square( polar_distance_threshold ) );
	Real const nonpolar_dis2_threshold( numeric::square( nonpolar_distance_threshold ) );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( !is_partner1[i] ) continue;
		Residue const & rsd1( pose.residue(i) );
		for ( Size j=1; j<= pose.total_residue(); ++j ) {
			if ( !is_partner2[j] ) continue;
			Residue const & rsd2( pose.residue(j) );
			if ( rsd1.nbr_atom_xyz().distance_squared( rsd2.nbr_atom_xyz() ) <=
					numeric::square( rsd1.nbr_radius() + rsd2.nbr_radius() + nonpolar_distance_threshold ) ) {
				bool found_contact( false );
				for ( Size ii=1; ii<= rsd1.nheavyatoms() && !found_contact; ++ii ) {
					bool const ii_polar( rsd1.atom_type(ii).is_donor() || rsd1.atom_type(ii).is_acceptor() );
					for ( Size jj=1; jj<= rsd2.nheavyatoms() && !found_contact; ++jj ) {
						bool const jj_polar( rsd2.atom_type(jj).is_donor() || rsd2.atom_type(jj).is_acceptor() );
						Real const dis2( rsd1.xyz( ii ).distance_squared( rsd2.xyz(jj) ) );
						if ( ( ii_polar && jj_polar && dis2 <= polar_dis2_threshold ) ||
								( !ii_polar && !jj_polar && dis2 <= nonpolar_dis2_threshold ) ) {
							found_contact = true;
						}
					}
				}
				if ( found_contact ) {
					is_interface[i] = is_interface[j] = true;
				}
			}
		}
	}
	return is_interface;

}


////
// returns is_interface, which includes positions in partner1 and partner2
//
bools
find_interface_by_heavyatom_contacts(
	Real const heavyatom_distance_threshold,
	bools const & is_partner1,
	bools const & is_partner2,
	Pose const & pose
)
{
	bools is_interface( pose.total_residue(), false );
	Real const dis2_threshold( numeric::square( heavyatom_distance_threshold ) );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( !is_partner1[i] ) continue;
		Residue const & rsd1( pose.residue(i) );
		for ( Size j=1; j<= pose.total_residue(); ++j ) {
			if ( !is_partner2[j] ) continue;
			Residue const & rsd2( pose.residue(j) );
			if ( rsd1.nbr_atom_xyz().distance_squared( rsd2.nbr_atom_xyz() ) <=
					numeric::square( rsd1.nbr_radius() + rsd2.nbr_radius() + heavyatom_distance_threshold ) ) {
				bool found_contact( false );
				for ( Size ii=1; ii<= rsd1.nheavyatoms() && !found_contact; ++ii ) {
					for ( Size jj=1; jj<= rsd2.nheavyatoms() && !found_contact; ++jj ) {
						if ( rsd1.xyz( ii ).distance_squared( rsd2.xyz(jj) ) <= dis2_threshold ) {
							found_contact = true;
						}
					}
				}
				if ( found_contact ) {
					is_interface[i] = is_interface[j] = true;
				}
			}
		}
	}
	return is_interface;

}




void
mutate_mutable_positions(
	Sizes const & dna_mutation_positions,
	Pose & pose
)
{
	for ( Sizes::const_iterator pos = dna_mutation_positions.begin(); pos != dna_mutation_positions.end();
			++pos ) {
		AA new_aa; new_aa = AA( first_DNA_aa + static_cast< int >( numeric::random::uniform() * 4 ) );
		devel::dna::make_base_pair_mutation( pose, *pos, new_aa );
	}
}


bool // returns TRUE on success
get_flexibility_from_frank_resfiles(
	Pose const & pose, // for pdb_info-based resfile reading
	string const & pdbid,
	string resfile_dir,
	string const & resfile_tag, // eg _cleanwat_0001_
	bools & is_flexible,
	bools & is_evalpos
)
{
	is_flexible.clear(); is_flexible.resize( pose.total_residue(), false );
	is_evalpos.clear(); is_evalpos.resize( pose.total_residue(), false );

	if ( resfile_dir[ resfile_dir.size()-1 ] != '/' ) resfile_dir.push_back('/');

	string predict_resfile( resfile_dir + pdbid + resfile_tag + "predict.resfile" );
	string eval_resfile( resfile_dir + pdbid + resfile_tag + "eval.resfile" );
	if ( !utility::file::file_exists( eval_resfile ) ) {
		eval_resfile = resfile_dir + pdbid + resfile_tag + ".resfile";
	}
	if ( !utility::file::file_exists( eval_resfile ) ) {
		cout << "unable to find resfiles: " << pdbid << ' ' << resfile_dir << ' ' << resfile_tag << endl;
		cerr << "unable to find resfiles: " << pdbid << ' ' << resfile_dir << ' ' << resfile_tag << endl;
		return false; // FAILURE

	} else {
		if ( !utility::file::file_exists( predict_resfile ) ) predict_resfile = eval_resfile;
		runtime_assert( utility::file::file_exists( predict_resfile ) &&
			utility::file::file_exists( eval_resfile ) );
		fill_bools_from_resfile( pose, predict_resfile, is_flexible );
		fill_bools_from_resfile( pose, eval_resfile, is_evalpos );
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			TRI.Trace << "resfile setup: " << I(4,i) << " is_flexible: " << is_flexible[i] <<
				" is_evalpos: " << is_evalpos[i] << endl;
		}
	}
	return true;
}



#endif
