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
#ifndef INCLUDED_apps_pilot_phil_chains_HH
#define INCLUDED_apps_pilot_phil_chains_HH

#include <apps/pilot/phil/phil.hh>
#include <core/pose/util.hh>

static basic::Tracer TR_CHAINS_HH( "apps.pilot.phil.chains_hh" );

///////////////////////////////////////////////////////////////////////////////
void
add_termini_at_protein_chainbreaks(
	Pose & pose,
	Real const distance_threshold = 2.0
)
{
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		if ( rsd.is_protein() ) {
			if ( !rsd.is_upper_terminus() ) {
				runtime_assert( pose.residue(i+1).is_protein() && pose.chain(i) == pose.chain(i+1) );
				// check for a chainbreak
				Real const dis( rsd.xyz("C").distance( pose.residue(i+1).xyz("N") ) );
				if ( dis > distance_threshold ) {
					TR_CHAINS_HH << "intra-chain protein chainbreak! " << i << ' ' << dis << ' ' <<
						pose.residue(i).name() << ' ' << pose.residue(i+1).name() << endl;
					add_upper_terminus_type_to_pose_residue( pose, i );
					add_lower_terminus_type_to_pose_residue( pose, i+1 );
					if ( pose.chain(i) == pose.chain(i+1 ) ) pose.conformation().insert_chain_ending(i);
				}
			}
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
// for example, so that cart_bonded won't score across there
//
// not adding chain endings in Conformation
//
//
void
add_cutpoints_at_protein_chainbreaks(
	Pose & pose,
	Real const distance_threshold = 2.0,
	string const tag=string("")
)
{
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		if ( rsd.is_protein() ) {
			if ( !rsd.is_upper_terminus() ) {
				runtime_assert( pose.residue(i+1).is_protein() && pose.chain(i) == pose.chain(i+1) );
				// check for a chainbreak
				Real const dis( rsd.xyz("C").distance( pose.residue(i+1).xyz("N") ) );
				if ( dis > distance_threshold ) {
					TR_CHAINS_HH << "intra-chain protein chainbreak! " << i << ' ' << dis << ' ' << tag << endl;
					add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, i );
					add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, i+1 );
					// if ( pose.chain(i) == pose.chain(i+1 ) ) pose.conformation().insert_chain_ending(i);
				}
			}
		}
	}
}



///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// note- also inserts chain-endings
//
void
add_termini_at_dna_chainbreaks(
	Pose & pose,
	Real const distance_threshold = 2.4 // optimal is ~1.6
)
{
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		if ( rsd.is_DNA() ) {
			if ( !rsd.is_lower_terminus() && ( i==1 || ( !pose.residue(i-1).is_DNA() ) ) ) {
				TR_CHAINS_HH << "missing lower terminus variant! " << i << ' ' << pose.total_residue() << endl;
				add_lower_terminus_type_to_pose_residue( pose, i );
			}
			if ( !rsd.is_upper_terminus() ) {
				if ( i==pose.total_residue() || !pose.residue(i+1).is_DNA() ) {
					TR_CHAINS_HH << "missing upper terminus variant! " << i << ' ' << pose.total_residue() << endl;
					add_upper_terminus_type_to_pose_residue( pose, i );
				} else {
					Real const dis( pose.residue(i).xyz("O3'").distance( pose.residue(i+1).xyz("P") ) );
					if ( dis > distance_threshold ) {
						TR_CHAINS_HH << "intra-chain dna chainbreak! " << i << ' ' << dis << ' ' << endl;
						add_upper_terminus_type_to_pose_residue( pose, i );
						add_lower_terminus_type_to_pose_residue( pose, i+1 );
						if ( pose.chain(i) == pose.chain(i+1 ) ) pose.conformation().insert_chain_ending(i);
					}
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector1< Size >
find_all_dna_chainbreaks( pose::Pose const & pose )
{
	Real const threshold( 2.5 ); // ideal is ~1.6

	vector1< Size > breaks;
	for ( Size i=1; i<= pose.total_residue()-1; ++i ) {
		conformation::Residue const & rsd1( pose.residue(i  ) );
		conformation::Residue const & rsd2( pose.residue(i+1) );
		if ( rsd1.is_DNA() && rsd2.is_DNA() ) {
			if ( rsd1.xyz("O3'").distance( rsd2.xyz("P") ) > threshold ) {
				TR_CHAINS_HH.Trace << "chainbreak: "<< i << ' ' << rsd1.xyz("O3'").distance( rsd2.xyz("P")) << std::endl;
				breaks.push_back( i );
			}
		}
	}
	return breaks;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
is_single_duplex(
	Pose const & pose,
	string const tag = ""
)
{
	Size const nres( pose.total_residue() );
	Size const nbp( nres/2 );


	scoring::dna::BasePartner const & partner( scoring::dna::retrieve_base_partner_from_pose( pose ) );

	bool retval( pose.conformation().num_chains() == 2 && pose.conformation().chain_end( 1 ) == nbp );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {

		if ( !pose.residue(i).is_DNA() || !partner[i] ) retval = false;

		if ( i <= nbp && partner[i] != nres - i + 1 ) retval = false;

		bool const is_lower( i == 1 || i == nbp+1 ), is_upper( i == nbp || i == nres );
		if (  is_lower && !pose.residue(i).is_lower_terminus() ) retval = false;
		if ( !is_lower &&  pose.residue(i).is_lower_terminus() ) retval = false;
		if (  is_upper && !pose.residue(i).is_upper_terminus() ) retval = false;
		if ( !is_upper &&  pose.residue(i).is_upper_terminus() ) retval = false;

		TR_CHAINS_HH.Trace << "is_single_duplex: retval " << retval << I(4,i) <<
			" chain " << pose.chain(i) <<
			" partner " << I(4,partner[i]) <<
			" is_lower " << is_lower << ' ' << pose.residue(i).is_lower_terminus() <<
			" is_upper " << is_upper << ' ' << pose.residue(i).is_upper_terminus() << ' ' << tag << endl;
	}
	return retval;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
get_chain_centroid_residue(
	Size const chain,
	Pose const & pose
)
{
	runtime_assert( chain>=1 && chain <= num_chains( pose ) );
	Size const cb( chain_begin(chain,pose)), ce( chain_end(chain,pose)), cl(ce-cb+1 );

	Vector centroid(0,0,0);
	for ( Size i=cb; i<= ce; ++i ) {
		centroid += pose.residue(i).xyz("CA");
	}
	centroid /= cl;

	Real mindis2( 1e6 );
	Size cenpos(0);
	for ( Size i=cb; i<= ce; ++i ) {
		Real const dis2( pose.residue(i).xyz("CA").distance_squared( centroid ) );
		if ( dis2 < mindis2 ) {
			mindis2 = dis2;
			cenpos = i;
		}
	}

	TR_CHAINS_HH.Trace << "get_chain_centroid_residue: chain: " << chain << " cenpos: " << cenpos <<
		" mindis: " << F(9,3,sqrt(mindis2)) << endl;

	runtime_assert( cenpos );
	runtime_assert( pose.chain(cenpos) == chain );

	return cenpos;


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
get_residues_centroid_residue(
	Sizes const poslist,
	Pose const & pose
)
{
	runtime_assert( !poslist.empty() );

	Vector centroid(0,0,0);
	foreach_( Size pos, poslist ) centroid += pose.residue(pos).xyz("CA");
	centroid /= poslist.size();

	Real mindis2( 1e6 );
	Size cenpos(0);
	foreach_( Size i, poslist ) {
		Real const dis2( pose.residue(i).xyz("CA").distance_squared( centroid ) );
		if ( dis2 < mindis2 ) {
			mindis2 = dis2;
			cenpos = i;
		}
	}

	TR_CHAINS_HH.Trace << "get_residues_centroid_residue: cenpos: " << cenpos <<
		" mindis: " << F(9,3,sqrt(mindis2)) << endl;

	runtime_assert( cenpos );
	runtime_assert( has_element( cenpos, poslist ) );

	return cenpos;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// remove everything except termini and disulfides
///
bools
remove_funny_variants(
	Pose & pose,
	string const filename = "unknown_filename"
)
{
	bools residue_changed( pose.total_residue(), false );

	// try removing variants
	for ( Size i=1; i<= pose.total_residue(); ++i ) {

		// remove all variants!!! eg  "PRO:pro_hydroxylated_case1" is causing us trouble
		for ( Size i_variant = FIRST_VARIANT; i_variant <= N_VARIANTS; ++i_variant ) {
			VariantType const vt = VariantType(i_variant);
			if (
				vt == DISULFIDE ||
				vt == LOWER_TERMINUS_VARIANT ||
				vt == UPPER_TERMINUS_VARIANT ||
				vt == CUTPOINT_LOWER ||
				vt == CUTPOINT_UPPER
			) {
				continue;
			}
			if ( pose.residue(i).has_variant_type( vt ) ) {
				TR_CHAINS_HH.Trace << "Try to remove variant " << vt << " at position " << i << ' ' <<
					pose.residue(i).name() << ' ' << filename << endl;
				//fflush( stdout );
				remove_variant_type_from_pose_residue( pose, vt, i );
				residue_changed[i] = true;
			}
		}
	}
	return residue_changed;
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif
