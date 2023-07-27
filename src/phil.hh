// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


#ifndef INCLUDED_apps_pilot_phil_phil_HH
#define INCLUDED_apps_pilot_phil_phil_HH


/// @file
/// @brief

#include <core/types.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
//#include <core/pack/task/residue_selector/ResidueSelector.hh>
#include <core/chemical/ResidueTypeFinder.hh>
// #include <core/chemical/ResidueTypeSelector.hh>
// #include <core/pack/pack_rotamers_envdep.hh>
#include <basic/options/util.hh>
#include <basic/basic.hh>
//#include <basic/options/util.hh>
//#include <basic/prof.hh>

#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh> // for add_variant_type_to_pose_residue, was in pose/util.hh
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
//#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/BasePartner.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>

#include <utility/vector1.functions.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/tools/make_map.hh>

#include <devel/blab/tools.hh>
#include <devel/blab/typedef.hh>
#include <devel/dna/util.hh> // dna_aa_from_oneletter_code

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <basic/Tracer.hh>

//// USING
using namespace core;
using namespace pose;
using namespace conformation;
using namespace chemical;
using namespace scoring;
//using import_pose::pose_from_pdb;
// using namespace options;
// namespace OK = OptionKeys;

// lazy:
using utility::vector1;
using utility::tools::make_vector1;
using ObjexxFCL::string_of;
using ObjexxFCL::lead_zero_string_of;
using ObjexxFCL::is_float;
using ObjexxFCL::is_int;
using ObjexxFCL::float_of; // dont use for very small numbers (below 1e-37 ish)
using ObjexxFCL::double_of;
using ObjexxFCL::int_of;
using ObjexxFCL::stripped; // trimmed copy of input string, leading AND trailing whitespace
using namespace basic::options;

using kinematics::FoldTree;
using kinematics::Stub;
using kinematics::RT;
using kinematics::Jump;
using kinematics::MoveMap;
using kinematics::MoveMapOP;
using scoring::dna::set_base_partner;
using scoring::dna::retrieve_base_partner_from_pose;

using namespace std;
using ObjexxFCL::format::F;
using ObjexxFCL::format::E;
using ObjexxFCL::format::A;
using ObjexxFCL::format::I;

using numeric::random::uniform;
using numeric::random::random_range;
using numeric::random::gaussian;
using numeric::random::random_permutation;

using core::chemical::ResidueTypeFinder;
// using core::chemical::ResidueTypeSelector;

//using numeric::random::random_element;
template< class T >
T const &
random_element( utility::vector1< T > const & v )
{
	return numeric::random::rg().random_element( v );
}

template< class T >
void
my_random_permutation( utility::vector1< T > & v )
{
	numeric::random::random_permutation( v, numeric::random::rg() );
}

template< class T >
std::string
string_join(
	utility::vector1< T > const & v,
	string const & sep,
	string const empty_retval = ""
)
{
	if ( v.empty() ) return empty_retval;

	string retval;
	for( Size i=1, i_end=v.size(); i<= i_end; ++i ) {
		retval += sep + string_of( v[i] );
	}
	return retval.substr(sep.size()); // trim off leading sep
}

// using std::string;
// using std::endl;

/// from devel/blab/tools.hh:
using devel::blab::split_to_vector1;
//using devel::blab::filebase;
using devel::blab::bools_from_sizes;
using devel::blab::sizes_from_bools;
/// from devel/blab/typedef.hh
using devel::blab::strings;
using devel::blab::Reals;
using devel::blab::Sizes;
using devel::blab::bools;


typedef numeric::xyzMatrix< Real > Matrix;
typedef utility::vector1< Vector > Vectors;


std::string
filebase( std::string const & file )
{
	size_t found = file.find_last_of("/\\");
	if ( found == std::string::npos ) return file;
	else return file.substr(found+1);
}

//std::string const CENTROID_DNA("centroid_dna");
using chemical::CENTROID_DNA;

//typedef utility::vector1< ResidueCOP > ResidueCOPs; // should in Residue.fwd.hh

/// @details  Helper function to remove the leading directory stuff
// inline
// std::string
// filebase( std::string const & file )
// {
//  size_t found = file.find_last_of("/\\");
//  if ( found == std::string::npos ) return file;
//  else return file.substr(found+1);
// }

/// @details  Like an assert, but always there
#define pbassert( tmp ) if ( !( tmp ) ) utility_exit_with_message("pbassert failed!");

/// put this here since it's going to change w/trunk, update code in one place
inline
Size
pose_resid_from_PDB( Pose const & pose, int const pdbpos, char const pdbchain )
{
	return pose.pdb_info()->pdb2pose( pdbchain, pdbpos );
}

/// @details  Useful helper for setting a standard foldtree
void
set_chain_end_fold_tree(
	Pose & pose
)
{
	Conformation const & conf( pose.conformation() );
	kinematics::FoldTree f( pose.size() );
	Size const root( conf.chain_end( 1 ) );
	for ( Size i=2; i<= conf.num_chains(); ++i ) {
		Size const anchor  ( conf.chain_end  ( i ) );
		Size const cutpoint( conf.chain_begin( i ) - 1 );
		f.new_jump( root, anchor, cutpoint );
	}
	f.reorder( root );
	pose.fold_tree( f );
}

// silly helper
bool
is_small( Real const x )
{
	return ( std::abs(x)<1e-3 );
}

/// @details  Split a string to a vector1
// utility::vector1< std::string >
// split_to_vector1( std::string const & s )
// {
//  utility::vector1< std::string > v;

//  std::istringstream l( s );
//  std::string tag;
//  l >> tag;
//  while ( !l.fail() ) {
//   v.push_back( tag );
//   l >> tag;
//  }
//  return v;
// }

// /// @details  Split a string to a vector1 using sep as a separator
// utility::vector1< std::string >
// split_to_vector1(
//          std::string s,
//          std::string const & sep
//          )
// {
//  utility::vector1< std::string > v;

//  Size pos( s.find( sep ) );
//  while ( pos != std::string::npos ) {
//   v.push_back( s.substr( 0, pos ) );
//   s.erase( 0, pos + sep.size() );
//   pos = s.find( sep );
//  }
//  assert( s.find( sep ) == std::string::npos );
//  v.push_back( s );
//  return v;
// }


// void
// read_lines_from_file(
//            string const & filename,
//            strings & lines
//            )
// {
//  PROF_START( util::GB_SETUP_FOR_PACKING );
//  utility::io::izstream data( filename );
//  string line;
//  while ( getline( data,line ) ) lines.push_back( line );
//  data.close();
//  PROF_STOP ( util::GB_SETUP_FOR_PACKING );
// }


/// helper
void
get_mean_sdev( Reals const & vals, Real & mean, Real & sdev )
{
	mean = sdev = 0.0;
	if ( vals.empty() ) return;

	for ( Reals::const_iterator val= vals.begin(); val != vals.end(); ++val ) mean += *val;
	mean /= vals.size();
	for ( Reals::const_iterator val= vals.begin(); val != vals.end(); ++val ) sdev += ( *val - mean ) * ( *val - mean );
	sdev = std::sqrt( sdev / vals.size() );
}

// map keys as a utility::vector1
template < typename T1, typename T2 >
utility::vector1< T1 >
get_keys( std::map< T1, T2 > const & m )
{
	utility::vector1< T1 > ks;
	for ( typename std::map< T1,T2 >::const_iterator it= m.begin(); it != m.end(); ++it ) ks.push_back( it->first );
	return ks;
}

// map keys as a std::vector
template < typename T1, typename T2 >
std::vector< T1 >
get_vkeys( std::map< T1, T2 > const & m )
{
	std::vector< T1 > ks;
	for ( typename std::map< T1,T2 >::const_iterator it= m.begin(); it != m.end(); ++it ) ks.push_back( it->first );
	return ks;
}

bool
has_substring( std::string const & full, std::string const & sub )
{
	return ( full.find( sub ) != std::string::npos );
}

bool
has_substring( std::string const & full, const char * sub )
{
	return ( full.find( sub ) != std::string::npos );
}

char
torsion2big_bin(
	Real const phi,
	Real const psi,
	Real const omega
)
{
	if ( std::abs( omega ) < 90 ) {
		return 'O'; // cis-omega
	} else if ( phi >= 0.0 ) {
		if ( -100 < psi && psi <= 100 ) {
			return 'G'; // alpha-L
		} else {
			return 'E'; // E
		}
	} else {
		if ( -125 < psi && psi <= 50 ) {
			return 'A'; // helical
		} else {
			return 'B'; // extended
		}
	}
	return 'X';
}

std::string
torsion2big_bin_string(
	Size const pos1,
	Size const pos2,
	Pose const & pose,
	bool const allow_nonprotein = false
)
{
	runtime_assert( pos1 <= pos2 );
	runtime_assert( pos2 <= pose.total_residue() );
	string bb;
	for ( Size i=pos1; i<= pos2; ++i ) {
		Residue const & irsd( pose.residue(i) );
		if ( !irsd.is_protein() ) {
			if ( allow_nonprotein ) {
				bb.push_back( 'X' );
			} else {
				utility_exit_with_message("torsion2big_bin_string:: nonprotein residue");
			}
		} else {
			Real phi( pose.phi(i) ), psi( pose.psi(i) ), omega( pose.omega(i) );
			if ( irsd.is_lower_terminus() ) { phi = -60.0; }                   // neutral value
			if ( irsd.is_upper_terminus() ) { psi = 150.0; omega = 179.99; }   // neutral values
			bb.push_back( torsion2big_bin( phi, psi, omega ) );
		}
	}
	return bb;

}


Size
get_jump_to_position(
	Size const seqpos,
	Pose const & pose
)
{
	kinematics::FoldTree const & f( pose.fold_tree() );
	for ( Size j=1; j<= f.num_jump(); ++j ) {
		if ( Size( f.downstream_jump_residue(j) ) == seqpos ) return j;
	}
	utility_exit_with_message("NO jump to "+string_of( seqpos ) );
	return 0;
}


Size
chain_begin( Size const ch, Pose const & pose )
{
	return pose.conformation().chain_begin(ch);
}

Size
chain_end( Size const ch, Pose const & pose )
{
	return pose.conformation().chain_end(ch);
}

Size
chain_begin( Pose const & pose, Size const ch )
{
	return pose.conformation().chain_begin(ch);
}

Size
chain_end( Pose const & pose, Size const ch )
{
	return pose.conformation().chain_end(ch);
}

Size
num_chains( Pose const & pose )
{
	return pose.conformation().num_chains();
}

/// @details  This function will make a sequence mutation while trying to preserve the variants
void
make_sequence_change(
	Size const seqpos,
	AA const & new_aa,
	pose::Pose & pose,
	Size const which_his_variant = 1
)
{
	conformation::Residue const & current_rsd( pose.residue( seqpos ) );
	if ( current_rsd.aa() == new_aa ) return; // already done

	// silly hacking
	static ResidueTypeSet const * fa_standard_rsd_set
		( & ( *ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) ) );
	ResidueTypeSet const * current_rsd_set( &(*pose.residue_type_set_for_pose()));
	// ResidueTypeSet const * current_rsd_set( &( *current_rsd.residue_type_set() ) );
	bool const current_rsd_is_fa_standard_rsd( ( fa_standard_rsd_set == current_rsd_set ) );

	ResidueTypeCOPs rsd_types
		( current_rsd_set->get_all_types_with_variants_aa( new_aa, current_rsd.type().variant_types() ) );
	// ResidueTypeCOPs rsd_types
	//  ( ResidueTypeSelector().set_aa( new_aa ).match_variants( current_rsd.type() ).select( current_rsd.residue_type_set() ) );

	Size rsd_types_index( 1 );
	std::string const errmsg
		( "make_sequence_change failed: new_aa= "+name_from_aa(new_aa)+" rsd_types.size()= "+string_of( rsd_types.size() ) );

	if ( new_aa == aa_his ) {
		if ( which_his_variant > rsd_types.size() ) utility_exit_with_message( errmsg );
		if ( rsd_types.size() != 2 && current_rsd_is_fa_standard_rsd ) utility_exit_with_message( errmsg );
		rsd_types_index = which_his_variant;
	} else if ( rsd_types.size() != 1 ) {
		utility_exit_with_message( errmsg );
	}

	conformation::ResidueOP new_rsd( ResidueFactory::create_residue( *(rsd_types[ rsd_types_index ] ),
		current_rsd, pose.conformation() ) );
	pose.replace_residue( seqpos, *new_rsd, false );
}

conformation::ResidueOP
get_vanilla_protein_residue( char const name1 )
{
	chemical::ResidueTypeSetCOP rsd_set( chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	chemical::ResidueTypeCOP new_rsd_type
		( ResidueTypeFinder( *rsd_set ).name1( name1 ).get_representative_type());
	runtime_assert( new_rsd_type->is_protein() );
	return conformation::ResidueFactory::create_residue( *new_rsd_type );
}

conformation::ResidueOP
get_vanilla_dna_residue( char const name1 )
{
	chemical::ResidueTypeSetCOP rsd_set( chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	AA const aa( devel::dna::dna_aa_from_oneletter_code( name1 ) );
	chemical::ResidueTypeCOP new_rsd_type( ResidueTypeFinder( *rsd_set ).aa( aa ).get_representative_type());
	runtime_assert( new_rsd_type->is_DNA() );
	return ResidueFactory::create_residue( *new_rsd_type );
}


static Real const init_phi( -150.0 ), init_psi( 150.0 ), init_omega( 180.0 );

bool
init_torsions_still_present(
	bools const & is_flexible,
	Pose const & pose
)
{
	using basic::subtract_degree_angles;

	bool done( true );
	for ( Size i=1; i<= pose.size(); ++i ) {
		if ( is_flexible[i] && pose.residue(i).is_protein() ) {
			Residue const & rsd( pose.residue(i) );
			if ( ( std::abs( subtract_degree_angles( pose.phi  (i), init_phi   ) ) < 1e-3 || rsd.is_lower_terminus() ) &&
					( std::abs( subtract_degree_angles( pose.psi  (i), init_psi   ) ) < 1e-3 || rsd.is_upper_terminus() ) &&
					( std::abs( subtract_degree_angles( pose.omega(i), init_omega ) ) < 1e-3 || rsd.is_upper_terminus() ) ) {
				done = false;
				break;
			}
		}
	}
	return !done;
}

void
pose_from_pdb( Pose & pose, string const & filename ) {
	ResidueTypeSetCOP rsd_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	import_pose::pose_from_file( pose, *rsd_set, filename , false /* read_fold_tree */, core::import_pose::PDB_file );
}


/// read a centroid pose
void
cenpose_from_pdb( Pose & pose, string const filename )
{
	ResidueTypeSetCOP rsd_set( ChemicalManager::get_instance()->residue_type_set( CENTROID ) );
	pose_from_file( pose, *rsd_set, filename , false, core::import_pose::PDB_file);
}

/// read a centroid pose
void
cendnapose_from_pdb( Pose & pose, string const filename )
{
	ResidueTypeSetCOP rsd_set( ChemicalManager::get_instance()->residue_type_set( CENTROID_DNA ) );
	pose_from_file( pose, *rsd_set, filename , false, core::import_pose::PDB_file);
}

/// useful vector1 function
///
// template < class T >
// bool
// has_element( vector1< T > const & v, T const & t )
// {
//  return ( std::find( v.begin(), v.end(), t ) != v.end() );
// }
// /// never can remember the order...
// template < class T >
// bool
// has_element( T const & t, vector1< T > const & v )
// {
//  return ( std::find( v.begin(), v.end(), t ) != v.end() );
// }
using devel::blab::has_element;

template < class T >
Size
vector1_index( vector1< T > const & v, T const & t )
{
	typename vector1< T >::const_iterator const it( std::find( v.begin(), v.end(), t ) );
	if ( it == v.end() ) return 0; // signal absence
	else return ( it - v.begin() ) + 1;
}

template < class T >
Size
vector1_index( T const & t, vector1< T > const & v )
{
	typename vector1< T >::const_iterator const it( std::find( v.begin(), v.end(), t ) );
	if ( it == v.end() ) return 0; // signal absence
	else return ( it - v.begin() ) + 1;
}

template < class T >
void
remove_element( T const & t, vector1< T > & v )
{
	typename vector1< T >::iterator it( std::find( v.begin(), v.end(), t ) );
	runtime_assert( it != v.end() );
	v.erase( it );
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ScoreFunctionOP
// create_env_indep_score_function( ScoreFunction const & old_scorefxn )
// {
//  ScoreFunctionOP new_scorefxn( old_scorefxn.clone() );

//  new_scorefxn->set_weight( hbond_sr_bb, ( old_scorefxn.get_weight(       hbond_sr_bb ) +
//                       old_scorefxn.get_weight(   env_hbond_sr_bb ) +
//                       old_scorefxn.get_weight( envno_hbond_sr_bb ) ) );
//  new_scorefxn->set_weight( hbond_lr_bb, ( old_scorefxn.get_weight(       hbond_lr_bb ) +
//                       old_scorefxn.get_weight(   env_hbond_lr_bb ) +
//                       old_scorefxn.get_weight( envno_hbond_lr_bb ) ) );
//  new_scorefxn->set_weight( hbond_bb_sc, ( old_scorefxn.get_weight(       hbond_bb_sc ) +
//                       old_scorefxn.get_weight(   env_hbond_bb_sc ) +
//                       old_scorefxn.get_weight( envno_hbond_bb_sc ) ) );
//  new_scorefxn->set_weight( hbond_sc   , ( old_scorefxn.get_weight(       hbond_sc ) +
//                       old_scorefxn.get_weight(   env_hbond_sc ) +
//                       old_scorefxn.get_weight( envno_hbond_sc ) ) );

//  new_scorefxn->set_weight(   env_hbond_sr_bb, 0.0 );
//  new_scorefxn->set_weight( envno_hbond_sr_bb, 0.0 );
//  new_scorefxn->set_weight(   env_hbond_lr_bb, 0.0 );
//  new_scorefxn->set_weight( envno_hbond_lr_bb, 0.0 );
//  new_scorefxn->set_weight(   env_hbond_bb_sc, 0.0 );
//  new_scorefxn->set_weight( envno_hbond_bb_sc, 0.0 );
//  new_scorefxn->set_weight(   env_hbond_sc, 0.0 );
//  new_scorefxn->set_weight( envno_hbond_sc, 0.0 );

//  /// Now for env-elec
//  Real const env_elec_weight2hack_elec_weight_factor( 0.5 ); // very rough guess...
//  runtime_assert( fabs( old_scorefxn.get_weight( env_elec_sc ) ) < 1e-3 );
//  new_scorefxn->set_weight( fa_elec, env_elec_weight2hack_elec_weight_factor * old_scorefxn.get_weight( env_elec ) );
//  new_scorefxn->set_weight( env_elec, 0.0 );

//  runtime_assert( !pack::score_function_has_local_context_dependent_methods( *new_scorefxn ) );

//  // TR.Trace << "create_env_indep_score_function: OLD" << endl;
//  // old_scorefxn.show( TR.Trace );
//  // TR.Trace << "create_env_indep_score_function: NEW" << endl;
//  // new_scorefxn->show( TR.Trace );

//  return new_scorefxn;
// }


core::chemical::AA
dna_aa_from_oneletter_code( char const c )
{
	using namespace chemical;
	switch ( c ) {
	case 'a' : return na_ade;
	case 'c' : return na_cyt;
	case 'g' : return na_gua;
	case 't' : return na_thy;
	default :
		utility_exit_with_message("unrecognized oneletter code for dna "+string(1,c) );
	}
	return aa_unk;
}

/// return a random element from a utility::vector1
template< class T >
vector1< T >
random_subset( Size const size, vector1< T > const & fullset )
{
	vector1< T > subset( fullset );

	while ( subset.size() > size ) {
		Size const pos( numeric::random::random_range( 0, subset.size()-1 ) ); // random_range returns int
		if ( pos >= subset.size() ) continue; // shouldnt happen
		subset.erase( subset.begin() + pos );
	}
	return subset;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details  Returns a random vector uniformly distributed on the unit sphere
Vector
random_unit_vector()
{
	// phi and theta both in radians
	Real const phi( 2 * numeric::constants::d::pi * numeric::random::rg().uniform() );
	Real const theta
		( std::acos(numeric::sin_cos_range( 1.0 - 2.0 * numeric::random::rg().uniform() ) ) );
	assert( 0 <=   phi &&   phi <= 2 * numeric::constants::d::pi );
	assert( 0 <= theta && theta <=     numeric::constants::d::pi );
	return core::kinematics::Stub().spherical( phi, theta, 1.0 );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Size
get_turnlen( string const & turn )
{
	if ( turn == "-" ) return 0;
	else return turn.size();
}






///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
void
translate_residue(
	Residue & rsd,
	Vector const & nbr_atom_xyz
)
{
	Vector const translate( nbr_atom_xyz - rsd.nbr_atom_xyz() );
	for ( Size i=1; i<= rsd.natoms(); ++i ) {
		rsd.set_xyz( i, rsd.xyz( i ) + translate );
	}
}

void
append_virtual_residue( Pose & pose )
{
	if ( pose.residue( pose.total_residue() ).name() == "VRT" ) return; // already done

	/// this stuff here is a little silly, but I don't want to change it now...
	Vector xyz_orig;
	{
		/// get the dnachains
		Sizes dnachains;
		for ( Size i=1; i<= num_chains(pose); ++i ) {
			if ( pose.residue(chain_begin(pose,i)).is_DNA() ) dnachains.push_back(i );
		}
		if ( dnachains.empty() ) {
			// no DNA
			xyz_orig = pose.residue( pose.total_residue()/2 ).nbr_atom_xyz();
		} else {
			//pbassert( dnachains.size() == 2 && dnachains[1]+1 == dnachains[2] ); // could relax?
			Size const dna_root( ( chain_begin( pose, dnachains.front() ) + chain_end( pose, dnachains.front() ) )/2 );
			xyz_orig = pose.residue( dna_root ).xyz("P");
		}
	}

	ResidueOP vrtrsd
		( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map( "VRT" ) ) );
	translate_residue( *vrtrsd, xyz_orig );
	pose.append_residue_by_jump( *vrtrsd, pose.total_residue()/2 );
	pose.conformation().insert_chain_ending( pose.total_residue()-1 );


	if ( scoring::dna::has_base_partner( pose ) ) {
		scoring::dna::retrieve_nonconst_base_partner_from_pose( pose ).resize( pose.total_residue() );
	}
}

void
set_simple_fold_tree_that_respects_cutpoint_variants( Pose & pose )
{
	FoldTree f( pose.total_residue() );
	for ( Size i=1; i< pose.total_residue(); ++i ) {
		if ( ( pose.residue( i  ).has_variant_type( chemical::CUTPOINT_LOWER ) &&
				pose.residue( i+1).has_variant_type( chemical::CUTPOINT_UPPER ) ) ||
				( pose.residue( i  ).is_upper_terminus() ) ||
				( pose.residue( i+1 ).name() == "VRT" ) ) {
			f.new_jump( i, i+1, i );
		}
	}
	f.reorder(1);
	pose.fold_tree(f);
}

///////////////////////////////////////////////////////////////////////////////////////
//
// returns 0 if no vrt residue
//
Size
get_vrt_pos( Pose const & pose )
{
	Size vrtpos(0);
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).name() == "VRT" ) {
			pbassert( !vrtpos );
			vrtpos = i;
		}
	}
	return vrtpos;
}
Size
get_first_vrt_pos( Pose const & pose )
{
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).name() == "VRT" ) return i;
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
create_subpose(
	bools const & subset,
	Pose const & full_pose,
	Pose & pose,
	Sizes & old2new
)
{
	old2new.clear();
	for ( Size i=1; i<= full_pose.total_residue(); ++i ) old2new.push_back( i );

	pose = full_pose;

	for ( Size i=full_pose.total_residue(); i>= 1; --i ) {
		if ( !subset[i] && ( i == 1 || subset[i-1] ) ) {
			// first residue of an excluded segment
			Size j(i);
			while ( j< pose.total_residue() && !subset[j+1] ) ++j;
			assert( j <= pose.total_residue() && !subset[j] );
			//TR.Trace << "delete: " << i << ' ' << j << endl;
			pose.conformation().delete_residue_range_slow( i, j );
			for ( Size k=i; k<=j; ++k ) old2new[k] = 0;
			Size const ndel( j-i+1 );
			for ( Size k=j+1; k<= full_pose.total_residue(); ++k ) if ( old2new[k] ) old2new[k] -= ndel;
		}
	}

	// debug old2new
	Size seen(0);
	for ( Size i=1; i<= full_pose.total_residue(); ++i ) {
		if ( old2new[i] ) {
			runtime_assert( old2new[i] == seen+1 );
			++seen;
		}
	}
	runtime_assert( seen == pose.total_residue() );

	// update base partner if necessary
	if ( scoring::dna::has_base_partner( full_pose ) ) {
		scoring::dna::BasePartner const full_partner( scoring::dna::retrieve_base_partner_from_pose( full_pose ) );
		for ( Size i=1; i<= full_pose.total_residue(); ++i ) {
			// sanity check on full_partner
			runtime_assert( full_partner[i] == 0 || full_partner[full_partner[i]] == i );
		}
		scoring::dna::BasePartner & partner( scoring::dna::retrieve_nonconst_base_partner_from_pose( pose ) );
		partner.resize( pose.total_residue() );
		for ( Size i=1; i<= pose.total_residue(); ++i ) partner[i] = 0;
		Sizes seen;
		for ( Size i=1; i<= full_pose.total_residue(); ++i ) {
			Size const j( old2new[i] );
			if ( j ) {
				runtime_assert( j<= pose.total_residue() );
				if ( full_partner[i] && old2new[ full_partner[i] ] ) {
					Size const jp( old2new[ full_partner[i] ] );
					if ( partner[j] ) {
						runtime_assert( partner[j] == jp );
					} else {
						partner[j] = jp;
					}
					if ( partner[jp] ) {
						runtime_assert( partner[jp] == j );
					} else {
						partner[jp] = j;
					}
				}
			}
		}
	}
}

////////////////////////////////////////////////
string
comma_separated_weighted_energies(
	EnergyMap const & energies,
	EnergyMap const & wts,
	bool const skip_zero_energies = false
)
{
	ostringstream out;
	for ( Size i=1; i<= n_score_types; ++i ) {
		ScoreType const st = ScoreType(i);
		if ( wts[st] == 0.0 ) continue;
		if ( skip_zero_energies && energies[st] == 0.0 ) continue;
		out << st << ':' << F(9,3,energies[st] * wts[st] )<< ',';
	}
	string retval( out.str() );
	if ( retval.empty() ) return retval;
	retval.erase( std::remove( retval.begin(), retval.end(),' '), retval.end() );
	retval.erase( retval.end()-1, retval.end() ); // get last ','
	return retval;

}

//
Real
get_scoreline_score_for_tag(
	string const & scoreline,
	string const & tag
)
{
	runtime_assert( scoreline.find(tag+": ") != string::npos );
	istringstream l( scoreline.substr( scoreline.find(tag+": ") ) );
	Real value;
	string tmp;
	l >> tmp >> value;
	runtime_assert( !l.fail() );
	return value;
}

#endif
