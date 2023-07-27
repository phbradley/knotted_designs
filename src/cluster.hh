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
#ifndef INCLUDED_apps_pilot_phil_cluster_HH
#define INCLUDED_apps_pilot_phil_cluster_HH

#include <apps/pilot/phil/phil.hh>
// #include <apps/pilot/phil/phil_options.hh>
// #include <apps/pilot/phil/phil_io.hh>
// #include <apps/pilot/phil/types.hh>

// #include <apps/pilot/phil/stub_transform.hh>
// //#include <apps/pilot/phil/stub_frag.hh>
// #include <apps/pilot/phil/symscores.hh>

// #include <devel/blab/loops/util.hh>
// #include <core/io/silent/SilentFileData.hh>
// #include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>

// #include <basic/datacache/BasicDataCache.hh>
// #include <basic/datacache/CacheableString.hh>
// #include <basic/datacache/CacheableStringFloatMap.hh>
// #include <core/pose/datacache/CacheableDataType.hh>

// //#include <core/id/AtomID.hh>
// #include <core/chemical/AtomType.hh>
// #include <core/scoring/packstat/compute_sasa.hh>
// #include <core/scoring/EnergyGraph.hh>
// #include <core/scoring/rms_util.hh>
// #include <core/scoring/sasa.hh>
// #include <core/scoring/constraints/AtomPairConstraint.hh>
// #include <core/scoring/dna/base_geometry.hh>
// #include <core/scoring/dna/BasePartner.hh>
// #include <core/scoring/ScoreFunctionFactory.hh>
// #include <core/scoring/symmetry/SymmetricScoreFunction.hh>
// #include <core/scoring/Energies.hh>
// #include <core/scoring/methods/EnergyMethodOptions.hh>
// #include <core/scoring/methods/EnergyMethodCreator.hh>
// #include <core/scoring/methods/WholeStructureEnergy.hh>

// #include <core/fragment/ConstantLengthFragSet.hh>

// #include <protocols/simple_moves/symmetry/SymMinMover.hh>
// #include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
// #include <protocols/simple_moves/symmetry/SymRotamerTrialsMover.hh>
// #include <protocols/relax/FastRelax.hh>
// #include <protocols/abinitio/ClassicAbinitio.hh>

// #include <protocols/sasa_scores/sasapack.hh>

// #include <devel/blab/classic_frags/TorsionFragment.hh>
// #include <devel/dna/util.hh>

// #include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>
// #include <protocols/moves/MonteCarlo.hh>
#include <core/pose/util.hh>
// #include <core/conformation/ResidueFactory.hh>
// #include <core/conformation/util.hh>
// #include <core/chemical/ChemicalManager.hh>
// #include <core/conformation/symmetry/SymDof.hh>
// #include <core/conformation/symmetry/SymmetryInfo.hh>
// #include <core/scoring/symmetry/SymmetricScoreFunction.hh>
// #include <core/pose/symmetry/util.hh>
// #include <core/chemical/VariantType.hh>

// #include <core/pack/task/PackerTask.hh>
// #include <core/pack/task/TaskFactory.hh>
// #include <core/pack/task/operation/TaskOperations.hh>
// #include <devel/init.hh>

// #include <numeric/random/random.hh>
// #include <numeric/random/random_permutation.hh>
#include <numeric/model_quality/rms.hh>
// #include <numeric/model_quality/maxsub.hh>

#include <basic/Tracer.hh>
// #include <basic/database/open.hh>

// #include <utility/io/izstream.hh>

// //////// option key includes
// #include <basic/options/keys/out.OptionKeys.gen.hh>
// #include <basic/options/keys/phil.OptionKeys.gen.hh>
// #include <basic/options/keys/run.OptionKeys.gen.hh>
// #include <basic/options/keys/dna.OptionKeys.gen.hh>
// #include <basic/options/keys/relax.OptionKeys.gen.hh>
// #include <basic/options/keys/abinitio.OptionKeys.gen.hh>

#include <fstream>

// static numeric::random::RandomGenerator RG(12321); // <- Magic number, do not change it!!!

static basic::Tracer TR_CLUSTER( "apps.pilot.phil.cluster_hh" );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
adjust_sse_turn_boundary(
	char const sse_bb,
	string & sse,
	string & turn,
	int & sse_end
)
{
	// slide junction backward
	while ( sse.size() && sse[ sse.size()-1 ] != sse_bb ) {
		turn.insert( turn.begin(), sse[ sse.size()-1 ] );
		sse.erase( sse.size()-1, 1 );
		--sse_end;
	}

	// slide junction forward
	while ( turn.size() && turn[0] == sse_bb ) {
		sse.push_back( turn[0] );
		turn.erase( 0, 1 );
		++sse_end;
	}
	if (  sse.size() ) { runtime_assert( sse[ sse.size()-1 ] == sse_bb ); }
	if ( turn.size() ) { runtime_assert( turn[0] != sse_bb ); }

}

/////
void
adjust_turn_sse_boundary(
	char const sse_bb,
	string & turn,
	string & sse,
	int & sse_begin
)
{
	// slide junction forward
	while ( sse.size() && sse[0] != sse_bb ) {
		turn.push_back( sse[0] );
		sse.erase( 0, 1 );
		++sse_begin;
	}
	// slide junction backward
	while ( turn.size() && turn[ turn.size()-1 ] == sse_bb ) {
		sse.insert( sse.begin(), turn[ turn.size()-1 ] );
		turn.erase( turn.size()-1, 1 );
		--sse_begin;
	}

	if (  sse.size() ) { runtime_assert( sse[0] == sse_bb ); }
	if ( turn.size() ) { runtime_assert( turn[ turn.size()-1 ] != sse_bb ); }

}


void
adjust_turn_sse_turn_boundaries(
	char const sse_bb,
	string & turn1,
	string & sse,
	string & turn2,
	int & sse_begin,
	int & sse_end
)
{
	while ( true ) { // have to loop
		string const oldturn1( turn1 ), oldsse( sse ), oldturn2( turn2 );


		// we have to handle a tricky situation where sse starts with sse_bb but then has a non sse_bb char
		while ( true ) {
			//TR.Trace << "adjust_turn_sse_turn_boundaries: loop " << turn1 << ' ' << sse << ' ' << turn2 << endl;
			if ( sse.empty() ) break;
			Size const split_pos( sse.size()/2 );
			string const half1( sse.substr( 0, split_pos ) ), half2( sse.substr( split_pos ) );
			runtime_assert( sse == half1 + half2 );
			if ( half1.empty() || half2.empty() ) break; // too confusing
			if ( half1[0] == sse_bb && half1 != string( half1.size(), sse_bb ) ) {
				turn1.push_back( sse[0] );
				sse.erase( 0,1 );
				++sse_begin;
				continue;
			}
			if ( half2[ half2.size()-1 ] == sse_bb && half2 != string( half2.size(), sse_bb ) ) {
				turn2.insert( turn2.begin(), sse_bb );
				sse.erase( sse.size()-1, 1 );
				--sse_end;
				continue;
			}
			break;
		}
		//TR.Trace << "adjust_turn_sse_turn_boundaries: postloop " << turn1 << ' ' << sse << ' ' << turn2 << endl;

		adjust_turn_sse_boundary( sse_bb, turn1, sse, sse_begin );
		adjust_sse_turn_boundary( sse_bb, sse, turn2, sse_end );

		runtime_assert( int( sse.size() ) == sse_end - sse_begin + 1 );

		//TR.Trace << "adjust_turn_sse_turn_boundaries: final " << turn1 << ' ' << sse << ' ' << turn2 << endl;

		if ( turn1 == oldturn1 && sse == oldsse && turn2 == oldturn2 ) break;
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class VectorsDistanceMetricNoSuper {
public:

	Real
	operator()( Vectors const & coords1, Vectors const & coords2 ) const
	{
		runtime_assert( coords1.size() == coords2.size() );
		if ( coords1.empty() ) return 0.0;
		Real rmsd(0.0);
		for ( Size i=1; i<= coords1.size(); ++i ) {
			rmsd += coords1[i].distance_squared( coords2[i] );
		}
		return sqrt( rmsd / coords1.size() );
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class VectorsDistanceMetricSuper {
public:

	Real
	operator()( Vectors const & coords1, Vectors const & coords2 ) const
	{
		runtime_assert( coords1.size() == coords2.size() );
		if ( coords1.empty() ) return 0.0;
		return numeric::model_quality::calc_rms( coords1, coords2 );
	}
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Real
get_cluster_average(
	string const & scoretag,
	Sizes const & members,
	vector1< map< string, Real > > const & all_scores
)
{
	Real avg(0);
	Size count(0);
	for ( Size i=1; i<= members.size(); ++i ) {
		if ( scoretag == "refold_rmsd" ) {
			bool const member_was_refolded( all_scores[ members[i] ].find( string("refolded") )->second > 0.5 );
			if ( member_was_refolded ) {
				avg += all_scores[ members[i] ].find( scoretag )->second;
				++count;
			}
		} else {
			avg += all_scores[ members[i] ].find( scoretag )->second;
			++count;
		}
	}
	if ( count ) avg/= count;
	return avg;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string
get_cluster_bbtag_component(
	Size const index,
	Sizes const & members,
	strings const & all_bbtags
)
{
	map< string, Size > count;
	for ( Size i=1; i<= members.size(); ++i ) {
		strings const l( split_to_vector1( all_bbtags[ members[i] ], "." ) );
		if ( l.size()<index ) {
			TR_CLUSTER.Trace << "get_cluster_bbtag_component: bad bbtag: "<< all_bbtags[members[i] ] << ' ' << index << endl;
		}
		++count[ l[index] ];
	}
	vector1< std::pair< Size, string > > l;
	for ( map< string, Size >::const_iterator it= count.begin(); it != count.end(); ++it ) {
		l.push_back( make_pair( it->second, it->first ) );
	}
	std::sort( l.begin(), l.end() );
	Size const topcount( l.back().first ), perc( ( 100 * topcount )/members.size() );
	string percstring("("+string_of(perc)+")"), toptag( l.back().second );
	if ( toptag.empty() ) toptag = "-";

	ostringstream out;
	out << ObjexxFCL::format::RJ(5,toptag) << ObjexxFCL::format::RJ(6,percstring);
	return out.str(); // length 11
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


#endif
