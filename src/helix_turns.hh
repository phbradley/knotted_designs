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
#ifndef INCLUDED_apps_pilot_phil_helix_turns_HH
#define INCLUDED_apps_pilot_phil_helix_turns_HH

#include <apps/pilot/phil/phil.hh>
#include <apps/pilot/phil/phil_options.hh>
#include <apps/pilot/phil/stub_frag.hh>
#include <core/optimization/powell.hh>
// #include <apps/pilot/phil/rms.hh>
// #include <apps/pilot/phil/sym_basic.hh>
// #include <apps/pilot/phil/types.hh>

// #include <numeric/xyzVector.io.hh>


static basic::Tracer TR_HELIX_TURNS( "apps.pilot.phil.helix_turns_hh" );


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
compute_helix_axis_and_midpoint(
	Vectors const & coords, // C-alpha coords
	Vector & axis,
	Vector & midpoint
)
{
	runtime_assert( !coords.empty() );
	Size const natoms( min( Size(4), coords.size() ) );
	Real const norm( 1.0 / Real(natoms) );
	Vector bxyz( 0,0,0 ), exyz( 0,0,0 );
	for ( Size i=0; i<natoms; ++i ) {
		bxyz += norm * coords[             1 + i ];
		exyz += norm * coords[ coords.size() - i ];
	}
	axis = ( exyz - bxyz ).normalized();
	midpoint = 0.5 * ( bxyz + exyz );

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// this is pretty slow...
//

void
get_helix_window_stubs(
	Vectors const & coords,
	Size const helixbegin,
	Stub & start_stub,
	Stub & stop_stub
)
{
	using namespace optimization;
	/// from helical_params_test
	Real const ideal_rise( 1.494 ), ideal_twist( numeric::conversions::radians( 99.236 ) ), ideal_radius( 2.292 );

	Size const helix_window( 7 );
	// fit ideal helix to the coords, return starting stub
	static HelicalParamsFitMultifunc func;
	func.set_target_coords( helixbegin, helixbegin+helix_window-1, coords );
	func.optimize_tilt( true );

	Multivec params( 5 );
	params[1] = ideal_rise;
	params[2] = ideal_twist;
	params[3] = 0;
	params[4] = 0;
	params[5] = ideal_radius;

	Real const tolerance( 1e-3 );
	Size iterations;
	Real final_func_value;
	optimization::powell( params, func, tolerance, iterations, final_func_value );

	// TR_HELIX_TURNS.Trace << "get_helix_window_stubs: func " << F(9,3,final_func_value) << " params " <<
	//   F(9,3,params[1]) << F(9,3,params[2]) << F(9,3,params[3]) << F(9,3,params[4]) << F(9,3,params[5]) << endl;

	func.get_helix_stubs_for_params( params, start_stub, stop_stub );

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
compute_turn_rmsds(
	Vectors const & coords,
	vector1< Vectors > const & pdb_coords,
	Real & avg_rmsd,
	Real & med_rmsd
)
{
	Real const min_rmsd( 1e-2 ); // dont include self-rmsds

	// compute rmsds to all turn coords
	Reals rmsds;
	for ( Size i=1; i<= pdb_coords.size(); ++i ) {
		runtime_assert( coords.size() == pdb_coords[i].size() );
		Real const rmsd( numeric::model_quality::calc_rms( coords, pdb_coords[i] ) );
		if ( rmsd>min_rmsd ) rmsds.push_back( rmsd );
	}
	std::sort( rmsds.begin(), rmsds.end() );

	med_rmsd =rmsds[ rmsds.size()/2 ];

	Size const n_outliers( rmsds.size()/10 ); // was 5
	for ( Size i=1; i<= n_outliers; ++i ) rmsds.pop_back();
	for ( Size i=1; i<= rmsds.size(); ++i ) avg_rmsd += rmsds[i];

	avg_rmsd/= rmsds.size();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
compute_turn_strain(
	string const & turn,
	Size const turnbegin,
	Vectors const & ca_coords,
	Real & avg_rmsd,
	Real & avg_rmsd_normed, // 0 to 1
	Real & med_rmsd
)
{
	static map< string, vector1< Vectors > > all_turn_coords;
	static map< string, Reals > all_turn_avg_rmsds;

	Size const turn_buffer( 7 );

	if ( all_turn_coords.empty() ) { //read info on bb coords
		string const pdb_coords_file( option[ my_options::pdb_coords_file ] );
		ifstream data( pdb_coords_file.c_str() );
		string line;
		while ( getline( data, line ) ) {
			strings const l( split_to_vector1( line ) );
			if ( l[1] == "turn_coords:" ) {
				string const turn( l[2] );
				if ( all_turn_coords.find( turn ) == all_turn_coords.end() ) all_turn_coords[ turn ];
				Vectors coords;
				Size const coordslen( int_of( l[6] ) );
				runtime_assert( coordslen == turn.size() + 2*turn_buffer );
				for ( Size i=1; i<= coordslen; ++i ) {
					coords.push_back( Vector( float_of( l[ 3*i+4 ] ), float_of( l[ 3*i+5 ] ), float_of( l[3*i+6 ] ) ) );
				}
				all_turn_coords.find( turn )->second.push_back( coords );
				//TR.Trace << "read turn: " << turn << ' ' << all_turn_coords.find( turn )->second.size() << endl;
			}
		}
		data.close();
	} // initialize the coordinates /////////////////////////////////////////////////////////////////////

	avg_rmsd = med_rmsd = 0;

	if ( all_turn_coords.count(turn) == 0 ) { // dont have data on this turn
		avg_rmsd = med_rmsd = 100.0;
		avg_rmsd_normed = 1.0;
		return;
	}

	if ( all_turn_avg_rmsds.count(turn) == 0 ) { // precompute similarities
		all_turn_avg_rmsds[turn];
		Reals & rmsds( all_turn_avg_rmsds.find(turn)->second ); rmsds.clear();
		vector1< Vectors > const & pdb_coords( all_turn_coords.find( turn )->second );
		TR_HELIX_TURNS.Trace << "precompute turn rmsds: " << turn << ' ' << pdb_coords.size() << endl;
		foreach_( Vectors const & coords, pdb_coords ) {
			Real avg_rmsd_tmp, med_rmsd_tmp;
			compute_turn_rmsds( coords, pdb_coords, avg_rmsd_tmp, med_rmsd_tmp );
			rmsds.push_back( avg_rmsd_tmp );
		}
		std::sort( rmsds.begin(), rmsds.end() );
	}


	Size const coordsbegin( turnbegin-turn_buffer ), coordsend( turnbegin+turn.size()+turn_buffer-1);
	runtime_assert( coordsend <= ca_coords.size() );

	Vectors coords;
	for ( Size i= coordsbegin; i<= coordsend; ++i ) coords.push_back( ca_coords[i] );

	compute_turn_rmsds( coords, all_turn_coords.find( turn )->second, avg_rmsd, med_rmsd );

	Size nbetter(0);
	Reals const & pdb_avg_rmsds( all_turn_avg_rmsds.find(turn)->second );
	foreach_( Real pdb_rmsd, pdb_avg_rmsds ) {
		if ( pdb_rmsd < avg_rmsd ) ++nbetter;
	}

	avg_rmsd_normed = Real( nbetter )/pdb_avg_rmsds.size();

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// backwards compatible
void
compute_turn_strain(
	string const & turn,
	Size const turnbegin,
	Vectors const & ca_coords,
	Real & avg_rmsd,
	Real & med_rmsd
)
{
	Real avg_rmsd_normed;
	compute_turn_strain( turn, turnbegin, ca_coords, avg_rmsd, avg_rmsd_normed, med_rmsd );
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
build_turn_library_v2(
	strings const & turns,
	map< string, StubFrags > & all_frags,
	bool const verbose = false,
	bool const clear_library_at_start = true
)
{
	static bool init( false );

	static map< string, vector1< TurnCoordsLine > > all_turn_coords;

	Size const turn_buffer( 7 ); // should be the same as helix_window

	if ( !init ) { //read info on bb coords
		init = true;
		string const pdb_coords_file( option[ my_options::pdb_coords_file ] );
		ifstream data( pdb_coords_file.c_str() );
		string line;
		while ( getline( data, line ) ) {
			strings const l( split_to_vector1( line ) );
			if ( l[1] == "turn_coords:" ) {
				TurnCoordsLine tcline;
				string const turn( l[2] );
				if ( all_turn_coords.find( turn ) == all_turn_coords.end() ) all_turn_coords[ turn ];
				//Vectors coords;
				Size const coordslen( int_of( l[6] ) );
				runtime_assert( coordslen == turn.size() + 2*turn_buffer );
				for ( Size i=1; i<= coordslen; ++i ) {
					tcline.coords.push_back( Vector( float_of( l[ 3*i+4 ] ), float_of( l[ 3*i+5 ] ), float_of( l[3*i+6 ] ) ) );
				}
				Size const torsions_start( 3*coordslen+8 );
				runtime_assert( l[ torsions_start-1 ] == "bb_torsions:" );
				for ( Size i=0; i< turn.size(); ++i ) { //read the bb torsions
					tcline.torsions.push_back( make_vector1( Real( float_of( l[ torsions_start+3*i   ] ) ),
						Real( float_of( l[ torsions_start+3*i+1 ] ) ),
						Real( float_of( l[ torsions_start+3*i+2 ] ) ) ) );
				}
				tcline.tag = filebase( l.back() ) + "_" + l[ l.size()-2 ] + "_" + l[ l.size()-1 ] +"_"+turn;
				tcline.sequence = l[4];
				runtime_assert( tcline.sequence.size() == turn.size() );
				all_turn_coords.find( turn )->second.push_back( tcline );
				//TR.Trace << "read turn: " << turn << ' ' << all_turn_coords.find( turn )->second.size() << endl;
			}

		}
		data.close();
	} // initialize the coordinates /////////////////////////////////////////////////////////////////////


	if ( clear_library_at_start ) all_frags.clear();
	for ( strings::const_iterator turn= turns.begin(); turn != turns.end(); ++turn ) {
		if ( all_turn_coords.find( *turn ) == all_turn_coords.end() ) {
			utility_exit_with_message("no PDB turn frags for turn "+(*turn));
		}
		all_frags[ *turn ];
		Size const turnlen( turn->size() );
		StubFrags & frags( all_frags.find( *turn )->second );
		//vector1< std::pair< string, Vectors > > const & pdb_turn_coords( all_turn_coords.find( *turn )->second );
		vector1< TurnCoordsLine > const & pdb_turn_coords( all_turn_coords.find( *turn )->second );
		Reals all_strains;
		for ( Size ii=1; ii<= pdb_turn_coords.size(); ++ii ) {
			TurnCoordsLine const & tcl( pdb_turn_coords[ii] );
			Vectors const & turn_coords( tcl.coords );
			runtime_assert( turn_coords.size() == turnlen + 2*turn_buffer );
			Size const stub1_helix_begin( 1 ), stub2_helix_begin( turn_buffer+turnlen+1 );
			Stub stub1, stub2, tmpstub;
			get_helix_window_stubs( turn_coords, stub1_helix_begin, tmpstub, stub1 ); // stop_stub of helix before turn
			get_helix_window_stubs( turn_coords, stub2_helix_begin, stub2, tmpstub ); // start_stub of helix after turn
			StubFrag f;
			f.rt = RT( stub1, stub2 );
			for ( Size j=1; j<= turnlen; ++j ) f.coords.push_back( stub1.global2local( turn_coords[turn_buffer+j] ));
			Real avg_rmsd, med_rmsd;
			compute_turn_strain( *turn, turn_buffer+1, turn_coords, avg_rmsd, med_rmsd );
			f.raw_strain = avg_rmsd;
			all_strains.push_back( f.raw_strain );
			f.torsions = tcl.torsions;
			f.sequence = tcl.sequence;
			f.id = tcl.tag;
			f.segtype = turn_type;
			frags.push_back( f );

			if ( verbose ) { //get some geometry params
				Real const axis_angle( degrees( acos( stub1.M.col_x().dot( stub2.M.col_x() ) ) ) );
				// dihedral
				Vector const p1( stub1.v - stub1.M.col_x() ), p2 ( stub1.v ), p3( stub2.v ), p4( stub2.v + stub2.M.col_x() );
				Real const axis_twist( numeric::dihedral_degrees( p1,p2,p3,p4) );
				TR_HELIX_TURNS.Trace << "turn_geom " << *turn << " axis_enddist: " << F(9,3,stub1.v.distance(stub2.v) ) <<
					" axis_angle: " << F(9,3,axis_angle) << " axis_twist: " << F(9,3,axis_twist) << endl;
			}
			if ( verbose ) { //get some geometry params, a different way
				Vector axis1, midpoint1, axis2, midpoint2;
				Vectors helix1_coords( turn_coords ), helix2_coords( turn_coords );
				helix1_coords.erase( helix1_coords.begin() + turn_buffer, helix1_coords.end() );
				helix2_coords.erase( helix2_coords.begin(), helix2_coords.begin() + turn_buffer + turnlen );
				runtime_assert( helix1_coords.size() == turn_buffer );
				runtime_assert( helix2_coords.size() == turn_buffer );
				compute_helix_axis_and_midpoint( helix1_coords, axis1, midpoint1 );
				compute_helix_axis_and_midpoint( helix2_coords, axis2, midpoint2 );
				Vector const lastca1( turn_coords[ turn_buffer] ),
					endpoint1( midpoint1 + axis1 * ( lastca1-midpoint1).dot( axis1 ) );
				Vector const firstca2( turn_coords[ turn_buffer+turnlen+1 ] ),
					startpoint2( midpoint2 + axis2 * ( firstca2-midpoint2).dot( axis2 ) );
				Real const axis_angle( degrees( acos( axis1.dot( axis2 ) ) ) );
				// dihedral
				Vector const p1( endpoint1 - axis1 ), p2 ( endpoint1 ), p3( startpoint2 ), p4( startpoint2 + axis2 );
				Real const axis_twist( numeric::dihedral_degrees( p1,p2,p3,p4) );
				TR_HELIX_TURNS.Trace << "turn_geom2 " << *turn <<
					" axis_enddist: " << F(9,3,endpoint1.distance(startpoint2) ) <<
					" axis_angle: " << F(9,3,axis_angle) << " axis_twist: " << F(9,3,axis_twist) << endl;
			}
			//if ( frags.size() >= max_frags ) break;
		} // loop over pdb turns of this type

		Size const nfrags( frags.size() );

		for ( Size i=1; i<= nfrags; ++i ) {
			StubFrag & f( frags[i] );
			Size nbetter(0);
			Real const epsilon( 1e-3 );
			for ( Reals::const_iterator s= all_strains.begin(); s!= all_strains.end(); ++s ) {
				if ( *s + epsilon < f.raw_strain ) ++nbetter;
			}

			f.norm_strain = Real( nbetter );
			// nbetter will range from 0 to nfrags-1
			if ( nfrags > 1 ) f.norm_strain /= ( nfrags-1 ); // now range in [0,1]
			//TR.Trace << "norm_strain_turn: " << *turn << F(9,3,f.raw_strain) << F(9,3,f.norm_strain) << endl;
		}


		TR_HELIX_TURNS.Trace << "build_turn_library:: Read " << frags.size() << " turns of type " << *turn << endl;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// assumes that each repeat goes: helix1, turn1, helix2, turn2
///
/// but we don't know which is the inner helix and which is the outer helix...
///
void
compute_helix_axis_angles(
	Pose const & pose,
	Size const base_repeat,
	Size const repeatlen,
	Size const helix1_len,
	Size const turn1_len,
	Size const helix2_len,
	//Size const turn2_len,
	Real & helix1_dist,
	Real & helix1_twist,
	Real & helix2_dist,
	Real & helix2_twist
)
{
	Size const repeatpos( repeatlen/2 );

	Vector center, t, n; // n is the symmetry axis unit vector
	Size const pos1( (base_repeat-1)*repeatlen + repeatpos ), pos2( base_repeat * repeatlen + repeatpos );
	Residue const & rsd1( pose.residue( pos1 ) ), &rsd2( pose.residue(pos2) );
	Stub const
		stub1( rsd1.xyz("N"), rsd1.xyz("CA"), rsd1.xyz("C") ),
		stub2( rsd2.xyz("N"), rsd2.xyz("CA"), rsd2.xyz("C") );

	Real theta;
	get_stub_transform_data( stub1, stub2, center, n, t, theta );

	Reals hdists, htwists;

	Size const base_repeat_offset( ( base_repeat-1)*repeatlen );

	for ( Size r=1; r<= 2; ++r ) {
		Size const hbegin( r == 1 ? base_repeat_offset + 1 : base_repeat_offset + helix1_len + turn1_len + 1 ),
			hend( r == 1 ? hbegin + helix1_len-1 : hbegin + helix2_len-1 );
		Vectors coords;
		string ss;
		string const bb( torsion2big_bin_string( hbegin, hend, pose ) );
		for ( Size i=hbegin; i<= hend; ++i ) {
			coords.push_back( pose.residue(i).xyz("CA") );
			ss.push_back( pose.secstruct(i) );
		}
		Vector axis, midpoint;
		bool const is_a_helix( ( std::count( ss.begin(), ss.end(), 'H' ) > std::count( ss.begin(), ss.end(), 'E' ) ) ||
			( std::count( bb.begin(), bb.end(), 'A' ) > std::count( bb.begin(), bb.end(), 'B' ) ) );

		if ( is_a_helix && coords.size() > 4 ) {
			compute_helix_axis_and_midpoint( coords, axis, midpoint );
		} else {
			midpoint = 0.5 * ( coords.front() + coords.back() );
			axis = coords.back() - coords.front();
		}

		/// distance from midpoint to the symmetry axis
		Vector radiusv( midpoint - center );
		radiusv -= n.dot( radiusv ) * n; // make normal to rotation axis vector
		hdists.push_back( radiusv.length() );

		/// angle of helix axis twist
		Vector const p1( midpoint + axis ), p2( midpoint ), p3( midpoint - radiusv ),
			p4( ( axis.dot(n)>0 ) ? p3 + n : p3 - n );
		runtime_assert( fabs( n.dot( (center - p3).normalized() ) )>0.99 ); // confirm p3 is on the symmetry axis

		htwists.push_back( numeric::dihedral_degrees( p1, p2, p3, p4 ) );
	}

	helix1_dist  = hdists [1];
	helix1_twist = htwists[1];

	helix2_dist  = hdists [2];
	helix2_twist = htwists[2];

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// assumes that each repeat goes: helix1, turn1, helix2, turn2
///
/// but we don't know which is the inner helix and which is the outer helix...
///
void
compute_helix_axis_angles(
	Vectors const & ca_coords,
	Size const base_repeat,
	Size const repeatlen,
	Size const helix1_len,
	Size const turn1_len,
	Size const helix2_len,
	//Size const turn2_len,
	Real & helix1_dist,
	Real & helix1_twist,
	Real & helix2_dist,
	Real & helix2_twist
)
{
	Size const repeatpos( repeatlen/2 );

	Vector center, t, n; // n is the symmetry axis unit vector
	Size const pos1( (base_repeat-1)*repeatlen + repeatpos ), pos2( base_repeat * repeatlen + repeatpos );
	//Residue const & rsd1( pose.residue( pos1 ) ), &rsd2( pose.residue(pos2) );
	Stub const
		stub1( ca_coords[pos1], ca_coords[pos1+1], ca_coords[pos1+2] ),
		stub2( ca_coords[pos2], ca_coords[pos2+1], ca_coords[pos2+2] );

	Real theta;
	get_stub_transform_data( stub1, stub2, center, n, t, theta );

	Reals hdists, htwists;

	Size const base_repeat_offset( ( base_repeat-1)*repeatlen );

	for ( Size r=1; r<= 2; ++r ) {
		Size const hbegin( r == 1 ? base_repeat_offset + 1 : base_repeat_offset + helix1_len + turn1_len + 1 ),
			hend( r == 1 ? hbegin + helix1_len-1 : hbegin + helix2_len-1 );
		Vectors coords;
		for ( Size i=hbegin; i<= hend; ++i ) coords.push_back( ca_coords[i] );
		Vector axis, midpoint;
		if ( coords.size() > 4 ) {
			compute_helix_axis_and_midpoint( coords, axis, midpoint );
		} else {
			midpoint = 0.5 * ( coords.front() + coords.back() );
			axis = coords.back() - coords.front();
		}

		/// distance from midpoint to the symmetry axis
		Vector radiusv( midpoint - center );
		radiusv -= n.dot( radiusv ) * n; // make normal to rotation axis vector
		hdists.push_back( radiusv.length() );

		/// angle of helix axis twist
		Vector const p1( midpoint + axis ), p2( midpoint ), p3( midpoint - radiusv ),
			p4( ( axis.dot(n)>0 ) ? p3 + n : p3 - n );
		runtime_assert( fabs( n.dot( (center - p3).normalized() ) )>0.99 ); // confirm p3 is on the symmetry axis

		htwists.push_back( numeric::dihedral_degrees( p1, p2, p3, p4 ) );
	}

	helix1_dist  = hdists [1];
	helix1_twist = htwists[1];

	helix2_dist  = hdists [2];
	helix2_twist = htwists[2];

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
