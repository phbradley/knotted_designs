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
#ifndef INCLUDED_apps_pilot_phil_helix_frags_HH
#define INCLUDED_apps_pilot_phil_helix_frags_HH

#include <apps/pilot/phil/phil.hh>
#include <apps/pilot/phil/phil_options.hh>
#include <apps/pilot/phil/rms.hh>
#include <apps/pilot/phil/sym_basic.hh>
#include <apps/pilot/phil/types.hh>

#include <utility/io/izstream.hh>

#include <numeric/xyzVector.io.hh>


static basic::Tracer TR_HELIX_FRAGS( "apps.pilot.phil.helix_frags_hh" );


namespace helix_frags {
Size const helix_window( 7 );
Size const strand_window( 4 );
Real const hh_contact_dis( 8.0 ), hh_contact_dis_for_library( 9.0 );
Real const ss_contact_dis( 7.0 ), ss_contact_dis_for_library( 8.0 );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Stub
get_helix_stub(
	Vectors const & coords,
	Size const helixbegin,
	Size const window = helix_frags::helix_window
)
{
	static bool init( false );
	static Vectors ideal_coords;
	static Vector ideal_p1, ideal_p2, ideal_p3;

	runtime_assert( window <= helix_frags::helix_window );

	//static Size const helix_window( helix_frags::helix_window );
	if ( !init ) {
		init = true;
		Pose pose;
		string const perfect_helix_file( option[ my_options::perfect_helix_file ] );
		pose_from_pdb( pose, perfect_helix_file );
		runtime_assert( pose.total_residue() > helix_frags::helix_window );
		Vectors tmpcoords;
		for ( Size i=1; i<= helix_frags::helix_window; ++i ) {
			ideal_coords.push_back( pose.residue(i  ).xyz("CA") );
			tmpcoords   .push_back( pose.residue(i+1).xyz("CA") );
		}

		/// define the central axis of the helix
		Stub const stub1( tmpcoords[1], tmpcoords[2], tmpcoords[3] );
		superimpose_coords( ideal_coords, tmpcoords );
		Stub const stub2( tmpcoords[1], tmpcoords[2], tmpcoords[3] );

		Vector center, axis, t;
		Real theta;
		Real const dev( get_stub_transform_data( stub1, stub2, center, axis, t, theta ) );
		TR_HELIX_FRAGS.Trace << "get_helix_stub: helix_axis: dev= " << F(9,3,dev) <<
			" trans= " << F( 9,3, t.dot(axis) ) <<
			" theta= " << F( 9,3, numeric::conversions::degrees( theta ) ) << endl;

		// project ideal_coords[1] onto the axis
		center += ( (ideal_coords[1]-center).dot(axis)*axis );
		runtime_assert( fabs( center.dot(axis) - ideal_coords[1].dot(axis) )<1e-3 );

		if ( axis.dot( ideal_coords[2]-ideal_coords[1] )<0 ) axis *= -1;

		// get local coords for some points that will define a helix stub
		Stub const stub( ideal_coords[1], ideal_coords[2], ideal_coords[3] );
		ideal_p1 = stub.global2local( center+axis );
		ideal_p2 = stub.global2local( center );
		//ideal_p3 = stub.global2local( ideal_coords[1] );
	}

	runtime_assert( coords.size() >= helixbegin + window-1 );
	runtime_assert( ideal_coords.size() >= window );

	/// superimpose ideal helix onto coords in the window
	Vectors idl_coords( ideal_coords ); idl_coords.resize( window );
	Vectors mod_coords;
	for ( Size i=1; i<= window; ++i ) mod_coords.push_back( coords[ helixbegin+i-1] );
	superimpose_coords( mod_coords, idl_coords );

	Real rmsd(0);
	for ( Size i=1; i<= window; ++i ) rmsd+= mod_coords[i].distance_squared( idl_coords[i] );
	rmsd = sqrt( rmsd/window );

	//TR.Trace << "get_helix_stub:  rmsd= " << F(9,3,rmsd) << endl;

	Stub const stub( idl_coords[1], idl_coords[2], idl_coords[3] );

	// the x-column of this stub will be parallel with the axis (p1-p2)
	// the y-column will point towards idl_coords[1]
	// the center of the stub will be at p1, which is along the axis, aligned with idl_coords+axis
	//   (ie, shifted fwd from idl_coords[1] position by axis(length=1))
	return kinematics::Stub( stub.local2global( ideal_p1 ), stub.local2global( ideal_p2 ), idl_coords[1] );

	// old way:
	// /// now compute a stub from idl_coords
	// Vector bxyz( 0,0,0 ), exyz( 0,0,0 );
	// for ( Size i=1; i<=4; ++i ) {
	//  bxyz += 0.25 * idl_coords[i];
	//  exyz += 0.25 * idl_coords[helix_window-i+1];
	// }
	// Vector const axis( ( exyz - bxyz ).normalized() );
	// return kinematics::Stub( bxyz+axis, bxyz, idl_coords[1] ); // pointing along the helix axis...
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Stub
get_strand_stub(
	Vectors const & coords,
	Size const strandbegin
)
{
	//Size const strand_window( 4 );
	Vector const & a( coords[ strandbegin ] ), &b( coords[ strandbegin+1 ]  ), &c( coords[ strandbegin+2 ] ),
		&d( coords[ strandbegin+3 ] );
	Vector const ab( 0.5*(a+b) ), cd( 0.5*(c+d) );//, axis( (cd-ab).normalized() );

	// centered at cd, x-axis running from ab to cd, y-axis
	// { // debugging
	//  Stub const stub( cd, ab, ab + (a-b) + (c-d) );
	//  TR_HELIX_FRAGS.Trace << "strand_stub: " <<
	//   " g2l(ab)= " << stub.global2local( ab ) <<
	//   " g2l(cd)= " << stub.global2local( cd ) <<
	//   " g2l(a)= " << stub.global2local( a ) <<
	//   " g2l(b)= " << stub.global2local( b ) <<
	//   " g2l(c)= " << stub.global2local( c ) <<
	//   " g2l(d)= " << stub.global2local( d ) << endl;
	// }
	return Stub( cd, ab, ab + (a-b) + (c-d) );

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct HelixHelixFrag {
	Vectors h1_coords;
	Vectors h2_coords;
	RT rt_12;
	RT rt_21;
	Size nsim;
	Real mindis2;
};

Real
helix_helix_frag_distance(
	HelixHelixFrag const & a,
	HelixHelixFrag const & b
)
{
	Real msd1(0.0), msd2(0.0), msd3(0.0), msd4(0.0);
	for ( Size i=1, i_e=a.h1_coords.size(); i<= i_e; ++i ) {
		msd1 += a.h1_coords[i].distance_squared( b.h1_coords[i] );
		msd2 += a.h2_coords[i].distance_squared( b.h2_coords[i] );
		msd3 += a.h1_coords[i].distance_squared( b.h2_coords[i] );
		msd4 += a.h2_coords[i].distance_squared( b.h1_coords[i] );
	}
	return sqrt( min( msd1, min( msd2, min( msd3, msd4 ) ) ) / a.h1_coords.size() );

}


class HelixHelixFragDistanceMetric {
public:
	Real
	operator() ( HelixHelixFrag const & a, HelixHelixFrag const & b ) const
	{
		return helix_helix_frag_distance( a, b );
	}

};

typedef vector1< HelixHelixFrag > HelixHelixFrags;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// bool // true if we found a contact
// get_helix_contact_frag(
//             Vectors const & h1_coords,
//             Vectors const & h2_coords,
//             HelixHelixFrag & f,
//             Real const contact_dis = 8.0
//             )
// {
//  static Size const helix_window( helix_frags_helix_window );
//  runtime_assert( h1_coords.size() >= helix_window && h2_coords.size() >= helix_window );
//  Real const contact_dis2( numeric::square( contact_dis ) );
//  Real dis2, mindis2( 1e6 );
//  Size ibest(0), jbest(0);
//  for ( Size i=1; i<= h1_coords.size(); ++i ) {
//   for ( Size j=1; j<= h2_coords.size(); ++j ) {
//    dis2 = h1_coords[i].distance_squared( h2_coords[j] );
//    if ( dis2 < mindis2 ) {
//     mindis2 = dis2;
//     ibest = i;
//     jbest = j;
//    }
//   }
//  }
//  if ( mindis2 < contact_dis2 ) {
//   //TR.Trace << "HH_contact " << F(9,3,sqrt( mindis2 ) ) << endl;
//   /// get stubs around these guys
//   Size const h1_start( min( max( 1, int(ibest)-3 ), int(h1_coords.size())-int(helix_window)+1 ) );
//   Size const h2_start( min( max( 1, int(jbest)-3 ), int(h2_coords.size())-int(helix_window)+1 ) );
//   Stub const h1_stub( get_helix_stub( h1_coords, h1_start ) ), h2_stub( get_helix_stub( h2_coords, h2_start ) );
//   f.h1_coords.clear(); f.h2_coords.clear();
//   for ( Size i=1; i<= helix_window; ++i ) {
//    f.h1_coords.push_back( h2_stub.global2local( h1_coords[ h1_start+i-1 ] ) );
//    f.h2_coords.push_back( h1_stub.global2local( h2_coords[ h2_start+i-1 ] ) );
//   }
//   f.rt_12 = RT( h1_stub, h2_stub );
//   f.rt_21 = RT( h2_stub, h1_stub );
//   f.mindis2 = mindis2;
//   return true;
//  } else {
//   return false;
//  }
// }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size // num rsd contacts: 0 if no contacts, positive otherwise
get_sse_contact_frag(
	Vectors const & h1_coords,
	Vectors const & h2_coords,
	SegmentType const & segtype,
	Real const contact_dis,
	HelixHelixFrag & f
)
{
	runtime_assert( segtype == alpha_type || segtype == beta_type );
	Size const window( segtype == alpha_type ? helix_frags::helix_window : helix_frags::strand_window );
	if ( h1_coords.size()<window || h2_coords.size()<window ) return false;
	Real const contact_dis2( numeric::square( contact_dis ) );
	Real dis2, mindis2( 1e6 );
	Size ibest(0), jbest(0), num_contacts(0);
	for ( Size i=1; i<= h1_coords.size(); ++i ) {
		for ( Size j=1; j<= h2_coords.size(); ++j ) {
			dis2 = h1_coords[i].distance_squared( h2_coords[j] );
			if ( dis2 < contact_dis2 ) {
				++num_contacts;
				if ( dis2 < mindis2 ) {
					mindis2 = dis2;
					ibest = i;
					jbest = j;
				}
			}
		}
	}
	if ( mindis2 < contact_dis2 ) {
		//TR.Trace << "HH_contact " << F(9,3,sqrt( mindis2 ) ) << endl;
		/// get stubs around these guys
		Size const h1_start( min( max( 1, int(ibest)-int(window)/2 ), int(h1_coords.size())-int(window)+1 ) );
		Size const h2_start( min( max( 1, int(jbest)-int(window)/2 ), int(h2_coords.size())-int(window)+1 ) );
		Stub const h1_stub( segtype == alpha_type ?
			get_helix_stub ( h1_coords, h1_start ) :
			get_strand_stub( h1_coords, h1_start ) );
		Stub const h2_stub( segtype == alpha_type ?
			get_helix_stub ( h2_coords, h2_start ) :
			get_strand_stub( h2_coords, h2_start ) );
		f.h1_coords.clear(); f.h2_coords.clear();
		for ( Size i=1; i<= window; ++i ) {
			f.h1_coords.push_back( h2_stub.global2local( h1_coords[ h1_start+i-1 ] ) );
			f.h2_coords.push_back( h1_stub.global2local( h2_coords[ h2_start+i-1 ] ) );
		}
		f.rt_12 = RT( h1_stub, h2_stub );
		f.rt_21 = RT( h2_stub, h1_stub );
		f.mindis2 = mindis2;
		return num_contacts;
	} else {
		return 0;
	}
}

bool
get_strand_contact_frag(
	Vectors const & h1_coords,
	Vectors const & h2_coords,
	HelixHelixFrag & f,
	Real const contact_dis = helix_frags::ss_contact_dis
)
{
	return get_sse_contact_frag( h1_coords, h2_coords, beta_type, contact_dis, f );
}


bool
get_helix_contact_frag(
	Vectors const & h1_coords,
	Vectors const & h2_coords,
	HelixHelixFrag & f,
	Real const contact_dis = helix_frags::hh_contact_dis
)
{
	return get_sse_contact_frag( h1_coords, h2_coords, alpha_type, contact_dis, f );
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool // true if we found a contact
get_helix_contact_coords(
	Size const helix1_begin,
	Size const helix1_end,
	Vectors const & h1_coords,
	Size const helix2_begin,
	Size const helix2_end,
	Vectors const & h2_coords,
	Vectors & h1_contact_coords,
	Vectors & h2_contact_coords
)
{
	Size const helix_window( helix_frags::helix_window );
	runtime_assert( helix1_end-helix1_begin+1 >= helix_window && helix2_end-helix2_begin+1 >= helix_window );
	Real const contact_dis2( 8*8 );
	Real dis2, mindis2( 1e6 );
	Size ibest(0), jbest(0);
	for ( Size i=helix1_begin; i<= helix1_end; ++i ) {
		for ( Size j=helix2_begin; j<= helix2_end; ++j ) {
			dis2 = h1_coords[i].distance_squared( h2_coords[j] );
			if ( dis2 < mindis2 ) {
				mindis2 = dis2;
				ibest = i;
				jbest = j;
			}
		}
	}
	if ( mindis2 < contact_dis2 ) {
		//TR.Trace << "HH_contact " << F(9,3,sqrt( mindis2 ) ) << endl;
		/// get stubs around these guys
		Size const h1_start( min( max( int(helix1_begin), int(ibest)-3 ), int(helix1_end)-int(helix_window)+1 ) );
		Size const h2_start( min( max( int(helix2_begin), int(jbest)-3 ), int(helix2_end)-int(helix_window)+1 ) );
		// IS THIS TOO SLOW? MAYBE COME UP WITH A FASTER GET_HELIX_STUB ROUTINE??
		Stub const h1_stub( get_helix_stub( h1_coords, h1_start ) ), h2_stub( get_helix_stub( h2_coords, h2_start ) );
		h1_contact_coords.clear(); h2_contact_coords.clear();
		for ( Size i=1; i<= helix_window; ++i ) {
			h1_contact_coords.push_back( h2_stub.global2local( h1_coords[ h1_start+i-1 ] ) );
			h2_contact_coords.push_back( h1_stub.global2local( h2_coords[ h2_start+i-1 ] ) );
		}
		return true;
	} else {
		return false;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool // true if we found a contact
get_helix_contact_coords(
	Vectors const & h1_coords,
	Vectors const & h2_coords,
	Vectors & h1_contact_coords,
	Vectors & h2_contact_coords
)
{
	return get_helix_contact_coords( 1, h1_coords.size(), h1_coords, 1, h2_coords.size(), h2_coords,
		h1_contact_coords, h2_contact_coords );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
build_helix_transform_library(
	HelixHelixFrags & all_frags,
	bool const include_sequence_neighbors = false
)
{
	static bool init( false );

	Size const helix_window( helix_frags::helix_window );

	static vector1< Vectors > all_helix_coords;
	static map< string, Sizes > pdb_helices;

	if ( !init ) { //read info on bb coords
		init = true;
		string const pdb_coords_file( option[ my_options::pdb_coords_file ] );
		ifstream data( pdb_coords_file.c_str() );
		string line;
		while ( getline( data, line ) ) {
			strings const l( split_to_vector1( line ) );
			if ( l[1] == "helix_coords:" || l[1] == "my_helix_coords:" ) {
				Size const helixlen( int_of( l[2] ) );
				Vectors coords;
				for ( Size i=1; i<= helixlen; ++i ) {
					coords.push_back( Vector( float_of( l[ 3*i+1 ] ), float_of( l[ 3*i+2 ] ), float_of( l[3*i+3 ] ) ) );
				}
				string const pdbchain( l[ l.size() ] +"_"+l[ l.size()-1] );
				all_helix_coords.push_back( coords );
				pdb_helices[ pdbchain ].push_back( all_helix_coords.size() );
			}
		}
		data.close();
		TR_HELIX_FRAGS.Trace << "Read " << all_helix_coords.size() << " helix coords from file " << pdb_coords_file << endl;
	} // initialize coords



	// look for helix-helix packing interactions
	HelixHelixFrag f;
	for ( map< string, Sizes >::const_iterator p= pdb_helices.begin(); p != pdb_helices.end(); ++p ) {
		Sizes const & helix_indices( p->second );
		if ( helix_indices.size()>=3 ) {
			for ( Size ii=1; ii<= helix_indices.size(); ++ii ) {
				Size const h1( helix_indices[ii] );
				Vectors const & h1_coords( all_helix_coords[h1] );
				if ( h1_coords.size() < helix_window ) continue;
				Size const jj_start( include_sequence_neighbors ? ii+1 : ii+2 );
				for ( Size jj=jj_start; jj<= helix_indices.size(); ++jj ) { // NOTE: not taking ii->ii+1 contacts
					// look for contacts between these helices
					Size const h2( helix_indices[jj] );
					Vectors const & h2_coords( all_helix_coords[h2] );
					if ( h2_coords.size() < helix_window ) continue;
					bool const found_contact
						( get_helix_contact_frag( h1_coords, h2_coords, f,
						helix_frags::hh_contact_dis_for_library ) ); // slightly large than dflt
					if ( found_contact ) {
						// TR_HELIX_FRAGS.Trace << "HH_contact_info mindis: " << F(9,3,sqrt(f.mindis2)) << ' ' <<
						//  all_frags.size()+1 << ' ' << ii << ' ' << jj << ' ' << p->first <<endl;
						all_frags.push_back( f );
					}
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
build_strand_transform_library(
	HelixHelixFrags & all_frags
)
{
	static bool init( false );

	Size const strand_window( helix_frags::strand_window );

	static vector1< Vectors > all_strand_coords;
	static map< string, Sizes > pdb_strands;

	if ( !init ) { //read info on bb coords
		init = true;
		string const pdb_coords_file( option[ my_options::pdb_coords_file ] );
		ifstream data( pdb_coords_file.c_str() );
		string line;
		while ( getline( data, line ) ) {
			strings const l( split_to_vector1( line ) );
			if ( l[1] == "strand_coords:" || l[1] == "my_strand_coords:" ) {
				Size const strandlen( int_of( l[2] ) );
				Vectors coords;
				for ( Size i=1; i<= strandlen; ++i ) {
					coords.push_back( Vector( float_of( l[ 3*i+1 ] ), float_of( l[ 3*i+2 ] ), float_of( l[3*i+3 ] ) ) );
				}
				string const pdbchain( l[ l.size() ] +"_"+l[ l.size()-1] );
				all_strand_coords.push_back( coords );
				pdb_strands[ pdbchain ].push_back( all_strand_coords.size() );
			}
		}
		data.close();
		TR_HELIX_FRAGS.Trace << "Read " << all_strand_coords.size() << " strand coords from file " << pdb_coords_file <<
			endl;
	} // initialize coords



	// look for helix-helix packing interactions
	HelixHelixFrag f;
	for ( map< string, Sizes >::const_iterator p= pdb_strands.begin(); p != pdb_strands.end(); ++p ) {
		Sizes const & strand_indices( p->second );
		if ( strand_indices.size()>=3 ) {
			for ( Size ii=1; ii<= strand_indices.size(); ++ii ) {
				Size const h1( strand_indices[ii] );
				Vectors const & h1_coords( all_strand_coords[h1] );
				if ( h1_coords.size() < strand_window ) continue;
				for ( Size jj=ii+2; jj<= strand_indices.size(); ++jj ) {
					// look for contacts between these strands
					Size const h2( strand_indices[jj] );
					Vectors const & h2_coords( all_strand_coords[h2] );
					if ( h2_coords.size() < strand_window ) continue;
					bool const found_contact
						( get_strand_contact_frag( h1_coords, h2_coords, f,
						helix_frags::ss_contact_dis_for_library ) ); // slightly large than dflt
					if ( found_contact ) {
						// TR_HELIX_FRAGS.Trace << "SS_contact_info mindis: " << F(9,3,sqrt(f.mindis2)) << ' ' <<
						//  all_frags.size()+1 << ' ' << ii << ' ' << jj << ' ' << p->first <<endl;
						all_frags.push_back( f );
					}
				}
			}
		}
	}
	if ( false ) {
		Size counter(0);
		foreach_ ( HelixHelixFrag const & f1, all_frags ) {
			++counter;
			Size nsim(0);
			Real mindis(1e6), dis;
			foreach_ ( HelixHelixFrag const & f2, all_frags ) {
				dis = helix_helix_frag_distance( f1, f2 );
				if ( dis<1e-3 ) continue; // same frag
				if ( dis <= 3 ) ++nsim;
				mindis = min( mindis, dis );
			}
			TR_HELIX_FRAGS.Trace << "build_strand_transform_library: nsim: " << I(4,nsim) << " mindis: " << F(9,3,mindis) <<
				I(6,counter) << I(6,all_frags.size()) << endl;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
load_helix_coords_from_file(
	string const filename,
	bool const clean,
	vector1< Vectors > & all_helix_coords
)
{
	runtime_assert( utility::file::file_exists( filename ) );

	utility::io::izstream data( filename );
	string line;
	while ( getline( data, line ) ) {
		strings const l( split_to_vector1( line ) );
		if ( ( clean && l[1] == "clean_helix_coords:" ) || ( !clean && l[1] == "helix_coords:" ) ) {
			Size const helixlen( int_of( l[2] ) );
			Vectors coords;
			for ( Size i=1; i<= helixlen; ++i ) {
				coords.push_back( Vector( float_of( l[ 3*i+1 ] ), float_of( l[ 3*i+2 ] ), float_of( l[3*i+3 ] ) ) );
			}
			all_helix_coords.push_back( coords );
		}
	}

	data.close();
	TR_HELIX_FRAGS.Trace << "Read " << all_helix_coords.size() << " helix coords from file " << filename <<
		" clean= " << clean << endl;

	runtime_assert( all_helix_coords.size() );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
compute_helix_strain(
	Size const helixbegin,
	Size const helixend,
	Vectors const & ca_coords,
	bool const use_clean_helix_coords,
	Real & avg_avg_rmsd,
	Real & avg_med_rmsd
)
{
	Size const helix_window( 7 );

	runtime_assert( helixend>= helixbegin);

	static vector1< Vectors > all_unclean_helix_coords, all_clean_helix_coords;

	vector1< Vectors > & all_helix_coords( use_clean_helix_coords ? all_clean_helix_coords : all_unclean_helix_coords );


	if ( all_helix_coords.empty() ) {
		load_helix_coords_from_file( option[ my_options::pdb_coords_file ], use_clean_helix_coords, all_helix_coords );
		numeric::random::random_permutation( all_helix_coords, numeric::random::rg() );
	}


	Size const helixlen( helixend - helixbegin + 1 );

	avg_avg_rmsd = avg_med_rmsd = 0;

	if ( helixlen < helix_window ) return;

	Size const big_helixlen( 20 );
	for ( Size i=1; i<= helixlen - helix_window+1; ++i ) {
		Vectors coords;
		for ( Size pos = helixbegin+i-1; pos<= helixbegin+i-1+helix_window-1; ++pos ) {
			coords.push_back( ca_coords[pos] );
		}
		Size const max_rmsds( 1000 ); // consider increasing??
		Reals rmsds;
		for ( Size ii=1; ii<= all_helix_coords.size(); ++ii ) {
			Size const pdb_helixlen( all_helix_coords[ii].size() );
			if ( pdb_helixlen >= i+helix_window-1 || pdb_helixlen >= big_helixlen ) {
				//cout << i << ' ' << ii << ' ' << pdb_helixlen << endl;
				Vectors pdb_coords;
				Size const pdb_start( ( pdb_helixlen >= i+helix_window-1 ) ? i : pdb_helixlen/2 - 4 );
				runtime_assert( pdb_helixlen >= pdb_start+helix_window-1 );
				for ( Size j=pdb_start; j<= pdb_start+helix_window-1; ++j ) pdb_coords.push_back( all_helix_coords[ii][j] );
				runtime_assert( coords.size() == pdb_coords.size() );
				rmsds.push_back( numeric::model_quality::calc_rms( coords, pdb_coords ) );
				if ( rmsds.size() >= max_rmsds ) break;
			}
		}

		std::sort( rmsds.begin(), rmsds.end() );

		Real const median_rmsd( rmsds[ rmsds.size()/2 ] );

		Size const n_outliers( rmsds.size()/15 );
		runtime_assert( n_outliers < rmsds.size() );
		for ( Size i=1; i<= n_outliers; ++i ) rmsds.pop_back();

		Real avg_rmsd( 0 );
		for ( Size i=1; i<= rmsds.size(); ++i ) avg_rmsd += rmsds[i];
		avg_rmsd/= rmsds.size();

		avg_avg_rmsd += avg_rmsd;
		avg_med_rmsd += median_rmsd;
	} // loop over helix windows
	Size const nwindows( helixlen - helix_window + 1 );
	avg_avg_rmsd /= nwindows;
	avg_med_rmsd /= nwindows;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
compute_helix_strain(
	Size const helixbegin,
	Size const helixend,
	Vectors const & ca_coords,
	Real & avg_avg_rmsd,
	Real & avg_med_rmsd
)
{
	bool const use_clean_helix_coords( false ); // old behavior
	compute_helix_strain( helixbegin, helixend, ca_coords, use_clean_helix_coords, avg_avg_rmsd, avg_med_rmsd );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
string
analyze_helix_coords(
	string const prefix,
	Size const helix_begin,
	Size const helix_end,
	Vectors const & coords,
	string const tag=string()
)
{
	Size const helix_len( helix_end-helix_begin+1), window(7); // maybe make smaller?

	//
	vector1< Stub > helix_stubs;
	for ( Size i=helix_begin; i<= helix_end-window+1; ++i ) {
		helix_stubs.push_back( get_helix_stub( coords, i, window ) );
	}

	// figure out total rise, rotation, etc
	Size const nstubs( helix_stubs.size() );
	runtime_assert( nstubs == helix_len-window+1 );

	Real total_twist(0), total_rise(0), total_bend(0), max_bend(0);
	for ( Size i=1; i<nstubs; ++i ) {
		// figure out rise, twist, bend
		// rise is just distance between the stub origins ('v'-vectors)
		// first take the bend out (align the x-axes)
		// then the twist is the rotation about the x-axis

		Stub stub1( helix_stubs[i] );
		Stub const & stub2( helix_stubs[i+1] );

		// how to rotate stub1 x-col to align with stub2 x-col?
		Vector const bend_axis( ( stub1.M.col_x().cross( stub2.M.col_x() ) ).normalized_any() );
		Real const bend_angle( std::acos( numeric::sin_cos_range( stub1.M.col_x().dot( stub2.M.col_x() ) ) ) );
		Matrix const R1( rotation_matrix( bend_axis, bend_angle ) );

		stub1.M = R1 * stub1.M;

		runtime_assert( stub1.M.col_x().dot( stub2.M.col_x() ) > 0.999 ); // should be parallel

		Vector const twist_axis( ( stub1.M.col_y().cross( stub2.M.col_y() ) ).normalized_any() );
		Real twist_angle( std::acos( numeric::sin_cos_range( stub1.M.col_y().dot( stub2.M.col_y() ) ) ) );
		//Matrix const R1( rotation_matrix( bend_axis, bend_angle ) );

		runtime_assert( twist_angle>-1e-3 );
		if ( twist_angle > 1e-3 ) {
			Real const dp( twist_axis.dot(stub1.M.col_x()) );
			runtime_assert( dp<-0.999 || dp>0.999 );
			if ( dp<0 ) twist_angle *= -1;
		}

		Real const rise( stub1.v.distance(stub2.v) );
		total_rise  += rise;
		total_twist += twist_angle;
		total_bend  += bend_angle;

		max_bend = max( max_bend, bend_angle );

		if ( tag.size() ) {
			TR_HELIX_FRAGS.Trace << "analyze_helix_coords: " << I(4,helix_begin+i-1 ) <<
				" rise: "  << F(9,3,rise ) <<
				" twist: " << F(9,3,degrees(twist_angle)) <<
				" bend: "  << F(9,3,degrees(bend_angle)) <<
				" bend_axis: "  << F(9,3,bend_axis.x()) << F(9,3,bend_axis.y()) << F(9,3,bend_axis.z()) << ' ' << tag << endl;
		}
	}

	Real avg_dihedral_1254(0);
	for ( Size i=helix_begin; i+4<=helix_end; ++i ) {
		// dihedral angle that captures something about twist
		Real const dihedral( numeric::dihedral_degrees( coords[i], coords[i+1], coords[i+4], coords[i+3] ) );
		avg_dihedral_1254 += dihedral;
		if ( tag.size() ) {
			TR_HELIX_FRAGS.Trace << "dihedral_1254: " << F(9,3,dihedral) << ' ' << i << ' ' << tag << endl;
		}
	}
	avg_dihedral_1254 /= ( helix_end-helix_begin-3 );

	ostringstream out;
	out <<
		prefix << "_avg_rise: "  << F(9,3,total_rise/(nstubs-1)) << ' ' <<
		prefix << "_avg_twist: " << F(9,3,degrees(total_twist/(nstubs-1))) << ' ' <<
		prefix << "_avg_bend: "  << F(9,3,degrees(total_bend/(nstubs-1))) << ' ' <<
		prefix << "_max_bend: "  << F(9,3,degrees(max_bend)) << ' ' <<
		prefix << "_avg_dihedral_1254: " << F(9,3,avg_dihedral_1254);

	return out.str();
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


#endif
