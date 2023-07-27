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
#ifndef INCLUDED_apps_pilot_phil_stub_frag_HH
#define INCLUDED_apps_pilot_phil_stub_frag_HH

#include <apps/pilot/phil/phil.hh>
#include <apps/pilot/phil/phil_options.hh>
//#include <apps/pilot/phil/phil_io.hh>
#include <apps/pilot/phil/types.hh>
#include <apps/pilot/phil/rms.hh>


#include <apps/pilot/phil/stub_transform.hh>
#include <core/optimization/Multifunc.hh>
#include <numeric/random/random_permutation.hh>
#include <numeric/model_quality/rms.hh>
#include <fstream>


static basic::Tracer TR_STUB_FRAG( "apps.pilot.phil.stub_frag_hh" );


Size
get_segment_window( SegmentType const & t )
{
	switch (t ) {
	case alpha_type :
		return 7;
	case beta_type :
		return 4;
	default :
		return 0;
	}
	return 0;
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct HelixParams {
	Real rise; // Angstroms ~1.5A
	Real twist; // radians ~99 degrees
	Real tilt; // radians -- small
	Real tilt_direction; // radians -- 0 to 2*pi
	Real tilt_precession; // radians -- small

	Real ca_distance; // not sure about this guy...
	HelixParams():
		rise(0),
		twist(0),
		tilt(0),
		tilt_direction(0),
		tilt_precession(0),
		ca_distance(0)
	{}

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct StubFrag {
	string id;
	RT rt;
	Vectors coords;
	Real raw_strain;
	Real norm_strain;
	HelixParams hparams; // ACK
	vector1< Reals > torsions;
	string sequence;
	SegmentType segtype;
};

typedef utility::vector1< StubFrag > StubFrags;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
Real dummy_tilt_direction;

Stub // returns the final stub
generate_helix_coords(
	Stub const & start_stub,
	HelixParams const & params,
	Size const nres,
	Vectors & coords,
	bool const fill_in_coords = true,
	bool const erase_coords_at_start = true,
	Real & new_tilt_direction = dummy_tilt_direction
)
{
	// start stub is the stub that places the first c-alpha
	if ( fill_in_coords && erase_coords_at_start ) coords.clear();

	// define the direction in which we tilt, at the start...
	Vector tilt_axis( cos( params.tilt_direction ) * start_stub.M.col_y() +
		sin( params.tilt_direction ) * start_stub.M.col_z() );

	Stub stub( start_stub );
	for ( Size i=1; i<= nres; ++i ) {
		Vector const local_helix_axis( stub.M.col_x() );
		runtime_assert( fabs( tilt_axis.dot( local_helix_axis ) ) < 1e-2 );
		if ( false ) { // hacking
			Vector const point( stub.v ), normal( local_helix_axis );
			TR_STUB_FRAG.Trace << "a_points: " << I(3,i) <<
				" P: " << F(9,3,point.x()) << F(9,3,point.y()) << F(9,3,point.z()) <<
				" N: " << F(9,3,normal.x()) << F(9,3,normal.y()) << F(9,3,normal.z()) <<
				" T: " << F(9,3,tilt_axis.x()) << F(9,3,tilt_axis.y()) << F(9,3,tilt_axis.z()) << endl;
			Real const new_tilt_direction( atan2( tilt_axis.dot( stub.M.col_z() ),
				tilt_axis.dot( stub.M.col_y() ) ) );
			// Real const delta_i( i == 1 ? 0.0 :
			//           basic::subtract_radian_angles( new_tilt_direction, params.tilt_direction )/(i-1));
			Real const expected_tilt_direction
				( basic::periodic_range( params.tilt_direction + (i-1)*(params.tilt_precession-params.twist),
				2*numeric::constants::d::pi ) );
			TR_STUB_FRAG.Trace << "new_tilt_dir: " << F(9,3,degrees(new_tilt_direction)) <<
				" old_tilt_dir: " << F(9,3,degrees(params.tilt_direction) ) <<
				" exp_tilt_dir: " << F(9,3,degrees(expected_tilt_direction) ) <<
				" exp_dev: " << F(9,3,degrees(basic::subtract_radian_angles(new_tilt_direction,expected_tilt_direction)) ) <<
				" twist: " << F(9,3,degrees(params.twist)) <<
				" tilt: " << F(9,3,degrees(params.tilt)) <<
				" tilt_prec: " << F(9,3,degrees(params.tilt_precession)) <<
				" i: " << i << ' ' << nres << endl;
		}
		/// make the first C-alpha
		if ( fill_in_coords ) coords.push_back( stub.v + params.ca_distance * stub.M.col_y() );
		/// translate to make the next stub
		stub.v = stub.v + params.rise * local_helix_axis;
		/// rotate about x by the helical twist
		stub.M = numeric::rotation_matrix( local_helix_axis, params.twist ) * stub.M;
		/// rotate in the plane defined by tilt_direction
		stub.M = numeric::rotation_matrix( tilt_axis, params.tilt ) * stub.M;
		/// tilt precession
		tilt_axis = numeric::rotation_matrix( local_helix_axis, params.tilt_precession ) * tilt_axis;
		/// make tilt_axis perpendicular to the new local_helix_axis
		tilt_axis -= stub.M.col_x() * stub.M.col_x().dot( tilt_axis );
		tilt_axis.normalize();
	}
	/// done

	runtime_assert( fabs( tilt_axis.dot( stub.M.col_x() ) ) < 1e-2 );

	new_tilt_direction = atan2( tilt_axis.dot( stub.M.col_z() ),
		tilt_axis.dot( stub.M.col_y() ) );
	return stub;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Stub // returns the final stub
generate_helix_stub(
	Stub const & start_stub,
	HelixParams const & params,
	Size const nres
)
{
	static Vectors coords;
	runtime_assert( coords.empty() );
	return generate_helix_coords( start_stub, params, nres, coords, false, false );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
helix_frag_from_helix_params_and_length(
	Size const helixlen,
	HelixParams const & hparams,
	StubFrag & helix_frag
)
{
	Stub const start_stub; // default
	Stub const stop_stub( generate_helix_coords( start_stub, hparams, helixlen, helix_frag.coords ) );
	helix_frag.rt = RT( start_stub, stop_stub );
	helix_frag.hparams = hparams;
	helix_frag.segtype = alpha_type;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// this is so dumb
void
create_helix_frag_from_helix_params(
	Size const helixlen,
	string const id,
	HelixParams const & hparams,
	StubFrag & f
)
{
	helix_frag_from_helix_params_and_length( helixlen, hparams, f );
	f.id = id;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
shorten_helix_frag_at_nterm(
	Size const ntrim,
	StubFrag & f
)
{
	//StubFrag f_old( f );
	Size const old_nres( f.coords.size() );
	runtime_assert( old_nres > ntrim );
	HelixParams hp( f.hparams );

	Vectors dummy_coords;
	Real new_tilt_direction;
	generate_helix_coords( Stub(), hp, ntrim, dummy_coords, false, false, new_tilt_direction );

	// Real const expected_new_tilt_direction
	//  ( basic::periodic_range( hp.tilt_direction + ntrim*(hp.tilt_precession-hp.twist),
	//               2*numeric::constants::d::pi ) );


	// TR_STUB_FRAG.Trace << "shorten_helix_frag_at_nterm: ntrim= " << ntrim <<
	//  " old_tilt_dir: " << F(9,3,degrees(hp.tilt_direction) ) <<
	//  " new_tilt_dir: " << F(9,3,degrees(new_tilt_direction) ) <<
	//  " exp_tilt_dir: " << F(9,3,degrees(expected_new_tilt_direction) ) <<
	//  " exp_dev: " << F(9,3,degrees(basic::subtract_radian_angles(new_tilt_direction,expected_new_tilt_direction))) <<
	//  endl;

	hp.tilt_direction = new_tilt_direction;

	helix_frag_from_helix_params_and_length( old_nres - ntrim, hp, f );

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// HelixParams dont change for this guy
//
void
shorten_helix_frag_at_cterm(
	Size const ntrim,
	StubFrag & f
)
{
	Size const old_nres( f.coords.size() );
	runtime_assert( old_nres > ntrim );

	HelixParams const hp( f.hparams );
	helix_frag_from_helix_params_and_length( old_nres - ntrim, hp, f );
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class HelicalParamsFitMultifunc : public optimization::Multifunc {
public:
	typedef optimization::Multivec Multivec;

	HelicalParamsFitMultifunc( Vectors const & target_coords ):
		target_coords_( target_coords ),
		optimize_tilt_( false ),
		optimize_tilt_precession_( false )
	{}

	HelicalParamsFitMultifunc():
		target_coords_(),
		optimize_tilt_( false ),
		optimize_tilt_precession_( false )
	{}

	virtual
	Real
	operator()( Multivec const & wts ) const;


	virtual
	void
	dfunc( Multivec const & phipsi, Multivec & dE_dphipsi ) const
	{
		dE_dphipsi.clear(); dE_dphipsi.resize( phipsi.size(), 0.0 );
	}

	void
	set_target_coords( Size const bpos, Size const epos, Vectors const & coords );

	void
	optimize_tilt( bool const setting ) { optimize_tilt_ = setting; }

	void
	optimize_tilt_precession( bool const setting ) { optimize_tilt_precession_ = setting; }

	void
	get_helix_stubs_for_params(
		Multivec const & params,
		Stub & start_stub,
		Stub & stop_stub
	) const;

private:
	Vectors target_coords_;
	bool optimize_tilt_;
	bool optimize_tilt_precession_;


};


void
HelicalParamsFitMultifunc::set_target_coords(
	Size const bpos,
	Size const epos,
	Vectors const & coords
)
{
	target_coords_.clear();
	for ( Size i=bpos; i<= epos; ++i ) target_coords_.push_back( coords[i] );
}



Real
HelicalParamsFitMultifunc::operator()( Multivec const & params ) const
{
	HelixParams hparams;
	runtime_assert( fabs( hparams.tilt_precession )<1e-3 ); // sanity ctor check

	if ( optimize_tilt_precession_ ) {
		runtime_assert( optimize_tilt_ );
		runtime_assert( params.size() == 6 );
		hparams.rise            = params[1]; // 1.5;
		hparams.twist           = params[2]; //radians( 99.0 );
		hparams.tilt            = params[3];
		hparams.tilt_direction  = params[4];
		hparams.tilt_precession = params[5];
		hparams.ca_distance     = params[6];
	} else if ( optimize_tilt_ ) {
		runtime_assert( params.size() == 5 );
		hparams.rise           = params[1]; // 1.5;
		hparams.twist          = params[2]; //radians( 99.0 );
		hparams.tilt           = params[3];
		hparams.tilt_direction = params[4];
		hparams.ca_distance    = params[5]; //1.5; // guess
	} else {
		runtime_assert( params.size() == 3 );
		hparams.rise = params[1]; // 1.5;
		hparams.twist = params[2]; //radians( 99.0 );
		hparams.tilt = 0;
		hparams.tilt_direction = 0;

		hparams.ca_distance = params[3]; //1.5; // guess
	}


	Vectors coords;
	Stub start_stub;
	generate_helix_coords( start_stub, hparams, target_coords_.size(), coords );

	static Real const pi4( 4*numeric::constants::d::pi );
	using numeric::square;
	Real tilt_penalty(0);
	//if ( hparams.tilt           < 0   ) tilt_penalty += 100 * square( hparams.tilt );
	if ( hparams.tilt_direction < -pi4 ) tilt_penalty += 100 * square( pi4 + hparams.tilt_direction );
	if ( hparams.tilt_direction >  pi4 ) tilt_penalty += 100 * square( hparams.tilt_direction - pi4 );
	if ( hparams.tilt_precession < -pi4 ) tilt_penalty += 100 * square( pi4 + hparams.tilt_precession );
	if ( hparams.tilt_precession >  pi4 ) tilt_penalty += 100 * square( hparams.tilt_precession - pi4 );

	return numeric::model_quality::calc_rms( target_coords_, coords ) + tilt_penalty;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
HelicalParamsFitMultifunc::get_helix_stubs_for_params(
	Multivec const & params,
	Stub & start_stub,
	Stub & stop_stub
) const
{
	HelixParams hparams;
	runtime_assert( fabs( hparams.tilt_precession )<1e-3 ); // sanity ctor check

	if ( optimize_tilt_precession_ ) {
		runtime_assert( optimize_tilt_ );
		runtime_assert( params.size() == 6 );
		hparams.rise            = params[1]; // 1.5;
		hparams.twist           = params[2]; //radians( 99.0 );
		hparams.tilt            = params[3];
		hparams.tilt_direction  = params[4];
		hparams.tilt_precession = params[5];
		hparams.ca_distance     = params[6];
		hparams.tilt_precession = 0;
	} else if ( optimize_tilt_ ) {
		runtime_assert( params.size() == 5 );
		hparams.rise           = params[1]; // 1.5;
		hparams.twist          = params[2]; //radians( 99.0 );
		hparams.tilt           = params[3];
		hparams.tilt_direction = params[4];
		hparams.ca_distance    = params[5];
	} else {
		runtime_assert( params.size() == 3 );
		hparams.rise = params[1]; // 1.5;
		hparams.twist = params[2]; //radians( 99.0 );
		hparams.tilt = 0;
		hparams.tilt_direction = 0;

		hparams.ca_distance = params[3]; //1.5; // guess
	}


	Vectors coords;
	Stub start_stub_before_fitting; // default location
	Stub const final_stub_before_fitting( generate_helix_coords( start_stub, hparams, target_coords_.size(), coords ) );

	Stub const coords_stub_before_fitting( Stub( coords[1], coords[2], coords[3] ) );

	Vector const test1a( start_stub_before_fitting.global2local( coords[3] ) ),
		test2a( final_stub_before_fitting.global2local( coords[4] ) );

	// now superimpose
	superimpose_coords( target_coords_, coords );

	Stub const coords_stub_after_fitting( Stub( coords[1], coords[2], coords[3] ) );

	// figure out how the stubs have changed...
	Matrix const R( coords_stub_after_fitting.M * coords_stub_before_fitting.M.transposed() );
	Vector const translation( coords_stub_after_fitting.v - R * coords_stub_before_fitting.v );

	start_stub.M = R * start_stub_before_fitting.M;
	start_stub.v = R * start_stub_before_fitting.v + translation;

	stop_stub.M = R * final_stub_before_fitting.M;
	stop_stub.v = R * final_stub_before_fitting.v + translation;

	{  // DEBUGGING
		Vector const test1b( start_stub.global2local( coords[3] ) ),
			test2b( stop_stub.global2local( coords[4] ) );
		runtime_assert( test1a.distance_squared( test1b )<1e-3 );
		runtime_assert( test2a.distance_squared( test2b )<1e-3 );
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct HelixParamsLine {
	HelixParams hparams;
	string tag;
	Size helixlen;
};
typedef vector1< HelixParamsLine > HelixParamsLines;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
build_helix_library_v2(
	Sizes const & helixlens,
	Size const max_frags,
	map< Size, StubFrags > & all_frags,
	Sizes const src_helixlens = Sizes()
)
{
	using numeric::conversions::radians;
	static bool init( false );

	//Size const helix_window( 7 );

	static HelixParamsLines all_helix_params_lines;

	if ( !init ) { //read info on helical params
		init = true;
		string const pdb_coords_file( option[ my_options::helix_params_file ] );
		ifstream data( pdb_coords_file.c_str() );
		runtime_assert( data.good() );
		string line;
		HelixParamsLine hpline;
		bool found_full_tilt_lines( false );
		while ( getline( data, line ) ) {
			strings const l( split_to_vector1( line ) );
			if ( l[1] == "final_params_helix_full_tilt" ) {
				found_full_tilt_lines = true;
				hpline.hparams.rise           =          float_of( l[ 7] );
				hpline.hparams.twist          = radians( float_of( l[ 9] ) );
				hpline.hparams.tilt           = radians( float_of( l[11] ) );
				hpline.hparams.tilt_direction = radians( float_of( l[13] ) );
				hpline.hparams.ca_distance    =          float_of( l[15] );
				hpline.tag = l.back();
				hpline.helixlen = int_of( l[3] );
				if ( src_helixlens.size() && !has_element( hpline.helixlen, src_helixlens ) ) continue;
				all_helix_params_lines.push_back( hpline );
			} else if ( l[1] == "final_params_helix_full_tilt_precession" ) {
				runtime_assert( !found_full_tilt_lines );
				hpline.hparams.rise            =          float_of( l[ 7] );
				hpline.hparams.twist           = radians( float_of( l[ 9] ) );
				hpline.hparams.tilt            = radians( float_of( l[11] ) );
				hpline.hparams.tilt_direction  = radians( float_of( l[13] ) );
				hpline.hparams.tilt_precession = radians( float_of( l[15] ) );
				hpline.hparams.ca_distance     =          float_of( l[17] );
				hpline.tag = l.back();
				hpline.helixlen = int_of( l[3] );
				if ( src_helixlens.size() && !has_element( hpline.helixlen, src_helixlens ) ) continue;
				all_helix_params_lines.push_back( hpline );
			}
		}
		data.close();
		numeric::random::random_permutation( all_helix_params_lines, numeric::random::rg() );
		TR_STUB_FRAG.Trace << "Read " << all_helix_params_lines.size() << " helix lines from file " << pdb_coords_file <<
			endl;
	} // initialize the coordinates /////////////////////////////////////////////////////////////////////

	all_frags.clear();
	Vectors coords;
	for ( Sizes::const_iterator hlen= helixlens.begin(); hlen != helixlens.end(); ++hlen ) {
		Size const helixlen( *hlen );
		if ( all_frags.count(helixlen) ) continue;
		coords.reserve( helixlen );
		all_frags[ helixlen ];
		StubFrags & frags( all_frags.find( helixlen )->second );
		Reals all_strains;
		for ( Size ii=1; ii<= all_helix_params_lines.size(); ++ii ) {
			HelixParamsLine const & hpline( all_helix_params_lines[ii] );
			// if user passes in src_helixlens then we don't worry about this check:
			if ( src_helixlens.empty() && ( hpline.helixlen > Real(helixlen)*1.25 ||
					hpline.helixlen < Real(helixlen)*0.8 ) ) continue;
			// reconstruct some coords
			StubFrag f;
			create_helix_frag_from_helix_params( helixlen, hpline.tag, hpline.hparams, f );
			frags.push_back( f );
			if ( frags.size() >= max_frags ) break;
		} // loop over pdb helices

		TR_STUB_FRAG.Trace << "build_helix_library_v2:: Created " << frags.size() << " helix frags of length " <<
			helixlen << endl;
	}
}

struct TurnCoordsLine {
	string tag;
	Vectors coords;
	vector1< Reals > torsions;
	string sequence;
};

typedef StubFrag HelixStubFrag;
typedef utility::vector1< HelixStubFrag > HelixStubFrags;



//////////////////////////////////////////////////////////////////////////////////////////////////

class HelixBarrelMultifunc : public optimization::Multifunc {
public:
	typedef optimization::Multivec Multivec;

	HelixBarrelMultifunc():
		nrepeat_(0),
		params_scale_(10)
	{}

	virtual
	Real
	operator()( Multivec const & wts ) const;


	virtual
	void
	dfunc( Multivec const & phipsi, Multivec & dE_dphipsi ) const
	{
		dE_dphipsi.clear(); dE_dphipsi.resize( phipsi.size(), 0.0 );
	}

	void
	set_nrepeat_and_frags_and_helix_lengths(
		Size const nrepeat,
		StubFrag const & h1,
		StubFrag const & h2,
		StubFrag const & t1,
		StubFrag const & t2,
		Size const h1len,
		Size const h2len
	)
	{
		nrepeat_ = nrepeat;
		t1_ = t1;
		t2_ = t2;
		h1len_ = h1len;
		h2len_ = h2len;
		params_from_helix_frags( h1, h2, target_params_ );
	}

	void
	params_from_helix_frags( StubFrag const & h1,  StubFrag const & h2, Multivec & params ) const;

	void
	params_from_helix_params( HelixParams const & h1params, HelixParams const & h2params, Multivec & params ) const;

	void
	helix_frags_from_params( Multivec const & params, StubFrag & h1,  StubFrag & h2 ) const;

private:

	void
	unpack_params( Multivec const & params ) const;

private:
	Size nrepeat_;
	StubFrag t1_, t2_;
	Size h1len_, h2len_;
	Multivec target_params_;
	mutable HelixParams h1params_, h2params_;
	Real params_scale_;

};




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
HelixBarrelMultifunc::unpack_params( Multivec const & params ) const
{
	h1params_.rise           = params[1]/params_scale_;
	h1params_.twist          = params[2]/params_scale_;
	h1params_.tilt           = params[3]/params_scale_;
	h1params_.tilt_direction = params[4]/params_scale_;
	h1params_.ca_distance    = params[5]/params_scale_;

	h2params_.rise           = params[6]/params_scale_;
	h2params_.twist          = params[7]/params_scale_;
	h2params_.tilt           = params[8]/params_scale_;
	h2params_.tilt_direction = params[9]/params_scale_;
	h2params_.ca_distance    = params[10]/params_scale_;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
HelixBarrelMultifunc::operator()( Multivec const & params ) const
{
	runtime_assert( nrepeat_ );

	// unpack the multivec params into the helical params arrays h1params_ and h2params_ (they are mutable)
	unpack_params( params );

	///
	Stub const start_stub;

	RT const h1_rt( start_stub, generate_helix_stub( start_stub, h1params_, h1len_ ) );
	RT const h2_rt( start_stub, generate_helix_stub( start_stub, h2params_, h2len_ ) );

	// TR.Trace << "h1_rt_trans: " <<
	//  F(9,3,h1_rt.get_translation().x()) <<
	//  F(9,3,h1_rt.get_translation().y()) <<
	//  F(9,3,h1_rt.get_translation().z()) << endl;

	// TR.Trace << "h2_rt_trans: " <<
	//  F(9,3,h2_rt.get_translation().x()) <<
	//  F(9,3,h2_rt.get_translation().y()) <<
	//  F(9,3,h2_rt.get_translation().z()) << endl;

	Matrix const
		&R1( h1_rt.get_rotation() ), &R2( t1_.rt.get_rotation() ),
		&R3( h2_rt.get_rotation() ), &R4( t2_.rt.get_rotation() );

	Vector const
		&t1( h1_rt.get_translation() ), &t2( t1_.rt.get_translation() ),
		&t3( h2_rt.get_translation() ), &t4( t2_.rt.get_translation() );

	/// what happens to start_stub.v as we go through a single repeat?
	Matrix Rtot( R1 );
	Vector ttot( t1 );

	ttot += Rtot * t2; Rtot.right_multiply_by( R2 );
	ttot += Rtot * t3; Rtot.right_multiply_by( R3 );
	ttot += Rtot * t4; Rtot.right_multiply_by( R4 );


	Real theta, theta_dev, rise;
	{
		Vector center, n, t;
		get_stub_transform_data( start_stub, Stub( Rtot, ttot ), center, n, t, theta );
		theta = numeric::conversions::degrees( theta );
		theta_dev = numeric::square( basic::subtract_degree_angles( theta, 360.0 / nrepeat_ ) );
		//theta_dev = fabs( basic::subtract_degree_angles( theta, 360.0 / nrepeat_ ) );
		rise = n.dot(t);
	}

	Vector vfinal( start_stub.v ); // is 0.0
	for ( Size n=1; n<= nrepeat_; ++n ) {
		vfinal += ttot;
		ttot = Rtot * ttot;
		//vfinal = Rtot * vfinal + ttot;
	}


	if ( false ) { // DEBUGGING
		Stub stub( start_stub );
		for ( Size n=1; n<= nrepeat_; ++n ) {
			stub = t2_.rt.make_jump( h2_rt.make_jump( t1_.rt.make_jump( h1_rt.make_jump( stub ) ) ) );
		}
		Vector recompute_vfinal( stub.v );
		TR_STUB_FRAG.Trace << "recompute_vfinal_distance: " << F(9,3,vfinal.distance( recompute_vfinal ) ) << endl;
		vfinal = recompute_vfinal;
	}

	// constraint term
	Real cstE( 0.0 );
	Real const cstwt( 1000.0 );
	static Reals const cst_weights( make_vector1(
			1.0, // rise              Angstroms
			1.0, // twist             radians
			1.0, // tilt              radians
			0.1, // tilt_direction    radians
			1.0, // ca_distance       Angstroms
			1.0, // rise
			1.0, // twist
			1.0, // tilt
			0.1, // tilt_direction
			1.0 // ca_distance
		));
	for ( Size i=1; i<= 10; ++i ) {
		cstE += cst_weights[i] * numeric::square( ( params[i] - target_params_[i] )/params_scale_ );
	}

	Real const vfinal_dis2( vfinal.distance_squared( start_stub.v ) );

	//Real const func( vfinal_dis2 + cstwt * cstE + theta_dev );
	Real const func( cstwt * cstE + theta_dev + numeric::square( nrepeat_ * rise ) );

	TR_STUB_FRAG.Trace << "barrel_func " << F(9,3,func) <<
		" vfdis: " << F(9,3,sqrt( vfinal_dis2 ) ) <<
		" theta: " << F(9,3,theta ) <<
		" rise: " << F(9,3,rise ) <<
		" wtdcstE: " << F(9,3,cstwt*cstE ) <<
		" theta_dev: " << F(9,3,theta_dev ) << endl;
	TR_STUB_FRAG.Trace << "barrel_params " << F(9,3,func);
	for ( Size i=1; i<= 10; ++i ) TR_STUB_FRAG.Trace << ' ' << F(6,3,params[i]/params_scale_);
	TR_STUB_FRAG.Trace << endl;

	return func;

}

////////////////////////////////////////////
void
HelixBarrelMultifunc::params_from_helix_frags(
	StubFrag const & h1,
	StubFrag const & h2,
	Multivec & params
) const
{
	params_from_helix_params( h1.hparams, h2.hparams, params );
}

////////////////////////////////////////////
////////////////////////////////////////////
void
HelixBarrelMultifunc::params_from_helix_params(
	HelixParams const & h1params,
	HelixParams const & h2params,
	Multivec & params
) const
{
	params.resize( 10 );
	params[1] = params_scale_ * h1params.rise;
	params[2] = params_scale_ * h1params.twist;
	params[3] = params_scale_ * h1params.tilt;
	params[4] = params_scale_ * h1params.tilt_direction;
	params[5] = params_scale_ * h1params.ca_distance;

	params[6] = params_scale_ * h2params.rise;
	params[7] = params_scale_ * h2params.twist;
	params[8] = params_scale_ * h2params.tilt;
	params[9] = params_scale_ * h2params.tilt_direction;
	params[10]= params_scale_ * h2params.ca_distance;

}

////////////////////////////////////////////
// note that we don't set StubFrag::segtype here
void
HelixBarrelMultifunc::helix_frags_from_params(
	Multivec const & params,
	StubFrag & h1,
	StubFrag & h2
) const
{
	unpack_params( params );

	Stub const start_stub; // default
	Stub const h1_stop_stub( generate_helix_coords( start_stub, h1params_, h1len_, h1.coords ) );
	h1.rt = RT( start_stub, h1_stop_stub );
	h1.hparams = h1params_; // set by unpack_params

	Stub const h2_stop_stub( generate_helix_coords( start_stub, h2params_, h2len_, h2.coords ) );
	h2.rt = RT( start_stub, h2_stop_stub );
	h2.hparams = h2params_; // set by unpack_params

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



#endif
