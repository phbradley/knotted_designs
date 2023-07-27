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
#ifndef INCLUDED_apps_pilot_phil_sym_centroid_HH
#define INCLUDED_apps_pilot_phil_sym_centroid_HH

#include <apps/pilot/phil/phil.hh>
#include <apps/pilot/phil/types.hh>

#include <apps/pilot/phil/sym_basic.hh>
#include <apps/pilot/phil/sym_setup.hh> // setup_unbound_frag_symminfo
#include <apps/pilot/phil/symscores.hh>

#include <devel/blab/loops/util.hh> // switch_pose_to_residue_type_set
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <protocols/relax/FastRelax.hh>
#include <devel/blab/classic_frags/TorsionFragment.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

static basic::Tracer TR_SYM_CENTROID( "apps.pilot.phil.sym_centroid_hh" );





///////////////////////////////////////////////////////////////////////////////////////
/// insert protein fragment
void
insert_protein_fragment(
	Size const fragsize,
	Sizes const & fragseq_poslist,
	//             bool const insert_frag_sequence,
	//            bool const freeze_anchorseq,
	devel::blab::classic_frags::FragLib const & fraglib,
	Pose & pose
)
{
	//bool const freeze_helixseq( option[ my_options::freeze_helixseq ] );
	Size const nbb( 3 );

	// hacking
	if ( pose::symmetry::is_symmetric( pose ) && !pose::symmetry::symmetry_info( pose )->contiguous_monomers() ) {
		// special case to handle here
		conformation::symmetry::SymmetryInfo const & symminfo( *pose::symmetry::symmetry_info( pose ) );

		//Size const nres_protein( chain_end( 1, pose ) );

		devel::blab::classic_frags::TorsionFragmentLibrary const & lib( fraglib.library( fragsize ) );

		Size fragpos(0), nn(0);
		while ( true ) {
			fragpos = int( lib.size() * uniform() + 1 );
			if ( lib[ fragpos ].size() ) {
				nn = int( uniform() * lib[ fragpos ].size()) + 1;
				break;
			}
		}

		runtime_assert( lib[ fragpos ][ nn ].nbb() == nbb );

		for ( Size k=1; k<= fragsize; ++k ) {
			Size seqpos( fragpos + k-1 ), basepos( symminfo.bb_follows( seqpos ) );
			if ( !basepos ) { runtime_assert( symminfo.bb_is_independent( seqpos ) ); basepos = seqpos; }
			else { runtime_assert( !symminfo.bb_is_independent( seqpos ) ); }
			pose.set_secstruct( basepos, lib[ fragpos ][ nn ].get_secstruct( k ) );
			for ( Size j=1; j<= nbb; ++j ) {
				pose.set_torsion( id::TorsionID( basepos, id::BB, j ), lib[ fragpos ][ nn ].get_torsion( k, j ) );
			}
		}
		return; /////////////////////////////// EARLY RETURN !! /////////////////////////////////////////////////
	}

	Size nrepeat, repeatlen, base_repeat;
	parse_symminfo( pose, nrepeat, repeatlen, base_repeat );


	//Size const nres_protein( nrepeat*repeatlen );//chain_end( 1, pose ) );
	//runtime_assert( nres_protein%repeatlen == 0 );

	devel::blab::classic_frags::TorsionFragmentLibrary const & lib( fraglib.library( fragsize ) );

	Size fragpos(0), nn(0);
	while ( true ) {
		fragpos = int( lib.size() * uniform() + 1 );
		if ( lib[ fragpos ].size() ) {
			nn = int( uniform() * lib[ fragpos ].size()) + 1;
			break;
		}
	}

	runtime_assert( lib[ fragpos ][ nn ].nbb() == nbb );

	Sizes repeats_for_inserting;
	if ( pose::symmetry::is_symmetric( pose ) ) repeats_for_inserting = make_vector1( base_repeat );
	else for ( Size i=1; i<= nrepeat; ++i ) repeats_for_inserting.push_back(i);

	for ( Sizes::const_iterator it= repeats_for_inserting.begin(); it != repeats_for_inserting.end(); ++it ) {
		Size const repeat( *it );
		for ( Size k=1; k<= fragsize; ++k ) {
			Size const rpos( (fragpos + (k-1) - 1)%repeatlen+1 ), pos( (repeat-1)*repeatlen + rpos );
			pose.set_secstruct( pos, lib[ fragpos ][ nn ].get_secstruct( k ) );
			for ( Size j=1; j<= nbb; ++j ) {
				pose.set_torsion( id::TorsionID( pos, id::BB, j ), lib[ fragpos ][ nn ].get_torsion( k, j ) );
			}
		}
		if ( fragseq_poslist.size() ) {
			for ( Size k=1; k<= fragsize; ++k ) {
				Size const rpos( (fragpos + (k-1) - 1)%repeatlen+1 ), pos( (repeat-1)*repeatlen + rpos );
				if ( std::find( fragseq_poslist.begin(), fragseq_poslist.end(), pos ) == fragseq_poslist.end() ) continue;
				//if ( freeze_anchorseq && pose.fold_tree().is_jump_point( pos ) ) continue;
				//if ( freeze_helixseq && pose.secstruct(pos) == 'H' ) continue;
				AA const new_aa( aa_from_oneletter_code( lib[ fragpos ][ nn ].get_sequence( k ) ) );
				if ( new_aa != pose.residue( pos ).aa() ) make_sequence_change( pos, new_aa, pose );
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////
void
make_frag_move(
	Size const fragsize,
	string const & tag,
	devel::blab::classic_frags::FragLib const & fraglib,
	Sizes const & fragseq_poslist,
	Pose & pose,
	protocols::moves::MonteCarlo & mc
)
{
	insert_protein_fragment( fragsize, fragseq_poslist, fraglib, pose ); // freeze_anchorseq
	//Real const score( mc.score_function()( pose ) );
	bool const mc_accept( mc.boltzmann( pose, tag ) );
	TR_SYM_CENTROID.Trace << "MCFRAG " << mc_accept << ' ' << " mc_temp: " << F(9,3,mc.temperature() ) <<
		" twistrise_E: " << F(9,3,pose.energies().total_energies()[ twistrise ] ) <<
		" twistrise_W: " << F(9,3,mc.score_function().get_weight( twistrise ) ) <<
		" ap_E: " << F(9,3,pose.energies().total_energies()[ atom_pair_constraint ] ) <<
		" ap_W: " << F(9,3,mc.score_function().get_weight( atom_pair_constraint ) ) <<
		// " cb_E: " << F(9,3,pose.energies().total_energies()[ chainbreak ] ) <<
		// " cb_E: " << F(9,3,mc.score_function().get_weight( chainbreak ) ) <<
		" last: " << F(9,3,mc.last_accepted_score()) <<
		" low: " << F(9,3,mc.lowest_score() ) << ' ' << tag << endl;
}

///////////////////////////////////////////////////////////////////////////////
Stub
get_beta_pairing_stub(
	Size const pos1,
	Size const pos2,
	Pose const & pose
)
{
	Residue const & rsd1( pose.residue(pos1) ), &rsd2( pose.residue(pos2) );
	Stub const
		stub1( rsd1.xyz("N"),rsd1.xyz("CA"),rsd1.xyz("C") ),
		stub2( rsd2.xyz("N"),rsd2.xyz("CA"),rsd2.xyz("C") );

	Vector center, symmetry_axis, t;
	Real theta;
	get_stub_transform_data( stub1, stub2, center, symmetry_axis, t, theta );

	runtime_assert( fabs( theta - numeric::constants::d::pi )<1e-3 || fabs( theta + numeric::constants::d::pi ) < 1e-3 );
	runtime_assert( fabs( t.dot( symmetry_axis )) < 1e-1 );

	Vector orig( 0.5*( rsd1.xyz("CA")+rsd2.xyz("CA") ) ), x( ( rsd2.xyz("CA") - rsd1.xyz("CA") ).normalized() ),
		y( symmetry_axis.cross(x));
	if ( y.dot( rsd2.xyz("C") - rsd1.xyz("C") )<0 ) y *= -1;
	// { // hacking
	//  Stub stub( orig, orig+x, orig+y );
	//  TR_SYM_CENTROID.Trace << "ORIG_STUB: " <<
	//   stub.M.col_x().dot( x ) << ' ' <<
	//   stub.M.col_y().dot( y ) << ' ' <<
	//   stub.M.col_z().dot( x.cross(y) ) << endl;
	// }
	return Stub( orig, orig+x, orig+y );
}




///////////////////////////////////////////////////////////////////////////////
// assumes a setup where strand1 is paired with strand-nrepeat_per_braid from braid2
void
calculate_braided_barrel_r_theta(
	Pose const & pose,
	Real & r,
	Real & theta
)
{


	Size nrepeat, repeatlen, base_repeat;
	parse_symminfo( pose, nrepeat, repeatlen, base_repeat );

	Size const nbraid(2);
	runtime_assert( nbraid == num_chains( pose ) - 1 );
	Size const nrepeat_per_braid( nrepeat/2 ), nres_protein( chain_end( nbraid, pose ) );
	Size const base_jump( base_repeat );
	conformation::symmetry::SymmetryInfoCOP symminfo( pose::symmetry::symmetry_info( pose ) );
	runtime_assert( symminfo->jump_is_independent( base_jump ) );
	Size const vrtpos1( pose.fold_tree().upstream_jump_residue(base_jump)),
		protpos1( pose.fold_tree().downstream_jump_residue(base_jump) );
	Size vrtpos2(0), protpos2(0);
	{ // look for the virtual that is paired with this one
		FoldTree const & f( pose.fold_tree() );
		for ( Size i=1; i<= f.num_jump(); ++i ) {
			if ( f.upstream_jump_residue(i) == vrtpos1 && f.downstream_jump_residue(i)>nres_protein ) {
				runtime_assert( !vrtpos2 );
				vrtpos2 = f.downstream_jump_residue(i);
			}
		}
		for ( Size i=1; i<= f.num_jump(); ++i ) {
			if ( f.upstream_jump_residue(i) == vrtpos2 ) {
				runtime_assert( !protpos2 );
				protpos2 = f.downstream_jump_residue(i);
				runtime_assert( protpos2<= nres_protein );
				runtime_assert( protpos2 > nrepeat_per_braid*repeatlen );
			}
		}
		runtime_assert( vrtpos2 && protpos2 );
	}

	// get central symmetry axis
	Vector central_symmetry_axis, central_symmetry_axis_center;
	{ // now construct a stub for this pairing
		Residue const & rsd1( pose.residue(1) ), &rsd2( pose.residue( repeatlen+1 ) );
		Stub const
			stub1( rsd1.xyz("N"),rsd1.xyz("CA"),rsd1.xyz("C") ),
			stub2( rsd2.xyz("N"),rsd2.xyz("CA"),rsd2.xyz("C") );

		Vector t;
		Real theta;
		get_stub_transform_data( stub1, stub2, central_symmetry_axis_center, central_symmetry_axis, t, theta );
	}

	Stub const beta_stub( get_beta_pairing_stub( protpos1, protpos2, pose ) );
	Vector const orig( beta_stub.v ), y( beta_stub.M.col_y() );
	Vector v( orig - central_symmetry_axis_center );
	v -= v.dot( central_symmetry_axis ) * central_symmetry_axis;
	Vector const proj( orig - v );
	runtime_assert( ( proj - central_symmetry_axis_center ).cross( central_symmetry_axis ).length() < 1e-3 );
	r = v.length();
	theta = numeric::dihedral_degrees( proj + central_symmetry_axis, proj, orig, orig + y );
}



///////////////////////////////////////////////////////////////////////////////
void
make_braided_barrel_r_theta_move( Real const r_sdev, Real const theta_sdev, Pose & pose )
{
	runtime_assert( false );
	Size nrepeat, repeatlen, base_repeat;
	parse_symminfo( pose, nrepeat, repeatlen, base_repeat );

	Size const base_jump( nrepeat + base_repeat );
	conformation::symmetry::SymmetryInfoCOP symminfo( pose::symmetry::symmetry_info( pose ) );

	runtime_assert( symminfo->jump_is_independent( base_jump ) );

	kinematics::Jump jump( pose.jump( base_jump ) );
	Real const delta_r    (     r_sdev * numeric::random::gaussian() );
	Real const delta_theta( theta_sdev * numeric::random::gaussian() );
	jump.set_rb_delta( 1, 1, delta_r  );
	jump.set_rb_delta( 4, 1, delta_theta  );
	jump.fold_in_rb_deltas();
	pose.set_jump( base_jump, jump );

}



///////////////////////////////////////////////////////////////////////////////
void
make_braided_barrel_beta_move(
	vector1< RTs > const & beta_transforms,
	Pose & pose
)
{
	static vector1< vector1< std::pair< Real, Real > > > r_thetas;
	if ( r_thetas.empty() ) {
		r_thetas.resize(2);
		string const r_thetas_file( option[ my_options::r_thetas_file ]() );
		ifstream data( r_thetas_file.c_str());
		runtime_assert( data.good() );
		string line, tag1, tag2, tag3;
		while ( getline( data, line ) ) {
			istringstream l( line );
			bool beta_bridged;
			Real r,theta;
			l >> tag1 >> beta_bridged >> tag2 >> r >> tag3 >> theta;
			if ( !l.fail() && tag1 == "beta_bridged" && tag2 == "r" && tag3 == "theta" ) {
				r_thetas[ 1+beta_bridged ].push_back( make_pair( r,theta ) );
			}
		}
		TR_SYM_CENTROID.Trace << "read " << r_thetas[1].size() << ' ' << r_thetas[2].size() << " from " <<
			option[ my_options::r_thetas_file ]() << endl;
		runtime_assert( !r_thetas[1].empty() && !r_thetas[2].empty() );
		data.close();
	}


	runtime_assert( !beta_transforms[ symmetrized_antiparallel_type ].empty() );
	runtime_assert( !beta_transforms[ symmetrized_antiparallel2_type ].empty() );
	using numeric::constants::d::pi;
	static Pose beta_pose;
	if ( beta_pose.total_residue() != 2 ) {
		//initialize the mini pose
		beta_pose.append_residue_by_bond( *get_vanilla_protein_residue('A') );
		beta_pose.append_residue_by_jump( *get_vanilla_protein_residue('A'), 1 );
	}

	// choose a random symmetrized antiparallel transform
	bool const beta_bridged( uniform()<0.5 );
	RT rt;
	{
		if ( beta_bridged ) rt = random_element( beta_transforms[ symmetrized_antiparallel_type ] );
		else                rt = random_element( beta_transforms[ symmetrized_antiparallel2_type ] );
		Residue const & rsd1( beta_pose.residue(1) ), &rsd2( beta_pose.residue(2) );
		id::StubID const stubid1( id::AtomID( rsd1.atom_index("N" ), 1 ),
			id::AtomID( rsd1.atom_index("CA"), 1 ),
			id::AtomID( rsd1.atom_index("C" ), 1 ) );
		id::StubID const stubid2( id::AtomID( rsd2.atom_index("N" ), 2 ),
			id::AtomID( rsd2.atom_index("CA"), 2 ),
			id::AtomID( rsd2.atom_index("C" ), 2 ) );
		beta_pose.conformation().set_stub_transform( stubid1, stubid2, rt );
	} // residue const & makes me nervous with conformational changes

	// now figure out what the symmetry axis is
	// Residue const & rsd1( beta_pose.residue(1) ), &rsd2( beta_pose.residue(2) );
	// Stub const
	//  stub1( rsd1.xyz("N"),rsd1.xyz("CA"),rsd1.xyz("C") ),
	//  stub2( rsd2.xyz("N"),rsd2.xyz("CA"),rsd2.xyz("C") );

	// Vector center, symmetry_axis, t;
	// Real theta;
	// get_stub_transform_data( stub1, stub2, center, symmetry_axis, t, theta );

	// runtime_assert( fabs( theta - pi )<1e-3 || fabs( theta + pi ) < 1e-3 );
	// runtime_assert( fabs( t.dot( symmetry_axis )) < 1e-3 );

	// Vector orig( 0.5*( rsd1.xyz("CA")+rsd2.xyz("CA") ) ), x( ( rsd2.xyz("CA") - rsd1.xyz("CA") ).normalized() ),
	//  y( symmetry_axis.cross(x));
	// if ( y.dot( rsd2.xyz("C") - rsd1.xyz("C") )<0 ) y *= -1;
	// vrt->set_xyz("ORIG", orig);
	// vrt->set_xyz("X", orig+x);
	// vrt->set_xyz("Y", orig+y);
	Stub const pairing_stub( get_beta_pairing_stub( 1, 2, beta_pose ) );
	RT const desired_rt( pairing_stub, Stub( beta_pose.residue(1).xyz("N"),
		beta_pose.residue(1).xyz("CA"),
		beta_pose.residue(1).xyz("C") ) );
	// which jump are we talking about?

	Size nrepeat, repeatlen, base_repeat;
	parse_symminfo( pose, nrepeat, repeatlen, base_repeat );

	Size const base_jump( base_repeat ), nres_protein( repeatlen * nrepeat ), nrepeat_per_braid( nrepeat/2 );
	conformation::symmetry::SymmetryInfoCOP symminfo( pose::symmetry::symmetry_info( pose ) );
	runtime_assert( symminfo->jump_is_independent( base_jump ) );
	Size const vrtpos( pose.fold_tree().upstream_jump_residue(base_jump)),
		protpos( pose.fold_tree().downstream_jump_residue(base_jump) );
	{ // set the jump
		Residue const & rsd1( pose.residue(vrtpos) ), &rsd2( pose.residue(protpos));
		id::StubID const stubid1( id::AtomID( rsd1.atom_index("ORIG" ), vrtpos ),
			id::AtomID( rsd1.atom_index("X"    ), vrtpos ),
			id::AtomID( rsd1.atom_index("Y"    ), vrtpos ) );
		id::StubID const stubid2( id::AtomID( rsd2.atom_index("N"    ), protpos ),
			id::AtomID( rsd2.atom_index("CA"   ), protpos ),
			id::AtomID( rsd2.atom_index("C"    ), protpos ) );
		pose.conformation().set_stub_transform( stubid1, stubid2, desired_rt );
		pose.set_jump( base_jump, pose.jump( base_jump ) );
	}

	{ // did it work?
		Size vrtpos2(0), protpos2(0);
		{ // look for the virtual that is paired with this one
			FoldTree const & f( pose.fold_tree() );
			for ( Size i=1; i<= f.num_jump(); ++i ) {
				if ( f.upstream_jump_residue(i) == vrtpos && f.downstream_jump_residue(i)>nres_protein ) {
					runtime_assert( !vrtpos2 );
					vrtpos2 = f.downstream_jump_residue(i);
				}
			}
			for ( Size i=1; i<= f.num_jump(); ++i ) {
				if ( f.upstream_jump_residue(i) == vrtpos2 ) {
					runtime_assert( !protpos2 );
					protpos2 = f.downstream_jump_residue(i);
					runtime_assert( protpos2<= nres_protein );
					runtime_assert( protpos2 > nrepeat_per_braid*repeatlen );
				}
			}
			runtime_assert( vrtpos2 && protpos2 );
		}
		Residue const & rsd1( pose.residue(protpos) ), &rsd2( pose.residue(protpos2) );
		RT const actual_rt( Stub( rsd1.xyz("N"), rsd1.xyz("CA"), rsd1.xyz("C") ),
			Stub( rsd2.xyz("N"), rsd2.xyz("CA"), rsd2.xyz("C") ) );
		Real const dev( rt.distance_squared( actual_rt ) );
		TR_SYM_CENTROID.Trace << "rt-dev " << dev << endl;
		runtime_assert( dev<1e-2 );
	}

	Real current_r, current_theta;
	calculate_braided_barrel_r_theta( pose, current_r, current_theta );

	std::pair< Real, Real > const desired_rtheta( random_element( r_thetas[ 1+beta_bridged ] ) );

	// try to modify these guys
	{
		Size const base_jump_from_center( nrepeat + base_repeat );

		runtime_assert( symminfo->jump_is_independent( base_jump_from_center ) );

		kinematics::Jump jump( pose.jump( base_jump_from_center ) );
		Real const delta_r    ( current_r - desired_rtheta.first );
		Real const delta_theta( basic::subtract_degree_angles( current_theta, desired_rtheta.second ) );
		jump.set_rb_delta( 1, 1, delta_r  );
		jump.set_rb_delta( 4, 1, delta_theta  );
		jump.fold_in_rb_deltas();
		pose.set_jump( base_jump_from_center, jump );

		Real new_r, new_theta;
		calculate_braided_barrel_r_theta( pose, new_r, new_theta );
		TR_SYM_CENTROID.Trace << "WORKED?? " <<
			current_r << ' ' << desired_rtheta.first << ' ' << new_r << ' ' <<
			current_theta << ' ' << desired_rtheta.second << ' ' << new_theta << ' ' <<
			basic::subtract_degree_angles( new_theta, desired_rtheta.second ) << ' ' <<
			endl;
		runtime_assert( fabs( basic::subtract_degree_angles( new_theta, desired_rtheta.second ) )<1e-3 );
	}
}


///////////////////////////////////////////////////////////////////////////////////////
void
get_antiparallel_braided_params(
	Pose const & pose,
	Real & d,
	Real & alpha,
	Real & rise
)
{
	Size const nbraid(2);
	Size nrepeat, repeatlen, base_repeat;
	parse_symminfo( pose, nrepeat, repeatlen, base_repeat );

	runtime_assert( nrepeat%nbraid==0 );
	Size const nrepeat_per_braid( nrepeat / nbraid );

	Size const //root( repeatlen * nrepeat_per_braid ), anchor( 2 * repeatlen * nrepeat_per_braid ),
		nres_braid( repeatlen * nrepeat_per_braid );
	//runtime_assert( pose.fold_tree().jump_exists( root, anchor ) );

	// figure out the symmetry axis
	// convention for this function is that rise should be positive
	// default convention from get_stub_transform_data is that theta is positive, which could lead to instability?
	Vector symmetry_axis,center;
	Real theta;
	{
		Vector t;
		Stub const
			istub( pose.residue(1).xyz("CA"),
			pose.residue(2).xyz("CA"),
			pose.residue(3).xyz("CA") ),
			jstub( pose.residue(repeatlen+1).xyz("CA"),
			pose.residue(repeatlen+2).xyz("CA"),
			pose.residue(repeatlen+3).xyz("CA") );
		get_stub_transform_data( istub, jstub, center, symmetry_axis, t, theta );
		rise = symmetry_axis.dot(t);
		if ( rise<0 ) {
			symmetry_axis *= -1;
			rise *= -1;
			theta *= -1;
		}
	}
	Vector symmetry_axis2,center2;
	{
		Real theta2, rise2;
		Vector t;
		Stub const
			istub( pose.residue(nres_braid+1).xyz("CA"),
			pose.residue(nres_braid+2).xyz("CA"),
			pose.residue(nres_braid+3).xyz("CA") ),
			jstub( pose.residue(nres_braid+repeatlen+1).xyz("CA"),
			pose.residue(nres_braid+repeatlen+2).xyz("CA"),
			pose.residue(nres_braid+repeatlen+3).xyz("CA") );
		get_stub_transform_data( istub, jstub, center, symmetry_axis2, t, theta2 );
		rise2 = symmetry_axis2.dot(t);
		if ( rise2<0 ) {
			symmetry_axis2 *= -1;
			rise2 *= -1;
			theta2 *= -1;
		}
		Real const dev( fabs( basic::subtract_radian_angles( theta, theta2 ) )+fabs( rise-rise2) );
		//TR_SYM_CENTROID.Trace << "symparams dev1: " << dev << endl;
		runtime_assert( dev<1e-1 );
	}

	if ( symmetry_axis2.distance( -1*symmetry_axis ) > 1e-1 ) {
		utility_exit_with_message("get_antiparallel_braided_params: pose is not braided: "+
			string_of( symmetry_axis2.distance( symmetry_axis ) ) );
	}

	// where are the centroids of the base repeats?
	Vectors centroids, projections;
	for ( Size r=0; r<2; ++r ) {
		Size const offset( nres_braid*r + (base_repeat-1)*repeatlen );
		Vector centroid(0,0,0);
		for ( Size i=1; i<= repeatlen; ++i ) centroid += pose.residue(offset+i).xyz("CA");
		centroid /= repeatlen;
		centroids.push_back( centroid );
		projections.push_back( center + symmetry_axis * symmetry_axis.dot( centroid - center ) ); // lies along the axis
	}

	// what are the current params?
	d = ( projections[2] - projections[1] ).dot( symmetry_axis );
	alpha = numeric::dihedral_radians( centroids[1], projections[1], projections[1] + symmetry_axis,
		projections[1] + symmetry_axis + ( centroids[2]-projections[2]));
}

///////////////////////////////////////////////////////////////////////////////////////
void
set_antiparallel_braided_pose(
	Real const d,
	Real const alpha,
	Pose & pose
)
{
	Size const nbraid(2);
	Size nrepeat, repeatlen, base_repeat;
	parse_symminfo( pose, nrepeat, repeatlen, base_repeat );

	runtime_assert( nrepeat%nbraid==0 );
	Size const nrepeat_per_braid( nrepeat / nbraid ), nres_braid( repeatlen * nrepeat_per_braid );
	runtime_assert( nbraid * nres_braid + 1 == pose.total_residue() );
	runtime_assert( pose.residue( pose.total_residue() ).name() == "VRT" );
	runtime_assert( num_chains( pose ) == nbraid+ 1 );

	// figure out the symmetry axis
	// convention for this function is that rise should be positive
	// default convention from get_stub_transform_data is that theta is positive, which could lead to instability?
	Vector symmetry_axis,center;
	Real theta, rise;
	{
		Vector t;
		Stub const
			istub( pose.residue(1).xyz("CA"),
			pose.residue(2).xyz("CA"),
			pose.residue(3).xyz("CA") ),
			jstub( pose.residue(repeatlen+1).xyz("CA"),
			pose.residue(repeatlen+2).xyz("CA"),
			pose.residue(repeatlen+3).xyz("CA") );
		get_stub_transform_data( istub, jstub, center, symmetry_axis, t, theta );
		rise = symmetry_axis.dot(t);
		if ( rise<0 ) {
			symmetry_axis *= -1;
			rise *= -1;
			theta *= -1;
		}
	}
	Vector symmetry_axis2,center2;
	{
		Real theta2, rise2;
		Vector t;
		Stub const
			istub( pose.residue(nres_braid+1).xyz("CA"),
			pose.residue(nres_braid+2).xyz("CA"),
			pose.residue(nres_braid+3).xyz("CA") ),
			jstub( pose.residue(nres_braid+repeatlen+1).xyz("CA"),
			pose.residue(nres_braid+repeatlen+2).xyz("CA"),
			pose.residue(nres_braid+repeatlen+3).xyz("CA") );
		get_stub_transform_data( istub, jstub, center2, symmetry_axis2, t, theta2 );
		rise2 = symmetry_axis2.dot(t);
		if ( rise2<0 ) {
			symmetry_axis2 *= -1;
			rise2 *= -1;
			theta2 *= -1;
		}
		Real const dev( fabs( basic::subtract_radian_angles( theta, theta2 ) )+fabs( rise-rise2) );
		//TR_SYM_CENTROID.Trace << "symparams dev2: " << dev << endl;
		runtime_assert( dev<1e-1 );
	}

	/// for moving the second braid we'll use set_stub_transform, which is not symmetrized yet
	Size const root( pose.total_residue() ), anchor1( nres_braid ), anchor2( 2 * nres_braid );
	runtime_assert( pose.fold_tree().jump_exists( root, anchor1 ) );
	runtime_assert( pose.fold_tree().jump_exists( root, anchor2 ) );
	runtime_assert( pose.fold_tree().root() == root ); // want to know that this rsd will stay fixed when we set jumps
	id::StubID const rootstubid( id::AtomID( pose.residue(root).atom_index("ORIG"), root ),
		id::AtomID( pose.residue(root).atom_index("X"   ), root ),
		id::AtomID( pose.residue(root).atom_index("Y"   ), root ) );

	id::StubID const anchorstubid( id::AtomID( pose.residue(anchor2).atom_index("N" ), anchor2 ),
		id::AtomID( pose.residue(anchor2).atom_index("CA"), anchor2 ),
		id::AtomID( pose.residue(anchor2).atom_index("C" ), anchor2 ) );
	// this must stay fixed while we fiddle with the pose:
	Stub const rootstub( pose.residue(root).xyz("ORIG"), pose.residue(root).xyz("X"), pose.residue(root).xyz("Y") );
	//Stub const anchor1stub(pose.residue(anchor1).xyz("N"),pose.residue(anchor1).xyz("CA"),pose.residue(anchor1).xyz("C"));

	// first, let's move the second braid so the symmetry axes are antiparallel (and colinear)
	// and center2 is transformed to center
	//if ( symmetry_axis2.distance( -1*symmetry_axis ) > 1e-6 ) {
	{
		Real const rotangle( std::acos( numeric::sin_cos_range( symmetry_axis2.dot( -1*symmetry_axis ) ) ) );
		Vector rotaxis( symmetry_axis2.cross( -1*symmetry_axis ) );
		rotaxis = rotaxis.normalized_any();
		numeric::xyzMatrix< Real > R( numeric::rotation_matrix( rotaxis, rotangle ) );
		runtime_assert( is_small( symmetry_axis.distance_squared( -1*R*symmetry_axis2 ) ) );
		Vector const v( center - R*center2 );
		// set the jump
		RT const desired_rt( rootstub, Stub( R * pose.residue(anchor2).xyz("N" ) + v,
			R * pose.residue(anchor2).xyz("CA") + v,
			R * pose.residue(anchor2).xyz("C" ) + v ) );
		// this routine is not symmetrized
		pose.conformation().set_stub_transform( rootstubid, anchorstubid, desired_rt );
	}

	{ // did it work?
		Vector symmetry_axis3,center3;
		Real theta3, rise3;
		Vector t;
		Stub const
			istub( pose.residue(nres_braid+1).xyz("CA"),
			pose.residue(nres_braid+2).xyz("CA"),
			pose.residue(nres_braid+3).xyz("CA") ),
			jstub( pose.residue(nres_braid+repeatlen+1).xyz("CA"),
			pose.residue(nres_braid+repeatlen+2).xyz("CA"),
			pose.residue(nres_braid+repeatlen+3).xyz("CA") );
		get_stub_transform_data( istub, jstub, center3, symmetry_axis3, t, theta3 );
		rise3 = symmetry_axis3.dot(t);
		if ( rise3<0 ) {
			symmetry_axis3 *= -1;
			rise3 *= -1;
			theta3 *= -1;
		}
		Real dev( fabs( basic::subtract_radian_angles( theta, theta3 ) )+fabs( rise-rise3) );
		//TR_SYM_CENTROID.Trace << "symparams dev3: " << dev << endl;
		runtime_assert( dev<1e-1 );
		dev = symmetry_axis3.distance( -1 * symmetry_axis );
		//TR_SYM_CENTROID.Trace << "axis dev3: " << dev << endl;
		runtime_assert( dev<1e-1 );
		dev = (center3-center).cross( symmetry_axis ).length();
		//TR_SYM_CENTROID.Trace << "center dev3: " << dev << endl;
		runtime_assert( dev<1e-1 );
	}


	// where are the centroids of the base repeats?
	Vectors centroids, projections;
	for ( Size r=0; r<2; ++r ) {
		Size const offset( nres_braid*r + (base_repeat-1)*repeatlen );
		Vector centroid(0,0,0);
		for ( Size i=1; i<= repeatlen; ++i ) centroid += pose.residue(offset+i).xyz("CA");
		centroid /= repeatlen;
		centroids.push_back( centroid );
		projections.push_back( center + symmetry_axis * symmetry_axis.dot( centroid - center ) ); // lies along the axis
	}

	// what are the current params?
	Real const current_d( ( projections[2] - projections[1] ).dot( symmetry_axis ) );
	Real const current_alpha( numeric::dihedral_radians( centroids[1], projections[1], projections[1] + symmetry_axis,
		projections[1] + symmetry_axis + ( centroids[2]-projections[2])));

	// so, we want to rotate braid two about the symmetry axis by (alpha - current_alpha)
	// and translate by (d-current_d)*symmetry_axis

	{
		numeric::xyzMatrix< Real > R( numeric::rotation_matrix( symmetry_axis, alpha - current_alpha ) );
		// want center to go to center + (d-current_d)*symmetry_axis
		Vector const v( center + (d-current_d)*symmetry_axis - R*center );
		// set the jump
		RT const desired_rt( rootstub, Stub( R * pose.residue(anchor2).xyz("N" ) + v,
			R * pose.residue(anchor2).xyz("CA") + v,
			R * pose.residue(anchor2).xyz("C" ) + v ) );
		// this routine is not symmetrized
		pose.conformation().set_stub_transform( rootstubid, anchorstubid, desired_rt );
	}

	// did it work ??
	{
		Size const offset( nres_braid + (base_repeat-1)*repeatlen );
		Vector new_centroid2(0,0,0);
		for ( Size i=1; i<= repeatlen; ++i ) new_centroid2 += pose.residue(offset+i).xyz("CA");
		new_centroid2 /= repeatlen;
		Vector const new_projection2( center + symmetry_axis * symmetry_axis.dot( new_centroid2 - center ) );

		Real const new_d( ( new_projection2 - projections[1] ).dot( symmetry_axis ) );
		Real const new_alpha( numeric::dihedral_radians( centroids[1], projections[1], projections[1] + symmetry_axis,
			projections[1] + symmetry_axis+ ( new_centroid2- new_projection2)));
		Real const d_dev( fabs( d - new_d ) ), alpha_dev( fabs( basic::subtract_radian_angles( alpha, new_alpha ) ) );
		//TR_SYM_CENTROID.Trace << "d_dev: " << d_dev << " alpha_dev: " << alpha_dev << endl;
		runtime_assert( d_dev+alpha_dev < 1e-1 );

		if ( false ) { // another test, slow
			Real d2, alpha2, rise3;
			get_antiparallel_braided_params( pose, d2, alpha2, rise3 );
			Real const d_dev( fabs( d - d2 ) ), alpha_dev( fabs( basic::subtract_radian_angles( alpha, alpha2 ) ) );
			//TR_SYM_CENTROID.Trace << "d_dev2: " << d_dev << " alpha_dev2: " << alpha_dev << endl;
			runtime_assert( d_dev+alpha_dev < 1e-1 );
		}
	}
}



///////////////////////////////////////////////////////////////////////////////////////
//
// after a fragment insertion, say, restore a braided configuration
//
void
set_parallel_braided_pose(
	Size const nbraid,
	Pose & pose
)
{
	Size nrepeat, repeatlen, base_repeat;
	parse_symminfo( pose, nrepeat, repeatlen, base_repeat );

	runtime_assert( nrepeat%nbraid == 0 );
	Size const nrepeat_per_braid( nrepeat/nbraid );
	Size const nres_braid( nrepeat_per_braid*repeatlen );

	// confirm fold tree is set up properly

	// figure out the symmetry axis
	Vector center, symmetry_axis, t;
	Real theta;
	{
		Stub const
			istub( pose.residue(1).xyz("CA"),
			pose.residue(2).xyz("CA"),
			pose.residue(3).xyz("CA") ),
			jstub( pose.residue(repeatlen+1).xyz("CA"),
			pose.residue(repeatlen+2).xyz("CA"),
			pose.residue(repeatlen+3).xyz("CA") );
		get_stub_transform_data( istub, jstub, center, symmetry_axis, t, theta );
	}

	// now we can figure out where the anchor backbones should go
	// rotation of (2pi/nbraid) about the symmetry axis
	numeric::xyzMatrix< Real > const R
		( numeric::rotation_matrix( symmetry_axis, ( 2 * numeric::constants::d::pi )/nbraid ) );
	Vector const v( center - R*center ); // so that x --> (Rx + v) preserves center

	Size const root( pose.total_residue() );
	Size braid_anchor(0);
	{
		FoldTree const & f( pose.fold_tree() );
		for ( Size i=1; i<= f.num_jump(); ++i ) {
			if ( f.upstream_jump_residue(i) == root ) {
				if ( braid_anchor ) {
					runtime_assert( ( f.downstream_jump_residue(i)-1 )%nres_braid + 1 == braid_anchor );
				} else {
					braid_anchor = ( f.downstream_jump_residue(i)-1 )%nres_braid + 1;
				}
			}
		}
		runtime_assert( braid_anchor );
	}
	runtime_assert( pose.residue( root ).name() == "VRT" );
	runtime_assert( pose.residue( root-1 ).is_protein() );
	Residue const & rootrsd( pose.residue(root) ), &braid1_anchor_rsd( pose.residue(braid_anchor) );
	Vector nxyz( braid1_anchor_rsd.xyz("N") ), caxyz( braid1_anchor_rsd.xyz("CA")), cxyz( braid1_anchor_rsd.xyz("C") );
	Stub const rootstub( rootrsd.xyz("ORIG"), rootrsd.xyz("X"), rootrsd.xyz("Y") );

	for ( Size i=1; i<= nbraid; ++i ) {
		Size const anchor( (i-1)*nres_braid + braid_anchor );
		runtime_assert( pose.fold_tree().jump_exists( anchor, root ) );
		RT const desired_rt( rootstub, Stub( nxyz, caxyz, cxyz ) );

		Residue const & rsd1( pose.residue(root) ), &rsd2( pose.residue(anchor) );
		id::StubID const stubid1( id::AtomID( rsd1.atom_index("ORIG" ), root ),
			id::AtomID( rsd1.atom_index("X"    ), root ),
			id::AtomID( rsd1.atom_index("Y"    ), root ) );

		id::StubID const stubid2( id::AtomID( rsd2.atom_index("N" ), anchor ),
			id::AtomID( rsd2.atom_index("CA"), anchor ),
			id::AtomID( rsd2.atom_index("C" ), anchor ) );

		// this routine is not symmetrized
		pose.conformation().set_stub_transform( stubid1, stubid2, desired_rt );
		//rotate the vectors
		nxyz  = R *  nxyz + v;
		caxyz = R * caxyz + v;
		cxyz  = R *  cxyz + v;
	}

}

///////////////////////////////////////////////////////////////////////////////////////
void
make_braided_frag_move(
	int const orientation,
	Size const nbraid,
	Size const fragsize,
	string const & tag,
	devel::blab::classic_frags::FragLib const & fraglib,
	Pose & pose,
	protocols::moves::MonteCarlo & mc
)
{
	Sizes const fragseq_poslist;
	Real d,alpha,rise;
	if ( orientation<0 ) get_antiparallel_braided_params( pose, d, alpha, rise );
	insert_protein_fragment( fragsize, fragseq_poslist, fraglib, pose ); // freeze_anchorseq
	if ( orientation>0 ) set_parallel_braided_pose( nbraid, pose );
	else if ( orientation<0 ) set_antiparallel_braided_pose( d, alpha, pose );
	// nothing to do if centroid_barrel (ie orientation=0)

	//Real const score( mc.score_function()( pose ) );
	bool const mc_accept( mc.boltzmann( pose, tag ) );
	TR_SYM_CENTROID.Trace << "MCFRAG " << mc_accept << ' ' << " mc_temp: " << F(9,3,mc.temperature() ) <<
		" twistrise_E: " << F(9,3,pose.energies().total_energies()[ twistrise ] ) <<
		" twistrise_W: " << F(9,3,mc.score_function().get_weight( twistrise ) ) <<
		" ap_E: " << F(9,3,pose.energies().total_energies()[ atom_pair_constraint ] ) <<
		" ap_W: " << F(9,3,mc.score_function().get_weight( atom_pair_constraint ) ) <<
		// " cb_E: " << F(9,3,pose.energies().total_energies()[ chainbreak ] ) <<
		// " cb_E: " << F(9,3,mc.score_function().get_weight( chainbreak ) ) <<
		" last: " << F(9,3,mc.last_accepted_score()) <<
		" low: " << F(9,3,mc.lowest_score() ) << ' ' << tag << endl;
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void
make_braid_move(
	int const orientation,
	Size const nbraid,
	string const & tag,
	Pose & pose,
	protocols::moves::MonteCarlo & mc
)
{
	if ( orientation>=0 ) return; // nothing to do
	runtime_assert( nbraid==2 );
	Real d,alpha,rise;
	get_antiparallel_braided_params( pose, d, alpha, rise );

	Real new_d( d + 0.25 * numeric::random::gaussian() * rise );
	while ( new_d >= rise/2 ) new_d -= rise;
	while ( new_d < -rise/2 ) new_d += rise;
	Real const new_alpha( alpha + numeric::random::gaussian() * numeric::conversions::radians( 30.0 ) );
	set_antiparallel_braided_pose( new_d, new_alpha, pose );

	bool const mc_accept( mc.boltzmann( pose, tag ) );
	TR_SYM_CENTROID.Trace << "MCBRAID " << mc_accept << ' ' << " mc_temp: " << F(9,3,mc.temperature() ) <<
		" twistrise_E: " << F(9,3,pose.energies().total_energies()[ twistrise ] ) <<
		" twistrise_W: " << F(9,3,mc.score_function().get_weight( twistrise ) ) <<
		" ap_E: " << F(9,3,pose.energies().total_energies()[ atom_pair_constraint ] ) <<
		" ap_W: " << F(9,3,mc.score_function().get_weight( atom_pair_constraint ) ) <<
		// " cb_E: " << F(9,3,pose.energies().total_energies()[ chainbreak ] ) <<
		// " cb_E: " << F(9,3,mc.score_function().get_weight( chainbreak ) ) <<
		" last: " << F(9,3,mc.last_accepted_score()) <<
		" low: " << F(9,3,mc.lowest_score() ) << ' ' << tag << endl;

}


///////////////////////////////////////////////////////////////////////////////////////

void
orient_trajectory_pose( Pose & pose )
{
	Size nrepeat, repeatlen, base_repeat;
	parse_symminfo( pose, nrepeat, repeatlen, base_repeat );

	pose::symmetry::make_asymmetric_pose( pose );

	Size const nres_protein( nrepeat * repeatlen );


	Stub const
		istub( pose.residue(1).xyz("CA"),
		pose.residue(2).xyz("CA"),
		pose.residue(3).xyz("CA") ),
		jstub( pose.residue(repeatlen+1).xyz("CA"),
		pose.residue(repeatlen+2).xyz("CA"),
		pose.residue(repeatlen+3).xyz("CA") );

	Vector center, symmetry_axis, t;
	Real theta;
	get_stub_transform_data( istub, jstub, center, symmetry_axis, t, theta );


	// transform so that symmetry axis is parallel to z-axis
	// projection of m_strand1_begin onto symmetry axis goes to (0,0,35)
	//
	Vector centroid(0,0,0);
	for ( Size i=1; i<= nres_protein; ++i ) centroid += pose.residue(i).xyz("CA");
	centroid/=nres_protein;
	Vector proj( centroid );
	proj = center + symmetry_axis * ( ( proj - center ).dot( symmetry_axis ) );
	runtime_assert( (proj-center).cross( symmetry_axis ).length()<1e-2 );

	Vector const dna_axis( 0,0,1 ), dna_center( 0,0,0 );

	Real const rotangle( std::acos( numeric::sin_cos_range( dna_axis.dot( symmetry_axis ) ) ) );
	Vector rotaxis( symmetry_axis.cross( dna_axis ) );
	rotaxis = rotaxis.normalized_any();
	numeric::xyzMatrix< Real > R( numeric::rotation_matrix( rotaxis, rotangle ) );
	runtime_assert( is_small( dna_axis.distance_squared( R*symmetry_axis ) ) );
	Vector const v( dna_center - R*proj );
	runtime_assert( is_small( dna_center.distance_squared( R*proj + v ) ) );

	pose.apply_transform_Rx_plus_v( R, v );

}


///////////////////////////////////////////////////////////////////////////////////////
void
simple_fold_abinitio(
	devel::blab::classic_frags::FragLib const & fraglib,
	Sizes const & fragseq_poslist, // empty if no seq changes
	Pose & pose,
	bool const use_twistrise_energy = false,
	Real const target_twist = 0.0,
	Real const target_rise = 0.0,
	Real const cycle_multiplier = 0.0,
	Real const twistrise_energy_weight = 1.0
)
{

	using namespace protocols::moves;
	using namespace devel::blab::classic_frags;

	Sizes const frag_sizes( fraglib.frag_sizes() );
	runtime_assert( frag_sizes.size() == 2 ); // not necessary
	Size const big_frag_size( max( frag_sizes ) ), small_frag_size( min( frag_sizes ) );

	bool const symmetry_was_on_at_start( conformation::symmetry::symmetry_is_on() );

	if (  symmetry_was_on_at_start && !pose::symmetry::is_symmetric( pose ) ) conformation::symmetry::turn_symmetry_off();
	if ( !symmetry_was_on_at_start &&  pose::symmetry::is_symmetric( pose ) ) conformation::symmetry::turn_symmetry_on();

	Size max_score0_cycles( 2000 );
	Size    score1_cycles( dry_run() ? 5 : 2000 );
	Size score2or5_cycles( dry_run() ? 5 : 2000 );
	Size    score3_cycles( dry_run() ? 5 : 4000 );
	Size score2or5_nloop( 10 );
	Size score3_nloop( 3 );

	if ( fragseq_poslist.size() ) { // not sure what the right scale factor is here...
		score1_cycles    *= 2;
		score2or5_cycles *= 2;
		score3_cycles    *= 2;
	}

	if ( cycle_multiplier>1e-3 ) {
		score1_cycles    = int( cycle_multiplier * Real( score1_cycles ) );
		score2or5_cycles = int( cycle_multiplier * Real( score2or5_cycles ) );
		score3_cycles    = int( cycle_multiplier * Real( score3_cycles ) );
	}

	ScoreFunctionOP score0( ScoreFunctionFactory::create_score_function( "score0.wts" ) );
	ScoreFunctionOP score1( ScoreFunctionFactory::create_score_function( "score1.wts" ) );
	ScoreFunctionOP score2( ScoreFunctionFactory::create_score_function( "score2.wts" ) );
	ScoreFunctionOP score3( ScoreFunctionFactory::create_score_function( "score3.wts" ) );
	ScoreFunctionOP score5( ScoreFunctionFactory::create_score_function( "score5.wts" ) );

	if ( pose::symmetry::is_symmetric( pose ) ) { runtime_assert( pose::symmetry::is_symmetric( *score0 ) ); }
	else { runtime_assert( !pose::symmetry::is_symmetric( *score0 ) ); }

	if ( use_twistrise_energy ) {
		TwistRiseEnergy twistrise_energy( target_twist, target_rise );
		score1->add_extra_method( twistrise, 0.1 * twistrise_energy_weight, twistrise_energy );
		score2->add_extra_method( twistrise, 0.3 * twistrise_energy_weight, twistrise_energy );
		score5->add_extra_method( twistrise, 0.2 * twistrise_energy_weight, twistrise_energy );
		score3->add_extra_method( twistrise, 0.5 * twistrise_energy_weight, twistrise_energy );
	}



	/// Monte Carlo setup
	Real const mc_temp( 2.0 );
	MonteCarloOP mc( new MonteCarlo( pose, *score0, mc_temp ) );


	/// score0
	mc->score_function( *score0 );
	for ( Size n=1; n<= max_score0_cycles; ++n ) {
		make_frag_move( big_frag_size, "stage0", fraglib, fragseq_poslist, pose, *mc );
		if ( !init_torsions_still_present( bools( pose.total_residue(), true ), pose ) ) {
			mc->reset( pose );
			break;
		}
	}
	mc->set_autotemp( true, mc_temp );
	mc->set_temperature( mc_temp );
	mc->reset( pose );

	// hacking
	static Size traj_run( 0 );
	++traj_run;
	string const traj_prefix("traj_"+lead_zero_string_of(traj_run,3));
	Size traj_counter(0);
	bool const make_trajectory( option[ my_options::make_trajectory ] );
	Pose trajectory_pose;
	protocols::moves::MCA const mc_status_for_pdb_output( MCA_accepted_score_beat_last );
	//protocols::moves::MCA const mc_status_for_pdb_output( MCA_accepted_score_beat_low );
	//bool const make_trajectory( false );


	{ /// score1 -- 9mers
		mc->score_function( *score1 );
		for ( Size n=1; n<= score1_cycles; ++n ) {
			make_frag_move( big_frag_size, "stage1", fraglib, fragseq_poslist, pose, *mc );
			if ( make_trajectory && mc->mc_accepted() == mc_status_for_pdb_output ) {
				++traj_counter;
				trajectory_pose = pose;
				orient_trajectory_pose( trajectory_pose );
				trajectory_pose.dump_pdb(traj_prefix+lead_zero_string_of(traj_counter,4)+".pdb");
			}
		}
		mc->recover_low( pose );
		mc->set_temperature( mc_temp );
	}

	{ /// alternating score2/score5 -- 9mers
		for ( Size m=1; m<= score2or5_nloop; ++m ) {
			if ( m%2==0 || m>7 ) { mc->score_function( *score2 ); }
			else { mc->score_function( *score5 ); }
			for ( Size n=1; n<= score2or5_cycles; ++n ) {
				make_frag_move( big_frag_size, "stage2or5", fraglib, fragseq_poslist, pose, *mc );
				if ( make_trajectory && mc->mc_accepted() == mc_status_for_pdb_output ) {
					++traj_counter;
					trajectory_pose = pose;
					orient_trajectory_pose( trajectory_pose );
					trajectory_pose.dump_pdb(traj_prefix+lead_zero_string_of(traj_counter,4)+".pdb");
				}
			}
			mc->recover_low( pose );
		}
		mc->set_temperature( mc_temp );
	}


	{ /// score3 -- 3mers
		for ( Size m=1; m<= score3_nloop; ++m ) {
			mc->score_function( *score3 );
			for ( Size n=1; n<= score3_cycles; ++n ) {
				make_frag_move( small_frag_size, "stage3", fraglib, fragseq_poslist, pose, *mc );
				if ( make_trajectory && mc->mc_accepted() == mc_status_for_pdb_output ) {
					++traj_counter;
					trajectory_pose = pose;
					orient_trajectory_pose( trajectory_pose );
					trajectory_pose.dump_pdb(traj_prefix+lead_zero_string_of(traj_counter,4)+".pdb");
				}
			}
			mc->recover_low( pose );
		}
	}

	mc->show_counters();
	mc->recover_low( pose );

	if ( symmetry_was_on_at_start ) conformation::symmetry::turn_symmetry_on();
	else conformation::symmetry::turn_symmetry_off();
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void
simple_fold_abinitio_peptide(
	devel::blab::classic_frags::FragLib const & fraglib,
	Sizes const & fragseq_poslist, // empty if no seq changes
	Pose & pose,
	bool const use_twistrise_energy = false,
	Real const target_twist = 0.0,
	Real const target_rise = 0.0,
	Real const cycle_multiplier = 0.0,
	Real const twistrise_energy_weight = 1.0
)
{

	using namespace protocols::moves;
	using namespace devel::blab::classic_frags;

	Sizes const frag_sizes( fraglib.frag_sizes() );
	runtime_assert( frag_sizes.size() == 2 ); // not necessary
	Size const big_frag_size( max( frag_sizes ) ), small_frag_size( min( frag_sizes ) );

	bool const symmetry_was_on_at_start( conformation::symmetry::symmetry_is_on() );

	if (  symmetry_was_on_at_start && !pose::symmetry::is_symmetric( pose ) ) conformation::symmetry::turn_symmetry_off();
	if ( !symmetry_was_on_at_start &&  pose::symmetry::is_symmetric( pose ) ) conformation::symmetry::turn_symmetry_on();

	Size max_score0_cycles( 2000 );
	Size    score1_cycles( dry_run() ? 5 : 2000 );
	Size score2or5_cycles( dry_run() ? 5 : 2000 );
	Size    score3_cycles( dry_run() ? 5 : 4000 );
	Size score2or5_nloop( 10 );
	Size score3_nloop( 3 );

	if ( fragseq_poslist.size() ) { // not sure what the right scale factor is here...
		score1_cycles    *= 2;
		score2or5_cycles *= 2;
		score3_cycles    *= 2;
	}

	if ( cycle_multiplier>1e-3 ) {
		score1_cycles    = int( cycle_multiplier * Real( score1_cycles ) );
		score2or5_cycles = int( cycle_multiplier * Real( score2or5_cycles ) );
		score3_cycles    = int( cycle_multiplier * Real( score3_cycles ) );
	}

	ScoreFunctionOP score0( ScoreFunctionFactory::create_score_function( "score0.wts" ) );
	ScoreFunctionOP score1( ScoreFunctionFactory::create_score_function( "score1.wts" ) );
	ScoreFunctionOP score2( ScoreFunctionFactory::create_score_function( "score2.wts" ) );
	ScoreFunctionOP score3( ScoreFunctionFactory::create_score_function( "score3.wts" ) );
	ScoreFunctionOP score5( ScoreFunctionFactory::create_score_function( "score5.wts" ) );

	if ( pose::symmetry::is_symmetric( pose ) ) { runtime_assert( pose::symmetry::is_symmetric( *score0 ) ); }
	else { runtime_assert( !pose::symmetry::is_symmetric( *score0 ) ); }

	if ( use_twistrise_energy ) {
		TwistRiseEnergy twistrise_energy( target_twist, target_rise );
		score1->add_extra_method( twistrise, 0.1 * twistrise_energy_weight, twistrise_energy );
		score2->add_extra_method( twistrise, 0.3 * twistrise_energy_weight, twistrise_energy );
		score5->add_extra_method( twistrise, 0.2 * twistrise_energy_weight, twistrise_energy );
		score3->add_extra_method( twistrise, 0.5 * twistrise_energy_weight, twistrise_energy );
	}


	/// backbone hbonds
	score1->set_weight( hbond_sr_bb, 0.5 );
	score1->set_weight( hbond_lr_bb, 0.5 );

	score2->set_weight( hbond_sr_bb, 0.5 );
	score2->set_weight( hbond_lr_bb, 0.5 );

	score5->set_weight( hbond_sr_bb, 1.5 );
	score5->set_weight( hbond_lr_bb, 1.5 );

	score3->set_weight( hbond_sr_bb, 1.0 );
	score3->set_weight( hbond_lr_bb, 1.0 );

	/// no compaction
	score1->set_weight( rg, 0.0 );
	score2->set_weight( rg, 0.0 );
	score5->set_weight( rg, 0.0 );
	score3->set_weight( rg, 0.0 );

	/// Monte Carlo setup
	Real const mc_temp( 2.0 );
	MonteCarloOP mc( new MonteCarlo( pose, *score0, mc_temp ) );


	/// score0
	mc->score_function( *score0 );
	for ( Size n=1; n<= max_score0_cycles; ++n ) {
		make_frag_move( big_frag_size, "stage0", fraglib, fragseq_poslist, pose, *mc );
		if ( !init_torsions_still_present( bools( pose.total_residue(), true ), pose ) ) {
			mc->reset( pose );
			break;
		}
	}
	mc->set_autotemp( true, mc_temp );
	mc->set_temperature( mc_temp );
	mc->reset( pose );

	// hacking
	static Size traj_run( 0 );
	++traj_run;
	string const traj_prefix("traj_"+lead_zero_string_of(traj_run,3));
	Size traj_counter(0);
	bool const make_trajectory( option[ my_options::make_trajectory ] );
	Pose trajectory_pose;
	protocols::moves::MCA const mc_status_for_pdb_output( MCA_accepted_score_beat_last );
	//protocols::moves::MCA const mc_status_for_pdb_output( MCA_accepted_score_beat_low );
	//bool const make_trajectory( false );


	{ /// score1 -- 9mers
		mc->score_function( *score1 );
		for ( Size n=1; n<= score1_cycles; ++n ) {
			make_frag_move( big_frag_size, "stage1", fraglib, fragseq_poslist, pose, *mc );
			if ( make_trajectory && mc->mc_accepted() == mc_status_for_pdb_output ) {
				++traj_counter;
				trajectory_pose = pose;
				orient_trajectory_pose( trajectory_pose );
				trajectory_pose.dump_pdb(traj_prefix+lead_zero_string_of(traj_counter,4)+".pdb");
			}
		}
		mc->recover_low( pose );
		mc->set_temperature( mc_temp );
	}

	{ /// alternating score2/score5 -- 9mers
		for ( Size m=1; m<= score2or5_nloop; ++m ) {
			if ( m%2==0 || m>7 ) { mc->score_function( *score2 ); }
			else { mc->score_function( *score5 ); }
			for ( Size n=1; n<= score2or5_cycles; ++n ) {
				make_frag_move( big_frag_size, "stage2or5", fraglib, fragseq_poslist, pose, *mc );
				if ( make_trajectory && mc->mc_accepted() == mc_status_for_pdb_output ) {
					++traj_counter;
					trajectory_pose = pose;
					orient_trajectory_pose( trajectory_pose );
					trajectory_pose.dump_pdb(traj_prefix+lead_zero_string_of(traj_counter,4)+".pdb");
				}
			}
			mc->recover_low( pose );
		}
		mc->set_temperature( mc_temp );
	}


	{ /// score3 -- 3mers
		for ( Size m=1; m<= score3_nloop; ++m ) {
			mc->score_function( *score3 );
			for ( Size n=1; n<= score3_cycles; ++n ) {
				make_frag_move( small_frag_size, "stage3", fraglib, fragseq_poslist, pose, *mc );
				if ( make_trajectory && mc->mc_accepted() == mc_status_for_pdb_output ) {
					++traj_counter;
					trajectory_pose = pose;
					orient_trajectory_pose( trajectory_pose );
					trajectory_pose.dump_pdb(traj_prefix+lead_zero_string_of(traj_counter,4)+".pdb");
				}
			}
			mc->recover_low( pose );
		}
	}

	mc->show_counters();
	mc->recover_low( pose );

	if ( symmetry_was_on_at_start ) conformation::symmetry::turn_symmetry_on();
	else conformation::symmetry::turn_symmetry_off();
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// this assumes that the jump points have been chosen so that the lower-numbered positions are
// donating hbonds to the residues before and after the upper-numbered positions
// (ie. the pleat at the beta-pairing is assumed to be lower-hbonded vs upper-hbonded)
//

void
make_beta_sheet_move(
	utility::vector1< utility::vector1< kinematics::RT > > const & beta_transforms,
	BetaPairingTypes const & jump_types,
	Pose & pose
)
{
	runtime_assert( pose::symmetry::is_symmetric( pose ) );
	conformation::symmetry::SymmetryInfo const & symminfo( *pose::symmetry::symmetry_info( pose ) );

	Size nres_protein( pose.total_residue() );
	while ( !pose.residue( nres_protein ).is_protein() ) --nres_protein;

	Sizes independent_jumps;
	kinematics::FoldTree const & f( pose.fold_tree() );
	for ( Size i=1; i<= f.num_jump(); ++i ) {
		if ( symminfo.jump_is_independent(i) &&
				f.upstream_jump_residue(i) <= nres_protein &&
				f.downstream_jump_residue(i) <= nres_protein ) { // intra-protein, independent jump
			independent_jumps.push_back( i );
		}
	}

	Size jumpno(0);
	{
		Size ntries(0);
		while ( true ) {
			++ntries;
			jumpno = random_element( independent_jumps );
			if ( jump_types[ jumpno ] != none_type ) break;
			runtime_assert( ntries<10000 );
		}
	}

	kinematics::RT const & rt( random_element( beta_transforms[ jump_types[ jumpno ] ] ) );

	Size const pos1( min( f.upstream_jump_residue( jumpno ), f.downstream_jump_residue( jumpno ) ) );
	Size const pos2( max( f.upstream_jump_residue( jumpno ), f.downstream_jump_residue( jumpno ) ) );

	Residue const & rsd1( pose.residue(pos1) ), &rsd2( pose.residue(pos2) );

	id::StubID const stubid1( id::AtomID( rsd1.atom_index("N" ), pos1 ),
		id::AtomID( rsd1.atom_index("CA"), pos1 ),
		id::AtomID( rsd1.atom_index("C" ), pos1 ) );
	id::StubID const stubid2( id::AtomID( rsd2.atom_index("N" ), pos2 ),
		id::AtomID( rsd2.atom_index("CA"), pos2 ),
		id::AtomID( rsd2.atom_index("C" ), pos2 ) );

	// this routine is not symmetrized
	pose.conformation().set_stub_transform( stubid1, stubid2, rt );

	// this will now symmetrize
	pose.set_jump( jumpno, pose.jump( jumpno ) );

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// this assumes that the jump points have been chosen so that the lower-numbered positions are
// donating hbonds to the residues before and after the upper-numbered positions
// (ie. the pleat at the beta-pairing is assumed to be lower-hbonded vs upper-hbonded)
//

// void
// make_beta_sheet_move(
//            utility::vector1< kinematics::RT > const & parallel_beta_transforms,
//            Pose & pose
//            )
// {
//  runtime_assert( pose::symmetry::is_symmetric( pose ) );
//  conformation::symmetry::SymmetryInfo const & symminfo( *pose::symmetry::symmetry_info( pose ) );

//  kinematics::RT const & rt( random_element( parallel_beta_transforms ) );

//  Size const nres_protein( pose.total_residue()-1 ); // ASSUMPTION

//  Sizes independent_jumps;
//  kinematics::FoldTree const & f( pose.fold_tree() );
//  for ( Size i=1; i<= f.num_jump(); ++i ) {
//   if ( symminfo.jump_is_independent(i) &&
//      f.upstream_jump_residue(i) <= (int) nres_protein &&
//      f.downstream_jump_residue(i) <= (int) nres_protein ) { // intra-protein, independent jump
//    independent_jumps.push_back( i );
//   }
//  }

//  Size const jumpno( random_element( independent_jumps ) );


//  Size const pos1( min( f.upstream_jump_residue( jumpno ), f.downstream_jump_residue( jumpno ) ) );
//  Size const pos2( max( f.upstream_jump_residue( jumpno ), f.downstream_jump_residue( jumpno ) ) );

//  Residue const & rsd1( pose.residue(pos1) ), &rsd2( pose.residue(pos2) );

//  id::StubID const stubid1( id::AtomID( rsd1.atom_index("N" ), pos1 ),
//               id::AtomID( rsd1.atom_index("CA"), pos1 ),
//               id::AtomID( rsd1.atom_index("C" ), pos1 ) );
//  id::StubID const stubid2( id::AtomID( rsd2.atom_index("N" ), pos2 ),
//               id::AtomID( rsd2.atom_index("CA"), pos2 ),
//               id::AtomID( rsd2.atom_index("C" ), pos2 ) );

//  // this routine is not symmetrized
//  pose.conformation().set_stub_transform( stubid1, stubid2, rt );

//  // this will now symmetrize
//  pose.set_jump( jumpno, pose.jump( jumpno ) );

// }


void
make_beta_move(
	utility::vector1< utility::vector1< kinematics::RT > > const & beta_transforms,
	BetaPairingTypes const & jump_types,
	Pose & pose,
	protocols::moves::MonteCarlo & mc
)
{
	make_beta_sheet_move( beta_transforms, jump_types, pose );

	string const tag( "beta" );
	bool const mc_accept( mc.boltzmann( pose, tag ) );

	TR_SYM_CENTROID.Trace << "MCBETA " << mc_accept << ' ' << " mc_temp: " << F(9,3,mc.temperature() ) <<
		" twistrise_E: " << F(9,3,pose.energies().total_energies()[ twistrise ] ) <<
		" twistrise_W: " << F(9,3,mc.score_function().get_weight( twistrise ) ) <<
		" ap_E: " << F(9,3,pose.energies().total_energies()[ atom_pair_constraint ] ) <<
		" ap_E: " << F(9,3,mc.score_function().get_weight( atom_pair_constraint ) ) <<
		" last: " << F(9,3,mc.last_accepted_score()) <<
		" low: " << F(9,3,mc.lowest_score() ) << ' ' << tag << endl;
}

void
make_braided_beta_move(
	int const orientation,
	Size const nbraid,
	utility::vector1< utility::vector1< kinematics::RT > > const & beta_transforms,
	BetaPairingTypes const & jump_types,
	Pose & pose,
	protocols::moves::MonteCarlo & mc
)
{
	if ( orientation == 0 ) { // silly hack to signal centroid_barrel
		make_braided_barrel_beta_move( beta_transforms, pose );
		// if ( uniform()<0.25 ) {
		//  make_braided_barrel_beta_move( beta_transforms, pose );
		// } else {
		//  make_braided_barrel_r_theta_move( 1.5 /*r_sdev*/, 30.0 /*theta_sdev*/, pose );
		// }
	} else {
		Real d,alpha,rise;
		if ( orientation<0 ) get_antiparallel_braided_params( pose, d, alpha, rise );
		make_beta_sheet_move( beta_transforms, jump_types, pose );
		if ( orientation>0 ) set_parallel_braided_pose( nbraid, pose );
		else set_antiparallel_braided_pose( d, alpha, pose );
	}
	string const tag( "beta" );
	bool const mc_accept( mc.boltzmann( pose, tag ) );

	TR_SYM_CENTROID.Trace << "MCBETA " << mc_accept << ' ' << " mc_temp: " << F(9,3,mc.temperature() ) <<
		" twistrise_E: " << F(9,3,pose.energies().total_energies()[ twistrise ] ) <<
		" twistrise_W: " << F(9,3,mc.score_function().get_weight( twistrise ) ) <<
		" ap_E: " << F(9,3,pose.energies().total_energies()[ atom_pair_constraint ] ) <<
		" ap_E: " << F(9,3,mc.score_function().get_weight( atom_pair_constraint ) ) <<
		" last: " << F(9,3,mc.last_accepted_score()) <<
		" low: " << F(9,3,mc.lowest_score() ) << ' ' << tag << endl;
}

// void
// make_beta_move(
//         utility::vector1< kinematics::RT > const & parallel_beta_transforms,
//         Pose & pose,
//         protocols::moves::MonteCarlo & mc
//         )
// {
//  make_beta_sheet_move( parallel_beta_transforms, pose );

//  string const tag( "beta" );
//  bool const mc_accept( mc.boltzmann( pose, tag ) );

//  TR_SYM_CENTROID.Trace << "MCBETA " << mc_accept << ' ' << " mc_temp: " << F(9,3,mc.temperature() ) <<
//   " twistrise_E: " << F(9,3,pose.energies().total_energies()[ twistrise ] ) <<
//   " twistrise_W: " << F(9,3,mc.score_function().get_weight( twistrise ) ) <<
//   " ap_E: " << F(9,3,pose.energies().total_energies()[ atom_pair_constraint ] ) <<
//   " ap_E: " << F(9,3,mc.score_function().get_weight( atom_pair_constraint ) ) <<
//   " last: " << F(9,3,mc.last_accepted_score()) <<
//   " low: " << F(9,3,mc.lowest_score() ) << ' ' << tag << endl;
// }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// does not initialize the jumps, so caller should already have set a reasonable/random transform
///////////////////////////////////////////////////////////////////////////////////////
void
simple_fold_abinitio_jumping(
	devel::blab::classic_frags::FragLib const & fraglib,
	vector1< vector1< kinematics::RT > > const & beta_transforms,
	BetaPairingTypes const & jump_types,
	Pose & pose,
	bool const use_twistrise_energy = false,
	Real const target_twist = 0.0,
	Real const target_rise = 0.0,
	Real const cycle_multiplier = 0.0
)
{
	bool const symmetry_was_on_at_start( conformation::symmetry::symmetry_is_on() );

	if (  symmetry_was_on_at_start && !pose::symmetry::is_symmetric( pose ) ) conformation::symmetry::turn_symmetry_off();
	if ( !symmetry_was_on_at_start &&  pose::symmetry::is_symmetric( pose ) ) conformation::symmetry::turn_symmetry_on();

	using namespace protocols::moves;
	using namespace devel::blab::classic_frags;

	Real const jump_move_frequency( 0.1 );


	Size max_score0_cycles( 2000 );
	Size    score1_cycles( dry_run() ? 5 : 2000 );
	Size score2or5_cycles( dry_run() ? 5 : 2000 );
	Size    score3_cycles( dry_run() ? 5 : 4000 );
	Size score2or5_nloop( 10 );
	Size score3_nloop( 3 );

	if ( cycle_multiplier>1e-3 ) {
		score1_cycles    = int( cycle_multiplier * Real( score1_cycles ) );
		score2or5_cycles = int( cycle_multiplier * Real( score2or5_cycles ) );
		score3_cycles    = int( cycle_multiplier * Real( score3_cycles ) );
	}


	ScoreFunctionOP score0( ScoreFunctionFactory::create_score_function( "score0.wts" ) );
	ScoreFunctionOP score1( ScoreFunctionFactory::create_score_function( "score1.wts" ) );
	ScoreFunctionOP score2( ScoreFunctionFactory::create_score_function( "score2.wts" ) );
	ScoreFunctionOP score3( ScoreFunctionFactory::create_score_function( "score3.wts" ) );
	ScoreFunctionOP score5( ScoreFunctionFactory::create_score_function( "score5.wts" ) );

	if ( pose::symmetry::is_symmetric( pose ) ) { runtime_assert( pose::symmetry::is_symmetric( *score0 ) ); }
	else { runtime_assert( !pose::symmetry::is_symmetric( *score0 ) ); }
	//  if ( pose::symmetry::is_symmetric( *score0 ) ) {
	//   /// must be the stupid commandline option
	//   /// need to convert to non-symmetric versions
	//   score0 = new ScoreFunction( *score0 );
	//   score1 = new ScoreFunction( *score1 );
	//   score2 = new ScoreFunction( *score2 );
	//   score3 = new ScoreFunction( *score3 );
	//   score5 = new ScoreFunction( *score5 );
	//   runtime_assert( !pose::symmetry::is_symmetric( *score0 ) );
	//  }
	// }


	if ( use_twistrise_energy ) {
		TwistRiseEnergy twistrise_energy( target_twist, target_rise );

		score1->add_extra_method( twistrise, 0.1 , twistrise_energy );
		score2->add_extra_method( twistrise, 0.3 , twistrise_energy );
		score5->add_extra_method( twistrise, 0.2 , twistrise_energy );
		score3->add_extra_method( twistrise, 0.5 , twistrise_energy );
	}

	// ramp the constraint weight?
	// this is being used to handle chainbreaks...
	score1->set_weight( atom_pair_constraint, 0.1 );
	score2->set_weight( atom_pair_constraint, 0.2 );
	score5->set_weight( atom_pair_constraint, 0.5 );
	score3->set_weight( atom_pair_constraint, 1.0 );


	/// Monte Carlo setup
	Real const mc_temp( 2.0 );
	MonteCarloOP mc( new MonteCarlo( pose, *score0, mc_temp ) );

	Sizes const fragseq_poslist; // empty

	/// score0
	mc->score_function( *score0 );
	for ( Size n=1; n<= max_score0_cycles; ++n ) {
		if ( uniform() < jump_move_frequency ) make_beta_move( beta_transforms, jump_types, pose, *mc );
		else make_frag_move( 9, "stage0", fraglib, fragseq_poslist, pose, *mc );
		if ( !init_torsions_still_present( bools( pose.total_residue(), true ), pose ) ) {
			mc->reset( pose );
			break;
		}
	}
	mc->set_autotemp( true, mc_temp );
	mc->set_temperature( mc_temp );
	mc->reset( pose );

	// hacking
	static Size traj_run( 0 );
	++traj_run;
	string const traj_prefix("traj_"+lead_zero_string_of(traj_run,3));
	Size traj_counter(0);
	bool const make_trajectory( option[ my_options::make_trajectory ] );



	{ /// score1 -- 9mers
		mc->score_function( *score1 );
		for ( Size n=1; n<= score1_cycles; ++n ) {
			if ( uniform() < jump_move_frequency ) {
				make_beta_move( beta_transforms, jump_types, pose, *mc );
			} else {
				make_frag_move( 9, "stage1", fraglib, fragseq_poslist, pose, *mc );
			}
			if ( make_trajectory && mc->mc_accepted() == protocols::moves::MCA_accepted_score_beat_low ) {
				++traj_counter;
				pose.dump_pdb(traj_prefix+lead_zero_string_of(traj_counter,4)+".pdb");
			}
		}
		mc->recover_low( pose );
		mc->set_temperature( mc_temp );
	}

	{ /// alternating score2/score5 -- 9mers
		for ( Size m=1; m<= score2or5_nloop; ++m ) {
			if ( m%2==0 || m>7 ) { mc->score_function( *score2 ); }
			else { mc->score_function( *score5 ); }
			for ( Size n=1; n<= score2or5_cycles; ++n ) {
				if ( uniform() < jump_move_frequency ) {
					make_beta_move( beta_transforms, jump_types, pose, *mc );
				} else {
					make_frag_move( 9, "stage2or5", fraglib, fragseq_poslist, pose, *mc );
				}
				if ( make_trajectory && mc->mc_accepted() == protocols::moves::MCA_accepted_score_beat_low ) {
					++traj_counter;
					pose.dump_pdb(traj_prefix+lead_zero_string_of(traj_counter,4)+".pdb");
				}
			}
			mc->recover_low( pose );
		}
		mc->set_temperature( mc_temp );
	}


	{ /// score3 -- 3mers
		for ( Size m=1; m<= score3_nloop; ++m ) {
			mc->score_function( *score3 );
			for ( Size n=1; n<= score3_cycles; ++n ) {
				if ( uniform() < jump_move_frequency ) {
					make_beta_move( beta_transforms, jump_types, pose, *mc );
				} else {
					make_frag_move( 3, "stage3", fraglib, fragseq_poslist, pose, *mc );
				}
				if ( make_trajectory && mc->mc_accepted() == protocols::moves::MCA_accepted_score_beat_low ) {
					++traj_counter;
					pose.dump_pdb(traj_prefix+lead_zero_string_of(traj_counter,4)+".pdb");
				}
			}
			mc->recover_low( pose );
		}
	}

	mc->show_counters();
	mc->recover_low( pose );

	if ( symmetry_was_on_at_start ) conformation::symmetry::turn_symmetry_on();
	else conformation::symmetry::turn_symmetry_off();
}


///////////////////////////////////////////////////////////////////////////////////////
void
simple_fold_abinitio_braided(
	int const orientation,
	Size const nbraid,
	BetaPairingTypes const & jump_types,
	vector1< vector1< kinematics::RT > > const & beta_transforms,
	devel::blab::classic_frags::FragLib const & fraglib,
	Pose & pose,
	bool const use_twistrise_energy = false,
	Real const target_twist = 0.0,
	Real const target_rise = 0.0
)
{
	using namespace protocols::moves;
	using namespace devel::blab::classic_frags;

	bool const symmetry_was_on_at_start( conformation::symmetry::symmetry_is_on() );
	bool const make_jump_moves( orientation == 0 || ( !jump_types.empty() && !beta_transforms.empty() ) );
	Real const jump_move_frequency( make_jump_moves ? 0.1 : -1 );

	if (  symmetry_was_on_at_start && !pose::symmetry::is_symmetric( pose ) ) conformation::symmetry::turn_symmetry_off();
	if ( !symmetry_was_on_at_start &&  pose::symmetry::is_symmetric( pose ) ) conformation::symmetry::turn_symmetry_on();

	Size max_score0_cycles( 2000 );
	Size    score1_cycles( dry_run() ? 5 : 2000 );
	Size score2or5_cycles( dry_run() ? 5 : 2000 );
	Size    score3_cycles( dry_run() ? 5 : 4000 );
	Size score2or5_nloop( 10 );
	Size score3_nloop( 3 );


	// if ( cycle_multiplier>1e-3 ) {
	//  score1_cycles    = int( cycle_multiplier * Real( score1_cycles ) );
	//  score2or5_cycles = int( cycle_multiplier * Real( score2or5_cycles ) );
	//  score3_cycles    = int( cycle_multiplier * Real( score3_cycles ) );
	// }

	ScoreFunctionOP score0( ScoreFunctionFactory::create_score_function( "score0.wts" ) );
	ScoreFunctionOP score1( ScoreFunctionFactory::create_score_function( "score1.wts" ) );
	ScoreFunctionOP score2( ScoreFunctionFactory::create_score_function( "score2.wts" ) );
	ScoreFunctionOP score3( ScoreFunctionFactory::create_score_function( "score3.wts" ) );
	ScoreFunctionOP score5( ScoreFunctionFactory::create_score_function( "score5.wts" ) );

	if ( pose::symmetry::is_symmetric( pose ) ) { runtime_assert( pose::symmetry::is_symmetric( *score0 ) ); }
	else { runtime_assert( !pose::symmetry::is_symmetric( *score0 ) ); }

	if ( use_twistrise_energy ) {
		TwistRiseEnergy twistrise_energy( target_twist, target_rise );
		score1->add_extra_method( twistrise, 0.1 , twistrise_energy );
		score2->add_extra_method( twistrise, 0.3 , twistrise_energy );
		score5->add_extra_method( twistrise, 0.2 , twistrise_energy );
		score3->add_extra_method( twistrise, 0.5 , twistrise_energy );
	}

	if ( make_jump_moves ) {
		// this is being used to handle chainbreaks...
		score1->set_weight( atom_pair_constraint, 0.1 );
		score2->set_weight( atom_pair_constraint, 0.2 );
		score5->set_weight( atom_pair_constraint, 0.5 );
		score3->set_weight( atom_pair_constraint, 1.0 );
	}

	/// Monte Carlo setup
	Real const mc_temp( 2.0 );
	MonteCarloOP mc( new MonteCarlo( pose, *score0, mc_temp ) );


	/// score0
	mc->score_function( *score0 );
	for ( Size n=1; n<= max_score0_cycles; ++n ) {
		if ( uniform() < jump_move_frequency ) {
			make_braided_beta_move( orientation, nbraid, beta_transforms, jump_types, pose, *mc );
		} else make_braided_frag_move( orientation, nbraid, 9, "stage0", fraglib, pose, *mc );
		// nothing happens if parallel:
		make_braid_move( orientation, nbraid, "stage0", pose, *mc );
		if ( !init_torsions_still_present( bools( pose.total_residue(), true ), pose ) ) {
			mc->reset( pose );
			break;
		}
	}
	mc->set_autotemp( true, mc_temp );
	mc->set_temperature( mc_temp );
	mc->reset( pose );

	// hacking
	static Size traj_run( 0 );
	++traj_run;
	string const traj_prefix("traj_"+lead_zero_string_of(traj_run,3));
	Size traj_counter(0);
	bool const make_trajectory( option[ my_options::make_trajectory ] );
	Pose trajectory_pose;
	protocols::moves::MCA const mc_status_for_pdb_output( MCA_accepted_score_beat_last );
	//protocols::moves::MCA const mc_status_for_pdb_output( MCA_accepted_score_beat_low );
	//bool const make_trajectory( false );


	{ /// score1 -- 9mers
		mc->score_function( *score1 );
		for ( Size n=1; n<= score1_cycles; ++n ) {
			if ( uniform() < jump_move_frequency ) {
				make_braided_beta_move( orientation, nbraid, beta_transforms, jump_types, pose, *mc );
			} else make_braided_frag_move( orientation, nbraid, 9, "stage1", fraglib, pose, *mc );
			if ( make_trajectory && mc->mc_accepted() == mc_status_for_pdb_output ) {
				++traj_counter;
				trajectory_pose = pose;
				orient_trajectory_pose( trajectory_pose );
				trajectory_pose.dump_pdb(traj_prefix+lead_zero_string_of(traj_counter,4)+".pdb");
			}
			make_braid_move( orientation, nbraid, "stage1", pose, *mc );
		}
		mc->recover_low( pose );
		mc->set_temperature( mc_temp );
	}

	{ /// alternating score2/score5 -- 9mers
		for ( Size m=1; m<= score2or5_nloop; ++m ) {
			if ( m%2==0 || m>7 ) { mc->score_function( *score2 ); }
			else { mc->score_function( *score5 ); }
			for ( Size n=1; n<= score2or5_cycles; ++n ) {
				if ( uniform() < jump_move_frequency ) {
					make_braided_beta_move( orientation, nbraid, beta_transforms, jump_types, pose, *mc );
				} else make_braided_frag_move( orientation, nbraid, 9, "stage2or5", fraglib, pose, *mc );
				if ( make_trajectory && mc->mc_accepted() == mc_status_for_pdb_output ) {
					++traj_counter;
					trajectory_pose = pose;
					orient_trajectory_pose( trajectory_pose );
					trajectory_pose.dump_pdb(traj_prefix+lead_zero_string_of(traj_counter,4)+".pdb");
				}
				make_braid_move( orientation, nbraid, "stage2or5", pose, *mc );
			}
			mc->recover_low( pose );
		}
		mc->set_temperature( mc_temp );
	}


	{ /// score3 -- 3mers
		for ( Size m=1; m<= score3_nloop; ++m ) {
			mc->score_function( *score3 );
			for ( Size n=1; n<= score3_cycles; ++n ) {
				if ( uniform() < jump_move_frequency ) {
					make_braided_beta_move( orientation, nbraid, beta_transforms, jump_types, pose, *mc );
				} else make_braided_frag_move( orientation, nbraid, 3, "stage3", fraglib, pose, *mc );
				if ( make_trajectory && mc->mc_accepted() == mc_status_for_pdb_output ) {
					++traj_counter;
					trajectory_pose = pose;
					orient_trajectory_pose( trajectory_pose );
					trajectory_pose.dump_pdb(traj_prefix+lead_zero_string_of(traj_counter,4)+".pdb");
				}
				make_braid_move( orientation, nbraid, "stage3", pose, *mc );
			}
			mc->recover_low( pose );
		}
	}

	mc->show_counters();
	mc->recover_low( pose );

	if ( symmetry_was_on_at_start ) conformation::symmetry::turn_symmetry_on();
	else conformation::symmetry::turn_symmetry_off();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// create a symmetric, protein-only pose with the same sequence
/// pick fragments using current ss and sequence
/// call the simple fold_abinitio routine above
///
Real // returns final score3 value
refold_repeat_pose(
	string const & repeatseq,
	Size nrepeat,
	Size const base_repeat,
	devel::blab::classic_frags::FragLib const & fraglib,
	Pose & pose,
	bool const abrelax = false,
	ScoreFunctionOP fa_scorefxn = 0,
	Size const abrelax_nrepeats_to_trim = 0,
	bool const use_asymmetric_centroid_scoring = false
)
{
	ScoreFunctionOP score3( ScoreFunctionFactory::create_score_function( "score3.wts" ) );
	runtime_assert( pose::symmetry::is_symmetric( *score3 ) );

	Size const repeatlen( repeatseq.size() );

	pose.clear();
	runtime_assert( !pose::symmetry::is_symmetric( pose ) );

	for ( Size i=1; i<= nrepeat; ++i ) {
		for ( Size j=0; j< repeatlen; ++j ) {
			pose.append_residue_by_bond( *get_vanilla_protein_residue( repeatseq[j] ), true ); // build_ideal_geometry
		}
	}


	add_lower_terminus_type_to_pose_residue( pose, 1 ); // not sure about this

	{ /// ADD A VIRTUAL RESIDUE AT THE END
		ResidueOP vrtrsd
			( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map( "VRTBB" ) ) );
		pose.append_residue_by_bond( *vrtrsd, true ); // since is polymer...

		kinematics::FoldTree f( pose.total_residue() );
		f.reorder( pose.total_residue() );
		pose.fold_tree( f );
	}

	Size nres_protein( repeatlen * nrepeat );
	pose.conformation().insert_chain_ending( nres_protein );
	runtime_assert( nres_protein == pose.total_residue()-1 );
	runtime_assert( num_chains( pose ) == 2 );
	runtime_assert( chain_end( 1, pose ) == nres_protein );

	for ( Size i=1; i<= nres_protein; ++i ) {
		pose.set_phi  ( i, init_phi   );
		pose.set_psi  ( i, init_psi   );
		pose.set_omega( i, init_omega );
		pose.set_secstruct( i, 'L' );
	}

	/// setup symminfo?
	conformation::symmetry::SymmetryInfo symminfo;
	setup_unbound_frag_symminfo( repeatlen, nrepeat, base_repeat, symminfo );


	// if ( !fraglib ) { // pick vall fragments
	//  runtime_assert( repeatseq.size() <= fragss_in.size() );
	//  runtime_assert( fragss_in.size()%repeatlen == 0 );
	//  runtime_assert( fragss_in.size() >= repeatlen );

	//  /// pick fragments
	//  Real const seq_weight( 1.0 ), ss_weight( 10.0 ); // ss trumps sequence...
	//  Sizes const fragsizes( make_vector1( 3, 9 ) );
	//  Size const nfrags( 200 );
	//  kinematics::MoveMap mm;
	//  mm.set_bb( true );
	//  mm.set_bb( pose.total_residue(), false );

	//  string fragss( fragss_in );
	//  if ( fragss.size() > nrepeat * repeatlen ) fragss.erase( nrepeat*repeatlen );
	//  while ( fragss.size() < nrepeat * repeatlen ) {
	//   Size const nrepeat_input_ss( fragss_in.size() / repeatlen );
	//   Size const reuse_repeat( 1+int(uniform()*nrepeat_input_ss) );
	//   fragss += fragss.substr( (reuse_repeat-1)*repeatlen, repeatlen );
	//   runtime_assert( fragss.size()%repeatlen == 0 );
	//  }
	//  runtime_assert( fragss.size() == nrepeat*repeatlen );
	//  TR_SYM_CENTROID.Trace << "refold_repeat_pose: pick_frags_for_ss: " << fragss << endl;

	//  //FragLib fraglib;
	//  fraglib = new devel::blab::classic_frags::FragLib();
	//  devel::blab::classic_frags::add_vall_fragments( fragsizes, nfrags, pose, mm, fragss, seq_weight, ss_weight,
	//                          *fraglib );
	// }

	/// switch to centroid
	devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID ); // was CENTROID_DNA

	/// now make symmetric
	if ( !use_asymmetric_centroid_scoring ) pose::symmetry::make_symmetric_pose( pose, symminfo );
	else {
		/// need this in order to preserve structural symmetry...
		symminfo_hack::repeatlen_ = repeatlen;
		symminfo_hack::nrepeat_ = nrepeat;
		symminfo_hack::base_repeat_ = base_repeat;
		runtime_assert( conformation::symmetry::symmetry_is_on() ); // hack
		conformation::symmetry::turn_symmetry_off();
		score3 = ScoreFunctionFactory::create_score_function( "score3.wts" );
		conformation::symmetry::turn_symmetry_on();
		runtime_assert( !pose::symmetry::is_symmetric( *score3 ) );
	}

	Sizes fragseq_poslist;
	simple_fold_abinitio( fraglib, fragseq_poslist, pose );

	Real final_score( (*score3)( pose ) );

	/// now switch back to fullatom, have to make asymmetric first:
	if ( !use_asymmetric_centroid_scoring ) {
		pose::symmetry::make_asymmetric_pose( pose );
	} else { // reset the hacky arrays
		symminfo_hack::repeatlen_   = 0;
		symminfo_hack::nrepeat_     = 0;
		symminfo_hack::base_repeat_ = 0;
	}

	devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
	pose::symmetry::make_symmetric_pose( pose, symminfo );

	if ( abrelax ) {
		/// do a fastrelax here
		if ( abrelax_nrepeats_to_trim ) {
			nrepeat -= abrelax_nrepeats_to_trim;
			nres_protein -= abrelax_nrepeats_to_trim*repeatlen;
			conformation::symmetry::SymmetryInfo symminfo2;
			setup_unbound_frag_symminfo( repeatlen, nrepeat, base_repeat, symminfo2 );
			pose::symmetry::make_asymmetric_pose( pose );
			Size const delete_begin( nrepeat*repeatlen + 1 ), delete_end( (nrepeat+abrelax_nrepeats_to_trim)*repeatlen );
			pose.conformation().delete_residue_range_slow( delete_begin, delete_end );
			pose::symmetry::make_symmetric_pose( pose, symminfo2 );
		}

		/// now we are symmetric, again, and have trimmed some repeats if desired
		protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

		MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
		movemap->set_bb (true);
		movemap->set_chi(true);
		movemap->set_bb ( pose.total_residue(), false );
		movemap->set_chi( pose.total_residue(), false );
		fastrelax.set_movemap( movemap );
		if ( !dry_run() ) fastrelax.apply( pose );

		final_score = (*fa_scorefxn)( pose );
	}

	return final_score;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// helper version that sets up vall-based fraglib first
//
Real // returns final score3 value
refold_repeat_pose(
	string const & repeatseq,
	string const & fragss_in,
	Size nrepeat,
	Size const base_repeat,
	Pose & pose,
	bool const abrelax = false,
	ScoreFunctionOP fa_scorefxn = 0,
	Size const abrelax_nrepeats_to_trim = 0,
	bool const use_asymmetric_centroid_scoring = false
)
{
	/// setup the fraglib using fragss_in
	devel::blab::classic_frags::FragLib fraglib;

	Size const repeatlen( repeatseq.size() );
	runtime_assert( repeatseq.size() <= fragss_in.size() );
	runtime_assert( fragss_in.size()%repeatlen == 0 );
	runtime_assert( fragss_in.size() >= repeatlen );

	/// pick fragments
	Real const seq_weight( 1.0 ), ss_weight( 10.0 ); // ss trumps sequence...
	Sizes const fragsizes( make_vector1( 3, 9 ) );
	Size const nfrags( 200 );
	kinematics::MoveMap mm;
	mm.set_bb( true );
	mm.set_bb( pose.total_residue(), false );

	string fragss( fragss_in );
	if ( fragss.size() > nrepeat * repeatlen ) fragss.erase( nrepeat*repeatlen );
	while ( fragss.size() < nrepeat * repeatlen ) {
		Size const nrepeat_input_ss( fragss_in.size() / repeatlen );
		Size const reuse_repeat( 1+int(uniform()*nrepeat_input_ss) );
		fragss += fragss.substr( (reuse_repeat-1)*repeatlen, repeatlen );
		runtime_assert( fragss.size()%repeatlen == 0 );
	}
	runtime_assert( fragss.size() == nrepeat*repeatlen );
	TR_SYM_CENTROID.Trace << "refold_repeat_pose: pick_frags_for_ss: " << fragss << endl;

	{
		Pose fragpose; // hacking
		for ( Size i=1; i<= nrepeat; ++i ) {
			for ( Size j=0; j< repeatlen; ++j ) {
				fragpose.append_residue_by_bond( *get_vanilla_protein_residue( repeatseq[j] ), true ); // build_ideal_geometry
			}
		}

		devel::blab::classic_frags::add_vall_fragments( fragsizes, nfrags, fragpose, mm, fragss, seq_weight, ss_weight,
			fraglib );
	}

	// now call other version
	return refold_repeat_pose( repeatseq, nrepeat, base_repeat, fraglib, pose,
		abrelax, fa_scorefxn, abrelax_nrepeats_to_trim,
		use_asymmetric_centroid_scoring );
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// helper version that sets up vall-based fraglib first
//
Real // returns final score3 value
refold_repeat_pose(
	string const & repeatseq,
	vector1< map< char, Real > > const & predss_in,
	Size nrepeat,
	Size const base_repeat,
	Pose & pose,
	bool const abrelax = false,
	ScoreFunctionOP fa_scorefxn = 0,
	Size const abrelax_nrepeats_to_trim = 0,
	bool const use_asymmetric_centroid_scoring = false
)
{
	/// setup the fraglib using fragss_in
	devel::blab::classic_frags::FragLib fraglib;

	Size const repeatlen( repeatseq.size() );
	runtime_assert( repeatseq.size() <= predss_in.size() );
	runtime_assert( predss_in.size()%repeatlen == 0 );
	runtime_assert( predss_in.size() >= repeatlen );

	/// pick fragments
	Real const seq_weight( 1.0 ), ss_weight( 10.0 ); // ss trumps sequence... maybe this is too high?
	Sizes const fragsizes( make_vector1( 3, 9 ) );
	Size const nfrags( 200 );
	kinematics::MoveMap mm;
	mm.set_bb( true );
	mm.set_bb( pose.total_residue(), false );

	vector1< map< char, Real > > predss( predss_in );
	if ( predss.size() > nrepeat * repeatlen ) predss.resize( nrepeat * repeatlen );
	while ( predss.size() < nrepeat * repeatlen ) {
		Size const nrepeat_input_ss( predss_in.size() / repeatlen );
		Size const reuse_repeat( 1+int(uniform()*nrepeat_input_ss) );
		for ( Size i=1; i<= repeatlen; ++i ) predss.push_back( predss_in[ (reuse_repeat-1)*repeatlen+i ] );
		runtime_assert( predss.size()%repeatlen == 0 );
	}
	runtime_assert( predss.size() == nrepeat*repeatlen );
	Real predss_conf(0.0);
	for ( Size i=1; i<= predss.size(); ++i ) predss_conf += max( predss[i]['H'], max( predss[i]['E'], predss[i]['L'] ));
	predss_conf /= predss.size();

	TR_SYM_CENTROID.Trace << "refold_repeat_pose: pick_frags_for_ss: predss_conf= " << F(9,3,predss_conf) << endl;

	{
		Pose fragpose; // hacking
		for ( Size i=1; i<= nrepeat; ++i ) {
			for ( Size j=0; j< repeatlen; ++j ) {
				fragpose.append_residue_by_bond( *get_vanilla_protein_residue( repeatseq[j] ), true ); // build_ideal_geometry
			}
		}

		devel::blab::classic_frags::add_vall_fragments( fragsizes, nfrags, fragpose, mm, predss, seq_weight, ss_weight,
			fraglib );
	}

	// now call other version
	return refold_repeat_pose( repeatseq, nrepeat, base_repeat, fraglib, pose,
		abrelax, fa_scorefxn, abrelax_nrepeats_to_trim,
		use_asymmetric_centroid_scoring );
}






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
