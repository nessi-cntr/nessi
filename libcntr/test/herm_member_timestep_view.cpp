#include "catch.hpp"
#include "cntr.hpp"

// TODO:
// Set matrixelement
// Incr
// Right/left multiply
// Interaction with herm_matrix_timestep 

#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
#define GREEN_TSTP_VIEW cntr::herm_matrix_timestep_view<double>
#define cfunction cntr::function<double>

TEST_CASE("Herm member timestep view","[Herm_member_timestep_view]"){
	int nt=50;
	int ntau=500;
	int size=2;
	int size2=5;
	double eps=1e-7;
	double beta=5.0;
	double h=0.01;
	// Square Green's function
	GREEN A=GREEN(nt,ntau,size,-1);
	cdmatrix a(size,size);
	a.setZero();
	a(0,0)=sqrt(2.0);
	a(0,1)=sqrt(2.0)*std::complex<double>(0.0,1.0);
	a(1,0)=sqrt(2.0)*std::complex<double>(0.0,-1.0);
	a(1,1)=-sqrt(2.0);
	cntr::green_from_H(A,0.0,a,beta,h,5,4,true);

	cfunction unity(nt,size);
	cdmatrix one(size,size);
	one.setZero();
	one(0,0)=1.0;
	one(0,1)=0.0;
	one(1,0)=0.0;
	one(1,1)=1.0;
	unity.set_constant(one);

	SECTION ("Set/Get herm_matrix_timestep_view"){
		// Check constructor from herm_matrix
		double err=0.0;
		for(int tstp=-1;tstp<nt;tstp++){
			GREEN_TSTP_VIEW Aview(tstp,A);
			err+=distance_norm2(tstp,Aview,A);
		}
		REQUIRE(err<eps);
		//Set/get herm_matrix_timestep
		err=0.0;
		for(int tstp=-1;tstp<nt;tstp++){
			GREEN_TSTP Atstp(tstp,ntau,size);
			A.get_timestep(tstp,Atstp);
			GREEN_TSTP_VIEW Aview(Atstp);
			err+=distance_norm2(tstp,Aview,Atstp);
		}
		REQUIRE(err<eps);
	}

	SECTION ("Increase and multiply"){
		double err=0.0;
		GREEN A2(nt,ntau,size);
		A2=A;
		// Test increase A+A=2A
		// TODO : Can't properly increase 
		for(int tstp=-1;tstp<nt;tstp++){
			GREEN_TSTP_VIEW Aview(tstp,A),Aview2(tstp,A);
			A2.right_multiply(tstp,unity,2.0);
			Aview2.incr_timestep(Aview,1.0);
			// err+=distance_norm2(tstp,Aview,A2);
		}
		REQUIRE(err<eps);

		// Test multiply 2A=2A
		err=0.0;
		A2=A;
		for(int tstp=-1;tstp<nt;tstp++){
			GREEN_TSTP_VIEW Aview(tstp,A);
			A2.right_multiply(tstp,unity,2.0);
			Aview.smul(2.0);
			// err+=distance_norm2(tstp,Aview,A2);
		}
		REQUIRE(err<eps);

	}
	// Check MPI broadcast
	#if CNTR_USE_MPI == 1
	// SECTION ("MPI - herm_matrix_timestep_view"){
	// 	MPI::Init(argc,argv);
    	// int numtasks=MPI::COMM_WORLD.Get_size();
    	// int taskid=MPI::COMM_WORLD.Get_rank();
    	// int mpi_imp=0;
    	// double err=0.0;
    	// GREEN A2(nt,ntau,size);
	// 	A2=A;
    	// for(int tstp=-1;tstp<nt;tstp++){
    	// 	GREEN_TSTP_VIEW Aview(tstp,A);
     	// 	Aview.MPI_reduce(mpi_imp);
     	// 	if(mpi_imp==taskid){
     	// 		A2.right_multiply(tstp,unity,numtasks);
     	// 		err+=distance_norm2(tstp,Aview,A2);
     	// 	}
    	// }
    	// REQUIRE(err<eps);
	// }
	#endif

	#if CNTR_USE_HDF5 == 1
	SECTION ("HDF5 - herm_matrix_timestep_view"){
		// double err=0.0;
		// for(int tstp=-1;tstp<nt;tstp++){
  //     		hid_t file_id = open_hdf5_file("herm_member_timestep_view.h5");
  //     		GREEN_TSTP_VIEW Aview(tstp,A);
  //     		GREEN_TSTP_VIEW Aview2(tstp,ntau,size,size,-1);

  //     		Aview.write_to_hdf5(file_id,"Gtstp");
  //     		close_hdf5_file(file_id);
  //     		// TODO: can't read hdf5 with herm_matrix_timestepview - should this even be possible ?
  //     		// Aview2.read_from_hdf5("herm_member_timestep_view.h5","Gtstp");
  //     		// close_hdf5_file(file_id);

  //     		err+=distance_norm2(tstp,Aview2,A);
  //     	}
  //     	REQUIRE(err<eps);
	}
	#endif
}
