#include "catch.hpp"
#include "cntr.hpp"
#include <iostream>
#include <fstream>
// TODO:
// Set matrixelement
// Incr
// Right/left multiply
// Interaction with herm_matrix_timestep 

#define cfunction cntr::function<double>
template <typename T>
T func_norm2(cntr::function<T> &f0, cntr::function<T> &f1)
{
	T err = 0.0;
	int nt = f0.nt_;
	assert(f0.nt_ == f1.nt_);
	int size1 = f0.size1_;
	int size2 = f0.size2_;
	int total_size = f0.total_size_;
	std::complex<T> diff;
	for(int idx = 0; idx < total_size; idx++)
	{
		diff = f0.data_[idx] - f1.data_[idx];
		err += diff.real() * diff.real() + diff.imag() * diff.imag();
	}
	err = sqrt(err);
	return err;
}

template <typename T>
T func_norm2(cntr::function<T> &f0, cdmatrix &f1)
{
	T err = 0.0;
	int nt = f0.nt_;
	int size1 = f0.size1_;
	int size2 = f0.size2_;
	assert(size1 == f1.rows());
	assert(size2 == f1.cols());
	cdmatrix val(size1,size2);
	for(int tstp=-1; tstp<=nt; tstp++) {
		f0.get_value(tstp,val);
		err += (val - f1).norm();
	}
	err = sqrt(err);
	return err;
}


template <typename T>
T func_norm2(cntr::function<T> &f0, std::complex<T> *f1)
{
	T err = 0.0;
	int nt = f0.nt_;
	for(int tstp = 0; tstp <= nt + 1; tstp++)
	{
		std::complex<T> diff = 0.0;
		diff = f0.data_[tstp] - f1[tstp];
		err += diff.real() * diff.real() +  diff.imag() * diff.imag();
	}
	err = sqrt(err);
	return err;
}

template <typename T>
T func_norm2(cntr::function<T> &f0, std::complex<T> constant)
{
	T err = 0.0;
	int nt = f0.nt_;
	for(int tstp = 0; tstp <= nt + 1; tstp++)
	{
		std::complex<T> diff = 0.0;
		diff = f0.data_[tstp] - constant;
		err += diff.real() * diff.real() +  diff.imag() * diff.imag();
	}
	err = sqrt(err);
	return err;
}




// Test get_set for cfunction
double setget(cfunction &A,cdmatrix &a){
	double toterr=0.0;
	cdmatrix tmp(A.size1_,A.size2_);
	for(int m=0;m<A.nt();m++){
		A.set_value(m,a);
		A.get_value(m,tmp);
		toterr+=(a-tmp).norm();
	}
	return toterr;
}



TEST_CASE("function","[function]"){
	int nt=50;
	int ntau=500;
	int size=2;
	int size2=5;
	double eps=1e-7;
	double eps2=1e-3; //Extrapolations within the bootstraping give worse accuracy
	double beta=5.0;
	double h=0.01;
	int kt=5;

	cfunction A(nt,size);

	SECTION ("Set/get for function"){
		// Square matrix
		cdmatrix a(size,size);
		a.setZero();
		a(0,0)=sqrt(2.0);
		a(0,1)=sqrt(2.0)*std::complex<double>(0.0,1.0);
		a(1,0)=sqrt(2.0)*std::complex<double>(0.0,-1.0);
		a(1,1)=-sqrt(2.0);

		// Rectangular matrix
		cdmatrix b(size,size2);
		b.setZero();
		for (int i=0;i<size;i++){
			for (int j=0;j<size2;j++){
				b(i,j)=i+j*2;
			}
		}


		cfunction B(nt,size,size2);
		// Writing/reading from squared function
		REQUIRE(setget(A,a)<eps);
		// Writing/reading from rectangular function
		REQUIRE(setget(B,b)<eps);
	}

	SECTION ("Extrapolation"){
		// Function
		cfunction B(nt,size);
		cdmatrix a(size,size),b(size,size),c(size,size);
		a.setZero();b.setZero();c.setZero();
		double t,toterr=0.0;
		for (int tstp=-1;tstp<=nt;tstp++){
			if(tstp==-1){
				t=0.0;
			}else{
				t=tstp*h;
			}
			c(0,0)=2.0*cos(t);
			c(0,1)=0.5*cos(t);
			c(1,0)=0.5*cos(t);
			c(1,1)=3.0*cos(t);

			A.set_value(tstp,c);
		}
		for(int tstp=0;tstp<=nt-1;tstp++){
			// Assume that function is set up to kt - first version of bootstrapping
			int n1=(tstp<=kt && tstp>=0 ? tstp : kt);
			for(int m=-1;m<=tstp;m++){
				A.get_value(m,a);
				B.set_value(m,a);
			}
			cntr::extrapolate_timestep(tstp,B,integration::I<double>(n1));
			A.get_value(tstp+1,a);
			B.get_value(tstp+1,b);
			toterr+=(a-b).norm();
		}
		REQUIRE(toterr<eps2);
	}

	SECTION ("Interpolation"){
		// Function
		cfunction B(nt,size);
		cdmatrix a0(size,size),a1(size,size),b0(size,size),b1(size,size),c(size,size);
		a0.setZero();a1.setZero();b0.setZero();b1.setZero();c.setZero();
		double t,toterr0=0.0,toterr1=0.0;
		for (int tstp=-1;tstp<=nt;tstp++){
			if(tstp==-1){
				t=0.0;
			}else{
				t=tstp*h;
			}
			c(0,0)=2.0*cos(t);
			c(0,1)=0.5*cos(t);
			c(1,0)=0.5*cos(t);
			c(1,1)=3.0*cos(t);

			A.set_value(tstp,c);
		}

		for(int tstp=1;tstp<nt;tstp++){
			// Assume that function is set up to kt - first version of bootstrapping
			int n1=(tstp<=kt && tstp>=0 ? tstp : kt);
			for(int m=-1;m<=tstp;m++){
				A.get_value(m,c);
				B.set_value(m,c);
			}
			// Test at time
			b0=cntr::interpolation(tstp,tstp,B,integration::I<double>(n1));
			A.get_value(tstp,a0);

			b1=cntr::interpolation(tstp,(tstp-0.5),B,integration::I<double>(n1));
			a1(0,0)=2.0*cos((tstp-0.5)*h);
			a1(0,1)=0.5*cos((tstp-0.5)*h);
			a1(1,0)=0.5*cos((tstp-0.5)*h);
			a1(1,1)=3.0*cos((tstp-0.5)*h);

			// std::cout << "---------------------------- " << std::endl;
			// std::cout << "tstp " << tstp << std::endl;
			// std::cout <<  a0 << std::endl;
			// std::cout <<  a1 << std::endl;

			// std::cout <<  b0 << std::endl;
			// std::cout <<  b1 << std::endl;

			// std::cout <<(a0-b0).norm() << std::endl;
			// std::cout <<(a1-b1).norm() << std::endl;

			// std::cout << "**************************** " << std::endl;

			// Test at time+dt/2
			toterr0+=(a0-b0).norm();
			toterr1+=(a1-b1).norm();
		}

		REQUIRE(toterr0<eps);
		REQUIRE(toterr1<eps2); //Problem is that for first steps low order interpolation makes an error of 10^-5
	}
	SECTION("new methods in function") {
		int nt = 1000;
		double err = 0.0;
		double eps = 1e-5;
		std::cout << "test incr and smul" << std::endl;
		{//incr and smul
			int size1 = 3;
			int size2 = 5;
			double alpha = 0.5;
			cntr::function<double> f0(nt, size1, size2), f1(nt, size1, size2);
			cdmatrix mat0(size1, size2), mat1(size1, size2);
			cdmatrix mat3(size1, size2);
			mat0 << 1.0, 2.0, 0.0, 4.0, 3.0, 2.1, 23, 12, 0.1, 0.0, 13, 231, 20.1, 2.1, 0.44;
			mat1 << -1.0, 0.5, 0.7, -2.1, 13.0, 2.9,-3,2,6.1,1.7, -3,-1,2.1,-2.1,8.44;

			mat3 = mat0 + alpha * mat1;
			f0.set_constant(mat0);
			f1.set_constant(mat1);
			f0.incr(f1, alpha);
			err = func_norm2(f0, mat3);
			REQUIRE(err < eps);

			mat3 = alpha * mat3;
			f0.smul(alpha);
			err = func_norm2(f0, mat3);
			REQUIRE(err < eps);
		}
		//std::cout << "testing multiplication" << std::endl;
		std::cout << "mult" << std::endl;
		{//multiply
			int size = 5;
			cntr::function<double> f0(nt, size), f1(nt, size);
			cdmatrix mat0(size, size), mat1(size, size);
			cdmatrix mat2;
			mat0 << 0.1,0.3,0.7,0.5,1.1,
			     0.1,0.3,0.7,0.5,1.1,
			     0.1,0.3,0.7,0.5,1.1,
			     0.1,0.3,0.7,0.5,1.1,
			     0.1,0.3,0.7,0.5,1.1;
			mat1 << 2, 7 ,4, 0.3, 2,
			     2, 7 ,4, 0.3, 2,
			     2, 7 ,4, 0.3, 2,
			     2, 7 ,4, 0.3, 2,
			     2, 7 ,4, 0.3, 2;

			mat2 = mat1 * mat0;
			f0.set_constant(mat0);
			f1.set_constant(mat1);
			f0.left_multiply(f1);
			err = func_norm2(f0, mat2);
			REQUIRE(err < eps);

			mat2 = mat0 * mat1;
			f0.set_constant(mat0);
			f1.set_constant(mat1);
			f0.right_multiply(f1);
			err = func_norm2(f0, mat2);
			REQUIRE(err < eps);

		}
		std::cout << "test matrixelements" << std::endl;
		{//matrixelements
			int size1 = 3;
			int size2 = 5;
			int i1, i2, j1, j2;
			double element1, element2;
			cntr::function<double> f0(nt, size1, size2), f1(nt, size1, size2), f2(nt, size1, size2), f3(nt, 1), f4(nt, 1);
			cdmatrix mat0(size1, size2), mat1(size1, size2), mattemp(size1, size2);
			mat0 << 1.0, 2.0, 0.0, 4.0, 3.0, 2.1, 23, 12, 0.1, 0.0, 13, 231, 20.1, 2.1, 0.44;
			mat1 << -1.0, 0.5, 0.7, -2.1, 13.0, 2.9, -3, 2, 6.1, 1.7, -3, -1, 2.1, -2.1, 8.44;
			//test
			f0.set_constant(mat0);
			f1.set_constant(mat1);
			for(i1 = 0; i1 < size1; i1++)
			{
				for(i2 = 0; i2 < size2; i2++)
				{
					for(j1 = 0; j1 < size1; j1++)
					{
						for(j2 = 0; j2 < size2; j2++)
						{
							//			mattemp = mat0;
							//std::cout << i1 << i2 << std::endl;
							//std::cout << j1 << j2 << std::endl;
							mat0(i1, i2) = mat1(j1, j2);
							f0.set_matrixelement(i1, i2, f1, j1, j2);
							err = func_norm2(f0, mat0);
							REQUIRE(err < eps);
							f2.set_constant(mat0);
							err = func_norm2(f0, f2);
							REQUIRE(err < eps);
						}
					}
				}
			}
			cdmatrix mat_1x1(1,1);
			for(i1 = 0; i1 < size1; i1++)
			{
				for(i2 = 0; i2 < size2; i2++)
				{
					f0.get_matrixelement(i1, i2, f3);
					mat_1x1(0,0) = mat0(i1, i2);
					f4.set_constant(mat_1x1);
					err = func_norm2(f3, mat_1x1);
					REQUIRE(err < eps);
					err = func_norm2(f3, f4);
					REQUIRE(err < eps);
				}
			}

		}

	}

#if CNTR_USE_HDF5 == 1
	SECTION("function hdf5 file read write") {
		int nt = 1000;
		double h = 0.01;
		double eps = 1e-10;
		{//test size = 1
			cntr::function<double> f0(nt, 1), f1;
			/*
			   double func[nt + 2];
			   for(int t = -1; t <= nt; t++)
			   {
			   double time = h * nt;
			   func[t + 1] = time * time * 0.5;
			   f0.set_value(t, func[t + 1]);
			   }
			   */
			cdmatrix constant=MatrixXcd::Zero(1,1);
			// std::complex<double> constant(0.0, 0.0);
			f1.set_constant(constant);
			f0.write_to_hdf5("test_func.h5", "data");
			f1.read_from_hdf5("test_func.h5", "data");
			double res = func_norm2(f0, f1);
			std::cout << "error is: " << res << std::endl;
			REQUIRE(res <= eps);
		}

	}
#endif

}

