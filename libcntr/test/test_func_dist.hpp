#include "cntr.hpp"
#include <cmath>
#include <complex>

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

template <typename T, class EigenMatrix>
T func_norm2(cntr::function<T> &f0, EigenMatrix &constant)
{
	T err = 0.0;
	int nt = f0.nt_;
	int size1 = f0.size1_;
	int size2 = f0.size2_;
	for(int tstp = -1; tstp <= nt; tstp++)
	{
		std::complex<T> diff = 0.0;
		for(int i = 0; i < size1; i++)
		{
			for(int j = 0; j < size2; j++)
			{
				diff = f0.ptr(tstp)[i * size2 + j] - constant(i,j);
//				std::cout << i << j << '\t' << f0.ptr(tstp)[i * size2 + j] << constant(i,j) << std::endl;
				err += diff.real() * diff.real() +  diff.imag() * diff.imag();
			}
		}
	}
	err = sqrt(err);
	return err;
}
