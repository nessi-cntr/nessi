#ifndef READ_INPUTFILE_H_
#define READ_INPUTFILE_H_

#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <string>
#include "string.h"

#include <eigen3/Eigen/Dense>

#ifndef CPLX
  #define CPLX cdouble
#endif

/////////////////////////////////////////////////////////////////////////////////////////////
// ROUTINES FOR INTERPRETING THE INPUT FILE (LIKE F0R LAYER)
// ... interpreting input file: 
#define INPUT_LINELEN 1024
void uncomment_line(char *line,int len){
   int i;
   int len1=strlen(line);
   for(i=0;i<len1;i++){
       if(line[i]=='#'){
	       line[i]='\0';
		   break;
	   }
   }
}
// read n consecutive double numbers from file:
void read_vector_from_file(char *file,std::vector<double> &data,int n){
    int j=-1;
	std::ifstream in;
	try{
	   data.resize(n,0);
	   in.open(file,std::ios::in);
	   if(!in.good()) throw(1);
	   for(j=0;j<n;j++){
			in >> data[j];
			if(!in.good()) throw(2);
	   }
	   in.close();
	}
	catch(int e){
		if(e==1){
			std::cerr <<  __FUNCTION__ << ": cannot open file " << file  << std::endl;
		}
		if(e==2){
			std::cerr << __FUNCTION__ << ": error reading from file " << file << std::endl;
			std::cerr <<  " at entry" << j << std::endl;
		}
		throw;
	}
}
// read first line in file containing the string flag
char* find_line_with_flag(char *file,const char *flag,char *line){
	std::ifstream in;
	char *pch=0;
	try{
		in.exceptions(std::ifstream::eofbit |std::ifstream::badbit | std::ifstream::failbit);
		in.open(file,std::ios::in);
		do{
			in.getline(line,INPUT_LINELEN);
			uncomment_line(line,INPUT_LINELEN);
			if((pch=strstr(line,flag))) break;
		}while(pch==NULL);	
		in.close();		
	}
	catch(...){
		std::cerr << "some error occurred in file " << file << std::endl;
		std::cerr << "while looking for line with " << flag << std::endl;
		throw;
	}
	return pch;
}
// split line into n words, separated bt characters in delim
// delim-characters in line are replaced by \0
int strtok(char *line,std::vector<char*> &words,const char *delim){
	char *pch;
    int n=0;
	words.resize(0);
	pch=strtok(line,delim);
	while(pch!=NULL){
	    n++;
		words.resize(n);
		words[n-1]=pch;
		pch=strtok(NULL,delim);
	}
	return n;
}
// find line containing consecutive fields "flag" "i1",
// and split rest of the line in words[0] ... words[n2-1]
template<typename T>
int find_line_with_fields(char *file,char *line,std::vector<char*> &words,const char *flag,T i1){
	std::ifstream in;
	char *pch=0;
	int n,l,n2=0;
	T i2;
	std::vector<char*> words1;
	try{
		in.exceptions(std::ifstream::eofbit |std::ifstream::badbit | std::ifstream::failbit);
		in.open(file,std::ios::in);
		do{
			in.getline(line,INPUT_LINELEN);
			uncomment_line(line,INPUT_LINELEN);
			if((pch=strstr(line,flag))){
				//std::cout << "1) found line matching line: " << line << std::endl;
				n=strtok(pch,words1," \t");
				//std::cout << "2) found " << n << " words in line" << std::endl;
				if(n>=2){
					std::stringstream(words1[1]) >> i2;
					//std::cout << "3) word1[1]= " << i2 << std::endl;
					if(i2==i1){
					    // found the line I am looking for:
						n2=n-2;
						words.resize(n2);
						for(l=0;l<n2;l++){
							words[l]=words1[l+2];
							//std::cout << "... word["<<l<<"]= " << words[l] << std::endl;
						}
						//std::cout << "4) return " << n2 << std::endl;
						break;
					}
				}
			}
		}while(1);	
		in.close();		
	}
	catch(...){
		std::cerr << "some error occurred in file " << file << std::endl;
		std::cerr << "while looking for line with two consecutive fields ... " << flag << " " << i1 << std::endl;
		throw;
	}
	return n2;
}
template<typename T1,typename T2>
int find_line_with_fields(char *file,char *line,std::vector<char*> &words,const char *flag,T1 i1,T2 j1){
	std::ifstream in;
	char *pch=0;
	int n,l,n2=0;
	T1 i2;
	T2 j2;
	std::vector<char*> words1;
	try{
		in.exceptions(std::ifstream::eofbit |std::ifstream::badbit | std::ifstream::failbit);
		in.open(file,std::ios::in);
		do{
			in.getline(line,INPUT_LINELEN);
			uncomment_line(line,INPUT_LINELEN);
			if((pch=strstr(line,flag))){
				//std::cout << "1) found line matching line: " << line << std::endl;
				n=strtok(pch,words1," \t");
				//std::cout << "2) found " << n << " words in line" << std::endl;
				if(n>=3){
					std::stringstream(words1[1]) >> i2;
					std::stringstream(words1[2]) >> j2;
					//std::cout << "3) word1[1]= " << i2 << std::endl;
					if(i2==i1 && j2==j1){
					    // found the line I am looking for:
						n2=n-3;
						words.resize(n2);
						for(l=0;l<n2;l++){
							words[l]=words1[l+3];
							//std::cout << "... word["<<l<<"]= " << words[l] << std::endl;
						}
						//std::cout << "4) return " << n2 << std::endl;
						break;
					}
				}
			}
		}while(1);	
		in.close();		
	}
	catch(...){
		std::cerr << "some error occurred in file " << file << std::endl;
		std::cerr << "while looking for line with three consecutive fields ... " << flag << " " << i1 << " " << j1 << std::endl;
		throw;
	}
	return n2;
}




// READ A PARAMETER FROM A SEQUENCE OF WORDS words[m],words[m+1],...
void read_param(std::vector<char*> &words,int m,double &data){
	try{
		if((int)words.size()<m+1) throw("read_param: too few words in line");			
		if(!(std::stringstream(words[m]) >> data)) 
			throw("read_param_tvector: cannot interpret arg");
	}
	catch(...){
		std::cerr << "error in read_param(double) " << std::endl;
		throw;
	}
}
void read_param(std::vector<char*> &words,int m,int &data){
	try{
		if((int)words.size()<m+1) throw("read_param: too few words in line");			
		if(!(std::stringstream(words[m]) >> data)) 
			throw("read_param_tvector: cannot interpret arg");
	}
	catch(...){
		std::cerr << "error in read_param(int) " << std::endl;
		throw;
	}
}

void read_param(std::vector<char*> &words,int m,bool &data){
	try{			
	  if(!(std::stringstream(words[m]) >> std::boolalpha >> data)) 
			throw("read_param_tvector: cannot interpret arg");
	}
	catch(...){
		std::cerr << "error in read_param(int) " << std::endl;
		throw;
	}
}

void read_param(std::vector<char*> &words,int m,CPLX &data){
  double real, imag;
      try{
		if((int)words.size()<m+2) throw("read_param: too few words in line");			
		//if(!(std::stringstream(words[m]) >> data.real())) 
		if(!(std::stringstream(words[m]) >> real)) 
			throw("read_param_tvector: cannot interpret arg");
		//if(!(std::stringstream(words[m+1]) >> data.imag())) 
		if(!(std::stringstream(words[m+1]) >> imag)) 
			throw("read_param_tvector: cannot interpret arg");
	}
	catch(...){
		std::cerr << "error in read_param(CPLX) " << std::endl;
		throw;
	}
      data.real(real);
      data.imag(imag);
}
void read_param(std::vector<char*> &words,int m,dvector &data,int dim){
	try{
	    data=dvector(dim);
		if((int)words.size()<m+dim) throw("read_param: too few words in line");			
		for(int i1=0;i1<dim;i1++){
			if(!(std::stringstream(words[m+i1]) >> data(i1))) 
				throw("read_param_tvector: cannot interpret arg");
		}
	}
	catch(...){
		std::cerr << "error in read_param(dvector) " << std::endl;
		throw;
	}
}

void read_param(std::vector<char*> &words,int m,std::vector<double> &data,int dim){
	try{
	    // std::vector<double> data(dim);
		if((int)words.size()<m+dim) throw("read_param: too few words in line");			
		for(int i1=0;i1<dim;i1++){
			if(!(std::stringstream(words[m+i1]) >> data[i1])) 
				throw("read_param_tvector: cannot interpret arg");
		}
	}
	catch(...){
		std::cerr << "error in read_param(dvector) " << std::endl;
		throw;
	}
}

// fine the first line containing "flag", and read the following number
template<typename T> void find_param(char *file,const char *flag,T &x){
	try{
		char line[INPUT_LINELEN],*pch;
		std::vector<char*> words;
		find_line_with_flag(file,flag,line);
		pch=strstr(line,flag);
		strtok(pch+strlen(flag),words," \t");
		read_param(words,0,x);
	}
	catch(...){
		std::cerr << "error in find_param(x) " << std::endl;
		throw;
	}
}

void find_param(char *file,const char *flag,dvector &x,int dim){
	try{

		char line[INPUT_LINELEN],*pch;
		std::vector<char*> words;
		find_line_with_flag(file,flag,line);
		pch=strstr(line,flag);
		strtok(pch+strlen(flag),words," \t");
		read_param(words,0,x,dim);
	}
	catch(...){
		std::cerr << "error in find_param(x) " << std::endl;
		throw;
	}
}

void find_param(char *file,const char *flag,std::vector<double> &x,int dim){
	try{
		x.resize(dim);
		char line[INPUT_LINELEN],*pch;
		std::vector<char*> words;
		find_line_with_flag(file,flag,line);
		pch=strstr(line,flag);
		strtok(pch+strlen(flag),words," \t");
		read_param(words,0,x,dim);
	}
	catch(...){
		std::cerr << "error in find_param(x) " << std::endl;
		throw;
	}
}

// READ a time-dependent function A SEQUENCE OF WORDS words[m],words[m+1],...
// if words[m] is of the form --FILENAME, read funtion from file
// otherwise, a time-independent value is read from words[m], words[m+1], ...
void read_param_tvector(std::vector<char*> &words,int m,std::vector<dvector> &data,int dim,int nt){
	char *pch;
	int i1,tstp;
	try{
	    data.resize(nt+2);
		for(tstp=-1;tstp<=nt;tstp++) data[tstp+1]=dvector(dim);
		// test whether the first word is a file
		if((int)words.size()<=m) throw("read_param_tvector: too few words in line");
		if((pch=strstr(words[m],"--"))){
			std::vector<double> tmp;
			read_vector_from_file(pch+2,tmp,(nt+2)*dim); // read as many numbers
			for(tstp=-1;tstp<=nt;tstp++){ 
				for(i1=0;i1<dim;i1++){ 
					data[tstp+1](i1)=tmp[(tstp+1)*dim+i1];
				}
			}
		}else{
			dvector tmp(dim);
			read_param(words,m,tmp,dim);
			for(tstp=-1;tstp<=nt;tstp++) data[tstp+1]=tmp;
		}
	}
	catch(...){
		std::cerr << "error in read_param_tvector(dvector) " << std::endl;
		throw;
	}
}
void read_param_tvector(std::vector<char*> &words,int m,std::vector<CPLX> &data,int nt){
	char *pch;
	int tstp;
	try{
	    data.resize(nt+2);
		if((int)words.size()<=m) throw("read_param_tvector: too few words in line");
		// test whether the first word is a file
		if((pch=strstr(words[m],"--"))){
			std::vector<double> tmp;
			read_vector_from_file(pch+2,tmp,(nt+2)*2); // read as many numbers
			for(tstp=-1;tstp<=nt;tstp++){ 
				data[tstp+1]=CPLX(tmp[(tstp+1)*2],tmp[(tstp+1)*2+1]);
			}
		}else{
			CPLX tmp;
			read_param(words,m,tmp);
			for(tstp=-1;tstp<=nt;tstp++) data[tstp+1]=tmp;
		}
	}
	catch(...){
		std::cerr << "error in read_param_tvector(CPLX) " << std::endl;
		throw;
	}
}
void read_param_tvector(std::vector<char*> &words,int m,std::vector<double> &data,int nt){
	char *pch;
	int tstp;
	try{
	    data.resize(nt+2);
		if((int)words.size()<=m) throw("read_param_tvector: too few words in line");
		// test whether the first word is a file
		if((pch=strstr(words[m],"--"))){
			std::vector<double> tmp;
			read_vector_from_file(pch+2,tmp,nt+2); // read as many numbers
			for(tstp=-1;tstp<=nt;tstp++){ 
				data[tstp+1]=tmp[tstp+1];
			}
		}else{
			double tmp;
			read_param(words,m,tmp);
			for(tstp=-1;tstp<=nt;tstp++) data[tstp+1]=tmp;
		}
	}
	catch(...){
		std::cerr << "error in read_param_tvector(double) " << std::endl;
		throw;
	}
}
///
void read_param_tvector(char *line,std::vector<dvector> &data,int dim,int nt){
	std::vector<char*> words;
    strtok(line,words," \t");
	read_param_tvector(words,0,data,dim,nt);
}
void read_param_tvector(char *line,std::vector<CPLX> &data,int nt){
	std::vector<char*> words;
    strtok(line,words," \t");
	read_param_tvector(words,0,data,nt);
}
void read_param_tvector(char *line,std::vector<double> &data,int nt){
	std::vector<char*> words;
    strtok(line,words," \t");
	read_param_tvector(words,0,data,nt);
}


void find_param_tvector(char *file,const char *flag,std::vector<dvector> &data,int dim,int nt){
	try{
		char line[INPUT_LINELEN],*pch;
		find_line_with_flag(file,flag,line);
		pch=strstr(line,flag);
		read_param_tvector(pch+strlen(flag),data,dim,nt);
	}
	catch(...){
		std::cerr << "error in find_param(x) " << std::endl;
		throw;
	}
}
void find_param_tvector(char *file,const char *flag,std::vector<CPLX> &data,int nt){
	try{
		char line[INPUT_LINELEN],*pch;
		find_line_with_flag(file,flag,line);
		pch=strstr(line,flag);
		read_param_tvector(pch+strlen(flag),data,nt);
	}
	catch(...){
		std::cerr << "error in find_param(x) " << std::endl;
		throw;
	}
}
void find_param_tvector(char *file,const char *flag,std::vector<double> &data,int nt){
	try{
		char line[INPUT_LINELEN],*pch;
		find_line_with_flag(file,flag,line);
		pch=strstr(line,flag);
		read_param_tvector(pch+strlen(flag),data,nt);
	}
	catch(...){
		std::cerr << "error in find_param(x) " << std::endl;
		throw;
	}
}

#undef CPLX

#endif


