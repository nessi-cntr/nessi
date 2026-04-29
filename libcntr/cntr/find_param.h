#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<sstream>
#include<boost/algorithm/string.hpp>

/// @private
/** \brief <b> Utility to split the input file lines into parameter name and value.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > Take a string specifying a parameter in a given format (input file: `__parametername==parametervalue`) and save the name and value as a vector of strings separately.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param & pv
* > [std::vector<std::string>] Vector of strings containing two elements, the name and value of the parameter
* @param line
* > [std::string] One line of the parameter inout file, which ios to be split.
*/
void split(std::vector<std::string> & pv, std::string line)
{
	line.erase(std::remove_if(line.begin(), line.end(),
				[](char &c) {
				return std::isspace<char>(c, std::locale::classic());
				}),
			line.end());

	pv.resize(2);
	auto peq = line.find("=");
	pv[0] = line.substr(0, peq);
	pv[1] = line.substr(peq + 1);
}

/** \brief <b> Function for reading out a double type parameter of a given input file.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > Parse an input file of given name (e.g. given by argv[1]) for a parameter (given in the input file in the format `__parametername==parametervalue`) and writing the value into a double.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param file_name
* > [std::string] Name of the input file.
* @param param
* > [std::string] Name of the parameter in the input file, in the format `__parametername`.
* @param &value
* > [double] Double to store the parameter value in.
*/
int find_param(std::string file_name, std::string param, double &value)
{
	std::string line, data;
	std::vector <std::string> pv;
	std::ifstream file(file_name.c_str());
	while(std::getline(file, line))
	{
		if (line.find_first_not_of("\t\n ") != std::string::npos)
		{
			split(pv, line);
			if (pv[0] == param) 
			{
				value = stod(pv[1]);
				return 0;
			}
		}
	}
	return 1;
}

/** \brief <b> Function for reading out an integer type parameter of a given input file.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > Parse an input file of given name (e.g. given by argv[1]) for a parameter (given in the input file in the format `__parametername==parametervalue`) and writing the value into an integer.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param file_name
* > [std::string] Name of the input file.
* @param param
* > [std::string] Name of the parameter in the input file, in the format `__parametername`.
* @param &value
* > [int] Integer to store the parameter value in.
*/
int find_param(std::string file_name, std::string param, int &value)
{
	std::string line, data;
	std::vector <std::string> pv;
	std::ifstream file(file_name.c_str());
	while(std::getline(file, line))
	{
		if (line.find_first_not_of("\t\n ") != std::string::npos)
		{
			split(pv, line);
			if (pv[0] == param) 
			{
				value = stoi(pv[1]);
				return 0;
			}
		}
	}
return 1;
}

/** \brief <b> Function for reading out a string type parameter of a given input file.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > Parse an input file of given name (e.g. given by argv[1]) for a parameter (given in the input file in the format `__parametername==parametervalue`) and writing the value into a string.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param file_name
* > [std::string] Name of the input file.
* @param param
* > [std::string] Name of the parameter in the input file, in the format `__parametername`.
* @param &value
* > [std::string] String to store the parameter value in.
*/
int find_param(std::string file_name, std::string param, std::string &value)
{
	std::string line, data;
	std::vector <std::string> pv;
	std::ifstream file(file_name.c_str());
	while(std::getline(file, line))
	{
		if (line.find_first_not_of("\t\n ") != std::string::npos)
		{
			split(pv, line);
			if (pv[0] == param) 
			{
				value = pv[1];
				return 0;
			}
		}
	}
return 1;
}
