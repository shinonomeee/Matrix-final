/*
 * @Description: Det.cpp
 * @version: 11.0
 * @Author: shinonome
 * @Date: 2021-03-11 21:28:32
 * @LastEditor: shinonome
 * @LastEditTime: 2021-07-25 19:10:00
 */

#include <cmath>
#include "classdef.h"
#define FILE_PATH \
	"/Users/chris/Desktop/线代Project/C++/15.0/15.0/matrix/tempfile/temp.txt"

template <class _Tp, size_t _Size_m, size_t _Size_n>
_Tp Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::det() {
	assert(_Squ_);
	if (rank < _Size_m)
		return 0;
	std::fstream finout(FILE_PATH, std::ios::in | std::ios::out | std::ios::trunc | std::ios::binary);
	if (!finout.is_open()) {
		std::cerr << "Something bad happened when trying to open the file!"
			  << std::endl;
		exit(EXIT_FAILURE);
	}
	return _Det(finout);
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
_Tp Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::_Det(std::fstream &finout) {
	array<int, _Size_m> arr, source;
	char ch;
	_Tp _Res = 0, temp = 1;
	for (int i = 0; i < _Size_m; i++)
		source[i] = i + 1;
	Full_Permutation(source, 0, _Size_m, finout);
	finout.seekg(std::ios::beg).clear();
	while (!finout.eof()) {
		for (tp i = 0; i < _Size_m; ++i) {
			finout >> ch;
			temp *= Matrix[i][ch - 49];
			arr[i] = ch - 48;
		}
		if (finout.eof())
			break;
		/*
     * std::cout << "arr = ";
     * for (auto i: arr)   std::cout << i << " ";
     * std::cout << std::endl;
     */
		_Res += pow(-1, Inverted_sequence(arr)) * temp;
		temp = 1;
	}
	finout.close();
	return _Res;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
int Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::Inverted_sequence(
	const array<int, _Size_m> &arr) {
	int _result = 0;
	for (tp i = 0; i < _Size_m; i++)
		for (tp j = i; j < _Size_m; j++)
			if (arr[i] > arr[j])
				++_result;
	return _result;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
void Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::Full_Permutation(
	array<int, _Size_m> &source,
	int begin,
	int end,
	std::fstream &finout) {
	int i;
	if (begin >= end)
		finout << source;
	else {
		for (i = begin; i < end; i++) {
			if (begin != i)
				swap(source[begin], source[i]);
			Full_Permutation(source, begin + 1, end, finout);
			if (begin != i)
				swap(source[begin], source[i]);
		}
	}
}
