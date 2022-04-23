/*
 * @Description: Other_Tran.cpp
 * @version: 10.0
 * @Author: shinonome
 * @Date: 2021-03-11 12:54:24
 * @LastEditor: shinonome
 * @LastEditTime: 2021-03-11 21:00:24
 */

#include "classdef.h"

template <class _Tp, size_t _Size_m, size_t _Size_n>
Linear_algebra::matrix<_Tp, _Size_n, _Size_m>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::Transpose() const {
	Linear_algebra::matrix<_Tp, _Size_n, _Size_m> _result;
	for (tp i = 0; i < _Size_m; i++)
		for (tp j = 0; j < _Size_n; j++)
			_result.Matrix[j][i] = Matrix[i][j];
	return _result;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
template <typename T>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::Col_Tran1(const tp Col,
							 const T coef) {
	assert(Col >= 0 && Col < _Size_n);
	for (tp i = 0; i < _Size_m; i++)
		this->Matrix[i][Col] *= coef;
	return *this;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::Col_Tran2(const tp Col_1,
							 const tp Col_2) {
	assert(Col_1 >= 0 && Col_1 < _Size_n && Col_2 >= 0 && Col_2 < _Size_n);
	for (tp i = 0; i < _Size_m; i++)
		swap(Matrix[i][Col_1], Matrix[i][Col_2]);
	return *this;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
template <typename T>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::Col_Tran3(const tp Col_1,
							 const T coef,
							 const tp Col_2) {
	assert(Col_1 >= 0 && Col_1 < _Size_n && Col_2 >= 0 && Col_2 < _Size_n);
	for (tp i = 0; i < _Size_m; i++)
		Matrix[i][Col_2] += coef * Matrix[i][Col_1];
	return *this;
}
