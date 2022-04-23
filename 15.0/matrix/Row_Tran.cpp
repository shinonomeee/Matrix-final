/*
 * @Description: Row_Tran.cpp
 * @version: 10.0
 * @Author: shinonome
 * @Date: 2021-03-09 00:15:51
 * @LastEditor: shinonome
 * @LastEditTime: 2021-03-09 00:43:36
 */

#include "classdef.h"

template <class _Tp, size_t _Size_m, size_t _Size_n>
template <typename T>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::Row_Tran1(const tp Row,
							 const T coef) {
	assert(Row >= 0 && Row < _Size_m);
	Matrix[Row] * coef;
	return *this;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::Row_Tran2(const tp Row_1,
							 const tp Row_2) {
	assert(Row_1 >= 0 && Row_1 < _Size_m && Row_2 >= 0 && Row_2 < _Size_m);
	Matrix[Row_1].swap(Matrix[Row_2]);
	return *this;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
template <typename T>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::Row_Tran3(const tp Row_1,
							 const T coef,
							 const tp Row_2) {
	assert(Row_1 >= 0 && Row_1 < _Size_m && Row_2 >= 0 && Row_2 < _Size_m);
	for (tp i = 0; i < _Size_n; i++)
		Matrix[Row_2][i] += coef * Matrix[Row_1][i];
	return *this;
}
