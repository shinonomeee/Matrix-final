/*
 * @Description: Basic.cpp
 * @version: 10.0
 * @Author: shinonome
 * @Date: 2021-03-08 23:34:13
 * @LastEditor: shinonome
 * @LastEditTime: 2021-03-09 18:50:00
 */

#include "classdef.h"

template <class _Tp, size_t _Size_m, size_t _Size_n>
void Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::_init_() {
	for (tp i = 0; i < _Size_m; i++)
		for (tp j = 0; j < _Size_n; j++)
			Matrix[i][j] = 0;
	rank = 0;
	_Squ_ = (_Size_m == _Size_n);
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
void Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::_init_(_Tp data) {
	fill(data);
	_Squ_ = (_Size_m == _Size_n);
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
template <typename T>
void Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::fill(T data) {
	for (tp i = 0; i < _Size_m; i++)
		Matrix[i].fill(data);
	if (data)
		rank = 1;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
void Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::swap(
	matrix<_Tp, _Size_m, _Size_n> &Value) {
	matrix<_Tp, _Size_m, _Size_n> Temp{ Value };
	Value = *this;
	*this = Temp;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
template <typename T>
void Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::swap(T &Value_1,
							 T &Value_2) {
	T temp{ Value_1 };
	Value_1 = Value_2;
	Value_2 = temp;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
_Tp Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::get_element_by_index(
	const tp pos_1,
	const tp pos_2) const {
	assert(pos_1 >= 0 && pos_1 < _Size_m && pos_2 >= 0 && pos_2 < _Size_n);
	return Matrix[pos_1][pos_2];
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
void Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::modify_element_by_index(
	const tp pos_1,
	const tp pos_2,
	const _Tp coef) {
	assert(pos_1 >= 0 && pos_1 < _Size_m && pos_2 >= 0 && pos_2 < _Size_n);
	Matrix[pos_1][pos_2] = coef;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
template <typename T>
T Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::max(const T &Value_1,
						     const T &Value_2) {
	if (Value_1 > Value_2)
		return Value_1;
	else
		return Value_2;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
template <typename T>
T Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::min(const T &Value_1,
						     const T &Value_2) {
	if (Value_1 > Value_2)
		return Value_2;
	else
		return Value_1;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
_Tp Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::max() {
	_Tp MAX = Matrix[0][0];
	for (tp i = 0; i < _Size_m; i++)
		for (tp j = 0; j < _Size_n; j++)
			if (MAX < Matrix[i][j])
				MAX = Matrix[i][j];
	return MAX;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
_Tp Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::min() {
	_Tp MIN = Matrix[0][0];
	for (tp i = 0; i < _Size_m; i++)
		for (tp j = 0; j < _Size_n; j++)
			if (MIN > Matrix[i][j])
				MIN = Matrix[i][j];
	return MIN;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
_Tp Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::Trace() {
	assert(_Squ_);
	_Tp sum = 0;
	for (tp i = 0; i < min(_Size_m, _Size_n); i++)
		sum += Matrix[i][i];
	return sum;
}
