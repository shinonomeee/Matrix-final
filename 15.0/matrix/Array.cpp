/*
 * @Description: 
 * @version: 
 * @Author: shinonome
 * @Date: 2021-09-24 20:13:39
 * @LastEditor: shinonome
 * @LastEditTime: 2022-04-23 18:44:56
 */

#include "classdef.h"

template <class _Tp, size_t _Size>
_Tp Linear_algebra::array<_Tp, _Size>::dot_product(
	const array<_Tp, _Size> &Value) const {
	_Tp _Res = 0;
	for (tp i = 0; i < _Size; i++)
		_Res += (*this)[i] * Value[i];
	return _Res;
}

template <class _Tp, size_t _Size>
_Tp Linear_algebra::array<_Tp, _Size>::Euclidean_norm() const {
	return pow(this->dot_product(*this), 0.5);
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
_Tp Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::self_dot_product(
	const tp Row) const {
	return Matrix[Row].dot_product(Matrix[Row]);
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
_Tp Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::dot_product(
	const tp Row,
	const array<_Tp, _Size_n> &Value) const {
	return Matrix[Row].dot_product(Value);
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
void Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::schmidt_orthogonalization() { }

template <class _Tp, size_t _Size_m, size_t _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::Schmidt_orthogonalization()
	const {
	matrix<_Tp, _Size_m, _Size_n> Temp;
	_Tp _temp = 0;
	for (tp i = 0; i < _Size_n; i++)
		Temp.modify_element_by_index(0, i, this->Matrix[0][i]);
	for (tp i = 1; i < _Size_m; i++) {
		for (tp j = 0; j < _Size_n; j++) {
			_temp = Temp.get_element_by_index(0, j);
			std::cout << _temp << std::endl;
			for (tp k = 0; k < i; k++) {
				_temp -=
					(Temp.dot_product(k, this->Matrix[i]) / Temp.self_dot_product(k)) *
					Temp.get_element_by_index(k, j);
			}
			std::cout << _temp << std::endl;
			Temp.modify_element_by_index(i, j, _temp);
		}
	}

	return Temp;
}
