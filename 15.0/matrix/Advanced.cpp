/*
 * @Description: Advanced.cpp
 * @version: 10.0
 * @Author: shinonome
 * @Date: 2021-03-10 20:00:21
 * @LastEditor: shinonome
 * @LastEditTime: 2021-03-11 19:39:17
 */

#include "classdef.h"

template <class _Tp, size_t _Size_m, size_t _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::Gaussian_elimination() {
	tp i, j, k, l;
	for (j = 0; j < _Size_m; j++) {
		l = 0;
		for (i = j; i < _Size_m; i++)
			for (k = 0; k < _Size_n; k++)
				if (!Matrix[i][k])
					l++;
		if (l == (_Size_m - j) * _Size_n)
			return *this;
		k = j;
		do {
			l = 0;
			if (!Matrix[j][k]) {
				for (i = j; i < _Size_m; i++) {
					if (Matrix[i][k]) {
						Row_Tran2(j, i);
						break;
					} else
						l++;
				}
			}
			if (l == _Size_m - j && l)
				k++;
		} while (l == _Size_m - j && l);
		Row_Tran1(j, 1 / Matrix[j][k]);
		for (i = j + 1; i < _Size_m; i++)
			Row_Tran3(j, -Matrix[i][k], i);
	}
	return *this;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::I_Tran() {
	this->Gaussian_elimination();
	for (tp i = _Size_m - 1; i > 0; i--) {
		for (tp j = i; j < _Size_n; j++) {
			if (Matrix[i][j]) {
				for (tp k = 0; k < i; k++)
					Row_Tran3(i, -Matrix[k][j], k);
				break;
			} else
				continue;
		}
	}
	return *this;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
Linear_algebra::tp
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::_rank_() const {
	tp i, j, k;
	matrix<_Tp, _Size_m, _Size_n> *Temp = new matrix<_Tp, _Size_m, _Size_n>;
	*Temp = *this;
	Temp->Gaussian_elimination();
	for (i = _Size_m - 1; i >= 0; i--) {
		k = 0;
		for (j = _Size_n - 1; j >= i; j--) {
			if (!Temp->get_element_by_index(i, j))
				k++;
		}
		if (k != _Size_n - i)
			return i + 1;
	}
	delete Temp;
	return 0;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>
Linear_algebra::matrix<_Tp, _Size_m, _Size_n>::Inverse() {
	assert(_Squ_ && rank == _Size_m);
	matrix<_Tp, _Size_m, 2 *_Size_n> *Temp =
		new matrix<_Tp, _Size_m, 2 * _Size_n>;
	for (tp i = 0; i < _Size_m; i++)
		for (tp j = 0; j < _Size_n; j++)
			Temp->modify_element_by_index(i, j, Matrix[i][j]);
	for (tp i = 0; i < _Size_m; i++)
		Temp->modify_element_by_index(i, _Size_n + i, 1);
	Temp->I_Tran();
	for (tp i = 0; i < _Size_m; i++) {
		for (tp j = 0; j < _Size_n; j++)
			Matrix[i][j] = Temp->get_element_by_index(i, _Size_n + j);
	}
	delete Temp;
	return *this;
}
