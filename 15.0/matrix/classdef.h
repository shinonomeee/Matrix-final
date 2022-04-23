/*
 * @Description: classdef.h
 * @version: 10.0
 * @Author: shinonome
 * @Date: 2021-03-08 22:04:26
 * @LastEditor: shinonome
 * @LastEditTime: 2021-03-11 19:39:31
 */

#ifndef Basedef_h
#define Basedef_h

#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

/* 宏定义了输出 matrix 和 array 的精度(PRECISION)和间隔(INTERVAL) */
#define PRECISION 2
#define INTERVAL 6

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * */

/* Linear_algebra */
/* namespace Linear_algebra begin */
namespace Linear_algebra {
typedef size_t tp;

/*
     以下定义了 Linear_algebra::array 类, 区别于 std::array
   在类中重载了基础运算符以及输入输出运算符
     * 注意：
     * 在 array 类中重载的运算符 + - * / 只对 Linear_algebra::array 生效
     * 而在类外声明的运算符 << >> 则对 std::array 和 Linear_algebra::array
   都生效
     * Linear_algebra::array 类不需要额外重载 = 运算符，由于 std::array
   类中已经重载了 = 运算符
     * 此外， Linear_algebra::array 也不需要额外重载 == != 运算符，由于
   std::array 类中已经重载了这两个运算符
     */

/* 在命名空间中定义模板类 */
template <class _Tp, size_t _Size>
struct array: public std::array<_Tp, _Size> {
	/* 向量点积 */
	_Tp dot_product(const array<_Tp, _Size> &Value) const;

	/* 求欧式范数 */
	_Tp Euclidean_norm() const;

	/* 含参构造函数 */
	// array(_Tp data) { this->fill(data); }

	/* 以下重载了两个 array 间的加减法运算 + - */
	array<_Tp, _Size> &operator+(const array<_Tp, _Size> &Value) {
		for (tp i = 0; i < _Size; i++)
			(*this)[i] += Value[i];
		return *this;
	}

	array<_Tp, _Size> &operator-(const array<_Tp, _Size> &Value) {
		for (tp i = 0; i < _Size; i++)
			(*this)[i] -= Value[i];
		return *this;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   * * * * * */

	/*
   * template <typename T>
   * _Tp &operator[](T _pos)
   * {
   *     assert(_pos >= 0 && _pos < _Size);
   *     return (*this)[_pos];
   * }
   */

	/* 以下重载了 array 对数乘和数除的运算符 */
	template <typename T> array<_Tp, _Size> &operator*(const T &Value) {
		for (tp i = 0; i < _Size; i++)
			(*this)[i] *= Value;
		return *this;
	}

	template <typename T> array<_Tp, _Size> &operator/(const T &Value) {
		for (tp i = 0; i < _Size; i++)
			(*this)[i] /= Value;
		return *this;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   * * * * * */
};
/* Linear_algebra::array 类定义结束 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * */

/* 重载了 输入 (istream) 和 输出 (ostream) 运算符 >> << */
template <class _Tp, size_t _Size>
std::ostream &operator<<(std::ostream &os,
			 const std::array<_Tp, _Size> &Value) {
	os << std::fixed << std::setprecision(PRECISION) << Value[0];
	for (tp i = 1; i < _Size; i++)
		os << std::fixed << std::setprecision(PRECISION) << std::setw(INTERVAL)
		   << Value[i];
	os << std::endl;
	return os;
}

template <class _Tp, size_t _Size>
std::istream &operator>>(std::istream &is, std::array<_Tp, _Size> &Value) {
	for (tp i = 0; i < _Size; i++)
		is >> Value[i];
	return is;
}

template <class _Tp, size_t _Size>
std::fstream &operator<<(std::fstream &fs,
			 const std::array<_Tp, _Size> &Value) {
	for (auto &ele: Value)
		fs << (char) (ele + 48);
	return fs;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * */

/*
 * 以下定义了容器类 matrix ，其中 matrix 容器为
 * Linear_algebra::array<Linear_algebra::array<_Tp, _Size_n>, _Size_m> 使用
 * array 作为容器的好处：有效避免了越界访问带来的问题 重载了矩阵的 + - * / == !=
 * 运算符，提升了其实用性
 */

template <class _Tp, size_t _Size_m, size_t _Size_n> class matrix {
	public:
	/* (_Size_m) * (_Size_n) 的二维数组 */
	typedef array<array<_Tp, _Size_n>, _Size_m> MAT;

	protected:
	MAT Matrix;
	tp rank;

	/* 初始化函数，将所有元素赋值为 0 */
	void _init_();
	void _init_(_Tp data);
	/* _Squ_ 用于判断是否为方阵 */
	bool _Squ_;

	FILE *fp;

	_Tp _Det(std::fstream &);
	int Inverted_sequence(const array<int, _Size_m> &);
	void Full_Permutation(array<int, _Size_m> &, int, int, std::fstream &);

	/* 求解矩阵的秩 */
	tp _rank_() const;

	/* 求矩阵自身的某些点积 */
	_Tp self_dot_product(const tp Row) const;
	_Tp dot_product(const tp Row, const array<_Tp, _Size_n> &Value) const;
	void schmidt_orthogonalization();

	public:
	/* 无参初始化函数 _init_() 会强制将所有元素赋值为 0 */
	matrix() {
		_init_();
	}

	/* 带参的构造函数 */
	matrix(_Tp data) {
		_init_(data);
	}

	void set_rank() {
		rank = _rank_();
	}
	tp get_rank() {
		return rank;
	}

	/*
   * template <typename T>
   * matrix(T data) { fill(data); }
   */

	/* 以下两个函数分别为：通过下标访问元素和通过下标修改元素 */
	_Tp get_element_by_index(const tp pos_1, const tp pos_2) const;
	void modify_element_by_index(const tp pos_1, const tp pos_2,
				     const _Tp coef);// pos_1, pos_2 为矩阵的下标 - 1

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   * * * * * */

	/* 转置矩阵，注意返回的是局部变量 */
	matrix<_Tp, _Size_n, _Size_m> Transpose() const;

	/*
  * 函数说明：
  * 在矩阵中，第 i 行第 j 列的元素一般被表示为 (i, j)
  * 而在二维数组中，其一般被表示为 array[i - 1][j - 1]
  ，类似的，在matrix容器中，其也被表示为 (i - 1, j - 1)
  *
  * template <typename T>
  * matrix<_Tp, _Size_m, _Size_n> Row_Tran1(const tp Row, const T coef); 将第
  Row 行乘系数 coef
  * matrix<_Tp, _Size_m, _Size_n> Row_Tran2(const tp Row_1, const tp Row_2);
  交换第 Row_1 - 1 行和第 Row_2 - 1 行
  * template <typename T>
  * matrix<_Tp, _Size_m, _Size_n> Row_Tran3(const tp Row_1,
                                            const T coef, tp const Row_2); 将第
  Row_1 - 1 行的 coef 倍加到第 Row_2 - 1 行
  * template <typename T>
  * matrix<_Tp, _Size_m, _Size_n> Col_Tran1(const tp Col, const T coef); 将第
  Col 列乘系数 coef
  * matrix<_Tp, _Size_m, _Size_n> Col_Tran2(const tp Col_1, const tp Col_2);
  交换第 Col_1 - 1 列和第 Col_2 - 1 列
  * template <typename T>
  * matrix<_Tp, _Size_m, _Size_n> Col_Tran3(const tp Col_1,
                                            const T coef, const tp Col_2); 将第
  Col_1 - 1 列的 coef 倍加到第 Col_2 - 1 列
  */

	template <typename T>
	matrix<_Tp, _Size_m, _Size_n> Row_Tran1(const tp Row, const T coef);
	matrix<_Tp, _Size_m, _Size_n> Row_Tran2(const tp Row_1, const tp Row_2);
	template <typename T>
	matrix<_Tp, _Size_m, _Size_n> Row_Tran3(const tp Row_1, const T coef, const tp Row_2);
	template <typename T>
	matrix<_Tp, _Size_m, _Size_n> Col_Tran1(const tp Col, const T coef);
	matrix<_Tp, _Size_m, _Size_n> Col_Tran2(const tp Col_1, const tp Col_2);
	template <typename T>
	matrix<_Tp, _Size_m, _Size_n> Col_Tran3(const tp Col_1, const T coef, const tp Col_2);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   * * * * * */

	/* 以 data填满整个矩阵 */
	template <typename T> void fill(T data);

	/*
   * 一些常用函数
   * 注意区分两个 swap 函数
   */

	void swap(matrix<_Tp, _Size_m, _Size_n> &Value);// 交换两个矩阵
	template <typename T> void swap(T &Value_1, T &Value_2);// 交换两个数
	template <typename T> T max(const T &Value_1, const T &Value_2);
	template <typename T> T min(const T &Value_1, const T &Value_2);
	_Tp max();
	_Tp min();

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   * * * * * */

	/* 求矩阵的迹 */
	_Tp Trace();

	/* 求行列式(其余三个函数放在了私有中) */
	_Tp det();

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   * * * * * */

	/*
   * 以下为最精华的算法：
   * 高斯消元
   * 单位阵变换
   */

	matrix<_Tp, _Size_m, _Size_n> Gaussian_elimination();// 高斯消元
	matrix<_Tp, _Size_m, _Size_n> I_Tran();// 将目标矩阵变为单位阵型

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   * * * * * */

	/* 求解逆矩阵 */
	matrix<_Tp, _Size_m, _Size_n> Inverse();

	/* 施密特正交化 返回一个局部变量 */
	matrix<_Tp, _Size_m, _Size_n> Schmidt_orthogonalization() const;

	/*
   * 以下重载了 matrix 类的一部分运算符
   * 基础运算符： + - * / (* / 为数乘和数除)
   * 赋值运算符： =
   */

	matrix<_Tp, _Size_m, _Size_n> &
	operator+(const matrix<_Tp, _Size_m, _Size_n> &Value) {
		for (tp i = 0; i < _Size_m; i++)
			this->Matrix[i] + Value.Matrix[i];
		return *this;
	}

	matrix<_Tp, _Size_m, _Size_n> &
	operator-(const matrix<_Tp, _Size_m, _Size_n> &Value) {
		for (tp i = 0; i < _Size_m; i++)
			this->Matrix[i] - Value.Matrix[i];
		return *this;
	}

	matrix<_Tp, _Size_m, _Size_n> &
	operator=(const matrix<_Tp, _Size_m, _Size_n> &Value) {
		for (tp i = 0; i < _Size_m; i++)
			this->Matrix[i] = Value.Matrix[i];
		return *this;
	}

	template <typename T>
	matrix<_Tp, _Size_m, _Size_n> &operator*=(const T &Value) {
		for (tp i = 0; i < _Size_m; i++)
			this->Matrix[i] * Value;
		return *this;
	}

	template <typename T>
	matrix<_Tp, _Size_m, _Size_n> &operator/(const T &Value) {
		for (tp i = 0; i < _Size_m; i++)
			this->Matrix[i] / Value;
		return *this;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   * * * * * */
};
/* matrix 类定义结束 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * */

/*
 * 以下重载了 matrix 类的 输入 (istream) 输出 (ostream) 运算符 >> <<
 * 可用于输入输出 Matrix
 */

template <class _Tp, size_t _Size_m, size_t _Size_n>
std::ostream &operator<<(std::ostream &os,
			 const matrix<_Tp, _Size_m, _Size_n> &Value) {
	for (tp i = 0; i < _Size_m; i++) {
		os << std::fixed << std::setprecision(PRECISION)
		   << Value.get_element_by_index(i, 0);
		for (tp j = 1; j < _Size_n; j++)
			os << std::fixed << std::setprecision(PRECISION) << std::setw(INTERVAL)
			   << Value.get_element_by_index(i, j);
		os << std::endl;
	}
	return os;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
std::istream &operator>>(std::istream &is,
			 matrix<_Tp, _Size_m, _Size_n> &Value) {
	_Tp temp;
	for (tp i = 0; i < _Size_m; i++)
		for (tp j = 0; j < _Size_n; j++) {
			is >> temp;
			Value.modify_element_by_index(i, j, temp);
		}
	Value.set_rank();
	return is;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * */

/* 以下重载了 matrix 的逻辑运算符 == 和 != */
template <class _Tp, size_t _Size_m, size_t _Size_n>
bool operator==(const matrix<_Tp, _Size_m, _Size_n> &Value_1,
		const matrix<_Tp, _Size_m, _Size_n> &Value_2) {
	for (tp i = 0; i < _Size_m; i++)
		for (tp j = 0; j < _Size_n; j++)
			if (Value_1.get_element_by_index(i, j) !=
			    Value_2.get_element_by_index(i, j))
				return false;
	return true;
}

template <class _Tp, size_t _Size_m, size_t _Size_n>
bool operator!=(const matrix<_Tp, _Size_m, _Size_n> &Value_1,
		const matrix<_Tp, _Size_m, _Size_n> &Value_2) {
	return !(Value_1 == Value_2);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * */

/* 矩阵乘法， 注意返回的是局部变量 */
template <class _Tp, size_t _Size_m, size_t _Size_n, size_t _Size_s>// 矩阵乘法
matrix<_Tp, _Size_m, _Size_s>
operator*(const matrix<_Tp, _Size_m, _Size_n> Value_1,
	  const matrix<_Tp, _Size_n, _Size_s> Value_2) {
	matrix<_Tp, _Size_m, _Size_s> Temp;
	for (tp i = 0; i < _Size_m; i++)
		for (tp j = 0; j < _Size_s; j++)
			for (tp k = 0; k < _Size_n; k++)
				/*
         * std::cout << "Temp[" << i << "][" << j << "]: " << Temp.Matrix[i][j]
         *           << "+=" << "Value_1[" << i <<"][" << k << "]: "
         * <<(Value_1.Matrix[k][i])
         *           << " * " << "Value_2[" << k <<"][" << j << "]: " <<
         * (Value_2.Matrix[j][k]) << std::endl;
         */
				Temp.modify_element_by_index(
					i,
					j,
					Temp.get_element_by_index(i, j) +
						Value_1.get_element_by_index(i, k) *
							Value_2.get_element_by_index(k, j));
	return Temp;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * */

}// namespace Linear_algebra
/* namespace Linear_algebra end */

/* Linear_algebra */

#endif /* Basedef_h */
