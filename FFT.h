/*****************************************************************************
* @brief :FFT
* 快速傅里叶变换相关算法
* https://zhuanlan.zhihu.com/p/399448009
* @author : acedtang
* @date : 2021/8/16
* @version : ver 1.0
* @inparam :
* @outparam :
*****************************************************************************/

#ifndef __SUN_FFT_H
#define __SUN_FFT_H

#include <algorithm>
#include <complex>
#include <array>
#include <cmath>
#include <cassert>
#include "Math/MathUtils.h"

namespace Sun {


	//返回2^0+2^1+...+2^(n) = 2^(n+1)-1
	template<int n>
	struct Proportional
	{
		static constexpr int value =  (1 << n) + Proportional<n - 1>().value;
	};
	template<>
	struct  Proportional<0>{
		static constexpr int value = 1;
	};

	template<class _T ,int _MAX_POWER = 8>
	class FFT
	{
	public:
		
		FFT() {
			for (int m = 0; m <= _MAX_POWER; ++m) {
				int n = (1 << m);
				_T cof = 2 * A_PI * (1.0 / (_T)(n));
				for (int k = 0; k < n; ++k) {
					_T tmp = cof* _T(k);
					W[n - 1 + k] = std::complex<_T>(cos(tmp), -sin(tmp));
				}
			}
			table.fill(0);
			for (int m = 1; m <= _MAX_POWER; ++m) {
				int n = (1 << m);
				for (int i = 0; i < n; ++i) {
					table[n - 1 + i] = reverseBit(i, m);
				}
			}
		}

		template<class _InputIter, class _OutIter>
		void run_fft(_InputIter start, _InputIter end, _OutIter ostart) {	
			size_t n = std::distance(start,end);		
			assert(isPowerOfTwo(n));

			//FFT将每两组合并为一组
			//总共多少组
			int times = n;
			//表示一组的长度
			int len =1 ;
			//W[K][N] = e^(-j*2pi* K/N);
			//初始化最底层的一组
			auto iter = ostart;
			for (int i = 0; i < n; ++i) {
				*iter = *(start + table[i + n-1]);
				++iter;
			}

			while (times > 1) {
				auto iter = ostart;
				//每两组作为一个整体循环
				int nt = times >> 1;
				for (int k = 0; k <nt; ++k) {
					//两组中的第一组
					auto iter2 = iter + len;
					for (int i = 0; i < len; ++i) {
						//第二组对应值乘以系数W[K][N]
						*(iter2) =   W[i + (len<<1) -1]* (*(iter2));
						//更新第一组X=G+H*W
						*iter = *iter + *iter2;
						++iter;
						++iter2;
					}
					//两组中的第二组
					for (int j = 0; j < len; ++j) {
						//X=G-W*H=G+W*H-W*H-W*H
						*(iter) = *(iter - len) - *iter - *iter;				
						++iter;
					}
				}
				len = len << 1;
				times = times >> 1;
			}

			return;
		}
			
	protected:
		//将一个数字按照m位进行翻转 (m=8为例: 01001011 => 11010010)
		inline int reverseBit(int k, int m) {
			int output = 0;
			//低位来自原数字的高位右移
			for (int i = 0; i < (m>>1); ++i) {
				output += (k & (1 << (m - i - 1))) >> (m - (i << 1) - 1);
			}
			//高位来自原数字的低位左移
			int t = (m&1)==0?1:0;
			for (int i = (m >> 1); i < m; ++i) {
				output += (k & (1 << (m - i - 1))) << t;
				t += 2;
			}
			return output;
		}

		//  N=2^M,0<=K<N, W[K+N-1] = e^(-j*2*pi* K/N);
		std::complex<_T> W[Proportional<_MAX_POWER>::value];
		//fft位反序映射表
		std::array<int, Proportional<_MAX_POWER>::value> table;
	};

	// just for test
	class DFT {
	public:
		void run_dft(float* data, int n ,std::complex<float>* output) {
			for (int k = 0; k < n; ++k) {
				output[k] = 0;
				for (int i = 0; i < n; ++i) {
					output[k] += std::complex<float>(cos(2 * A_PI * i * k / n), -sin(2 * A_PI * i * k / n)) * data[i];
				}
			}
		}
	};
}

#endif