// FFT.cpp : Defines the entry pounsigned for the console application.
//

#include "stdafx.h"
#include <cmath>
#include <iostream>
#include <iomanip>

template <unsigned B>	struct POWER2		{ enum { value = 2 * POWER2<B-1>::value, }; };
template <>			struct POWER2<1>	{ enum { value = 2, }; };
template <unsigned B>	struct POWER4		{ enum { value = 4 * POWER4<B-1>::value, }; };
template <>			struct POWER4<1>	{ enum { value = 4, }; };
template <unsigned B>	struct POWER8		{ enum { value = 8 * POWER8<B-1>::value, }; };
template <>			struct POWER8<1>	{ enum { value = 8, }; };

template <typename T>
class Complex {
public:
	Complex() : real(T(0)), imag(T(0)) {}
	Complex(const T& r, const T& i) : real(T(r)), imag(T(i)) {}
	~Complex() {}

	Complex<T> operator+ (const Complex<T>& c) const {
		return Complex<T>(real + c.real, imag + c.imag);
	}
	Complex<T> operator- (const Complex<T>& c) const {
		return Complex<T>(real - c.real, imag - c.imag);
	}
	Complex<T> operator* (const Complex<T>& c) const {
		return Complex<T>(real*c.real-imag*c.imag, imag*c.real+real*c.imag);
	}

	T real;
	T imag;
};

#define BIT 3
#define N (POWER2<BIT>::value)

typedef float R;
typedef Complex<R> C;

unsigned BitMirror[N];

static const R PI = R(3.14159265358979323846);

void init() {
	unsigned mirror = 0;
	for (unsigned i = 0; i < N; i++) {
		BitMirror[i] = mirror;
		unsigned mask = N;
		while (mirror & (mask >>= 1)) mirror &= ~mask;
		mirror |= mask;
	}
}

void BitReverse(const unsigned& nn, C* const data) {
	for (unsigned i = 0; i < nn; i++) {
		unsigned mirror = BitMirror[i];
		if (mirror > unsigned(i)) {
			C temp = data[i];
			data[i] = data[mirror];
			data[mirror] = temp;
		}
	}
}

C W(const unsigned& k) {
	const R theta = R(-(2 * PI / N)*k);
	return C(cos(theta), sin(theta));
}

C EI(const R& k) {
	return C(cos(k), sin(k));
}

void DFT(const unsigned& nn, const C* const input, C* const output) {
	for (unsigned i = 0; i < nn; i++) {
		output[i] = C();
		for (unsigned j = 0; j < nn; j++) {
			output[i] = output[i] + input[j] * W(i*j);
		}
	}
}

void Radix2(const unsigned& nn, const C* const input, C* const output) {

}
void Radix2(const unsigned& nn, C* const data) {
	for (unsigned step = 1; step < (nn); step <<= 1) {
		const unsigned jump = step << 1;
		for (unsigned factorId = 0; factorId < step; factorId++) {
			const C factor = EI(-(PI / R(step))* factorId);
			for (unsigned a = factorId; a<(nn); a += jump) {
				const unsigned b = a + step;
				const C wb = factor*data[b];
				data[b] = data[a] - wb;
				data[a] = data[a] + wb;
			}
		}
	}
	//const R pi = R(-3.14159265358979323846);
	////   Iteration through dyads, quadruples, octads and so on...
	//for (unsigned int Step = 1; Step < N; Step <<= 1) {
	//	//   Jump to the next entry of the same transform factor
	//	const unsigned int Jump = Step << 1;
	//	//   Angle increment
	//	const R delta = pi / R(Step);
	//	//   Auxiliary sin(delta / 2)
	//	const R Sine = sin(R(delta * .5));
	//	//   Multiplier for trigonometric recurrence
	//	const C Multiplier(R(-2. * Sine * Sine), sin(delta));
	//	//   Start value for transform factor, fi = 0
	//	C Factor(1., 0);
	//	//   Iteration through groups of different transform factor
	//	for (unsigned int Group = 0; Group < Step; ++Group) {
	//		//   Iteration within group 
	//		for (unsigned int Pair = Group; Pair < N; Pair += Jump) {
	//			//   Match position
	//			const unsigned int Match = Pair + Step;
	//			//   Second term of two-point transform
	//			const C Product = Factor * data[Match];
	//			//   Transform for fi + pi
	//			data[Match] = data[Pair] - Product;
	//			//   Transform for fi
	//			data[Pair] = data[Pair] + Product;
	//		}
	//		//   Successive transform factor via trigonometric recurrence
	//		Factor = Multiplier * Factor + Factor;
	//	}
	//}
}

using namespace std;
unsigned _tmain(unsigned argc, _TCHAR* argv[])
{
	C x[N];
	C X[N];
	init();

	cout << " N : " << N << endl;
	for (unsigned i = 0; i < N; i++) {
		x[i] = C(sin(R(i*2*PI/N)), 0);
	}
	DFT(N, x, X);
	BitReverse(N, x);
	Radix2(N, x);
	cout << " real " << "\t" << " imag " << endl;
	for (unsigned i = 0; i < N; i++) {
		cout << fixed << x[i].real - X[i].real << "\t" << x[i].imag - X[i].imag << endl;
		//cout << fixed << x[i].real << "\t" << x[i].imag << "\t" << fixed << X[i].real << "\t" << X[i].imag << endl;
	}

	return 0;
}

