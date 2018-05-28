#pragma once
#include "Matrix.h"
namespace mat {
	template<int N>
	using Rvec = mat::Matrix<float, N, 1>;

	namespace Detail {
		template<typename T>
		T constexpr sqrtNewtonRaphson(T x, T lo, T hi) {
			const T mid = T{ 0.5 } * (lo + x / lo);
			return (lo == hi) ? lo : sqrtNewtonRaphson(x, mid, lo);
		}
		constexpr int isqrt_helper(int x, int lo, int hi) {
			const int mid = (lo + hi + 1) / 2;
			return (lo == hi) ? lo : ( (x / mid < mid) ? isqrt_helper(x, lo, mid - 1) : isqrt_helper(x, mid, hi) );
		}
	}
	template<typename T>
	T constexpr sqrt(T x) {
		return Detail::sqrtNewtonRaphson(x, x, T{ 0 });
	}
	constexpr int isqrt(int x) {
		return Detail::isqrt_helper(x, 0, x / 2 + 1);
	}
	constexpr bool isSquare(int n) {
		return n == (isqrt(n)*isqrt(n));
	}
	template<typename T>
	constexpr T abs(T x) {
		return (x > T{ 0 }) ? x : -x;
	}
	template<typename T>
	constexpr T power(T x, int n) {
		return n == 0 ? T{ 1 } : x * power(x, n - 1);
	}
	constexpr long long int factorial(int n) {
		return n <= 1 ? 1 : (n * factorial(n - 1));
	}
	template<typename T>
	constexpr T TaylorSin(T x, int N) {
		return -(2 * (N % 2) - 1) * (power(x, 2 * N + 1) / factorial(2 * N + 1));
	}
	template<typename T>
	constexpr T TaylorCosine(T x, int N) {
		return -(2 * (N % 2) - 1) * (power(x, 2 * N) / factorial(2 * N));
	}
	template<typename T>
	constexpr T sine(T x) {
		T res{};
		for (size_t i = 0; i < 10; i++)	{
			res +=TaylorSin(x, i);
		}
		return res;
	}
	template<typename T>
	constexpr T cosine(T x) {
		T res{};
		for (size_t i = 0; i < 10; i++) {
			res += TaylorCosine(x, i);
		}
		return res;
	}
	template<typename T>
	constexpr T tangent(T x) {
		return sine(x) / cosine(x);
	}
	template<typename T, int N>
	constexpr T dot(const Matrix<T, N, 1> lhs, const Matrix<T, N, 1> rhs) {
		T res{};
		for (size_t i = 0; i < N; i++) {
			res += lhs.value[i] * rhs.value[i];
		}
		return res;
	}
	template<typename T, int N>
	constexpr T norm(const Matrix<T, N, 1>& x) {
		return sqrt(dot(x, x));
	}
	constexpr float distance(const Matrix<float, 2, 1> line1, Matrix<float, 2, 1> line2, const Matrix<float, 2, 1> point) {
		auto p1 = line1;
		auto p2 = line2;

		auto x = Matrix<float, 2, 1>{ 1.0f, 0.0f };
		auto y = Matrix<float, 2, 1>{ 0.0f, 1.0f };

		auto x0 = dot(point, x);
		auto y0 = dot(point, y);
		auto x1 = dot(p1, x);
		auto y1 = dot(p1, y);
		auto x2 = dot(p2, x);
		auto y2 = dot(p2, y);
		auto rel = p2 - p1;
		return ((y2 - y1)*x0 - (x2 - x1)*y0 + x2 * y1 - y2 * x1) / norm(rel);
	}
	template <typename T>
	constexpr Matrix<T, 3, 3> crossProductMatrix(const Matrix<T, 3, 1> a) {
		std::array<T, 9> tmp{
			0, a(3, 1),  a(2, 1),
			a(3, 1), 0, a(1, 1),
			a(2, 1), a(1, 1), 0
		};
		return AzMat::Matrix<T, 3, 3> {tmp};
	}

	template<typename T>
	constexpr Matrix<T, 4, 4> perspective(const T fovy, const T aspect, const T zNear, const T zFar) {
		T const tanHalfFovy = tangent(fovy / T{ 2 });
		auto i = T{ 1 } / (aspect * tanHalfFovy);
		auto j = T{ 1 } / (tanHalfFovy);
		auto l = -(zFar + zNear) / (zFar - zNear);
		auto m = -(T{ 2 } *zFar * zNear) / (zFar - zNear);
		return Matrix<T, 4, 4> {
			i, T{ 0 }, T{ 0 }, T{ 0 },
				T{ 0 }, j, T{ 0 }, T{ 0 },
				T{ 0 }, T{ 0 }, l, m,
				T{ 0 }, T{ 0 }, -T{ 1 }, T{ 0 }
		};
	}
	template <typename T, int M, int N>
	constexpr Matrix<T, N, M> transpose(const Matrix<T, M, N>& a) {
		std::array<T, N*M> tmp{};
		for (size_t i = 1; i <= M; i++) {
			for (size_t j = 1; j <= N; j++) {
				tmp[(j - 1)*M + (i - 1)] = a(i, j);
			}
		}
		return Matrix<T, N, M>{tmp};
	}
	template<typename T>
	constexpr Matrix<T, 4, 4> translationMatrix(const Matrix<T, 3> r) {
		return Matrix<float, 4, 4> {
			1.f, 0.f, 0.f, r(1),
				0.f, 1.f, 0.f, r(2),
				0.f, 0.f, 1.f, r(3),
				0.f, 0.f, 0.f, 1.f
		};
	}
	//special orthogonal group in R2 is the rotation matrices given by some radian scalar
	constexpr Matrix<float, 2, 2> SO2(const float theta) {
		return {
			cosine(theta), -sine(theta),
			sine(theta), cosine(theta)
		};
	}
	constexpr Matrix<float, 2, 2> SO2(R2Vec c) {
		return {
			c(1,1), -c(2,1),
			c(2,1), c(1,1)
		};
	}
	constexpr Matrix<float, 3, 3> rotateX(const float theta) {
		return {
			1.0f, 0.0f, 0.0f,
			0.0f, mat::cosine(theta), -mat::sine(theta),
			0.0f, mat::sine(theta), mat::cosine(theta) };
	}
	constexpr Matrix<float, 3, 3> rotateY(const float theta) {
		return {
			mat::cosine(theta), 0.0f, mat::sine(theta),
			0.0f, 1.0f, 0.0f,
			-mat::sine(theta), 0.0f, mat::cosine(theta) };
	}
	constexpr Matrix<float, 3, 3> rotateZ(const float theta) {
		return {
			mat::cosine(theta), -mat::sine(theta), 0.0f,
			mat::sine(theta), mat::cosine(theta), 0.0f,
			0.0f, 0.0f, 1.0f };
	}
	template<int N>
	std::array<R2Vec, N> rootsOfUnity(float base) {
		std::array<R2Vec, N> point{};
		auto first = R2Vec{ base, 0.0f };
		float Theta = (2.f*3.14159f) / N;
		for (size_t i = 0; i < N; i++) {
			first = SO<2>(Theta)*first;
			point[i] = first;
		}
		return point;
	}
	std::array<Rvec<3>, 36> cubePoints(mat::Matrix<float, 3, 3> M = mat::make_Iden<float, 3>(), Rvec<3> b = {0.f,0.f,0.f}) {
		auto zero = Rvec<3>{ 0.f, 0.f, 0.f };
		auto x = Rvec<3>{ 1.f, 0.f, 0.f };
		auto y = Rvec<3>{ 0.f, 1.f, 0.f };
		auto z = Rvec<3>{ 0.f, 0.f, 1.f };
		std::array<Rvec<3>, 3> basis{x, y, z};
		std::transform(basis.begin(), basis.end(), basis.begin(),
			[&](Rvec<3> x) { return M * x + b; });
		return {
			// Back Face
			zero, x, x+y,
			x+y, y, zero,

			z, x+z, x+y+z,
			Rvec<3>{ 1.f, 1.f, 1.f }, Rvec<3>{ 0.f, 1.f, 1.f },	Rvec<3>{ 0.f, 0.f, 1.f },

			Rvec<3>{ 0.f, 1.f, 1.f }, Rvec<3>{ 0.f, 1.f, 0.f }, Rvec<3>{ 0.f, 0.f, 0.f },
			Rvec<3>{ 0.f, 0.f, 0.f }, Rvec<3>{ 0.f, 0.f, 1.f }, Rvec<3>{ 0.f, 1.f, 1.f },

			Rvec<3>{ 1.f, 1.f, 1.f }, Rvec<3>{ 1.f, 1.f, 0.f }, Rvec<3>{ 1.f, 0.f, 0.f },
			Rvec<3>{ 1.f, 0.f, 0.f }, Rvec<3>{ 1.f, 0.f, 1.f }, Rvec<3>{ 1.f, 1.f, 1.f },

			Rvec<3>{ 0.f, 0.f, 0.f }, Rvec<3>{ 1.f, 0.f, 0.f }, Rvec<3>{ 1.f, 0.f, 1.f },
			Rvec<3>{ 1.f, 0.f, 1.f }, Rvec<3>{ 0.f, 0.f, 1.f }, Rvec<3>{ 0.f, 0.f, 0.f },

			Rvec<3>{ 0.f, 1.f, 0.f }, Rvec<3>{ 1.f, 1.f, 0.f }, Rvec<3>{ 1.f, 1.f, 1.f },
			Rvec<3>{ 1.f, 1.f, 1.f }, Rvec<3>{ 0.f, 1.f, 1.f }, Rvec<3>{ 0.f, 1.f, 0.f }
		};
	};

}