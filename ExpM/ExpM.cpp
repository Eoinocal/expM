
#include "stdafx.h"

// Matrix inversion routine.
// Uses lu_factorize and lu_substitute in uBLAS to invert a matrix
template<typename T>
boost::optional<ub::matrix<T>> invert(const ub::matrix<T>& input)
{
	// create a working copy of the input
	ub::matrix<T> A(input);
	// create a permutation matrix for the LU-factorization
	ub::permutation_matrix<std::size_t> pm(A.size1());

	// perform LU-factorization
	size_t res = lu_factorize(A, pm);
	if (res != 0) return {};

	// create identity matrix of "inverse"
	ub::matrix<T> inverse(ub::identity_matrix<T>(A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return inverse;
}

template<typename T>
ub::matrix<T> expm_taylor(const ub::matrix<T>& m, int iter = 50)
{
	ub::matrix<T> exp_m{ub::identity_matrix<T>{m.size1()}};

	auto tmp{m};
	T j = 1;

	exp_m += tmp;

	for (int i = 0; i < iter; ++i)
	{
		tmp = prod(tmp, m);
		tmp /= ++j;

		exp_m += tmp;
	}

	return exp_m;
}

template<typename T>
ub::matrix<T> expm_taylor(vex::Context& ctx, const ub::matrix<T>& A, int iter = 50)
{
	vcl::matrix<T> vcl_A(A.size1(), A.size2());
	vcl::copy(A, vcl_A);

	auto exp_A = expm_taylor(ctx, vcl_A, iter);

	ub::matrix<T> out_A{A.size1(), A.size2()};
	vcl::copy(exp_A, out_A);

	return out_A;
}

template<typename T>
vcl::matrix<T> expm_taylor(vex::Context& ctx, const vcl::matrix<T>& vcl_A, int iter = 50)
{
	vcl::matrix<T> tmp_A(vcl_A);
	vcl::matrix<T> exp_A(vcl::identity_matrix<T>(vcl_A.size1()));
	T j = 1;

	exp_A += vcl_A;

	for (int i = 0; i < iter; ++i)
	{
		tmp_A = vcl::linalg::prod(tmp_A, vcl_A);
		tmp_A /= ++j;

		exp_A += tmp_A;
	}
	
	return exp_A;
}

template<typename T>
ub::matrix<T> expm_pade_pq(const ub::matrix<T>& M, int p, int q)
{
	auto N = [&]() -> ub::matrix<T>
	{
		ub::matrix<T> coeff{ub::identity_matrix<T>{M.size1()}};
		ub::matrix<T> N_m{ub::identity_matrix<T>{M.size1()}};

		for (int k = 1; k < p; ++k)
		{
			coeff = static_cast<T>(p - k + 1) / ((p + q - k + 1) * k) * prod(coeff, M);
			N_m = N_m + coeff;
		}

		return N_m;
	};

	auto D = [&]() -> ub::matrix<T>
	{
		ub::matrix<T> coeff{ub::identity_matrix<T>{M.size1()}};
		ub::matrix<T> D_m{ub::identity_matrix<T>{M.size1()}};

		for (int k = 1; k < q; ++k)
		{
			coeff = static_cast<T>(q - k + 1) / ((p + q - k + 1) * (k)) * prod(coeff, -M);
			D_m = D_m + coeff;
		}

		return D_m;
	};

	return prod(*invert(D()), N());
}

template<typename T>
ub::matrix<T> expm_pade_ss(const ub::matrix<T>& M, int q, int j)
{
	ub::matrix<T> t{M / pow(T{2}, j)};

	auto eM = expm_pade_pq(t, q, q);

	for (int i = 0; i < j; ++i)
		eM = ub::prod(eM, eM);

	return eM;
}

template<typename T, typename S>
ub::matrix<T> expm_ode(const ub::matrix<T>& M, S&& stepper)
{
	using state_type = ub::matrix<T>;
	ub::matrix<T> eM{ub::identity_matrix<T>{M.size1()}};

	auto steps = ode::integrate_adaptive(stepper,
		[&](const state_type& X, state_type& dXdt, double t)
		{
			dXdt = ub::prod(M, X);
		},
		eM,
		0., 1., 0.01);
	
	return eM;
}

template<typename T>
ub::matrix<T> load(std::istream& src)
{
	unsigned w, h;

	src >> w;
	src >> h;

	ub::matrix<T> m{w, h};

	for (unsigned i = 0; i < w; ++i)
		for (unsigned j = 0; j < h; ++j)
		{
			src >> m.at_element(i, j);
		}

	return m;
}

template<typename F>
void compare(vex::Context ctx, F&& exp)
{
	using namespace std::string_literals;

	vex::profiler<> prof(ctx);

	for (const auto& name : {"mat1"s, "mat2"s})
	{
		std::cout << "Testing " << name << ": ";

		const ub::matrix<double> m = load<double>(std::ifstream(name + ".txt"));
		const ub::matrix<double> ans = load<double>(std::ifstream(name + "_ans.txt"));

		prof.tic_cl(name);
		const ub::matrix<double> exp_m = exp(m);
		prof.toc(name);

		std::cout << ub::norm_inf(exp_m - ans) << "\n";
	}

	std::cout << "Testing rand: ";

	std::uniform_real_distribution<double> distribution(-0.5, 0.5);
	std::default_random_engine generator;

	ub::matrix<double> m(256, 256);

	for (auto i = m.begin1(), e = m.end1(); i != e; ++i)
		std::generate(i.begin(), i.end(), [&]() { return distribution(generator) / m.size1(); });

	const ub::matrix<double> ans{m};

	prof.tic_cl("rand");
	const ub::matrix<double> exp_m = exp(m);
	prof.toc("rand");

	std::cout << ub::norm_inf(exp_m - ans) << "\n";
	
	std::cout << prof << std::endl;
	std::cout << std::endl;
}

int main()
{
	vex::Context ctx(vex::Filter::GPU && vex::Filter::DoublePrecision);
	if (!ctx) throw std::runtime_error("No GPUs with double precision found");
	std::cout << ctx << std::endl;

	compare(ctx, [](const ub::matrix<double>& m) { return expm_taylor(m, 20); });

	compare(ctx, [&](const ub::matrix<double>& m) { return expm_taylor(ctx, m, 20); });

	compare(ctx, [](const ub::matrix<double>& m) { return expm_pade_pq(m, 10, 10); });

	compare(ctx, [](const ub::matrix<double>& m) { return expm_pade_ss(m, 5, 8); });

	compare(ctx, [](const ub::matrix<double>& m) { return expm_ode(m, ode::make_controlled<ode::runge_kutta_cash_karp54<ub::matrix<double>>>(1.e-11, 1.e-11)); });

	system("pause");
	return 0;
}
