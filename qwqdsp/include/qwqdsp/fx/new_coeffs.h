#pragma once
#include <array>
#include <complex>

template<typename Sample>
struct FastCoeffs {
	static constexpr size_t complexCount = 5;
	static constexpr size_t realCount = 1;
	static constexpr size_t filterOrder = 11;
	static constexpr Sample fpass = Sample(20000.000000);
	static constexpr Sample fstop = Sample(24168.597105);
	static constexpr std::array<std::complex<Sample>, complexCount> complexPoles{{
		{Sample(-5561.982558545558), Sample(7721.564144482831)},
		{Sample(-3936.7227375705806), Sample(13650.197801811304)},
		{Sample(-2348.139919173165), Sample(17360.271619482097)},
		{Sample(-1177.5927594523114), Sample(19350.80475283362)},
		{Sample(-351.83634005467775), Sample(20192.238514865294)}
	}};
	static constexpr std::array<Sample, realCount> realPoles{{
		Sample(-6297.9970566057245)
	}};
	// Coeffs for direct bandlimited synthesis of a polynomial-segment waveform
	static constexpr std::array<std::complex<Sample>, complexCount> complexCoeffsDirect{{
		{Sample(-16437.98665421065), Sample(-7224.765894589879)},
		{Sample(7790.757942544076), Sample(9526.642450543113)},
		{Sample(-840.5346697125525), Sample(-6786.810488747387)},
		{Sample(-1524.1587677384439), Sample(2562.99591951713)},
		{Sample(754.6036883901486), Sample(-311.82754418710147)}
	}};
	static constexpr std::array<Sample, realCount> realCoeffsDirect{{
		Sample(10260.028875848706)
	}};
};
	
template<typename Sample>
struct BestCoeffs {
	static constexpr size_t complexCount = 12;
	static constexpr size_t realCount = 1;
	static constexpr size_t filterOrder = 25;
	static constexpr Sample fpass = Sample(20000.000000);
	static constexpr Sample fstop = Sample(21235.713920);
	static constexpr std::array<std::complex<Sample>, complexCount> complexPoles{{
		{Sample(-3341.949128109254), Sample(4268.074767165956)},
		{Sample(-2946.293299786584), Sample(8172.223758516495)},
		{Sample(-2415.1493874294415), Sample(11465.74957710688)},
		{Sample(-1865.2641743625645), Sample(14062.862883541336)},
		{Sample(-1375.4750370263282), Sample(16005.955390059582)},
		{Sample(-979.3741564693673), Sample(17403.81522252918)},
		{Sample(-678.4379295586048), Sample(18379.376540414632)},
		{Sample(-459.61279109499947), Sample(19050.194942808404)},
		{Sample(-296.1237905055368), Sample(19474.679698907876)},
		{Sample(-205.63167845077965), Sample(19797.36066586945)},
		{Sample(-76.3965030491654), Sample(19914.70412883792)},
		{Sample(-42.01130754433552), Sample(20027.196942402177)}
	}};
	static constexpr std::array<Sample, realCount> realPoles{{
		Sample(-3489.664242296751)
	}};
	// Coeffs for direct bandlimited synthesis of a polynomial-segment waveform
	static constexpr std::array<std::complex<Sample>, complexCount> complexCoeffsDirect{{
		{Sample(-12537.983371715141), Sample(-3551.1891483546524)},
		{Sample(9723.747737188487), Sample(6081.318808686397)},
		{Sample(-6156.288314097386), Sample(-7065.61069094887)},
		{Sample(2842.2177910952064), Sample(6615.118256303697)},
		{Sample(-404.2768278045403), Sample(-5249.949733520729)},
		{Sample(-996.7696878599824), Sample(3564.860744018354)},
		{Sample(1509.4005775326725), Sample(-2004.4340405566234)},
		{Sample(-1447.0499676056818), Sample(806.5497753985597)},
		{Sample(1006.4705263731172), Sample(-104.42463656860839)},
		{Sample(-420.402419817584), Sample(-276.7771081538272)},
		{Sample(102.48085265065643), Sample(279.18827194919294)},
		{Sample(-28.793807131041394), Sample(-55.92936818990831)}
	}};
	static constexpr std::array<Sample, realCount> realCoeffsDirect{{
		Sample(6807.244796746716)
	}};
};
	