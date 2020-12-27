#include "mini.h"

#ifdef AETHER_MINUIT
#include <Minuit2/MnMigrad.h>
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPlot.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MinosError.h>
#include <Minuit2/MnHesse.h>
#include <Minuit2/MnContours.h>
#endif

namespace aether {


#ifdef AETHER_MINUIT

/**
 * Minuit2 callback function object
 */
class MinimizeContextFCN : public ROOT::Minuit2::FCNBase
{
	MinimizeContext& ctx;
public:

	MinimizeContextFCN(MinimizeContext& ctx)
		:ctx{ctx}
	{
	}

	double operator()(const std::vector<double>& params) const override
	{
		return ctx(params);
	}

	double Up() const override
	{
		return ctx.delta_chi_squared();
	}
};




MinimizerResult
minimize_locally_Minuit2(const std::vector<double>& x0, MinimizeContext& context,std::ostream& os)
{
	MinimizeContextFCN fcn(context);

	// https://root.cern.ch/root/htmldoc/guides/minuit2/Minuit2.html#interpretation-of-parameter-errors
	// http://seal.cern.ch/documents/minuit/mnusersguide.pdf
	// http://seal.cern.ch/documents/minuit/mntutorial.pdf
	// http://seal.web.cern.ch/seal/documents/minuit/mnerror.pdf

	ROOT::Minuit2::MnUserParameters params;
	for (size_t i=0;i<x0.size();i++) {
		switch (i) {
		case 0:
			params.Add("v", x0.at(0),0.1);
			break;
		case 1:
			params.Add("α",x0.at(1),0.1);
			break;
		case 2:
			params.Add("δ",x0.at(2),0.1);
			break;
		default:
			std::string s("p");
			s.append(std::to_string(i));
			params.Add(s,x0.at(i),0.1);
		}
	}


	ROOT::Minuit2::MnMigrad migrad(fcn,params,2);
	if (context.fixed_ad()) {
		migrad.Fix(1);
		migrad.Fix(2);
	}
	ROOT::Minuit2::FunctionMinimum fm = migrad();	
	if (fm.IsValid()) {
		ROOT::Minuit2::MnHesse hesse;
		hesse(fcn,fm);
	}
	os << fm;

	MinimizerResult result;
	result.valid = fm.IsValid();
	result.y = fm.Fval();

	if (fm.Parameters().IsValid()) {
		unsigned int k=0;
		for (unsigned int i=0;i<x0.size();i++) {
			if (migrad.Parameter(i).IsFixed()) {
				result.x.push_back(x0.at(i));
			}
			else {
				result.x.push_back(fm.Parameters().Vec()[k]);
				k++;
			}
		}
	}

	if (fm.IsValid()) {
		ROOT::Minuit2::MnMinos minos(fcn,fm);
		for (unsigned int i=0;i<x0.size();i++) {
			if (migrad.Parameter(i).IsFixed()) {
				result.u.push_back(0.0);
			}
			else {
				auto error = minos.Minos(i);
				os << error;
				auto e = std::max(std::abs(error.Lower()),std::abs(error.Upper()));
				result.u.push_back(e);
			}			
		}
	}
	else {
		if (fm.Error().HesseFailed()) {
			os << "warning: Hesse matrix failed\n";
		}
		if (!fm.Error().IsAccurate()) {
			os << "warning: errors are not accurate\n";
		}

		if (fm.Error().IsAvailable()) {
			for (unsigned int i=0;i<fm.Error().Matrix().Nrow();i++) {
				result.u.push_back(std::sqrt(context.delta_chi_squared() * fm.Error().Matrix()(i,i)));
			}
		}
	}


	return result;
}

#endif

}// aether
