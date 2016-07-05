#ifndef LIBGCM_SUPERDUPERMODEL_HPP
#define LIBGCM_SUPERDUPERMODEL_HPP

#include <lib/rheology/models/Model.hpp>
#include <lib/rheology/models/ElasticModel.hpp>

namespace gcm {


class SuperDuperModel : public ElasticModel<3> {
public:
	typedef ContinualDamageOde                     InternalOde;
	typedef IdealPlasticFlowCorrector              Corrector;

	typedef typename InternalOde::Variables       OdeVariables;
};


}


#endif // LIBGCM_SUPERDUPERMODEL_HPP
