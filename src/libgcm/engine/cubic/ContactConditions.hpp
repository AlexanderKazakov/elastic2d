#ifndef LIBGCM_CUBIC_CONTACTCONDITIONS_HPP
#define LIBGCM_CUBIC_CONTACTCONDITIONS_HPP

#include <libgcm/engine/cubic/DefaultMesh.hpp>
#include <libgcm/rheology/materials/materials.hpp>


namespace gcm {
namespace cubic {


/**
 * Applying contact conditions for cubic meshes.
 * The approach is ghost (fixture) nodes on borders.
 * Setting appropriate values in ghost nodes
 * from contact neighbors before calculation,
 * we can calculate contacts in the same manner as inner nodes.
 */
template<typename TGrid>
class AbstractContactCopier {
public:
	typedef typename TGrid::PartIterator PartIterator;
	
	AbstractContactCopier(const PartIterator& boxA_, const PartIterator& boxB_) :
			boxA(boxA_), boxB(boxB_) {
		boxA.assertBoundsValid();
		boxB.assertBoundsValid();
		assert_eq(boxA.size(), boxB.size());
	}
	
	
	/**
	 * Apply contact conditions by appropriate values copying
	 * from mesh b to mesh a (b is immutable)
	 */
	virtual void apply(
			AbstractMesh<TGrid>& a, const AbstractMesh<TGrid>& b) = 0;
	
	
protected:
	PartIterator boxA;
	PartIterator boxB;
	
};


template<typename TGrid,
         typename MeshA, typename MeshB>
struct ContactCopier : public AbstractContactCopier<TGrid> {
	typedef AbstractContactCopier<TGrid>    Base;
	typedef typename Base::PartIterator     PartIterator;
	
	ContactCopier(const PartIterator& boxA_, const PartIterator& boxB_) :
			Base(boxA_, boxB_) { }
	
	virtual void apply(
			AbstractMesh<TGrid>& a, const AbstractMesh<TGrid>& b) override {
		const MeshB& meshB = dynamic_cast<const MeshB&>(b);
		      MeshA& meshA = dynamic_cast<      MeshA&>(a);
		
		PartIterator iterB = this->boxB.begin();
		for (PartIterator iterA  = this->boxA.begin();
		                  iterA != this->boxA.end(); ++iterA) {
			meshA._pde(iterA) = meshB.pde(iterB);
			++iterB;
		}
		assert_true(iterB == iterB.end());
	}
	
};


} // namespace cubic
} // namespace gcm


#endif // LIBGCM_CUBIC_CONTACTCONDITIONS_HPP
