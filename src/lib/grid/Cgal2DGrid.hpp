#ifndef LIBGCM_CGAL2DGRID_HPP
#define LIBGCM_CGAL2DGRID_HPP

#include <mpi.h>

#include <lib/grid/AbstractGrid.hpp>
#include <lib/linal/linal.hpp>
#include <lib/util/storage/Cgal2DTriangulationStorage.hpp>
#include <lib/util/snapshot/VtkCgal2DSnapshotter.hpp>


namespace gcm {
	/**
	 * 2D movable unstructured triangle grid by CGAL library
	 */
	class Cgal2DGrid : public AbstractGrid {
	public:
		struct Iterator {
			size_t iter = 0;
			Iterator(size_t value) : iter(value) { };
			const Iterator& operator*() { return *this; };
			Iterator& operator++() {
				iter++;
				return (*this);
			};
		};
		typedef Iterator ForwardIterator;
		ForwardIterator begin() const { return 0; };
		ForwardIterator end() const { return sizeOfRealNodes(); };
		typedef Iterator VtkIterator;
		VtkIterator vtkBegin() const { return begin(); };
		VtkIterator vtkEnd() const { return end(); };

	protected:
		/**
		 * @param it begin() <= iterator < end()
		 * @return index in std::vector
		 */
		size_t getIndex(const Iterator& it) const {
			return it.iter;
		};
		/** Coordinates, cells, etc ... */
		Cgal2DTriangulationStorage mesh;

	public:
		/** Read-only access to real coordinates */
		const linal::Vector3 coords(const Iterator& it) const {
			return mesh.getCoordinates(it.iter);
		};

		size_t sizeOfRealNodes() const {
			return mesh.verticesSize();
		};
		size_t sizeOfAllNodes() const {
			return sizeOfRealNodes();
		};

	protected:

		virtual void initializeImpl(const Task &task) override;
		virtual void beforeStageImpl() override { };
		virtual void afterStageImpl() override { };
		virtual void beforeStepImpl() override { };
		virtual void afterStepImpl() override { };
		virtual void recalculateMinimalSpatialStep() override { /* TODO for movable grid */ };

		virtual void applyInitialConditions(const Task& task) = 0;
		virtual void recalculateMaximalLambda() = 0;

		virtual void initializeImplImpl(const Task& task) = 0;

		USE_AND_INIT_LOGGER("gcm.Cgal2DGrid");
	};
}

#endif // LIBGCM_CGAL2DGRID_HPP
