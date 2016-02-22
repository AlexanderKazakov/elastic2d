#ifndef LIBGCM_CUBICGRID_HPP
#define LIBGCM_CUBICGRID_HPP

#include <lib/mesh/grid/StructuredGrid.hpp>
#include <lib/linal/linal.hpp>

namespace gcm {
	/**
	 * Non-movable structured rectangular grid
	 * @tparam TModel rheology model
	 */
	class CubicGrid : public StructuredGrid {
	public:
		virtual ~CubicGrid() { };
		struct Iterator : public linal::VectorInt<3> {
			linal::VectorInt<3> bounds = {0, 0, 0};
			Iterator(const linal::VectorInt<3> indices_, const linal::VectorInt<3> bounds_) {
				(*this) = indices_;
				bounds = bounds_;
			};
			const Iterator& operator*() { return *this; };
			using linal::VectorInt<3>::Matrix;
			using linal::VectorInt<3>::operator=;
		protected:
			int increment(const int i) {
				(*this)(i) = ((*this)(i) + 1) % this->bounds(i);
				return (*this)(i);
			};
		};

		/** slow-X fast-Z memory access efficient iterator */
		struct ForwardIterator : public Iterator {
			using Iterator::Iterator;
			Iterator& operator++() {
				!increment(2) && !increment(1) && ((*this)(0)++);
				return (*this);
			};
		};
		ForwardIterator begin() const { return ForwardIterator({0, 0, 0}, sizes); };
		ForwardIterator end() const { return ForwardIterator({sizes(0), 0, 0}, sizes); };

		/** slow-Z fast-X iterator */
		struct VtkIterator : public Iterator {
			using Iterator::Iterator;
			VtkIterator& operator++() {
				!increment(0) && !increment(1) && ((*this)(2)++);
				return (*this);
			};
		};
		VtkIterator vtkBegin() const { return VtkIterator({0, 0, 0}, sizes); };
		VtkIterator vtkEnd() const { return VtkIterator({0, 0, sizes(2)}, sizes); };

		/** Read-only access to real coordinates */
		const linal::Vector3 coords(const int x, const int y, const int z) const {
			return startR + linal::plainMultiply(linal::VectorInt<3>({x, y, z}), h);
		};
		const linal::Vector3 coords(const Iterator& it) const {
			return startR + linal::plainMultiply(it, h);
		};

		size_t sizeOfRealNodes() const {
			return (size_t) linal::directProduct(sizes);
		};
		size_t sizeOfAllNodes() const {
			return (size_t)linal::directProduct(sizes + 2 * accuracyOrder * linal::VectorInt<3>({1, 1, 1}));
		};

		linal::VectorInt<3> getSizes() const { return sizes; };
		int getAccuracyOrder() const { return accuracyOrder; };
		linal::Vector<3> getH() const { return h; };

	protected:
		/**
		 * @param it begin() <= iterator < end()
		 * @return index in std::vector
		 */
		size_t getIndex(const Iterator& it) const {
			return getIndex(it(0), it(1), it(2));
		};
		/**
		 * @param x x index < sizes(0)
		 * @param y y index < sizes(1)
		 * @param z z index < sizes(2)
		 * @return index in std::vector
		 */
		size_t getIndex(const int x, const int y, const int z) const {
			return (size_t)
					((2 * accuracyOrder + sizes(2)) * (2 * accuracyOrder + sizes(1)) * (x + accuracyOrder)
			       + (2 * accuracyOrder + sizes(2)) * (y + accuracyOrder)
			       + (z + accuracyOrder));
		};

		int accuracyOrder = 0; // order of accuracy of spatial interpolation

		linal::VectorInt<3> sizes = {0, 0, 0}; // numbers of nodes along each direction (on this core)
		linal::Vector<3> startR = {0, 0, 0}; // global coordinates of the first real node of the grid
		linal::Vector<3> h = {0, 0, 0}; // spatial steps along each direction

		virtual void initializeImpl(const Task &task) override;

		virtual void recalculateMaximalLambda() = 0;
		virtual void applyInitialConditions(const Task& task) = 0;

		virtual void initializeImplImpl(const Task& task) = 0;

		USE_AND_INIT_LOGGER("gcm.CubicGrid");
	};
}

#endif // LIBGCM_CUBICGRID_HPP
