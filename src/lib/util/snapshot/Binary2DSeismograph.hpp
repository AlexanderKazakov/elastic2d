#ifndef LIBGCM_BINARY2DSEISMOGRAPH_HPP
#define LIBGCM_BINARY2DSEISMOGRAPH_HPP

#include <lib/util/Logging.hpp>
#include <lib/util/snapshot/Snapshotter.hpp>

namespace gcm {
	/**
	 * For 2D binary seismography for inverse problem
	 */
	template<class TMesh>
	class Binary2DSeismograph : public Snapshotter {
	public:
		typedef typename TMesh::Model            Model;
		typedef typename Model::PdeVariables     PdeVariables;
		typedef GetSetter<PdeVariables>          GETSETTER;
		typedef typename GETSETTER::Getter       Getter;
		
		const std::string FILE_EXTENSION = std::string("bin");
		const std::string FOLDER_NAME    = std::string("2dseismo");

	protected:
		virtual void initializeImpl(const Task& task) override;
		virtual void beforeStatementImpl(const Statement&) override;
		virtual void snapshotImpl(const AbstractGrid* grid, const int step) override;
		virtual void afterStatement() override;

	private:
		USE_AND_INIT_LOGGER("gcm.Binary2DSeismograph")
		std::vector<precision> surface; // values storage on current time step to write
		Getter valuesGetter;
		size_t sizeY;
		real hY;
		
		std::ofstream fileStream;

		void writeHeadOfTable();
	};
}

#endif // LIBGCM_BINARY2DSEISMOGRAPH_HPP
