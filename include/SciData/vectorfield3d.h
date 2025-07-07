#ifndef VECTORFIELD3D_H
#define VECTORFIELD3D_H

#include <SciData/datagrid3d.h>
#include <Utilities/boost_math_utilities.h>

/**
 * @author Andrea Salvadori
 */
namespace CFF  { namespace SciData
{
	/**
	 * @brief Uniform 3D grid of 3D vectors.
	 *		  See the documentation of the CFF::SciData::DataGrid3D class for details.
	 */
	typedef DataGrid3D<CFF::Utilities::Vector3d> VectorField3D;

	typedef std::shared_ptr<VectorField3D> VectorField3DSP;
	typedef std::shared_ptr<const VectorField3D> ConstVectorField3DSP;
} }

#endif // VECTORFIELD3D_H
