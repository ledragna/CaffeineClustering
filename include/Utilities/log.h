#ifndef LOG_H
#define LOG_H

#include <iostream>
#include <string>

/**
 * @author Andrea Salvadori
 */
namespace CFF { namespace Utilities
{

/**
 * @brief throwAndPrintError
 *
 * @param errorMessage  The message to print to the error console.
 *                      It will also be passed as parameter to the exception.
 *
 * @param ExceptionType The type of the exception to throw.
 *                      Must have a constructor accepting an std::string as parameter.
 */
template<typename ExceptionType>
void throwAndPrintError(const std::string& errorMessage)
{
	std::cerr << errorMessage << std::endl;
	throw ExceptionType(errorMessage);
}

} }

#endif // LOG_H
