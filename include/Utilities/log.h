#ifndef LOG_H
#define LOG_H

#include <QString>
#include <iostream>

/**
 * @author Andrea Salvadori
 */
namespace SNS { namespace Utilities
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
void throwAndPrintError(const QString& errorMessage)
{
	std::cerr << errorMessage.toStdString() << std::endl;
	throw ExceptionType(errorMessage.toStdString());
}

} }

#endif // LOG_H
