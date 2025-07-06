#ifndef ITERABLE_PRIORITY_QUEUE_H
#define ITERABLE_PRIORITY_QUEUE_H

#include <queue>

/**
 * @author Andrea Salvadori
 */
namespace SNS { namespace Utilities
{

	/**
	 * @brief	Like std::priority_queue, but iterable.
	 *			Use it at your own risk!!!
	 */
	template<typename T,
			 typename Container = std::vector<T>,
			 typename Compare = std::less<typename Container::value_type>>
	class MyPriorityQueue: public std::priority_queue<T,Container,Compare>
	{
	public:
		// 'c' is the protected field related to the wrapped container of std::priority_queue
		typename Container::const_iterator begin() const { return std::priority_queue<T,Container,Compare>::c.cbegin(); }
		typename Container::const_iterator end() const { return std::priority_queue<T,Container,Compare>::c.cend(); }
		typename Container::iterator begin() { return std::priority_queue<T,Container,Compare>::c.begin(); }
		typename Container::iterator end()	{ return std::priority_queue<T,Container,Compare>::c.end(); }
	};

} }

#endif // ITERABLE_PRIORITY_QUEUE_H
