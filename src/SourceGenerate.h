/*

*/
#ifndef SOURCEGENERATE_HPP_
#define SOURCEGENERATE_HPP_

#include <random>
#include <aff3ct.hpp>
using namespace aff3ct::module;

		template <typename B = int>
		class sourcegenerate : public Source_random<B>
		{
		private:
			
		public:
			sourcegenerate(const int K, const int seed = 0);
		protected:
			
		};


#endif /*SOURCEGENERATE_HPP_*/
