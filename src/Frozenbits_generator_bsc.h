/*!
 * \file
 * \brief Class tools::Frozenbits_generator_5G.
 */
#ifndef FROZENBITS_GENERATOR_BSC_HPP_
#define FROZENBITS_GENERATOR_BSC_HPP_

#include <aff3ct.hpp>

namespace aff3ct
{
	namespace tools
	{
		class Frozenbits_generator_bsc : public Frozenbits_generator_file
		{
		private:
			const int m;
			const int N_max = 65536;

		public:
			Frozenbits_generator_bsc(const int K, const int N);

			~Frozenbits_generator_bsc();

			virtual Frozenbits_generator_bsc* clone() const;

		private:
			void evaluate();//store best_chn from file

		};
	}
}

#endif /* FROZENBITS_GENERATOR_BSC_HPP_ */
