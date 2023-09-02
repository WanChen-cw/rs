/*

*/
#ifndef SIMULATIONPOLAR_HPP_
#define SIMULATIONPOLAR_HPP_
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <aff3ct.hpp>
#include "Frozenbits_generator_bsc.h"
using namespace aff3ct;
void testpolar()
{
	//parament
	int   N = 65536;     // codeword size
	int   K = 50036;     // number of information bits
	int   fe = 100;     // number of frame errors
	int   seed = 0;     // PRNG seed for the AWGN channel
	float ber_min = 0.02f; // minimum  value
	float ber_max = 0.02f; // maximum  value
	float ber_step = 0.02f; // SNR step
	float R;                   // code rate (R=K/N)
	R = (float)K / (float)N;

	//module
	std::unique_ptr<module::Source_random<>>				source;
	std::unique_ptr<module::Encoder_polar<>>					encoder;
	std::unique_ptr<module::Modem_OOK_BSC<>>				modem;
	std::unique_ptr<module::Channel_binary_symmetric<>>		channel;
	std::unique_ptr<module::Decoder_polar_SC_naive<>>				decoder;
	std::unique_ptr<module::Monitor_BFER<>>					monitor;

	// module init 
	std::vector<bool> fb ;
	fb = std::vector<bool  >(N);
	tools::Frozenbits_generator_bsc fbg(K, N);
	fbg.generate(fb);
	//std::string filename = "65536.txt"; // 替换为实际的文件路径
	//std::ifstream inputFile(filename); // 创建文件输入流对象
	/*if (!inputFile.is_open()) {
		std::cout << "Failed to open file: " << filename << std::endl;
	}
	float value;
	while (inputFile >> value) {
		fb.push_back(value);
	}
	inputFile.close();*/

	int L = 16;
//std::vector<bool>&frozen_bits
	source = std::unique_ptr<module::Source_random				<>>(new module::Source_random			<>(K));
	encoder = std::unique_ptr<module::Encoder_polar				<>>(new module::Encoder_polar		<>(K, N, fb));
	modem = std::unique_ptr<module::Modem_OOK_BSC				<>>(new module::Modem_OOK_BSC           <>( N));
	channel = std::unique_ptr<module::Channel_binary_symmetric  <>>(new module::Channel_binary_symmetric<>( N));
	decoder = std::unique_ptr<module::Decoder_polar_SC_naive			<>>(new module::Decoder_polar_SC_naive			<>(K,N,fb));
	monitor = std::unique_ptr<module::Monitor_BFER				<>>(new module::Monitor_BFER			<>( K, fe));
	channel->set_seed(seed);

	//buffer init 
	std::vector<int  > ref_bits;
	std::vector<int  > enc_bits;
	std::vector<float> symbols;
	std::vector<float> ep;
	std::vector<float> noisy_symbols;
	std::vector<float> LLRs;
	std::vector<int  > dec_bits;
	ref_bits = std::vector<int  >(K);
	enc_bits = std::vector<int  >( N);
	symbols = std::vector<float>( N);
	ep = std::vector<float>(1);
	noisy_symbols = std::vector<float>(N);
	LLRs = std::vector<float>( N);
	dec_bits = std::vector<int  >(K);

	//tool using 
	std::unique_ptr<tools::Interleaver_core_random<>>         itl;     // a interleaver noise type
	std::unique_ptr<tools::Event_probability<>>               noise;     // a sigma noise type
	std::vector<std::unique_ptr<tools::Reporter>> reporters; // list of reporters dispayed in the terminal
	std::unique_ptr<tools::Terminal_std>          terminal;  // manage the output text in the terminal

	itl = std::unique_ptr<tools::Interleaver_core_random<>>(new tools::Interleaver_core_random<>(8 * N));
	// create a sigma noise type
	noise = std::unique_ptr<tools::Event_probability<>>(new tools::Event_probability<>());
	reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_noise<>(*noise)));
	// report the bit/frame error rates
	reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_BFER<>(*monitor)));
	// report the simulation throughputs
	reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_throughput<>(*monitor)));
	// create a terminal that will display the collected data from the reporters
	terminal = std::unique_ptr<tools::Terminal_std>(new tools::Terminal_std(reporters));
	// display the legend in the terminal
	terminal->legend();


	// loop over the ep
	for (auto ber = ber_max; ber >= ber_min; ber -= ber_step)
	{
		// compute the current sigma for the channel noise
		std::fill(ep.begin(), ep.end(), ber);
		noise->set_value(ber);
		// display the performance (BER and FER) in real time (in a separate thread)
		terminal->start_temp_report();

		// run the simulation chain
		while (!monitor->fe_limit_achieved() && !terminal->is_interrupt())
		{
			source->generate(ref_bits);
			encoder->encode(ref_bits, enc_bits);
			modem->modulate(enc_bits, symbols);
			channel->add_noise(ep, symbols, noisy_symbols);
			modem->demodulate(ep, noisy_symbols, LLRs);
			decoder->decode_siho(LLRs, dec_bits);
			monitor->check_errors(dec_bits, ref_bits);
		}

		// display the performance (BER and FER) in the terminal
		terminal->final_report();
		// reset the monitor for the next SNR
		monitor->reset();
		terminal->reset();
		// if user pressed Ctrl+c twice, exit the SNRs loop
		if (terminal->is_over()) break;
	}

}
#endif /*SIMULATIONPOLKAR_HPP_*/


