#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <fstream>
#include <aff3ct.hpp>
#include "framegenerate.h"
using namespace aff3ct;
using namespace aff3ct::module;
void testrs()
{
	int Frame_size = 500000;
	int FrameLength = 2226;//byte
	int datelength = 2213;	//byte
	int frame_itl_number = 120;
	int size_itl;	//bit
	int size_rs_N = 255;
	int size_rs_T = 16;
	int size_rs_K = size_rs_N - 2 * size_rs_T;
	int fe = 100;
	int seed=0;
	std::vector<float> ep;
	std::string filename = "data.txt"; // 替换为实际的文件路径
	std::ifstream inputFile(filename); 
	if (!inputFile.is_open()) {
		std::cout << "Failed to open file: " << filename << std::endl;
	}
	float value;
	while (inputFile >> value) {
		ep.push_back(value);
	}
	inputFile.close();
	int m=(int)std::ceil(std::log2(size_rs_N));
	//buffer
	std::vector<int  > ref_bits = std::vector<int  >(m * FrameLength* frame_itl_number);
	int number_rs = (FrameLength * frame_itl_number + size_rs_K - 1) / size_rs_K;
	std::vector<int  > enc_bits = std::vector<int  >(m * size_rs_N * number_rs);
	size_itl = m * size_rs_N * number_rs;
	std::vector<int  > itl_bits = std::vector<int  >(size_itl);
	std::vector<float> LLRs = std::vector<float>(size_itl);
	std::vector<float> itl_LLRs = std::vector<float>(size_itl);
	std::vector<int  > dec_bits = std::vector<int  >(m * FrameLength * frame_itl_number);

	//tool using 
	std::vector<std::unique_ptr<tools::Reporter>> reporters; // list of reporters dispayed in the terminal
	std::unique_ptr<tools::Terminal_std>          terminal;  // manage the output text in the terminal

	//source
	std::unique_ptr<framegenerate<>> source = std::unique_ptr<framegenerate <>>(new framegenerate <>(FrameLength, datelength));
	uint32_t frame_id = 3;
	uint32_t feedback_frame_id = 7;
	bool isReliableFrame = 0;

	//RS 
	tools::RS_polynomial_generator GF_poly(next_power_of_2(size_rs_N) - 1, size_rs_T);
	std::unique_ptr<module::Encoder<>>					encoder = std::unique_ptr<module::Encoder<>>(new module::Encoder_RS<>(size_rs_K, size_rs_N, GF_poly));

	//interleaver
	std::unique_ptr<tools::Interleaver_core<>> itl_core = std::unique_ptr<tools::Interleaver_core <>>(new tools::Interleaver_core_random<>(size_itl));
	std::unique_ptr<module::Interleaver<>>  itl1 = std::unique_ptr<module::Interleaver <>>(new module::Interleaver<>(*itl_core));

	//deinterleave
	std::unique_ptr<module::Interleaver<float>>  itl2 = std::unique_ptr<module::Interleaver <float>>(new module::Interleaver<float>(*itl_core));

	//decode
	std::unique_ptr<module::Decoder_RS_std<>> decoder = std::unique_ptr<module::Decoder_RS_std<>>(new module::Decoder_RS_std<>(size_rs_K, size_rs_N, GF_poly));

	//monitor
	std::unique_ptr<module::Monitor_BFER<>>	monitor = std::unique_ptr<module::Monitor_BFER	<>>(new module::Monitor_BFER<>(8 * FrameLength, fe));
	
	//tools
	// report the bit/frame error rates
	reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_BFER<>(*monitor)));
	// report the simulation throughputs
	reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_throughput<>(*monitor)));
	// create a terminal that will display the collected data from the reporters
	terminal = std::unique_ptr<tools::Terminal_std>(new tools::Terminal_std(reporters));
	// display the legend in the terminal
	terminal->legend();
	//monitor using
	//Modile Channel
	int subvector_size = 2000; // 可以根据需要调整子向量的大小,每多少bit信道特征变化
	std::unique_ptr<module::Modem<>>	modem1 = std::unique_ptr<module::Modem<>>(new module::Modem_OOK_BSC   <>(subvector_size));
	std::unique_ptr<module::Channel<>>	channel = std::unique_ptr<module::Channel<>>(new module::Channel_binary_symmetric<>(subvector_size));
	channel->set_seed(seed);
	int remaining_itl_bits = itl_bits.size() % subvector_size;
	std::unique_ptr<module::Modem<>>	modem2;
	std::unique_ptr<module::Channel<>>	channe2;
	if (remaining_itl_bits > 0) {
		modem2 = std::unique_ptr<module::Modem<>>(new module::Modem_OOK_BSC   <>(remaining_itl_bits));
		channe2 = std::unique_ptr<module::Channel<>>(new module::Channel_binary_symmetric<>(remaining_itl_bits));
		channe2->set_seed(seed);
	}
	std::unique_ptr<module::CRC<>>crc = std::unique_ptr<module::CRC_polynomial<>>(new module::CRC_polynomial<>((FrameLength - 4) * 8, "32-GZIP"));
	terminal->start_temp_report();
	// run the simulation chain
	while (!monitor->fe_limit_achieved() && !terminal->is_interrupt())
	{
		for (int i = 0; i < frame_itl_number; ++i) {
		int startIdx = i * m * FrameLength;
		int endIdx = (i + 1) *m * FrameLength;
		std::vector<int> segment= std::vector<int  >(m * FrameLength);
		source->generate(segment, frame_id, feedback_frame_id, isReliableFrame);
		for (size_t j = 0; j < segment.size(); ++j) {
			ref_bits[i * segment.size() + j] = segment[j];
		}
		}

		for (int i = 0; i < number_rs; ++i) {
			int startIdx = i * size_rs_K * 8;
			int endIdx = std::min((i + 1) * size_rs_K * 8, static_cast<int>(ref_bits.size()));
			std::vector<int> segment(ref_bits.begin() + startIdx, ref_bits.begin() + endIdx);
			if (endIdx < (i + 1) * size_rs_K * 8) {
				int zerosToAdd = (i + 1) * size_rs_K * 8 - endIdx;
				segment.insert(segment.end(), zerosToAdd, 0);
			}
			std::vector<int> enc_segment(8 * size_rs_N);
			encoder->encode(segment, enc_segment);
			for (size_t j = 0; j < enc_segment.size(); ++j) {
				enc_bits[i * size_rs_N * 8 + j] = enc_segment[j];
				//std::cout << "Element " << startIdx + j << ": " << enc_segment[j] << std::endl;
			}
		}
		itl1->interleave(enc_bits, itl_bits);
		//itl_bits = enc_bits;
		int num_subvectors = itl_bits.size() / subvector_size;
		int current_ep_index = 0;
		for (int i = 0; i < num_subvectors; i++) {
			
			// Extract a subvector containing subvector_size symbols
			std::vector<int> sub_itl_bits(itl_bits.begin() + i * subvector_size, itl_bits.begin() + (i + 1) * subvector_size);
			std::vector<float> sub_symbols = std::vector<float>(subvector_size);
			std::vector<float> sub_noisy_symbols = std::vector<float>(subvector_size);
			std::vector<float> sub_LLRs = std::vector<float>(subvector_size);
			// Get the current ep value from the ep vector
			float current_ep_value = ep[current_ep_index];
			// Call channel->add_noise with the current ep value and the subvector
			std::vector<float> current_ep(1, current_ep_value);
			modem1->modulate(sub_itl_bits, sub_symbols);
			channel->add_noise(current_ep, sub_symbols, sub_noisy_symbols);
			modem1->demodulate(current_ep, sub_noisy_symbols, sub_LLRs);
			// Copy the noisy subvector back to the main noisy_symbols vector
			std::copy(sub_LLRs.begin(), sub_LLRs.end(), LLRs.begin() + i * subvector_size);

			// Move to the next ep value in the ep vector
			current_ep_index = (current_ep_index + 1) % ep.size(); // 循环使用ep值
		}

		if (remaining_itl_bits > 0) {

		
			std::vector<int> remaining_sub_itl_bits(itl_bits.end() - remaining_itl_bits, itl_bits.end());
			// Extract a subvector containing subvector_size symbols
			std::vector<float> remaining_symbols = std::vector<float>(remaining_itl_bits);
			std::vector<float> remaining_noisy_symbols = std::vector<float>(remaining_itl_bits);
			std::vector<float> remaining_LLRs = std::vector<float>(remaining_itl_bits);
			// Get the current ep value from the ep vector
			float current_ep_value = ep[current_ep_index];
			// Call channel->add_noise with the current ep value and the remaining subvector
			std::vector<float> current_ep(1, current_ep_value); // 创建包含一个ep值的容器
			modem2->modulate(remaining_sub_itl_bits, remaining_symbols);
			channe2->add_noise(current_ep, remaining_symbols, remaining_noisy_symbols);
			modem2->demodulate(current_ep, remaining_noisy_symbols, remaining_LLRs);
			// Copy the noisy remaining subvector back to the main noisy_symbols vector
			std::copy(remaining_LLRs.begin(), remaining_LLRs.end(), LLRs.end() - remaining_itl_bits);
		}


		itl2->deinterleave(LLRs, itl_LLRs);
		//itl_LLRs = LLRs;
		for (int i = 0; i < number_rs; ++i) {
			int startIdx = i * size_rs_N * 8;
			int endIdx = std::min((i + 1) * size_rs_N * 8, static_cast<int>(itl_LLRs.size()));
			std::vector<float> segment(itl_LLRs.begin() + startIdx, itl_LLRs.begin() + endIdx);
			std::vector<int> dec_segment(8 * size_rs_K);
			decoder->decode_siho(segment, dec_segment);
			for (size_t j = 0; j < dec_segment.size(); ++j) {
				if (i * size_rs_K * 8 + j >= dec_bits.size())
					break;
				else
				dec_bits[i * size_rs_K * 8 + j] = dec_segment[j];
			}

		}

		for (int i = 0; i < frame_itl_number; ++i) {
			int startIdx = i * m * FrameLength;
			int endIdx = (i + 1) * m * FrameLength;
			//std::vector<int> segment = std::vector<int  >(m * (FrameLength));
			std::vector<int> segmentdec = std::vector<int  >(m *(FrameLength));
			for (size_t j = 0; j < segmentdec.size(); ++j) {
				//segment[j] = ref_bits[i * segment.size() + j];
				segmentdec[j] = dec_bits[i * segmentdec.size() + j];
			}
			//monitor->check_errors(dec_bits, ref_bits);
			std::vector<int> segment = std::vector<int  >(m * (FrameLength-4));
			std::vector<int> segmentcrc = std::vector<int  >(m * (FrameLength));
			std::copy(std::begin(segmentdec), std::end(segmentdec)-32, std::begin(segment));
			crc->build(segment, segmentcrc);
			monitor->check_errors(segmentcrc, segmentdec);
		}
		
	}
	// display the performance (BER and FER) in the terminal
	terminal->final_report();
	// reset the monitor for the next SNR
	monitor->reset();
	terminal->reset();
}