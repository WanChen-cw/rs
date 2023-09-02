#include "framegenerate.h"
#pragma comment(lib, "ws2_32.lib")


template <typename T>
void read_File(const std::string& filename, std::vector<T>& ep) {
	std::ifstream inputFile(filename);
	if (!inputFile.is_open()) {
		std::cout << "Failed to open channel parameter file: " << filename << std::endl;
		return;
	}

	float value;
	while (inputFile >> value) {
		ep.push_back(value);
	}

	inputFile.close();
}

int main() {
	//parameters init
	int Frame_size = 1000;		//number of  frames      (byte)
	int FrameLength = 2226;		//length of  frame 
	int datelength = 2213;		//length of  data
	int frame_itl_number = 128;	//Number of interleaved frames
	int size_rs_N = 255;		//parameters of RS
	int size_rs_T = 16;
	int size_rs_K = size_rs_N - 2 * size_rs_T;
	int vary_chl_bit = 1000;	//Every x bit channel condition change
	int fe = 100;				//Number of errored frames (simulation stop condition) (not used)
	int seed = 0;				//random number seed

	int window_size = 2560;
	//channel parameters ep reading
	std::vector<float> ep;
	std::string channel_filename = "../conf/channel/channelfile.txt";
	read_File(channel_filename, ep);

	//buffer
	int m = (int)std::ceil(std::log2(size_rs_N));				//Each byte of data requires m bits of binary representation
	int bit_length_source= m * FrameLength * frame_itl_number;	//The bit length of each transmission
	int number_rs = (FrameLength * frame_itl_number + size_rs_K - 1) / size_rs_K;//The number of rs codes required
	int bit_length_transmission = m * size_rs_N * number_rs;	//
	std::vector<int  > ref_bits =	std::vector<int >	(bit_length_source);
	std::vector<int  > enc_bits =	std::vector<int >	(bit_length_transmission);
	std::vector<int  > itl_bits =	std::vector<int  >	(bit_length_transmission);
	std::vector<float> LLRs		=	std::vector<float>	(bit_length_transmission);
	std::vector<float> itl_LLRs =	std::vector<float>	(bit_length_transmission);
	std::vector<int  > dec_bits =	std::vector<int  >	(m * FrameLength * frame_itl_number);
	std::vector<int>	sub_itl_bits;
	std::vector<float>	sub_symbols = std::vector<float>(vary_chl_bit);
	std::vector<float>	sub_noisy_symbols = std::vector<float>(vary_chl_bit);
	std::vector<float>	sub_LLRs = std::vector<float>(vary_chl_bit);
	std::vector<float>	current_ep(1);
	std::vector<float>	current_ep1(1,0.0034);

	//module init
	std::unique_ptr<framegenerate<>>				source = std::unique_ptr<framegenerate <>>(new framegenerate <>(FrameLength, datelength));
	tools::RS_polynomial_generator					GF_poly(next_power_of_2(size_rs_N) - 1, size_rs_T);
	std::unique_ptr<module::CRC<>>					crc = std::unique_ptr<module::CRC_polynomial<>>(new module::CRC_polynomial<>((FrameLength - 4) * 8, "32-GZIP"));
	std::unique_ptr<module::Encoder<>>				encoder = std::unique_ptr<module::Encoder<>>(new module::Encoder_RS<>(size_rs_K, size_rs_N, GF_poly));
	std::unique_ptr<tools::Interleaver_core<>>		itl_core = std::unique_ptr<tools::Interleaver_core <>>(new tools::Interleaver_core_random<>(bit_length_transmission));
	std::unique_ptr<module::Interleaver<>>			itl1 = std::unique_ptr<module::Interleaver <>>(new module::Interleaver<>(*itl_core));
	std::unique_ptr<module::Modem<>>				modem1 = std::unique_ptr<module::Modem<>>(new module::Modem_OOK_BSC   <>(vary_chl_bit));
	std::unique_ptr<module::Channel<>>				channel = std::unique_ptr<module::Channel<>>(new module::Channel_binary_symmetric<>(vary_chl_bit));
	channel->set_seed(seed);
	int remaining_itl_bits = itl_bits.size() % vary_chl_bit;
	std::unique_ptr<module::Modem<>>	modem2;
	std::unique_ptr<module::Channel<>>	channe2;
	if (remaining_itl_bits > 0) {
			modem2 = std::unique_ptr<module::Modem<>>(new module::Modem_OOK_BSC   <>(remaining_itl_bits));
			channe2 = std::unique_ptr<module::Channel<>>(new module::Channel_binary_symmetric<>(remaining_itl_bits));
			channe2->set_seed(seed);
	}
	std::unique_ptr<module::Interleaver<float>>		itl2 = std::unique_ptr<module::Interleaver <float>>(new module::Interleaver<float>(*itl_core));
	std::unique_ptr<module::Decoder_RS_std<>>		decoder = std::unique_ptr<module::Decoder_RS_std<>>(new module::Decoder_RS_std<>(size_rs_K, size_rs_N, GF_poly));
	std::unique_ptr<module::Monitor_BFER<>>			monitor = std::unique_ptr<module::Monitor_BFER	<>>(new module::Monitor_BFER<>(8 * FrameLength, fe));
	//tools
	std::vector<std::unique_ptr<tools::Reporter>>		reporters; // list of reporters dispayed in the terminal
	std::unique_ptr<tools::Terminal_std>				terminal;  // manage the output text in the terminal
	reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_BFER<>(*monitor)));			// report the bit/frame error rates
	reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_throughput<>(*monitor)));		// report the simulation throughputs
	//terminal = std::unique_ptr<tools::Terminal_std>(new tools::Terminal_std(reporters));					// create a terminal that will display the collected data from the reporters
	//terminal->legend();																						// display the legend in the terminal
	//terminal->start_temp_report();
	//-------------------------------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------------------------------- 
	// socket通信
	WSADATA wsaData;
	WSAStartup(MAKEWORD(2, 2), &wsaData);
	SOCKET clientSocket;
	sockaddr_in serverAddr;
	clientSocket = socket(AF_INET, SOCK_STREAM, 0);
	memset(&serverAddr, 0, sizeof(serverAddr));
	serverAddr.sin_family = AF_INET;
	serverAddr.sin_port = htons(12345); // Server port
	serverAddr.sin_addr.s_addr = inet_addr("127.0.0.1"); // Server IP address
	connect(clientSocket, (struct sockaddr*)&serverAddr, sizeof(serverAddr));
	//-------------------------------------------------------------------------------------------------------------------
	std::vector<char>	llrbyteData(LLRs.size() * sizeof(float));			//将llr数据浮点数向量转换为字节流便于发送
	std::vector<uint32_t>feedback_frame_id = std::vector<uint32_t>(1, 0);	//接收到的最后反馈帧号
	uint32_t frame_id_cur = 0;												//当前发送帧号
	uint32_t frame_id_next = 1;												//已发送的最大传输帧号id加一
	uint32_t feedback_frame_id_last = 0;									//上一个反馈帧号
	bool isReliableFrame = 0;
	int error = 0;															//当前是否错误
	std::vector<char> frame_id_recv;										//反馈
	//------------------------------------------------------------------------------------------------------------------- 
	//sumulation
	std::cout << "client! " << std::endl;
	std::ofstream outFile("../conf/result/frane_id.txt"); // 文件用于写入
	while (feedback_frame_id[0] < Frame_size) {
		for (int i = 0; i < frame_itl_number; ++i) {
			if (error == 0){
			frame_id_cur = frame_id_next;
			frame_id_next++;
			}
		else {
			if (feedback_frame_id[0] != feedback_frame_id_last){
				frame_id_cur = frame_id_next;
				frame_id_next++;
				feedback_frame_id_last = feedback_frame_id[0];
				error = 0;
			}
		}
		if ((frame_id_cur - feedback_frame_id[0]) > window_size){
			error = 1;
			frame_id_next = frame_id_cur;
			frame_id_cur = feedback_frame_id[0];
		}
		if (outFile.is_open()) { // 确保文件成功打开
			outFile << frame_id_cur<< std::endl;
		}
		int startIdx = i * m * FrameLength;
		int endIdx = (i + 1) * m * FrameLength;
		std::vector<int> segment = std::vector<int  >(m * FrameLength);
		source->generate(segment, frame_id_cur++, 0, isReliableFrame);
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
		int num_subvectors = itl_bits.size() / vary_chl_bit;
		int current_ep_index = 0;
		for (int i = 0; i < num_subvectors; i++) {
		sub_itl_bits.assign(itl_bits.begin() + i * vary_chl_bit, itl_bits.begin() + (i + 1) * vary_chl_bit);
		current_ep[0]= ep[current_ep_index];
		modem1->modulate(sub_itl_bits, sub_symbols);
		channel->add_noise(current_ep, sub_symbols, sub_noisy_symbols);
		modem1->demodulate(current_ep1, sub_noisy_symbols, sub_LLRs);
		// Copy the noisy subvector back to the main noisy_symbols vector
		std::copy(sub_LLRs.begin(), sub_LLRs.end(), LLRs.begin() + i * vary_chl_bit);
		// Move to the next ep value in the ep vector
		current_ep_index = (current_ep_index + 1) % ep.size(); // 循环使用ep值
		}
		if (remaining_itl_bits > 0) {
		std::vector<int>	remaining_sub_itl_bits(itl_bits.end() - remaining_itl_bits, itl_bits.end());
		std::vector<float>	remaining_symbols = std::vector<float>(remaining_itl_bits);
		std::vector<float>	remaining_noisy_symbols = std::vector<float>(remaining_itl_bits);
		std::vector<float>	remaining_LLRs = std::vector<float>(remaining_itl_bits);
		modem2->modulate(remaining_sub_itl_bits, remaining_symbols);
		channe2->add_noise(current_ep, remaining_symbols, remaining_noisy_symbols);
		modem2->demodulate(current_ep1, remaining_noisy_symbols, remaining_LLRs);
		// Copy the noisy remaining subvector back to the main noisy_symbols vector
		std::copy(remaining_LLRs.begin(), remaining_LLRs.end(), LLRs.end() - remaining_itl_bits);
	}
	//-------------------------------------------------------------------------------------------------------------------
		//（此次仿真发送解除调制数据）
		memcpy(llrbyteData.data(), LLRs.data(), llrbyteData.size());
		send(clientSocket, llrbyteData.data(), llrbyteData.size(), 0);
	//-------------------------------------------------------------------------------------------------------------------
		frame_id_recv.resize(frame_itl_number * sizeof(uint32_t));
		recv(clientSocket, frame_id_recv.data(), frame_id_recv.size(), 0);
		memcpy(feedback_frame_id.data(), frame_id_recv.data(), frame_id_recv.size());
		std::cout << std::endl << "feedback_frame_id: " << feedback_frame_id[0];
	}
	//-------------------------------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------------------------------- 
	outFile.close();			// 关闭文件frame_id
	//terminal->final_report();	// display the performance (BER and FER) in the terminal
	monitor->reset();			// reset the monitor for the next SNR
	//terminal->reset();
	closesocket(clientSocket);	//断开连接
	WSACleanup();
	return 0;
}


//#include "simulationRS.hpp"
//#include "simulationpolar.hpp"
//int main(int argc, char** argv)
//{
//	//testrs();
//	testpolar();
//}
